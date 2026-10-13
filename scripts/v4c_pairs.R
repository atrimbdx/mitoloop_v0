#!/usr/bin/env Rscript
# v4c_pairs.R - are the loops at mitochondrial promoters the same in two
# metabolic states of one cell type?
#
# Purpose
#   For each experiment of a paired design (two states of one cell type, plus =
#   higher mitochondrial content), compares the peakC loops of the two states
#   at mitochondrial and at control promoters, against the technical
#   variability measured within each state.
#
# Inputs
#   --libraries   libraries_pairs.tsv from `v4c_pairs.py export`; columns
#                 label, experiment, state (and genome, resolution,
#                 depth_control, depth_mito, n_viewpoints, not used here)
#   --export      directory of <label>.tsv.gz from `v4c_pairs.py export`:
#                 columns name, set (mito | control), chrom, pos, then one
#                 column of raw counts per offset in bp
#
# Outputs (in --out)
#   pairs_viewpoints.tsv  one row per experiment and viewpoint: experiment,
#                         name, set, chrom, pos, within_minus, within_plus,
#                         between_1, between_2 (Jaccard of the peak sets),
#                         n_peaks_minus, n_peaks_plus, within, between
#   pairs_peaks.tsv       one row per union peak: name, summit, distance,
#                         called_minus, called_plus, count_minus, count_plus,
#                         log2_ratio, log2_within_minus, log2_within_plus,
#                         log2_between_1, log2_between_2, experiment, set,
#                         chrom, vp_pos
#   pairs_summary.tsv     one row per experiment and set: n_viewpoints,
#                         concordance_within, concordance_between,
#                         concordance_loss, n_peaks, shared,
#                         amplitude_between, amplitude_within, noise_q95,
#                         excess, direction
#   pairs_meta.tsv        one row per experiment, the statistics of both sets
#                         side by side and their differences mito minus
#                         control (d_concordance_loss, d_excess, d_direction)
#   pairs.svg             concordance and excess, one point per experiment
#
# Method, for one experiment
#   Replicates of a state are pooled. The deeper state is then thinned to the
#   depth of the other by binomial thinning (each contact kept with
#   probability f, the ratio of depths). With --thin-on control (default) the
#   depths are taken on the control viewpoints, so that mitochondrial
#   viewpoints that genuinely gain contacts keep them; `all` uses every
#   viewpoint, `none` skips the thinning.
#
#   Each pooled state is split in two halves by a binomial draw with
#   probability 1/2, which gives two independent halves of a Poisson count.
#   That gives four profiles per viewpoint, all of the same depth:
#
#       A1 A2   halves of mito-      B1 B2   halves of mito+
#
#   within  = A1 vs A2 and B1 vs B2     (same state, technical difference)
#   between = A1 vs B1 and A2 vs B2     (the two states)
#
#   Peaks are called in each half separately (peakC single.analysis, one
#   window size --wsize for everything) and two peaks match when their
#   summits lie within --tol. The comparison is then
#
#       concordance   Jaccard of the peak sets, within against between
#       excess        share of loops whose |log2| contact ratio between the
#                     states exceeds the --excess-quantile quantile of the
#                     |log2| taken between two halves of one state. With no
#                     state-dependent loop it returns 1 - --excess-quantile
#                     in expectation.
#       direction     mean log2(plus / minus) at the union peaks
#
#   Because the four profiles carry the same number of contacts, `within` is
#   the technical floor of `between` and no depth correction is needed after
#   it. A state-specific loop shows up as between < within, and as an excess
#   above 1 - --excess-quantile. A median |log2| would not detect a minority
#   of changing loops, which leaves it untouched; hence the excess statistic.
#
#   The same is computed on control promoters. Whether mitochondrial
#   promoters behave differently is read on the difference of differences,
#   mito minus control, and combined over the experiments with a sign test.
#   One binomial split is drawn per run; --seed sets the random draws of the
#   thinning, the split and the --n-vp subsample, and repeating the run with
#   other seeds gives replicate splits.
#
# Usage
#   Rscript scripts/v4c_pairs.R compare \
#       --libraries results/export_pairs/libraries_pairs.tsv \
#       --export results/export_pairs --wsize 3 --out results/pairs_out \
#       [--cores 8]
#   Rscript scripts/v4c_pairs.R --help     (all options and defaults)
#
# References
#   Geeven et al. 2018. peakC: a flexible, non-parametric peak calling package
#     for 4C and Capture-C data. Nucleic Acids Research.
#   Jaccard 1912. The distribution of the flora in the alpine zone. New
#     Phytologist. (Jaccard index)
#   Packages: peakC (github deWitLab/peakC), isotone, data.table, parallel.

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})
if (requireNamespace("isotone", quietly = TRUE))
  suppressPackageStartupMessages(library(isotone))

# ---------------------------------------------------------------------------
# arguments

defaults <- list(cores = 4L, wsize = 3, dist_min = 50e3, dist_max = 1e6,
                 qwr = 1, qwd = 1.5, tol = 20e3, out = "pairs_out", seed = 1L,
                 min_bins = 20L, n_vp = NA_real_, min_count = 5,
                 plus_state = "plus", minus_state = "minus",
                 thin_on = "control", excess_quantile = 0.95,
                 libraries = NA_character_, export = NA_character_)
numeric_opts <- c("cores", "wsize", "dist_min", "dist_max", "qwr", "qwd",
                  "tol", "seed", "min_bins", "n_vp", "min_count",
                  "excess_quantile")

usage <- function() {
  cat(
"usage: Rscript v4c_pairs.R compare --libraries FILE --export DIR [options]

required
  --libraries        libraries_pairs.tsv written by `v4c_pairs.py export`
  --export           directory of the <label>.tsv.gz tables

options (default in brackets)
  --out              output directory [pairs_out]
  --wsize            peakC window size, in bins [3]
  --qwd, --qwr       peakC quantile thresholds on the difference and ratio
                     to the background [1.5, 1]
  --dist-min         peaks closer to the viewpoint are ignored by peakC, bp
                     [50000]
  --dist-max         peaks farther from the viewpoint are dropped, bp [1e6]
  --min-bins         viewpoints with fewer finite bins (half on each side)
                     are skipped [20]
  --tol              two summits match within this distance, bp [20000]
  --min-count        peaks with fewer than 2 x this many contacts (both
                     states together) are left out of excess, amplitude and
                     direction [5]
  --excess-quantile  quantile of the within-state |log2| used as the noise
                     threshold of `excess` [0.95]
  --thin-on          viewpoints whose depth sets the thinning factor:
                     control, all or none (no thinning) [control]
  --plus-state       state label of the higher mitochondrial content [plus]
  --minus-state      state label of the other state [minus]
  --n-vp             analyse a random subset of this many viewpoints per
                     experiment [all]
  --seed             seed of the thinning, the binomial split and --n-vp [1]
  --cores            parallel processes (mclapply) [4]
")
}

parse_args <- function(a) {
  opt <- c(list(cmd = a[1]), defaults)
  i <- 2
  while (i <= length(a)) {
    key <- gsub("-", "_", sub("^--", "", a[i]))
    if (!key %in% names(defaults)) stop("unknown option: ", a[i])
    if (i + 1 > length(a)) stop("option ", a[i], " needs a value")
    opt[[key]] <- a[i + 1]
    i <- i + 2
  }
  for (k in numeric_opts) opt[[k]] <- as.numeric(opt[[k]])
  if (!opt$thin_on %in% c("control", "all", "none"))
    stop("--thin-on must be control, all or none")
  if (!(opt$excess_quantile > 0 && opt$excess_quantile < 1))
    stop("--excess-quantile must lie between 0 and 1")
  for (k in c("libraries", "export"))
    if (is.na(opt[[k]])) stop("--", k, " is required")
  opt
}

# ---------------------------------------------------------------------------
# reading

read_export <- function(dir, labels) {
  out <- list()
  for (a in labels) {
    f <- file.path(dir, paste0(a, ".tsv.gz"))
    if (!file.exists(f)) { message("missing ", f); next }
    d <- fread(f, sep = "\t", na.strings = "NA", header = TRUE)
    m <- as.matrix(d[, -(1:4)])
    storage.mode(m) <- "double"
    rownames(m) <- d$name
    out[[a]] <- list(m = m, set = setNames(d$set, d$name),
                     chrom = setNames(d$chrom, d$name),
                     pos = setNames(d$pos, d$name),
                     offsets = as.numeric(colnames(d)[-(1:4)]))
  }
  out
}

# pooled counts of one state: replicates summed, NA where every replicate is NA
pool_state <- function(ex, labels) {
  ms <- lapply(labels, function(a) ex[[a]]$m)
  ms <- Filter(Negate(is.null), ms)
  if (length(ms) == 0) return(NULL)
  nm <- Reduce(intersect, lapply(ms, rownames))
  if (length(nm) == 0) return(NULL)
  acc <- matrix(0, length(nm), ncol(ms[[1]]), dimnames = list(nm, NULL))
  seen <- matrix(0L, length(nm), ncol(ms[[1]]))
  for (m in ms) {
    v <- m[nm, , drop = FALSE]
    ok <- is.finite(v)
    acc[ok] <- acc[ok] + v[ok]
    seen <- seen + ok
  }
  acc[seen == 0] <- NA_real_
  acc
}

# ---------------------------------------------------------------------------
# peakC

prep <- function(vec, offsets, vp, res, min_bins) {
  if (is.null(vec) || all(!is.finite(vec))) return(NULL)
  centre <- floor(vp / res) * res
  pos <- centre + offsets + res / 2
  keep <- is.finite(vec)
  if (sum(keep & pos < vp) < min_bins / 2 || sum(keep & pos > vp) < min_bins / 2)
    return(NULL)
  tot <- sum(vec[keep])
  if (tot <= 0) return(NULL)
  cbind(pos[keep], 1e6 * vec[keep] / tot)
}

run_peakc <- function(profile, vp, opt) {
  if (is.null(profile)) return(NULL)
  r <- tryCatch({
    x <- peakC::single.analysis(profile, vp.pos = vp, wSize = opt$wsize,
                                qWd = opt$qwd, qWr = opt$qwr,
                                minDist = opt$dist_min)
    list(peak = x$peak, pos = x$ratios[, 1], ratio = x$ratios[, 2])
  }, error = function(e) NULL)
  if (is.null(r)) return(NULL)
  r$peak <- r$peak[abs(r$peak - vp) <= opt$dist_max]
  r
}

# runs of consecutive significant bins -> one summit each
regions <- function(peak, pos, ratio, res) {
  if (length(peak) == 0) return(NULL)
  p <- sort(unique(peak))
  grp <- cumsum(c(1, diff(p) > res * 1.01))
  rbindlist(lapply(split(p, grp), function(x) {
    rr <- ratio[match(x, pos)]
    if (all(!is.finite(rr))) rr <- rep(1, length(x))
    data.table(summit = x[which.max(rr)], summit_ratio = max(rr, na.rm = TRUE))
  }))
}

# Jaccard between two summit sets, two summits matching within tol
jaccard <- function(a, b, tol) {
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  if (length(a) == 0 || length(b) == 0) return(0)
  used <- rep(FALSE, length(b))
  hit <- 0L
  for (x in a) {
    d <- abs(b - x)
    d[used] <- Inf
    j <- which.min(d)
    if (is.finite(d[j]) && d[j] <= tol) { used[j] <- TRUE; hit <- hit + 1L }
  }
  hit / (length(a) + length(b) - hit)
}

# contacts in the window around a summit, from a raw count vector
window_count <- function(vec, offsets, vp, summit, res, wsize) {
  centre <- floor(vp / res) * res
  pos <- centre + offsets + res / 2
  k <- which(abs(pos - summit) <= (wsize + 0.5) * res)
  if (length(k) == 0) return(NA_real_)
  v <- vec[k]
  if (all(!is.finite(v))) return(NA_real_)
  sum(v, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# one viewpoint of one experiment

do_viewpoint <- function(name, halves, full, offsets, vp, res, opt) {
  # halves: list(A1, A2, B1, B2) of raw count vectors, equal depth
  # full:   list(A, B) of raw count vectors, already depth matched
  pk <- lapply(halves, function(v) {
    r <- run_peakc(prep(v, offsets, vp, res, opt$min_bins), vp, opt)
    if (is.null(r)) return(NULL)
    regions(r$peak, r$pos, r$ratio, res)
  })
  s <- lapply(pk, function(x) if (is.null(x)) numeric(0) else x$summit)
  names(s) <- c("A1", "A2", "B1", "B2")

  jac <- c(within_minus  = jaccard(s$A1, s$A2, opt$tol),
           within_plus   = jaccard(s$B1, s$B2, opt$tol),
           between_1     = jaccard(s$A1, s$B1, opt$tol),
           between_2     = jaccard(s$A2, s$B2, opt$tol))

  # union of the summits seen in any half, merged within tol
  all_s <- sort(unlist(s, use.names = FALSE))
  peaks <- NULL
  if (length(all_s)) {
    grp <- cumsum(c(1, diff(all_s) > opt$tol))
    uni <- vapply(split(all_s, grp), function(x) x[ceiling(length(x) / 2)],
                  numeric(1))
    called_minus <- vapply(uni, function(u)
      any(abs(c(s$A1, s$A2) - u) <= opt$tol), logical(1))
    called_plus <- vapply(uni, function(u)
      any(abs(c(s$B1, s$B2) - u) <= opt$tol), logical(1))
    wc <- function(v) vapply(uni, function(u)
      window_count(v, offsets, vp, u, res, opt$wsize), numeric(1))
    cm <- wc(full$A); cp <- wc(full$B)
    # The four halves carry the same number of contacts, so a ratio taken
    # between two halves of one state and a ratio taken between the states
    # are equally precise. Comparing a half-depth technical ratio with a
    # full-depth ratio between states would make the technical floor look
    # larger than it is.
    a1 <- wc(halves[[1]]); a2 <- wc(halves[[2]])
    b1 <- wc(halves[[3]]); b2 <- wc(halves[[4]])
    peaks <- data.table(
      name = name, summit = as.numeric(uni), distance = as.numeric(uni) - vp,
      called_minus = called_minus, called_plus = called_plus,
      count_minus = cm, count_plus = cp,
      log2_ratio = log2((cp + 1) / (cm + 1)),
      log2_within_minus = log2((a2 + 1) / (a1 + 1)),
      log2_within_plus = log2((b2 + 1) / (b1 + 1)),
      log2_between_1 = log2((b1 + 1) / (a1 + 1)),
      log2_between_2 = log2((b2 + 1) / (a2 + 1)))
  }
  list(jac = jac, peaks = peaks,
       n_peaks_minus = length(unique(round(c(s$A1, s$A2) / opt$tol))),
       n_peaks_plus = length(unique(round(c(s$B1, s$B2) / opt$tol))))
}

# ---------------------------------------------------------------------------
# one experiment

do_experiment <- function(exp_name, lib, dir, opt) {
  labs <- lib[experiment == exp_name]
  minus <- labs[state == opt$minus_state, label]
  plus <- labs[state == opt$plus_state, label]
  message(sprintf("\n%s: %d mito- and %d mito+ libraries", exp_name,
                  length(minus), length(plus)))
  ex <- read_export(dir, c(minus, plus))
  ex <- Filter(Negate(is.null), ex)
  if (!all(c(minus, plus) %in% names(ex))) {
    minus <- intersect(minus, names(ex)); plus <- intersect(plus, names(ex))
  }
  if (length(minus) == 0 || length(plus) == 0) {
    message("  a state has no readable table, experiment skipped"); return(NULL)
  }
  same_bins <- vapply(ex, function(e) identical(e$offsets, ex[[1]]$offsets),
                      logical(1))
  if (!all(same_bins)) {
    # pooling and thinning add counts bin by bin, so every library of an
    # experiment must share one resolution and window
    message("  libraries of this experiment differ in resolution or window, ",
            "experiment skipped"); return(NULL)
  }
  res <- diff(ex[[1]]$offsets[1:2])
  offsets <- ex[[1]]$offsets
  setinfo <- ex[[1]]$set
  posv <- ex[[1]]$pos

  A <- pool_state(ex, minus)          # mito-
  B <- pool_state(ex, plus)           # mito+
  if (is.null(A) || is.null(B)) { message("  pooling failed"); return(NULL) }
  nm <- intersect(rownames(A), rownames(B))
  A <- A[nm, , drop = FALSE]; B <- B[nm, , drop = FALSE]

  # Depth is taken on the control viewpoints by default, so that a change at
  # mitochondrial viewpoints does not enter the thinning factor.
  ref <- if (opt$thin_on == "all") nm else nm[setinfo[nm] == "control"]
  dA <- sum(A[ref, ], na.rm = TRUE)
  dB <- sum(B[ref, ], na.rm = TRUE)
  if (opt$thin_on == "none") {
    fA <- 1; fB <- 1
  } else {
    if (dA <= 0 || dB <= 0) {
      message("  no contact on the ", opt$thin_on, " viewpoints"); return(NULL)
    }
    fA <- min(1, dB / dA); fB <- min(1, dA / dB)
  }
  if (opt$thin_on == "none")
    message(sprintf("  control-viewpoint contacts: mito- %.3g, mito+ %.3g; no thinning",
                    dA, dB))
  else
    message(sprintf("  %s-viewpoint contacts: mito- %.3g, mito+ %.3g; thinning %.3f / %.3f",
                    opt$thin_on, dA, dB, fA, fB))

  thin <- function(M, f) {
    if (f >= 1) return(M)
    ok <- is.finite(M)
    M[ok] <- rbinom(sum(ok), round(M[ok]), f)
    M
  }
  half <- function(M) {
    ok <- is.finite(M)
    H <- M
    H[ok] <- rbinom(sum(ok), round(M[ok]), 0.5)
    list(H, M - H)
  }
  set.seed(opt$seed)
  A <- thin(A, fA); B <- thin(B, fB)
  hA <- half(A); hB <- half(B)

  keep <- nm[rowSums(is.finite(A), na.rm = TRUE) > 0 &
               rowSums(is.finite(B), na.rm = TRUE) > 0]
  if (is.finite(opt$n_vp) && opt$n_vp < length(keep))
    keep <- sort(sample(keep, opt$n_vp))
  message(sprintf("  %d viewpoints (%d mito, %d control)", length(keep),
                  sum(setinfo[keep] == "mito"), sum(setinfo[keep] == "control")))

  work <- function(name) {
    vp <- posv[[name]]
    do_viewpoint(name,
                 list(hA[[1]][name, ], hA[[2]][name, ],
                      hB[[1]][name, ], hB[[2]][name, ]),
                 list(A = A[name, ], B = B[name, ]),
                 offsets, vp, res, opt)
  }
  out <- if (opt$cores > 1)
    mclapply(keep, work, mc.cores = max(1L, as.integer(opt$cores)))
  else lapply(keep, work)
  names(out) <- keep
  bad <- vapply(out, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  out <- out[!bad]
  if (length(out) == 0) { message("  no viewpoint returned"); return(NULL) }

  vpt <- rbindlist(lapply(names(out), function(n) {
    j <- out[[n]]$jac
    data.table(experiment = exp_name, name = n, set = setinfo[[n]],
               chrom = ex[[1]]$chrom[[n]], pos = posv[[n]],
               within_minus = j[["within_minus"]], within_plus = j[["within_plus"]],
               between_1 = j[["between_1"]], between_2 = j[["between_2"]],
               n_peaks_minus = out[[n]]$n_peaks_minus,
               n_peaks_plus = out[[n]]$n_peaks_plus)
  }))
  vpt[, within := rowMeans(cbind(within_minus, within_plus), na.rm = TRUE)]
  vpt[, between := rowMeans(cbind(between_1, between_2), na.rm = TRUE)]

  pk <- rbindlist(lapply(names(out), function(n) {
    p <- out[[n]]$peaks
    if (is.null(p)) return(NULL)
    p[, `:=`(experiment = exp_name, set = setinfo[[n]],
             chrom = ex[[1]]$chrom[[n]], vp_pos = posv[[n]])]
    p
  }), fill = TRUE)

  list(viewpoints = vpt, peaks = pk, resolution = res,
       thin = c(minus = fA, plus = fB))
}

# ---------------------------------------------------------------------------

summarise <- function(vpt, pk, opt) {
  a <- vpt[is.finite(within) & is.finite(between),
           .(n_viewpoints = .N,
             concordance_within = mean(within, na.rm = TRUE),
             concordance_between = mean(between, na.rm = TRUE)),
           by = .(experiment, set)]
  a[, concordance_loss := concordance_within - concordance_between]
  if (!is.null(pk) && nrow(pk)) {
    ok <- pk[is.finite(log2_ratio) &
               (count_minus + count_plus) >= 2 * opt$min_count]
    # `excess` is the share of loops whose change between the states is
    # larger than the technical noise of the same experiment allows: the
    # threshold is the --excess-quantile quantile of |log2| taken between two
    # halves of one state, so a set with no state-dependent loop returns
    # 1 - --excess-quantile. A median would not do here, since a minority of
    # changing loops leaves it untouched.
    b <- ok[, {
      wi <- abs(c(log2_within_minus, log2_within_plus))
      be <- abs(c(log2_between_1, log2_between_2))
      thr <- quantile(wi, opt$excess_quantile, na.rm = TRUE, names = FALSE)
      .(n_peaks = .N,
        shared = mean(called_minus & called_plus, na.rm = TRUE),
        amplitude_between = median(be, na.rm = TRUE),
        amplitude_within = median(wi, na.rm = TRUE),
        noise_q95 = thr,
        excess = mean(be > thr, na.rm = TRUE),
        direction = mean(log2_ratio, na.rm = TRUE))
    }, by = .(experiment, set)]
    a <- merge(a, b, by = c("experiment", "set"), all.x = TRUE)
  }
  a[order(experiment, set)]
}

meta <- function(s) {
  vars <- intersect(c("concordance_loss", "excess", "direction",
                      "n_viewpoints", "n_peaks"), names(s))
  w <- dcast(s, experiment ~ set, value.var = vars)
  need <- c("concordance_loss_mito", "concordance_loss_control")
  if (!all(need %in% names(w))) return(NULL)
  w[, d_concordance_loss := concordance_loss_mito - concordance_loss_control]
  if (all(c("excess_mito", "excess_control") %in% names(w)))
    w[, d_excess := excess_mito - excess_control]
  if (all(c("direction_mito", "direction_control") %in% names(w)))
    w[, d_direction := direction_mito - direction_control]
  w[]
}

sign_test <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) return(list(n = 0, mean = NA_real_, p = NA_real_, pos = NA_integer_))
  k <- sum(x > 0)
  list(n = n, mean = mean(x), pos = k,
       p = if (n > 0) binom.test(k, n, 0.5)$p.value else NA_real_)
}

figure <- function(s, path, floor_excess = 0.05) {
  ok <- s[is.finite(concordance_within) & is.finite(concordance_between)]
  if (!nrow(ok)) return(invisible(NULL))
  svg(path, width = 9, height = 4.2)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(6, 4.2, 2.4, 1), mgp = c(2.6, 0.7, 0))
  cols <- c(mito = "#B2182B", control = "#2166AC")
  exps <- unique(ok$experiment)
  xs <- seq_along(exps)
  yl <- range(c(ok$concordance_within, ok$concordance_between), na.rm = TRUE)
  pad <- max(0.02, diff(yl) * 0.25)
  plot(NA, xlim = c(0.5, length(exps) + 0.5), ylim = c(yl[1] - pad, yl[2] + pad),
       xaxt = "n", xlab = "", ylab = "peak concordance (Jaccard)",
       main = "concordance", cex.main = 1)
  axis(1, at = xs, labels = exps, las = 2, cex.axis = 0.8)
  for (st in names(cols)) {
    d <- ok[set == st]
    if (!nrow(d)) next
    i <- match(d$experiment, exps) + ifelse(st == "mito", -0.08, 0.08)
    segments(i, d$concordance_within, i, d$concordance_between,
             col = cols[[st]])
    points(i, d$concordance_within, col = cols[[st]], pch = 1)
    points(i, d$concordance_between, col = cols[[st]], pch = 19)
  }
  legend("bottomleft", c("mito", "control", "within state", "between states"),
         col = c(cols[["mito"]], cols[["control"]], "grey30", "grey30"),
         pch = c(19, 19, 1, 19), bty = "n", cex = 0.75)
  if ("excess" %in% names(ok)) {
    yl <- c(0, max(0.2, max(ok$excess, na.rm = TRUE) * 1.15))
    plot(NA, xlim = c(0.5, length(exps) + 0.5), ylim = yl, xaxt = "n",
         xlab = "", ylab = "share of loops changing",
         main = "change beyond noise", cex.main = 1)
    axis(1, at = xs, labels = exps, las = 2, cex.axis = 0.8)
    abline(h = floor_excess, col = "grey60", lty = 2)
    for (st in names(cols)) {
      d <- ok[set == st]
      if (!nrow(d)) next
      i <- match(d$experiment, exps) + ifelse(st == "mito", -0.08, 0.08)
      points(i, d$excess, col = cols[[st]], pch = 19)
      segments(i, floor_excess, i, d$excess, col = cols[[st]])
    }
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------

cmd_compare <- function(opt) {
  if (!requireNamespace("peakC", quietly = TRUE))
    stop("peakC is missing: remotes::install_github('deWitLab/peakC')")
  lib <- fread(opt$libraries, sep = "\t", colClasses = list(character = "label"))
  for (k in c("experiment", "state", "label")) lib[[k]] <- as.character(lib[[k]])
  exps <- lib[, .(n_states = uniqueN(state)), by = experiment][n_states == 2,
                                                               experiment]
  skipped <- setdiff(unique(lib$experiment), exps)
  if (length(skipped))
    message("experiments without both states, skipped: ",
            paste(skipped, collapse = ", "))
  if (!length(exps)) stop("no experiment has both states")
  dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

  res <- lapply(exps, function(e) do_experiment(e, lib, opt$export, opt))
  names(res) <- exps
  res <- Filter(Negate(is.null), res)
  if (!length(res)) stop("no experiment produced a result")

  vpt <- rbindlist(lapply(res, `[[`, "viewpoints"), fill = TRUE)
  pk <- rbindlist(lapply(res, `[[`, "peaks"), fill = TRUE)
  s <- summarise(vpt, pk, opt)

  fwrite(vpt, file.path(opt$out, "pairs_viewpoints.tsv"), sep = "\t")
  if (nrow(pk)) fwrite(pk, file.path(opt$out, "pairs_peaks.tsv"), sep = "\t")
  fwrite(s, file.path(opt$out, "pairs_summary.tsv"), sep = "\t")

  cat("\nper experiment and promoter set\n")
  print(s)

  m <- meta(s)
  if (!is.null(m)) {
    fwrite(m, file.path(opt$out, "pairs_meta.tsv"), sep = "\t")
    cat("\nmito minus control, one value per experiment\n")
    cols <- intersect(c("experiment", "d_concordance_loss", "d_excess",
                        "d_direction"), names(m))
    print(m[, ..cols])
    for (v in setdiff(cols, "experiment")) {
      t <- sign_test(m[[v]])
      cat(sprintf("  %-20s mean %+.4f over %d experiments, %d positive, sign test p = %.3f\n",
                  v, t$mean, t$n, t$pos, t$p))
    }
    cat("\nA positive d_concordance_loss means the loops of mitochondrial\n",
        "promoters differ between the two states more than those of control\n",
        "promoters, beyond the technical floor measured within each state.\n",
        sep = "")
  }
  figure(s, file.path(opt$out, "pairs.svg"), 1 - opt$excess_quantile)
  cat("\nWritten to ", opt$out, "/\n", sep = "")
}

# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || args[1] == "help" || any(args %in% c("-h", "--help"))) {
  usage()
  quit(status = if (length(args)) 0 else 1)
}
opt <- parse_args(args)
if (identical(opt$cmd, "compare")) cmd_compare(opt) else
  stop("unknown subcommand: ", opt$cmd)
