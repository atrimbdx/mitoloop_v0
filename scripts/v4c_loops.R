#!/usr/bin/env Rscript
# v4c_loops.R - peakC loops on virtual-4C profiles, reciprocity, and loop
# strength against mtDNA copy number
#
# Purpose
#   R side of the mitoloop_v0 loop pipeline (v4c_pipeline.py runs the steps
#   before and after). Calls peaks with peakC on the virtual-4C profiles of
#   mitochondrial promoters, term by term, keeps the peaks confirmed from the
#   partner anchor, merges them into loops, and relates loop strength to the
#   mtDNA copy number of the matching GTEx tissue.
#
# Subcommands
#   scan        Window size for peakC, chosen on reproducibility between
#               libraries of one term, as in Geeven et al. 2018: peaks are
#               called separately in each library (peakC single.analysis) and
#               compared by the Jaccard index at the level of peaks (runs of
#               consecutive significant bins), two peaks matching when their
#               summits lie within --tol. The smallest window whose median
#               Jaccard reaches --plateau times the best one is proposed.
#   call        peakC on the mitochondrial viewpoints, term by term:
#               combined.analysis (rank product over replicates) when the term
#               has two libraries or more, single.analysis when it has one.
#               Writes the peaks and the partner anchors to extract.
#   reciprocal  peakC on the profile of each partner anchor, in the same term.
#               A mitochondrial peak is kept when the partner's profile has a
#               peak at the mitochondrial TSS (the reciprocity criterion of
#               peakHiC). Reciprocal peaks of all terms are merged into loops.
#   deseq       Counts at each loop in every library, DESeq2 with term as
#               factor, a normalisation factor per library and loop (contacts
#               of the viewpoint in that library times the share of contacts
#               that library puts at the loop's distance, taken on the control
#               viewpoints), and, for every pair of terms, the log2 fold change
#               shrunk with ashr. Input of `v4c_pipeline.py relate`.
#   assoc       Loop strength regressed on mtDNA copy number, one observation
#               per library. For each loop, log2(observed / expected contact),
#               the expected being the one used by `deseq`, is fitted on log2
#               mtCN and log10 depth with limma (vooma weights;
#               duplicateCorrelation for the correlation between libraries of
#               one GTEx tissue, the level at which mtCN is measured). Global
#               tests: mean t^2 (association in either direction) and mean t
#               (direction), with a null built by moving mtCN between GTEx
#               tissues. With --sets, the same permutations test whether the
#               loops of a gene set follow mtCN more than the other loops.
#               With --reference, the loops of those promoters are the
#               reference: every other gene, and every set, is compared with
#               them on the change in contact per doubling of mtCN, each loop
#               against reference loops of similar contact level.
#   tissue      The same normalised loop values as `assoc`, averaged into one
#               loop strength per library and per GTEx tissue, against mtCN.
#
# Inputs
#   --libraries  libraries.tsv of `v4c_pipeline.py select`: accession, term,
#                gtex, mtcn, depth, log10_depth, kept; only kept == 1 is used.
#   --export     folder of `v4c_pipeline.py export`: <accession>.tsv.gz with
#                name, chrom, pos, then one column per offset in bp; decay.tsv
#                (accession, one column per offset) for deseq, assoc, tissue.
#                For reciprocal, the export of the partner set.
#   --out        working folder: call writes into it, reciprocal reads
#                peaks_mito.tsv from it, deseq, assoc and tissue read
#                loops.tsv from it.
#   --sets       optional TSV with columns source, set, gene (assoc).
#   --reference  optional file of gene symbols, one per line (assoc, tissue).
#
# Outputs (in --out)
#   scan        scan.tsv: window, window_kb, median_j, q1, q3, peaks_per_lib
#   call        peaks_mito.tsv: term, gene, chrom, vp_pos, n_libs, start, end,
#               summit, summit_ratio; partners.tsv: name, chrom, pos, set
#   reciprocal  peaks_mito.tsv with a column reciprocal (0/1); loops.tsv:
#               gene, chrom, vp_pos, summit, start, end, n_terms, called (the
#               ENCODE terms in which the loop is reciprocal, |-separated),
#               loop_id (gene:summit)
#   deseq       deseq_pairs.tsv.gz: term_a, term_b, loop_id, gene, log2fc,
#               lfcSE, pvalue, padj, baseMean, called_a, called_b;
#               presence.tsv: gene, one 0/1 column per term
#   assoc       assoc_loops.tsv (one row per loop: log2fc_per_doubling, t,
#               p_limma, fdr_limma, p_perm, fdr_perm, ...), assoc_global.tsv
#               (statistic, observed, null_median, null_q95, p_perm),
#               assoc_null.tsv (the permutation null), assoc_sets.tsv and
#               assoc_genes.tsv when sets or a reference are given,
#               assoc_volcano.svg
#   tissue      tissue_strength_libraries.tsv (accession, term, gtex, mtcn,
#               depth, strength), tissue_strength.tsv (gtex, libraries, mtcn,
#               strength, sd), tissue_strength.svg
#
# Method
#   Profiles are converted to reads per million of the profile and given to
#   peakC with bin centres as positions; peaks closer than --dist-min to the
#   viewpoint are not tested and peaks beyond --dist-max are dropped. A loop's
#   count is the sum over --wsize bins centred on its summit. Its expected
#   count is the viewpoint's contacts in the --dist-min to --dist-max band
#   times the share of contacts that library puts at that distance on the
#   control viewpoints, so that differences in distance decay between
#   libraries do not read as differences in looping. In assoc the GLS fit is
#   written out so that each permutation costs two matrix-vector products;
#   it reproduces limma::lmFit with the same block, correlation and weights,
#   and the script checks the agreement. Permutations move mtCN between GTEx
#   tissues, the libraries of a tissue moving together, because mtCN is
#   measured per tissue.
#
# Usage
#   Rscript v4c_loops.R scan --libraries libraries.tsv --export export/mito
#   Rscript v4c_loops.R call --libraries libraries.tsv --export export/mito \
#       --wsize 3 --out loops
#   Rscript v4c_loops.R reciprocal --libraries libraries.tsv \
#       --export export/partners --wsize 3 --out loops
#   Rscript v4c_loops.R assoc --libraries libraries.tsv --export export/mito \
#       --wsize 3 --out loops [--sets gene_sets.tsv] [--reference controls.txt]
#   Rscript v4c_loops.R tissue --libraries libraries.tsv --export export/mito \
#       --wsize 3 --out loops
#   Rscript v4c_loops.R deseq --libraries libraries.tsv --export export/mito \
#       --wsize 3 --out loops
#   Rscript v4c_loops.R --help          (all options)
#
#   With --reference, the reference genes need loops: export the mito and
#   control sets together (v4c_pipeline.py export --set mito,control) and run
#   call, reciprocal, assoc and tissue on that export.
#
#   --wsize must be identical in call, reciprocal, deseq, assoc and tissue:
#   they all work on the same window around a loop summit.
#
# Packages: data.table, parallel, peakC (GitHub deWitLab/peakC), isotone,
#   caTools, limma (assoc); DESeq2 and ashr (deseq only).
#
# References
#   Geeven et al. 2018. peakC: a flexible, non-parametric peak calling
#       package for 4C and Capture-C data. Nucleic Acids Research.
#   Bianchi et al. 2019. Detailed regulatory interaction map of the human
#       heart facilitates gene discovery for cardiovascular disease. bioRxiv.
#       (peakHiC: virtual 4C peak calling on Hi-C, loops called from both
#       anchors)
#   Breitling et al. 2004. Rank products: a simple, yet powerful, new method
#       to detect differentially regulated genes in replicated microarray
#       experiments. FEBS Letters. (combined.analysis of peakC)
#   de Leeuw et al. 2009. Isotone optimization in R: pool-adjacent-violators
#       algorithm (PAVA) and active set methods. Journal of Statistical
#       Software. (isotone, used by peakC)
#   Durand et al. 2016. Juicer provides a one-click system for analyzing
#       loop-resolution Hi-C experiments. Cell Systems. (.hic, hicstraw)
#   Abdennur and Mirny 2020. Cooler: scalable storage for Hi-C data and other
#       genomically labeled arrays. Bioinformatics.
#   ENCODE Project Consortium et al. 2020. Expanded encyclopaedias of DNA
#       elements in the human and mouse genomes. Nature. (intact Hi-C data)
#   GTEx Consortium 2020. The GTEx Consortium atlas of genetic regulatory
#       effects across human tissues. Science.
#   Ritchie et al. 2015. limma powers differential expression analyses for
#       RNA-sequencing and microarray studies. Nucleic Acids Research.
#   Law et al. 2014. voom: precision weights unlock linear model analysis
#       tools for RNA-seq read counts. Genome Biology. (vooma)
#   Smyth et al. 2005. Use of within-array replicate spots for assessing
#       differential expression in microarray experiments. Bioinformatics.
#       (duplicateCorrelation)
#   Love et al. 2014. Moderated estimation of fold change and dispersion for
#       RNA-seq data with DESeq2. Genome Biology.
#   Stephens 2017. False discovery rates: a new deal. Biostatistics. (ashr)
#   Sen 1968. Estimates of the regression coefficient based on Kendall's tau.
#       Journal of the American Statistical Association. (Theil-Sen line)
#   Jaccard 1912. The distribution of the flora in the alpine zone. New
#       Phytologist. (Jaccard index of peak sets, scan)
#   Phipson and Smyth 2010. Permutation P-values should never be zero.
#       Statistical Applications in Genetics and Molecular Biology.
#       (permutation p-values, assoc)
#   caTools, R package (running statistics, used by peakC).

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})
# peakC calls isotone itself; loading it here only silences its message
if (requireNamespace("isotone", quietly = TRUE))
  suppressPackageStartupMessages(library(isotone))

# ---------------------------------------------------------------------------
# arguments

die <- function(...) {
  message("error: ", ...)
  quit(save = "no", status = 1)
}

COMMANDS <- c("scan", "call", "reciprocal", "deseq", "assoc", "tissue")
NEED_WSIZE <- c("call", "reciprocal", "deseq", "assoc", "tissue")

# name, default, type (chr, num, choice), allowed values, subcommands, help
OPTIONS <- list(
  list("libraries", NULL, "chr", NULL, "all",
       "libraries.tsv written by v4c_pipeline.py select (required)"),
  list("export", NULL, "chr", NULL, "all",
       "folder of tables written by v4c_pipeline.py export (required)"),
  list("out", "loops", "chr", NULL, "all", "working and output folder"),
  list("wsize", NA_real_, "num", NULL, NEED_WSIZE,
       "peakC window size in bins, no default: run `scan` and use the window it proposes"),
  list("cores", 4, "num", NULL, "all", "parallel processes (mclapply)"),
  list("dist_min", 50e3, "num", NULL, "all",
       "bp; peakC minDist, and lower edge of the distance band of the expected count"),
  list("dist_max", 1e6, "num", NULL, "all",
       "bp; peaks beyond are dropped, upper edge of the distance band"),
  list("alpha_fdr", 0.1, "num", NULL, c("call", "reciprocal"),
       "peakC alphaFDR of combined.analysis (two libraries or more)"),
  list("qwr", 1, "num", NULL, c("scan", "call", "reciprocal"),
       "peakC qWr, ratio threshold over the background"),
  list("qwd", 1.5, "num", NULL, c("scan", "call", "reciprocal"),
       "peakC qWd of single.analysis (one library)"),
  list("min_bins", 20, "num", NULL, c("scan", "call", "reciprocal"),
       "minimum finite bins in a profile (half on each side of the viewpoint)"),
  list("single", "yes", "choice", c("yes", "no"), c("call", "reciprocal"),
       "call peaks in terms with one library (single.analysis)"),
  list("per_term", "all", "choice", c("all", "one"), c("call", "reciprocal"),
       "libraries calling the peaks of a term: all, or one (the deepest)"),
  list("thin", NA_real_, "num", NULL, c("call", "reciprocal"),
       "thin every calling library to this many contacts (binomial); NA: no thinning"),
  list("tol", 20e3, "num", NULL, c("scan", "reciprocal"),
       "bp; scan: summit matching tolerance; reciprocal: a partner peak within tol + half a bin of the TSS confirms the loop"),
  list("wsizes", "1,3,5,7,9,11", "chr", NULL, "scan",
       "comma-separated window sizes (bins) to compare"),
  list("plateau", 0.9, "num", NULL, "scan",
       "smallest window whose median Jaccard reaches plateau x the best"),
  list("n_vp", 300, "num", NULL, "scan", "viewpoints sampled"),
  list("max_pairs", 3, "num", NULL, "scan", "library pairs per term"),
  list("seed", 1, "num", NULL, "all", "random seed"),
  list("deseq_min_libs", 2, "num", NULL, "deseq",
       "loops measured in fewer libraries are dropped"),
  list("deseq_min_count", 10, "num", NULL, "deseq",
       "loops with fewer contacts summed over libraries are dropped"),
  list("min_terms", 1, "num", NULL, c("assoc", "tissue"),
       "keep the loops called in at least this many terms"),
  list("min_count", 5, "num", NULL, c("assoc", "tissue"),
       "keep the loops whose median count over libraries reaches this"),
  list("min_loops", 20, "num", NULL, c("assoc", "tissue"),
       "stop when fewer loops remain"),
  list("pseudocount", 0.5, "num", NULL, c("assoc", "tissue"),
       "added to counts before log2(observed / expected)"),
  list("adjust_depth", "yes", "choice", c("yes", "no"), c("assoc", "tissue"),
       "log10 depth as covariate (assoc), removed before averaging (tissue)"),
  list("perms", 2000, "num", NULL, "assoc",
       "permutations of mtCN between GTEx tissues"),
  list("fdr", 0.1, "num", NULL, "assoc", "FDR for loops, sets and genes"),
  list("sets", NULL, "chr", NULL, "assoc",
       "TSV source, set, gene: gene sets tested against the other loops"),
  list("min_set", 10, "num", NULL, "assoc", "minimum loops in a tested set"),
  list("reference", NULL, "chr", NULL, c("assoc", "tissue"),
       "gene symbols whose loops are the reference (assoc) or left out (tissue)"),
  list("min_ref_loops", 50, "num", NULL, "assoc",
       "stop when the reference has fewer loops"),
  list("ref_strata", 10, "num", NULL, "assoc",
       "strata of contact level for the comparison with the reference")
)
OPT_NAMES <- vapply(OPTIONS, `[[`, "", 1)

usage <- function(cmd = NULL) {
  cat("Usage: Rscript v4c_loops.R <command> [--option value ...]\n",
      "Commands: ", paste(COMMANDS, collapse = ", "), "\n",
      "Options may be written --dist-min or --dist_min, and --key=value.\n",
      "Unknown options are an error.\n\n", sep = "")
  for (o in OPTIONS) {
    where <- if (identical(o[[5]], "all")) "all" else paste(o[[5]], collapse = ",")
    if (!is.null(cmd) && !identical(o[[5]], "all") && !(cmd %in% o[[5]])) next
    def <- if (is.null(o[[2]])) "none" else if (is.na(o[[2]])) "NA" else
      format(o[[2]], scientific = FALSE)
    vals <- if (!is.null(o[[4]])) paste0(" {", paste(o[[4]], collapse = ","), "}") else ""
    cat(sprintf("  --%-16s %s%s\n  %18s default %s; used by: %s\n",
                gsub("_", "-", o[[1]]), o[[6]], vals, "", def, where))
  }
}

parse_args <- function(a) {
  if (!length(a) || a[1] %in% c("-h", "--help", "help")) {
    usage(); quit(save = "no", status = if (length(a)) 0 else 1)
  }
  cmd <- a[1]
  if (!(cmd %in% COMMANDS))
    die("unknown command '", cmd, "'; one of ", paste(COMMANDS, collapse = ", "),
        " (see --help)")
  opt <- setNames(lapply(OPTIONS, `[[`, 2), OPT_NAMES)
  opt <- c(list(cmd = cmd), opt[!vapply(opt, is.null, logical(1))])
  i <- 2
  while (i <= length(a)) {
    tok <- a[i]
    if (tok %in% c("-h", "--help")) { usage(cmd); quit(save = "no", status = 0) }
    if (!startsWith(tok, "--")) die("unexpected argument '", tok, "'")
    tok <- sub("^--", "", tok)
    if (grepl("=", tok, fixed = TRUE)) {
      key <- sub("=.*$", "", tok); val <- sub("^[^=]*=", "", tok); step <- 1
    } else {
      key <- tok
      if (i + 1 > length(a) || startsWith(a[i + 1], "--"))
        die("option --", key, " needs a value")
      val <- a[i + 1]; step <- 2
    }
    key <- gsub("-", "_", key)
    if (!(key %in% OPT_NAMES))
      die("unknown option --", gsub("_", "-", key), " (see Rscript v4c_loops.R --help)")
    spec <- OPTIONS[[match(key, OPT_NAMES)]]
    if (spec[[3]] == "num") {
      v <- suppressWarnings(as.numeric(val))
      if (is.na(v) && toupper(val) != "NA")
        die("option --", gsub("_", "-", key), " needs a number, got '", val, "'")
      val <- v
    } else if (spec[[3]] == "choice" && !(val %in% spec[[4]])) {
      die("option --", gsub("_", "-", key), " takes ",
          paste(spec[[4]], collapse = " or "), ", got '", val, "'")
    }
    opt[[key]] <- val
    i <- i + step
  }
  for (k in c("libraries", "export"))
    if (is.null(opt[[k]])) die("--", k, " is required (see --help)")
  if (cmd %in% NEED_WSIZE && !is.finite(opt$wsize))
    die("--wsize is required for ", cmd, ": the peakC window size in bins ",
        "(an odd number). Choose it with\n  Rscript v4c_loops.R scan ",
        "--libraries ", opt$libraries, " --export <export of the mito set>\n",
        "which proposes the smallest window whose median Jaccard between ",
        "replicate libraries reaches --plateau times the best one, and use the ",
        "same value in call, reciprocal, deseq, assoc and tissue.")
  opt
}

# ---------------------------------------------------------------------------
# reading

read_libraries <- function(path) {
  lib <- fread(path, sep = "\t", colClasses = list(character = "accession"))
  lib <- lib[kept == 1]
  lib[, term := as.character(term)]
  lib
}

# which libraries call the peaks of a term: all of them, or only the deepest
call_libs <- function(lib, opt) {
  if (identical(opt$per_term, "one"))
    lib[order(-depth), .(accession = accession[1]), by = term]
  else
    lib[, .(accession = accession), by = term]
}


read_export <- function(dir, accs) {
  out <- list()
  for (a in accs) {
    f <- file.path(dir, paste0(a, ".tsv.gz"))
    if (!file.exists(f)) { message("missing ", f); next }
    d <- fread(f, sep = "\t", na.strings = "NA", header = TRUE)
    m <- as.matrix(d[, -(1:3)])
    storage.mode(m) <- "double"
    rownames(m) <- d$name
    out[[a]] <- list(m = m, chrom = setNames(d$chrom, d$name),
                     pos = setNames(d$pos, d$name),
                     offsets = as.numeric(colnames(d)[-(1:3)]))
  }
  if (!length(out)) die("no table of a kept library found in ", dir)
  out
}

# optional: bring every calling library to the same depth
thin_factors <- function(lib, use, opt) {
  d <- setNames(lib$depth, lib$accession)[use$accession]
  if (!is.finite(opt$thin)) return(setNames(as.list(rep(1, length(d))),
                                            names(d)))
  setNames(as.list(pmin(1, opt$thin / d)), names(d))
}


# one library profile as peakC wants it: bin centre, reads per million
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

# peakC on the replicate profiles of one viewpoint in one term
run_peakc <- function(profiles, vp, opt) {
  profiles <- Filter(Negate(is.null), profiles)
  k <- length(profiles)
  if (k == 0) return(NULL)
  if (k == 1 && opt$single != "yes") return(NULL)
  r <- tryCatch({
    if (k >= 2) {
      d <- profiles; d$num.exp <- k
      x <- peakC::combined.analysis(d, num.exp = k, vp.pos = vp,
                                    wSize = opt$wsize, alphaFDR = opt$alpha_fdr,
                                    qWr = opt$qwr, minDist = opt$dist_min)
      list(peak = x$peak, pos = x$dbR[, 1], ratio = x$ratio, k = k)
    } else {
      x <- peakC::single.analysis(profiles[[1]], vp.pos = vp,
                                  wSize = opt$wsize, qWd = opt$qwd,
                                  qWr = opt$qwr, minDist = opt$dist_min)
      list(peak = x$peak, pos = x$ratios[, 1], ratio = x$ratios[, 2], k = 1)
    }
  }, error = function(e) NULL)
  if (is.null(r)) return(NULL)
  r$peak <- r$peak[abs(r$peak - vp) <= opt$dist_max]
  r
}

# positions as plain integers in names (never scientific notation)
plain_int <- function(x) format(x, scientific = FALSE, trim = TRUE)

regions <- function(peak, pos, ratio, res) {
  if (length(peak) == 0) return(NULL)
  p <- sort(unique(peak))
  grp <- cumsum(c(1, diff(p) > res * 1.01))
  rbindlist(lapply(split(p, grp), function(x) {
    rr <- ratio[match(x, pos)]
    data.table(start = min(x), end = max(x), summit = x[which.max(rr)],
               summit_ratio = max(rr, na.rm = TRUE))
  }))
}

term_profiles <- function(ex, accs, name, vp, opt) {
  lapply(accs, function(a) {
    e <- ex[[a]]
    if (is.null(e) || !(name %in% rownames(e$m))) return(NULL)
    v <- e$m[name, ]
    f <- if (is.null(opt$thin_f)) 1 else opt$thin_f[[a]]
    if (is.finite(f) && f < 1) {          # same depth for every library
      ok <- is.finite(v)
      v[ok] <- rbinom(sum(ok), round(v[ok]), f)
    }
    prep(v, e$offsets, vp, opt$res, opt$min_bins)
  })
}

# ---------------------------------------------------------------------------
# scan

cmd_scan <- function(opt) {
  lib <- read_libraries(opt$libraries)
  ex <- read_export(opt$export, lib$accession)
  opt$res <- diff(ex[[1]]$offsets[1:2])
  set.seed(opt$seed)
  reps <- lib[, .(accs = list(accession)), by = term][lengths(accs) >= 2]
  if (nrow(reps) == 0) die("no term with two kept libraries")
  pairs <- rbindlist(lapply(seq_len(nrow(reps)), function(i) {
    cmb <- combn(reps$accs[[i]], 2)
    cmb <- cmb[, seq_len(min(ncol(cmb), opt$max_pairs)), drop = FALSE]
    data.table(term = reps$term[i], a = cmb[1, ], b = cmb[2, ])
  }))
  genes <- rownames(ex[[1]]$m)
  genes <- sample(genes, min(length(genes), opt$n_vp))
  ws <- as.integer(strsplit(opt$wsizes, ",")[[1]])
  cat(sprintf("%d terms with replicates, %d library pairs, %d viewpoints\n",
              nrow(reps), nrow(pairs), length(genes)))
  cat(sprintf("%6s %10s %12s %12s\n", "window", "median J", "IQR", "peaks/lib"))
  tab <- list()
  for (w in ws) {
    opt$wsize <- w
    one <- mclapply(seq_len(nrow(pairs)), function(i) {
      inter <- 0; uni <- 0; np <- 0
      for (g in genes) {
        vp <- ex[[pairs$a[i]]]$pos[[g]]
        pa <- term_profiles(ex, pairs$a[i], g, vp, opt)[[1]]
        pb <- term_profiles(ex, pairs$b[i], g, vp, opt)[[1]]
        if (is.null(pa) || is.null(pb)) next
        opt1 <- opt; opt1$single <- "yes"
        xa <- run_peakc(list(pa), vp, opt1)
        xb <- run_peakc(list(pb), vp, opt1)
        sa <- if (is.null(xa)) numeric(0) else regions(xa$peak, xa$pos, xa$ratio, opt$res)$summit
        sb <- if (is.null(xb)) numeric(0) else regions(xb$peak, xb$pos, xb$ratio, opt$res)$summit
        # peaks compared by their summits with a fixed tolerance in kb, so
        # that wider windows, which mark wider stretches, are not favoured
        ma <- sum(vapply(sa, function(x) any(abs(sb - x) <= opt$tol), logical(1)))
        mb <- sum(vapply(sb, function(x) any(abs(sa - x) <= opt$tol), logical(1)))
        m <- (ma + mb) / 2
        inter <- inter + m
        uni <- uni + length(sa) + length(sb) - m
        np <- np + length(sa) + length(sb)
      }
      c(j = if (uni > 0) inter / uni else NA, peaks = np / 2)
    }, mc.cores = opt$cores)
    m <- do.call(rbind, one)
    q <- quantile(m[, "j"], c(.25, .5, .75), na.rm = TRUE)
    tab[[length(tab) + 1]] <- data.table(window = w, window_kb = w * opt$res / 1e3,
                                         median_j = q[2], q1 = q[1], q3 = q[3],
                                         peaks_per_lib = median(m[, "peaks"]))
    cat(sprintf("%4dkb %10.3f %5.3f-%-6.3f %12.0f\n", w * opt$res / 1e3, q[2],
                q[1], q[3], median(m[, "peaks"])))
  }
  tab <- rbindlist(tab)
  best <- max(tab$median_j, na.rm = TRUE)
  pick <- tab[median_j >= opt$plateau * best][1]
  cat(sprintf("\nproposed window: %d bins (%g kb), median Jaccard %.3f, best %.3f\n",
              pick$window, pick$window_kb, pick$median_j, best))
  dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
  fwrite(tab, file.path(opt$out, "scan.tsv"), sep = "\t")
}

# ---------------------------------------------------------------------------
# call

cmd_call <- function(opt) {
  set.seed(opt$seed)                       # thinning draws, if any
  lib <- read_libraries(opt$libraries)
  ex <- read_export(opt$export, lib$accession)
  opt$res <- diff(ex[[1]]$offsets[1:2])
  use <- call_libs(lib, opt)
  opt$thin_f <- thin_factors(lib, use, opt)
  genes <- rownames(ex[[1]]$m)
  terms <- sort(unique(lib$term))
  cat(sprintf("%d terms, %d viewpoints, window %d bins, %s\n", length(terms),
              length(genes), opt$wsize,
              if (identical(opt$per_term, "one"))
                "one library per term (the deepest)" else
                "every library of the term"))
  if (is.finite(opt$thin)) cat(sprintf(
    "  thinned to %.3g contacts; factors %.2f to %.2f\n", opt$thin,
    min(unlist(opt$thin_f)), max(unlist(opt$thin_f))))
  res <- mclapply(genes, function(g) {
    e <- ex[[which(sapply(ex, function(x) g %in% rownames(x$m)))[1]]]
    vp <- e$pos[[g]]; ch <- e$chrom[[g]]
    rbindlist(lapply(terms, function(t) {
      r <- run_peakc(term_profiles(ex, use[term == t]$accession, g, vp, opt),
                     vp, opt)
      if (is.null(r)) return(NULL)
      rg <- regions(r$peak, r$pos, r$ratio, opt$res)
      if (is.null(rg)) return(NULL)
      cbind(data.table(term = t, gene = g, chrom = ch, vp_pos = vp,
                       n_libs = r$k), rg)
    }))
  }, mc.cores = opt$cores)
  peaks <- rbindlist(res)
  dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
  if (nrow(peaks) == 0) die("no peak called in any term; nothing written")
  fwrite(peaks, file.path(opt$out, "peaks_mito.tsv"), sep = "\t")
  partners <- unique(peaks[, .(chrom, pos = summit)])
  partners[, name := paste0("P_", chrom, "_", plain_int(pos))]
  partners[, set := "partner"]
  fwrite(partners[, .(name, chrom, pos, set)],
         file.path(opt$out, "partners.tsv"), sep = "\t")
  cat(sprintf("%d peak regions, %d partner anchors to extract\n",
              nrow(peaks), nrow(partners)))
  print(peaks[, .(regions = .N, viewpoints = uniqueN(gene)), by = .(term, n_libs)])
}

# ---------------------------------------------------------------------------
# reciprocal

cmd_reciprocal <- function(opt) {
  set.seed(opt$seed)                       # thinning draws, if any
  lib <- read_libraries(opt$libraries)
  ex <- read_export(opt$export, lib$accession)
  opt$res <- diff(ex[[1]]$offsets[1:2])
  use <- call_libs(lib, opt)
  opt$thin_f <- thin_factors(lib, use, opt)
  peaks <- fread(file.path(opt$out, "peaks_mito.tsv"), sep = "\t")
  peaks[, partner := paste0("P_", chrom, "_", plain_int(summit))]
  keys <- unique(peaks[, .(term, partner, summit)])
  cat(sprintf("%d partner profiles to call\n", nrow(keys)))
  calls <- mclapply(seq_len(nrow(keys)), function(i) {
    r <- run_peakc(term_profiles(ex, use[term == keys$term[i]]$accession,
                                 keys$partner[i], keys$summit[i], opt),
                   keys$summit[i], opt)
    if (is.null(r)) numeric(0) else r$peak
  }, mc.cores = opt$cores)
  names(calls) <- paste(keys$term, keys$partner)
  peaks[, reciprocal := mapply(function(t, p, vp) {
    pk <- calls[[paste(t, p)]]
    as.integer(length(pk) > 0 && any(abs(pk - vp) <= opt$tol + opt$res / 2))
  }, term, partner, vp_pos)]
  fwrite(peaks, file.path(opt$out, "peaks_mito.tsv"), sep = "\t")
  cat(sprintf("reciprocal: %d of %d peak regions (%.0f%%)\n",
              sum(peaks$reciprocal), nrow(peaks), 100 * mean(peaks$reciprocal)))

  rc <- peaks[reciprocal == 1]
  if (nrow(rc) == 0) {
    fwrite(data.table(gene = character(), chrom = character(),
                      vp_pos = numeric(), summit = numeric(), start = numeric(),
                      end = numeric(), n_terms = integer(), called = character(),
                      loop_id = character()),
           file.path(opt$out, "loops.tsv"), sep = "\t")
    cat("no reciprocal peak: loops.tsv written with its header only\n")
    return(invisible(NULL))
  }
  loops <- rbindlist(lapply(split(rc, rc$gene), function(x) {
    bins <- sort(unique(unlist(mapply(function(s, e) seq(s, e, by = opt$res),
                                      x$start, x$end))))
    comp <- cumsum(c(1, diff(bins) > opt$res * 1.01))
    rbindlist(lapply(split(bins, comp), function(b) {
      inside <- x[summit >= min(b) & summit <= max(b)]
      tab <- inside[, .(n = .N, r = max(summit_ratio)), by = summit][order(-n, -r)]
      called <- unique(x[start <= max(b) & end >= min(b)]$term)
      data.table(gene = x$gene[1], chrom = x$chrom[1], vp_pos = x$vp_pos[1],
                 summit = tab$summit[1], start = min(b), end = max(b),
                 n_terms = length(called),
                 called = paste(called, collapse = "|"))
    }))
  }))
  loops[, loop_id := paste0(gene, ":", summit)]
  fwrite(loops, file.path(opt$out, "loops.tsv"), sep = "\t")
  cat(sprintf("%d loops over %d viewpoints; called in a median of %g terms\n",
              nrow(loops), uniqueN(loops$gene), median(loops$n_terms)))
}

# ---------------------------------------------------------------------------
# deseq

# Counts at each loop in every library, and the expected count used to
# normalise them: contacts of the viewpoint in that library times the share of
# contacts the library puts at the loop's distance, taken on the control
# viewpoints (decay.tsv). w = 1 where the library has data for the loop.
loop_matrix <- function(opt, lib, ex, loops) {
  offs <- ex[[1]]$offsets
  band <- abs(offs) >= opt$dist_min & abs(offs) <= opt$dist_max
  decay <- fread(file.path(opt$export, "decay.tsv"), sep = "\t", header = TRUE,
                 colClasses = list(character = "accession"))
  D <- as.matrix(decay[, -1]); rownames(D) <- decay$accession
  frac <- D / rowSums(D[, band, drop = FALSE])
  h <- floor(opt$wsize / 2)
  L <- nrow(loops); N <- nrow(lib)
  cnt <- matrix(0, L, N); nf <- matrix(NA_real_, L, N); w <- matrix(0, L, N)
  for (j in seq_len(N)) {
    e <- ex[[lib$accession[j]]]
    for (i in seq_len(L)) {
      g <- loops$gene[i]
      if (!(g %in% rownames(e$m))) next
      v <- e$m[g, ]
      centre <- floor(loops$vp_pos[i] / opt$res) * opt$res
      k <- match(loops$summit[i] - opt$res / 2 - centre, offs)
      if (is.na(k)) next
      cols <- (k - h):(k + h); cols <- cols[cols >= 1 & cols <= length(offs)]
      if (all(!is.finite(v[cols]))) next
      tot <- sum(v[band], na.rm = TRUE)
      fw <- sum(frac[lib$accession[j], cols])
      if (tot <= 0 || fw <= 0) next
      cnt[i, j] <- sum(v[cols], na.rm = TRUE)
      nf[i, j] <- tot * fw
      w[i, j] <- 1
    }
  }
  list(cnt = cnt, nf = nf, w = w)
}

cmd_deseq <- function(opt) {
  miss <- c("DESeq2", "ashr")[!vapply(c("DESeq2", "ashr"), requireNamespace,
                                      logical(1), quietly = TRUE)]
  if (length(miss))
    die("deseq needs the Bioconductor package(s) ", paste(miss, collapse = ", "),
        ", not installed. Install them with BiocManager::install(c(\"DESeq2\", ",
        "\"ashr\")), or use assoc and tissue, which need limma only.")
  suppressPackageStartupMessages(library(DESeq2))
  lib <- read_libraries(opt$libraries)
  ex <- read_export(opt$export, lib$accession)
  lib <- lib[accession %in% names(ex)]
  opt$res <- diff(ex[[1]]$offsets[1:2])
  loops <- fread(file.path(opt$out, "loops.tsv"), sep = "\t")
  lm <- loop_matrix(opt, lib, ex, loops)
  cnt <- lm$cnt; nf <- lm$nf; w <- lm$w; N <- nrow(lib)
  keep <- rowSums(w) >= opt$deseq_min_libs & rowSums(cnt) >= opt$deseq_min_count
  cnt <- cnt[keep, , drop = FALSE]; nf <- nf[keep, , drop = FALSE]
  w <- w[keep, , drop = FALSE]; loops <- loops[keep]
  gm <- apply(nf, 1, function(x) exp(mean(log(x[is.finite(x)]))))
  nf <- nf / gm
  nf[!is.finite(nf)] <- 1
  cat(sprintf("%d loops, %d libraries, %d terms\n", nrow(cnt), N,
              uniqueN(lib$term)))

  terms <- sort(unique(lib$term))
  code <- setNames(sprintf("T%02d", seq_along(terms)), terms)
  cd <- data.frame(term = factor(code[lib$term], levels = code),
                   row.names = lib$accession)
  storage.mode(cnt) <- "integer"
  dimnames(cnt) <- list(loops$loop_id, lib$accession)
  dds <- DESeqDataSetFromMatrix(cnt, cd, ~ term)
  normalizationFactors(dds) <- nf
  if (any(w == 0)) assays(dds, withDimnames = FALSE)[["weights"]] <- w
  dds <- tryCatch(DESeq(dds, quiet = TRUE),
                  error = function(e) DESeq(dds, fitType = "local", quiet = TRUE))
  co <- attr(dispersionFunction(dds), "coefficients")
  if (!is.null(co))
    cat(sprintf("dispersion trend: asymptotic %.4f (CV %.2f at high counts), extra-Poisson %.3f\n",
                co[1], sqrt(co[1]), co[2]))

  present <- sapply(terms, function(t) {
    idx <- lib$term == t
    rowSums(w[, idx, drop = FALSE]) > 0
  })
  called <- lapply(strsplit(loops$called, "|", fixed = TRUE), unique)
  pr <- t(combn(terms, 2))
  cat(sprintf("%d pairs of terms\n", nrow(pr)))
  out <- mclapply(seq_len(nrow(pr)), function(i) {
    a <- pr[i, 1]; b <- pr[i, 2]
    ok <- present[, a] & present[, b]
    if (!any(ok)) return(NULL)
    ct <- c("term", code[[a]], code[[b]])
    r <- results(dds, contrast = ct)
    s <- tryCatch(lfcShrink(dds, contrast = ct, res = r, type = "ashr",
                            quiet = TRUE), error = function(e) NULL)
    if (is.null(s)) return(NULL)
    data.table(term_a = a, term_b = b, loop_id = loops$loop_id,
               gene = loops$gene, log2fc = s$log2FoldChange, lfcSE = s$lfcSE,
               pvalue = r$pvalue, padj = r$padj, baseMean = r$baseMean,
               called_a = as.integer(sapply(called, function(x) a %in% x)),
               called_b = as.integer(sapply(called, function(x) b %in% x)))[ok]
  }, mc.cores = opt$cores)
  out <- rbindlist(out)
  fwrite(out, file.path(opt$out, "deseq_pairs.tsv.gz"), sep = "\t")

  # viewpoint present in a term: at least one kept library has its profile
  genes <- rownames(ex[[1]]$m)
  pres <- sapply(terms, function(t) {
    acc <- lib[term == t]$accession
    as.integer(Reduce(`|`, lapply(acc, function(a)
      genes %in% rownames(ex[[a]]$m)[rowSums(is.finite(ex[[a]]$m)) > 0])))
  })
  fwrite(data.table(gene = genes, pres), file.path(opt$out, "presence.tsv"),
         sep = "\t")
  cat(sprintf("%d loop x pair rows written\n", nrow(out)))
}

# ---------------------------------------------------------------------------
# assoc: loop strength regressed on mtCN, one observation per library

# weighted GLS for every loop at once, with a common within-block correlation
# rho and fixed observation weights: the same model as limma::lmFit(block,
# correlation, weights), written for speed because it is refitted at every
# permutation. Returns the coefficient of column `col`, its t, and the residual
# variance, for the moderation that follows.
gls_fit <- function(E, W, X, block, rho, col) {
  n <- ncol(E); p <- ncol(X)
  V <- outer(block, block, "==") * rho; diag(V) <- 1
  out_b <- out_s2 <- out_se <- numeric(nrow(E))
  for (i in seq_len(nrow(E))) {
    s <- 1 / sqrt(W[i, ])
    Vi <- V * outer(s, s)
    R <- chol(Vi)
    Xs <- backsolve(R, X, transpose = TRUE)
    ys <- backsolve(R, E[i, ], transpose = TRUE)
    qx <- qr(Xs)
    b <- qr.coef(qx, ys)
    r <- ys - Xs %*% b
    s2 <- sum(r^2) / (n - p)
    XtXi <- chol2inv(qr.R(qx))
    out_b[i] <- b[col]
    out_s2[i] <- s2
    out_se[i] <- sqrt(XtXi[col, col])      # unscaled standard error
  }
  list(b = out_b, s2 = out_s2, useful = out_se, df = n - p)
}

# The same GLS, prepared once for permutations of the mtCN column.
# Only that column changes between permutations, and it is a relabelling of
# the G values held by the GTEx tissues: x = A v, with A the library-by-tissue
# indicator. Per loop, with R the Cholesky factor of the loop's covariance and
# Q an orthonormal basis of the whitened fixed columns (intercept, depth),
# P = (I - QQ') R^-T A is fixed; by the Frisch-Waugh theorem the mtCN
# coefficient is (c'v) / (v'Kv) with c = P'y_res and K = P'P, and the residual
# variance follows from ||y_res||^2. Each permutation is then two
# matrix-vector products over all loops: the result is identical to gls_fit,
# for a small fraction of the time.
gls_prepare <- function(E, W, Z, A, block, rho) {
  n <- ncol(E); L <- nrow(E); G <- ncol(A); q <- ncol(Z)
  V <- outer(block, block, "==") * rho; diag(V) <- 1
  C <- matrix(0, L, G); K <- matrix(0, L, G * G); yy <- numeric(L)
  for (i in seq_len(L)) {
    s <- 1 / sqrt(W[i, ])
    R <- chol(V * outer(s, s))
    Zs <- backsolve(R, Z, transpose = TRUE)
    As <- backsolve(R, A, transpose = TRUE)
    ys <- backsolve(R, E[i, ], transpose = TRUE)
    Q <- qr.Q(qr(Zs))
    P <- As - Q %*% crossprod(Q, As)
    yr <- ys - Q %*% crossprod(Q, ys)
    C[i, ] <- crossprod(P, yr)
    K[i, ] <- as.vector(crossprod(P))
    yy[i] <- sum(yr^2)
  }
  list(C = C, K = K, yy = yy, df = n - q - 1)
}

gls_fast <- function(prep, v) {
  xy <- as.vector(prep$C %*% v)
  xx <- as.vector(prep$K %*% as.vector(outer(v, v)))
  b <- xy / xx
  s2 <- pmax(prep$yy - b^2 * xx, 1e-12) / prep$df
  list(b = b, s2 = s2, useful = 1 / sqrt(xx), df = prep$df)
}

# empirical Bayes moderation of the variances (limma::squeezeVar), then t
moderated_t <- function(g) {
  sv <- limma::squeezeVar(g$s2, df = g$df)
  t <- g$b / (g$useful * sqrt(sv$var.post))
  attr(t, "df") <- g$df + min(sv$df.prior, 1e6)
  t
}

read_sets <- function(path, genes, min_loops) {
  s <- fread(path, sep = "\t")
  s <- s[gene %in% genes]
  s[, key := paste(source, set, sep = "::")]
  idx <- split(s$gene, s$key)
  idx[vapply(idx, function(x) sum(genes %in% x) >= min_loops, logical(1))]
}

cmd_assoc <- function(opt) {
  suppressPackageStartupMessages(library(limma))
  set.seed(opt$seed)
  lib <- read_libraries(opt$libraries)
  ex <- read_export(opt$export, lib$accession)
  lib <- lib[accession %in% names(ex)]
  opt$res <- diff(ex[[1]]$offsets[1:2])
  loops <- fread(file.path(opt$out, "loops.tsv"), sep = "\t")
  if (opt$min_terms > 1) {
    n0 <- nrow(loops)
    loops <- loops[n_terms >= opt$min_terms]
    cat(sprintf("%d of %d loops called in at least %d terms\n", nrow(loops),
                n0, opt$min_terms))
  }
  lm <- loop_matrix(opt, lib, ex, loops)

  # loops measured in every library, with enough contacts to be measured
  complete <- rowSums(lm$w) == ncol(lm$w)
  enough <- apply(lm$cnt, 1, median) >= opt$min_count
  keep <- complete & enough
  cat(sprintf("%d loops; %d measured in every library, %d of them with a median count >= %g\n",
              nrow(loops), sum(complete), sum(keep), opt$min_count))
  if (sum(keep) < opt$min_loops)
    die(sprintf("%d loops left, fewer than --min-loops %g", sum(keep), opt$min_loops))
  cnt <- lm$cnt[keep, , drop = FALSE]; nf <- lm$nf[keep, , drop = FALSE]
  loops <- loops[keep]
  gm <- apply(nf, 1, function(x) exp(mean(log(x))))
  nf <- nf / gm

  # log2 of observed over expected contact at the loop
  E <- log2((cnt + opt$pseudocount) / nf)
  dimnames(E) <- list(loops$loop_id, lib$accession)

  # reference loops (the control promoters), when given: the other loops are
  # compared with them, gene by gene and set by set
  ref_idx <- rep(FALSE, nrow(loops))
  if (!is.null(opt$reference)) {
    ref_genes <- unique(scan(opt$reference, what = "", quiet = TRUE))
    ref_idx <- loops$gene %in% ref_genes
    cat(sprintf("reference: %d loops at %d promoters of %s; %d loops tested against them\n",
                sum(ref_idx), uniqueN(loops$gene[ref_idx]), opt$reference,
                sum(!ref_idx)))
    if (sum(ref_idx) < opt$min_ref_loops)
      die(sprintf("%d reference loops, fewer than --min-ref-loops %g",
                  sum(ref_idx), opt$min_ref_loops))
  }
  test_idx <- if (any(ref_idx)) !ref_idx else rep(TRUE, nrow(loops))

  mt_lib <- log2(lib$mtcn)
  dep <- log10(lib$depth)
  block <- as.integer(factor(lib$gtex))
  gtex <- sort(unique(lib$gtex))
  mt_gtex <- sapply(gtex, function(g) median(log2(lib[gtex == g]$mtcn)))
  gi <- match(lib$gtex, gtex)
  one <- rep(1, nrow(lib))
  covs <- if (identical(opt$adjust_depth, "no")) NULL else cbind(depth = dep)
  design <- function(x) cbind(intercept = one, mtcn = x, covs)
  X0 <- cbind(intercept = one, covs)

  # observation weights from the mean-variance trend and the correlation
  # between libraries of one GTEx tissue, both estimated without mtCN in the
  # model, so that they stay the same under every permutation
  v <- vooma(E, X0)
  cf <- duplicateCorrelation(v$E, X0, block = block, weights = v$weights)
  v <- vooma(E, X0, block = block, correlation = cf$consensus)
  cf <- duplicateCorrelation(v$E, X0, block = block, weights = v$weights)
  rho <- cf$consensus
  W <- v$weights
  cat(sprintf("%d libraries, %d terms, %d GTEx tissues (%d distinct mtCN values)\n",
              nrow(lib), uniqueN(lib$term), length(gtex),
              length(unique(mt_gtex))))
  cat(sprintf("correlation between libraries of one GTEx tissue: %.2f\n", rho))

  A <- sapply(gtex, function(g) as.numeric(lib$gtex == g))
  if (max(abs(A %*% mt_gtex - mt_lib)) > 1e-9)
    die("mtCN differs between libraries of one GTEx tissue")
  t0 <- Sys.time()
  prep <- gls_prepare(E, W, X0, A, block, rho)
  cat(sprintf("model prepared for %d loops in %.0f s\n", nrow(E),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  obs_fit <- gls_fast(prep, mt_gtex)
  t_obs <- moderated_t(obs_fit)
  df_t <- attr(t_obs, "df")
  # limma's own fit, for the reported per-loop p-values
  fit <- eBayes(lmFit(v, design(mt_lib), block = block, correlation = rho))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  if (max(abs(tt$t - t_obs)) > 1e-4 * max(1, max(abs(t_obs))))
    warning("internal GLS and limma disagree; check the input")

  # sets and genes, as index vectors into the loops
  sets <- if (!is.null(opt$sets)) read_sets(opt$sets, loops$gene[test_idx], opt$min_set)
          else list()
  set_idx <- lapply(sets, function(s) which(test_idx & loops$gene %in% s))
  gene_idx <- if (any(ref_idx)) split(which(test_idx), loops$gene[test_idx]) else list()
  # Genes and sets are compared on the change in contact per doubling of mtCN
  # (the coefficient), not on t: t grows with how well a loop is measured, so
  # a well-covered loop would look more associated than an average reference
  # loop for the same change. Each loop is also compared with reference loops
  # of similar contact level (--ref-strata quantile strata of the reference
  # loops' mean log2 O/E, deciles by default),
  # in case the shared change differs between strong and weak loops.
  lvl <- rowMeans(E)
  n_str <- as.integer(opt$ref_strata)
  if (any(ref_idx)) {
    br <- quantile(lvl[ref_idx], seq(0, 1, length.out = n_str + 1))
    br[1] <- -Inf; br[length(br)] <- Inf
    stratum <- as.integer(cut(lvl, unique(br), include.lowest = TRUE))
  }
  comp <- function(b, idx) {
    if (!length(idx)) return(numeric(0))
    if (any(ref_idx)) {
      rm <- tapply(b[ref_idx], factor(stratum[ref_idx],
                                       levels = sort(unique(stratum))), mean)
      adj <- b - rm[as.character(stratum)]
      vapply(idx, function(i) mean(adj[i], na.rm = TRUE), numeric(1))
    } else {
      vapply(idx, function(i) mean(b[i]) - mean(b[setdiff(which(test_idx), i)]),
             numeric(1))
    }
  }

  n_hits <- function(t) sum(p.adjust(2 * pt(-abs(t), df = df_t), "BH") <= opt$fdr)
  # mean t^2 once the k loops with the largest |t| are set aside, k being the
  # number of loops that pass the FDR in the observed data. Computed the same
  # way in every permutation, it asks whether the association survives
  # without the loops that are individually significant.
  global <- function(t, k) {
    t <- t[test_idx]
    t2 <- sort(t^2, decreasing = TRUE)
    c(mean_t2 = mean(t2), mean_t = mean(t), n_fdr = n_hits(t),
      mean_t2_trim = mean(t2[-seq_len(k)]))
  }
  k_trim <- max(n_hits(t_obs[test_idx]), 1L)
  g_obs <- global(t_obs, k_trim)
  s_obs <- comp(obs_fit$b, set_idx)
  q_obs <- comp(obs_fit$b, gene_idx)

  # permutation: mtCN moved between GTEx tissues, the libraries of a tissue
  # moving together; depth, the shared shift and the correlation structure
  # stay in place
  B <- opt$perms
  chunks <- split(seq_len(B), rep_len(seq_len(opt$cores), B))
  seeds <- sample.int(1e9, length(chunks))
  res <- mclapply(seq_along(chunks), function(ci) {
    set.seed(seeds[ci])
    n <- length(chunks[[ci]])
    exceed <- numeric(length(t_obs))
    G <- matrix(NA_real_, n, length(g_obs))
    S <- matrix(NA_real_, n, length(s_obs))
    Q <- matrix(NA_real_, n, length(q_obs))
    for (k in seq_len(n)) {
      fp <- gls_fast(prep, mt_gtex[sample.int(length(gtex))])
      tp <- moderated_t(fp)
      exceed <- exceed + (abs(tp) >= abs(t_obs))
      G[k, ] <- global(tp, k_trim)
      if (length(s_obs)) S[k, ] <- comp(fp$b, set_idx)
      if (length(q_obs)) Q[k, ] <- comp(fp$b, gene_idx)
    }
    list(exceed = exceed, G = G, S = S, Q = Q)
  }, mc.cores = opt$cores)
  exceed <- Reduce(`+`, lapply(res, `[[`, "exceed"))
  G <- do.call(rbind, lapply(res, `[[`, "G")); colnames(G) <- names(g_obs)
  p_loop <- (1 + exceed) / (B + 1)

  two_sided <- function(obs, null) {
    c0 <- mean(null)
    (1 + sum(abs(null - c0) >= abs(obs - c0))) / (length(null) + 1)
  }
  one_sided <- function(obs, null) (1 + sum(null >= obs)) / (length(null) + 1)
  # the permutation p cannot go below 1/(B+1); a z against the permutation
  # null (its own mean and spread) gives a continuous p for the FDR over many
  # genes
  z_of <- function(obs, null) (obs - mean(null)) / sd(null)
  glob <- data.table(
    statistic = names(g_obs), observed = g_obs,
    null_median = apply(G, 2, median),
    null_q95 = apply(G, 2, quantile, 0.95),
    p_perm = c(one_sided(g_obs[1], G[, 1]), two_sided(g_obs[2], G[, 2]),
               one_sided(g_obs[3], G[, 3]), one_sided(g_obs[4], G[, 4])))

  dist <- loops$summit - loops$vp_pos
  tab <- data.table(loop_id = loops$loop_id, gene = loops$gene,
                    reference = as.integer(ref_idx),
                    chrom = loops$chrom, vp_pos = loops$vp_pos,
                    summit = loops$summit, distance = dist,
                    n_terms = loops$n_terms,
                    mean_log2_oe = rowMeans(E),
                    log2fc_per_doubling = tt$logFC, t = tt$t,
                    p_limma = tt$P.Value, fdr_limma = tt$adj.P.Val,
                    p_perm = p_loop, fdr_perm = p.adjust(p_loop, "BH"))
  dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
  fwrite(tab[order(p_limma)], file.path(opt$out, "assoc_loops.tsv"), sep = "\t")
  fwrite(glob, file.path(opt$out, "assoc_global.tsv"), sep = "\t")
  fwrite(as.data.table(G), file.path(opt$out, "assoc_null.tsv"), sep = "\t")

  cat(sprintf("\nglobal tests on the %d tested loops, %d permutations of mtCN between GTEx tissues\n",
              sum(test_idx), B))
  lab <- c(mean_t2 = "association, either direction (mean t^2)",
           mean_t = "direction (mean t; > 0 = loops stronger with mtCN)",
           n_fdr = sprintf("loops at FDR <= %g (limma)", opt$fdr),
           mean_t2_trim = sprintf("association without the %d top loops", k_trim))
  for (i in seq_len(nrow(glob)))
    cat(sprintf("  %-52s observed %8.3f   null median %8.3f   p = %.4f\n",
                lab[glob$statistic[i]], glob$observed[i], glob$null_median[i],
                glob$p_perm[i]))
  if (any(ref_idx))
    cat(sprintf("  mean t: tested loops %+.2f, reference loops %+.2f\n",
                mean(t_obs[test_idx]), mean(t_obs[ref_idx])))

  summarise <- function(M, obs, idx, names_) {
    if (!length(obs)) return(NULL)
    M <- do.call(rbind, M)
    z <- vapply(seq_along(obs), function(i) z_of(obs[i], M[, i]), numeric(1))
    out <- data.table(name = names_,
                      loops = lengths(idx),
                      viewpoints = vapply(idx, function(i) uniqueN(loops$gene[i]),
                                          integer(1)),
                      mean_t = vapply(idx, function(i) mean(t_obs[i]), numeric(1)),
                      mean_log2fc_per_doubling = vapply(idx, function(i)
                        mean(tt$logFC[i]), numeric(1)),
                      difference_log2fc = obs, z = z,
                      p_perm = vapply(seq_along(obs), function(i)
                        two_sided(obs[i], M[, i]), numeric(1)),
                      p_z = 2 * pnorm(-abs(z)))
    out[, fdr := p.adjust(p_z, "BH")]
    out[, direction := ifelse(fdr > opt$fdr, "as the reference",
                        ifelse(difference_log2fc > 0, "follows mtCN more",
                               "follows mtCN less"))]
    out[order(p_z)]
  }
  ref_label <- if (any(ref_idx)) "the reference loops" else "the other loops"
  if (length(s_obs)) {
    st <- summarise(lapply(res, `[[`, "S"), s_obs, set_idx, names(sets))
    setnames(st, "name", "set")
    fwrite(st, file.path(opt$out, "assoc_sets.tsv"), sep = "\t")
    cat(sprintf("\n%d gene sets with >= %d loops, each against %s; %d at FDR <= %g\n",
                nrow(st), opt$min_set, ref_label, sum(st$fdr <= opt$fdr), opt$fdr))
    print(head(st[, .(set, loops, mean_log2fc_per_doubling, difference_log2fc,
                      p_z, fdr, direction)], 12))
  }
  if (length(q_obs)) {
    gt <- summarise(lapply(res, `[[`, "Q"), q_obs, gene_idx, names(gene_idx))
    setnames(gt, "name", "gene")
    fwrite(gt, file.path(opt$out, "assoc_genes.tsv"), sep = "\t")
    cat(sprintf("\n%d genes, each against %s: %d follow mtCN more, %d less (FDR <= %g)\n",
                nrow(gt), ref_label, sum(gt$direction == "follows mtCN more"),
                sum(gt$direction == "follows mtCN less"), opt$fdr))
    print(head(gt[, .(gene, loops, mean_log2fc_per_doubling, difference_log2fc,
                      p_z, fdr, direction)], 12))
  }
  draw_assoc(file.path(opt$out, "assoc_volcano.svg"), tab[test_idx][order(p_limma)],
             glob, opt$fdr, B)
  cat(sprintf("\nWritten to %s/\n", opt$out))
}

draw_assoc <- function(path, tab, glob, fdr, B) {
  W <- 720; H <- 520; pad <- 80; pw <- W - 2 * pad; ph <- H - 2 * pad
  x <- tab$log2fc_per_doubling; y <- -log10(tab$p_limma)
  xr <- max(abs(x)) * 1.05; yr <- max(y, 1) * 1.05
  sx <- function(v) pad + (v + xr) / (2 * xr) * pw
  sy <- function(v) pad + ph - v / yr * ph
  hit <- tab$fdr_limma <= fdr
  s <- c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" font-family="Helvetica,Arial" font-size="11">', W, H, W, H),
         sprintf('<rect width="%d" height="%d" fill="white"/>', W, H),
         sprintf('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#999"/>', pad, pad, pw, ph),
         sprintf('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#bbb"/>', sx(0), pad, sx(0), pad + ph))
  for (i in order(hit)) {
    s <- c(s, sprintf('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" fill-opacity="%.2f"><title>%s %.3f p=%.2g</title></circle>',
                      sx(x[i]), sy(y[i]), if (hit[i]) 3.2 else 2.2,
                      if (hit[i]) "#c44" else "#888", if (hit[i]) 0.85 else 0.35,
                      tab$loop_id[i], x[i], tab$p_limma[i]))
  }
  for (v in pretty(c(-xr, xr), 5)) if (abs(v) <= xr)
    s <- c(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle">%g</text>', sx(v), pad + ph + 16, v))
  for (v in pretty(c(0, yr), 5)) if (v <= yr)
    s <- c(s, sprintf('<text x="%d" y="%.1f" text-anchor="end">%g</text>', pad - 6, sy(v) + 4, v))
  g1 <- glob[statistic == "mean_t2"]; g2 <- glob[statistic == "mean_t"]
  s <- c(s,
    sprintf('<text x="%d" y="%d" font-size="12">association of loop strength with mtCN: p = %.3f; direction: mean t %+.2f, p = %.3f</text>',
            pad, pad - 26, g1$p_perm, g2$observed, g2$p_perm),
    sprintf('<text x="%d" y="%d" fill="#666">one point per loop; red: FDR &lt;= %g; global p from %d permutations of mtCN between GTEx tissues</text>',
            pad, pad - 10, fdr, B),
    sprintf('<text x="%.0f" y="%d" text-anchor="middle">change in loop contact per doubling of mtDNA copy number (log2)</text>', pad + pw / 2, H - 24),
    sprintf('<text x="24" y="%.0f" text-anchor="middle" transform="rotate(-90 24 %.0f)">-log10 p (limma, libraries of one GTEx tissue as blocks)</text>', pad + ph / 2, pad + ph / 2),
    "</svg>")
  writeLines(s, path)
}


# ---------------------------------------------------------------------------
# tissue: loop strength against mtDNA copy number, one point per tissue
#
# Same loops, same normalisation and same depth adjustment as `assoc`, summed
# into one number per library instead of one test per loop. Each loop is
# centred on its own mean across libraries, so a loop that is strong
# everywhere does not dominate the average; what is left is how far that library
# sits from the catalogue's usual level. Libraries of one tissue are then
# averaged, since mtCN is measured per GTEx tissue and not per library.

theil_sen <- function(x, y) {
  n <- length(x)
  if (n < 3) return(c(NA_real_, NA_real_))
  sl <- c()
  for (i in 1:(n - 1)) for (j in (i + 1):n)
    if (x[j] != x[i]) sl <- c(sl, (y[j] - y[i]) / (x[j] - x[i]))
  if (!length(sl)) return(c(NA_real_, NA_real_))
  b <- median(sl)
  c(median(y - b * x), b)
}

spearman <- function(x, y) {
  if (length(x) < 3) return(NA_real_)
  suppressWarnings(cor(rank(x), rank(y)))
}

cmd_tissue <- function(opt) {
  set.seed(opt$seed)
  lib <- read_libraries(opt$libraries)
  ex <- read_export(opt$export, lib$accession)
  lib <- lib[accession %in% names(ex)]
  opt$res <- diff(ex[[1]]$offsets[1:2])
  loops <- fread(file.path(opt$out, "loops.tsv"), sep = "\t")
  if (opt$min_terms > 1) loops <- loops[n_terms >= opt$min_terms]
  lm <- loop_matrix(opt, lib, ex, loops)
  keep <- rowSums(lm$w) == ncol(lm$w) &
    apply(lm$cnt, 1, median) >= opt$min_count
  if (sum(keep) < opt$min_loops)
    die(sprintf("%d loops left, fewer than --min-loops %g", sum(keep), opt$min_loops))
  cnt <- lm$cnt[keep, , drop = FALSE]; nf <- lm$nf[keep, , drop = FALSE]
  loops <- loops[keep]
  if (!is.null(opt$reference)) {
    ref <- unique(scan(opt$reference, what = "", quiet = TRUE))
    drop <- loops$gene %in% ref
    cat(sprintf("%d loops at reference promoters left out\n", sum(drop)))
    cnt <- cnt[!drop, , drop = FALSE]; nf <- nf[!drop, , drop = FALSE]
    loops <- loops[!drop]
  }
  cat(sprintf("%d loops at %d promoters, %d libraries\n", nrow(loops),
              uniqueN(loops$gene), nrow(lib)))

  gm <- apply(nf, 1, function(x) exp(mean(log(x))))
  E <- log2((cnt + opt$pseudocount) / (nf / gm))
  E <- E - rowMeans(E)                      # each loop on its own scale

  if (identical(opt$adjust_depth, "yes")) {
    d <- lib$log10_depth
    d <- d - mean(d)
    E <- E - outer(as.vector((E %*% d) / sum(d^2)), d)
    cat("depth removed from every loop before averaging\n")
  }
  strength <- colMeans(E)

  per_lib <- data.table(accession = lib$accession, term = lib$term,
                        gtex = lib$gtex, mtcn = lib$mtcn,
                        depth = lib$depth, strength = strength)
  per_tis <- per_lib[, .(libraries = .N, mtcn = mtcn[1],
                         strength = mean(strength),
                         sd = if (.N > 1) sd(strength) else NA_real_),
                     by = gtex][order(mtcn)]

  rho <- spearman(log2(per_tis$mtcn), per_tis$strength)
  ts <- theil_sen(log2(per_tis$mtcn), per_tis$strength)
  cat(sprintf("\n%d tissues, Spearman rho = %+.3f, Theil-Sen slope = %+.4f log2 per doubling of mtCN\n",
              nrow(per_tis), rho, ts[2]))
  print(per_tis)

  dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
  fwrite(per_lib, file.path(opt$out, "tissue_strength_libraries.tsv"), sep = "\t")
  fwrite(per_tis, file.path(opt$out, "tissue_strength.tsv"), sep = "\t")
  svg_tissue(per_tis, per_lib, rho, ts, nrow(loops),
             file.path(opt$out, "tissue_strength.svg"))
  cat(sprintf("\nWritten to %s/\n", opt$out))
}

svg_tissue <- function(tis, libs, rho, ts, n_loops, path) {
  W <- 760; H <- 560; pad <- 86
  pw <- W - pad - 30; ph <- H - pad - 70
  lx <- log2(libs$mtcn); ly <- libs$strength
  tx <- log2(tis$mtcn); ty <- tis$strength
  xr <- range(c(lx, tx), na.rm = TRUE); yr <- range(c(ly, ty), na.rm = TRUE)
  xr <- xr + c(-1, 1) * diff(xr) * 0.06
  yr <- yr + c(-1, 1) * diff(yr) * 0.10
  px <- function(v) pad + (v - xr[1]) / diff(xr) * pw
  py <- function(v) pad + ph - (v - yr[1]) / diff(yr) * ph
  s <- c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" font-family="Helvetica,Arial,sans-serif" font-size="12">', W, H),
         sprintf('<rect width="%d" height="%d" fill="white"/>', W, H))
  ticks <- pretty(xr, 6)
  for (t in ticks) if (t >= xr[1] && t <= xr[2]) {
    s <- c(s, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#e6e6e6"/>', px(t), pad, px(t), pad + ph),
           sprintf('<text x="%.1f" y="%.1f" text-anchor="middle" fill="#444">%.0f</text>', px(t), pad + ph + 18, 2^t))
  }
  for (t in pretty(yr, 5)) if (t >= yr[1] && t <= yr[2]) {
    s <- c(s, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#e6e6e6"/>', pad, py(t), pad + pw, py(t)),
           sprintf('<text x="%.1f" y="%.1f" text-anchor="end" fill="#444">%+.2f</text>', pad - 8, py(t) + 4, t))
  }
  for (i in seq_along(lx))
    s <- c(s, sprintf('<circle cx="%.1f" cy="%.1f" r="2.4" fill="#bbbbbb"/>', px(lx[i]), py(ly[i])))
  if (is.finite(ts[2])) {
    x0 <- xr[1]; x1 <- xr[2]
    s <- c(s, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#444" stroke-width="1.5" stroke-dasharray="5,4"/>',
                      px(x0), py(ts[1] + ts[2] * x0), px(x1), py(ts[1] + ts[2] * x1)))
  }
  for (i in seq_along(tx)) {
    s <- c(s, sprintf('<circle cx="%.1f" cy="%.1f" r="5.5" fill="#B2182B" stroke="white" stroke-width="1.5"/>', px(tx[i]), py(ty[i])),
           sprintf('<text x="%.1f" y="%.1f" font-size="10" fill="#333">%s</text>', px(tx[i]) + 8, py(ty[i]) + 3.5, tis$gtex[i]))
  }
  s <- c(s,
    sprintf('<line x1="%d" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#333"/>', pad, pad + ph, pad + pw, pad + ph),
    sprintf('<line x1="%d" y1="%d" x2="%d" y2="%.1f" stroke="#333"/>', pad, pad, pad, pad + ph),
    sprintf('<text x="%d" y="%d" font-size="13">loop strength at mitochondrial promoters against mitochondrial DNA copy number</text>', pad, pad - 44),
    sprintf('<text x="%d" y="%d">%d loops, one red point per GTEx tissue, grey points the individual libraries</text>', pad, pad - 26, n_loops),
    sprintf('<text x="%d" y="%d" fill="#666">Spearman rho = %+.2f; Theil-Sen slope %+.3f log2 per doubling</text>', pad, pad - 10, rho, ts[2]),
    sprintf('<text x="%.0f" y="%d" text-anchor="middle">mitochondrial DNA copies per cell (log scale)</text>', pad + pw / 2, H - 24),
    sprintf('<text x="22" y="%.0f" text-anchor="middle" transform="rotate(-90 22 %.0f)">mean loop strength, log2 over the catalogue average</text>', pad + ph / 2, pad + ph / 2),
    "</svg>")
  writeLines(s, path)
}

# ---------------------------------------------------------------------------

opt <- parse_args(commandArgs(trailingOnly = TRUE))
switch(opt$cmd,
       scan = cmd_scan(opt),
       call = cmd_call(opt),
       reciprocal = cmd_reciprocal(opt),
       deseq = cmd_deseq(opt),
       assoc = cmd_assoc(opt),
       tissue = cmd_tissue(opt))
