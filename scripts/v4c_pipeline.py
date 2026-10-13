#!/usr/bin/env python3
"""
v4c_pipeline.py - virtual 4C profiles from Hi-C libraries, library selection,
export for peakC, and a pairwise comparison of terms against mtDNA copy number.

Purpose
-------
Python side of the mitoloop_v0 loop pipeline. For every Hi-C library of a
library table, extract one virtual-4C profile per viewpoint (the TSS of a
mitochondrial or control gene, or a partner anchor), keep the libraries deep
enough to call loops, and write the profiles as tables for v4c_loops.R, which
calls loops with peakC, filters them by reciprocity and relates their strength
to mtDNA copy number. The `relate` subcommand works on the outputs of
`v4c_loops.R deseq`.

Subcommands
-----------
  extract  Raw virtual-4C profiles: contacts between an anchor (TSS +/-
           --anchor-half, or a partner bin) and every bin within +/- --window,
           observed counts without normalisation, summed over the anchor
           bins. One npz per library; --index runs one library (for job
           arrays). Anchors come from the GTF (mito and control genes) or from
           a TSV (the partner anchors written by `v4c_loops.R call`).
  select   Library depth (contacts over the --depth-set viewpoints) and the
           depth filter: libraries below median - --k-sd SD of log10 depth are
           dropped, but never below --min-keep libraries; if the cutoff would
           leave fewer, the --min-keep deepest are kept.
  export   One table per kept library for R: a viewpoint per row, a bin per
           column, for the viewpoint set(s) of --set (mito, partner, control,
           or a comma-separated list such as mito,control). When the exported sets include mito or control viewpoints,
           also decay.tsv: per library, the counts summed over the control
           viewpoints at each distance, used by v4c_loops.R for the
           distance-matched expected count.
  relate   Pairwise route between terms: from the DESeq2 results of
           `v4c_loops.R deseq`, the change in looping between every two terms,
           against their difference in mtDNA copy number.

Inputs
------
  --tissues   library table (TSV), one row per Hi-C library, with columns
              accession, term_name (ENCODE biosample term), gtex_tissue
              (matching GTEx tissue), mtcn (mtDNA copy number of that GTEx
              tissue) and at least one location column. Location columns are
              tried in the order of --location-columns (default s3_url,
              download_url, path); a value is used when it is an http(s) URL
              or an existing local file. An optional `format` column (hic,
              cool, mcool) overrides the format guessed from the extension.
              Column names are options (--accession-column, --term-column,
              --gtex-column, --mtcn-column).
  --gtf       GENCODE GTF (gzipped or plain), with --genes (mitochondrial gene
              symbols, one per line) and --controls (control symbols). The TSS
              is that of the MANE Select transcript when tagged, otherwise
              that of the gene record.
  --anchors   instead of the GTF, a TSV with name, chrom, pos[, set] (the
              partners.tsv written by `v4c_loops.R call`).
  Hi-C        .hic (local or URL, read with hicstraw), .cool or .mcool (local,
              read with cooler; raw counts in the `count` column).

Outputs
-------
  extract  <outdir>/<accession>.npz with keys counts (viewpoint x bin, raw
           counts, NaN where unavailable or below --min-raw), raw_total,
           names, chrom, pos, set (mito, control or partner), offsets (bp),
           accession, term, gtex, mtcn, resolution, anchor_half, via (the
           location column used).
  select   libraries.tsv: accession, term, gtex, mtcn, depth, log10_depth,
           kept (1 or 0).
  export   <outdir>/<accession>.tsv.gz: name, chrom, pos, then one column per
           offset in bp; decay.tsv: accession, then one column per offset.
  relate   relate_stats.tsv (statistic, rho, p_perm, rho_power_partial,
           rho_vs_power), relate_pairs.tsv (one row per pair of terms) and
           relate_size.svg.

Method
------
Profiles are raw observed counts: peakC normalises each profile itself and
the count models of v4c_loops.R need raw counts. Neighbouring viewpoints
within --group-span are read with one query: the anchor interval of the
group (first TSS - --anchor-half to last TSS + --anchor-half) against the
window. hicstraw returns the intra-chromosomal records whose bin starts fall
in the two intervals, in either orientation. A cooler stores the upper
triangle only, so the cooler reader queries the whole window against itself
and keeps the same records, which makes the two formats give identical
counts on identical contacts. The anchor band itself carries no usable value
and is left as NaN.

`relate` sums, for every pair of terms, the absolute shrunken log2 fold
changes over the loops, the number of significant changes (--fdr, --min-fold)
and their signed sum oriented towards the term of higher mtDNA copy number,
each per shared viewpoint, and correlates them (Spearman) with the difference
in mtDNA copy number. p-values come from permutations of mtDNA copy number
between GTEx tissues, the terms of one tissue moving together. A partial
Spearman on the statistical power of the pair (log harmonic mean of depths)
is reported beside it.

Usage
-----
    python3 v4c_pipeline.py extract --tissues tissues.tsv \\
        --gtf data/gencode.gtf.gz --genes genes.txt --controls controls.txt \\
        --outdir raw [--index N]
    python3 v4c_pipeline.py select --tissues tissues.tsv --indir raw \\
        --out libraries.tsv
    python3 v4c_pipeline.py export --indir raw --libraries libraries.tsv \\
        --set mito --outdir export/mito
    Rscript v4c_loops.R scan --libraries libraries.tsv --export export/mito
    Rscript v4c_loops.R call --libraries libraries.tsv --export export/mito \\
        --wsize 3 --out loops
    python3 v4c_pipeline.py extract --tissues tissues.tsv \\
        --anchors loops/partners.tsv --outdir raw_partners [--index N]
    python3 v4c_pipeline.py export --indir raw_partners \\
        --libraries libraries.tsv --set partner --outdir export/partners
    Rscript v4c_loops.R reciprocal --libraries libraries.tsv \\
        --export export/partners --wsize 3 --out loops
    Rscript v4c_loops.R assoc --libraries libraries.tsv \\
        --export export/mito --wsize 3 --out loops
    # pairwise route
    Rscript v4c_loops.R deseq --libraries libraries.tsv \\
        --export export/mito --wsize 3 --out loops
    python3 v4c_pipeline.py relate --libraries libraries.tsv --loops loops \\
        --outdir relate

References
----------
  Durand et al. 2016. Juicer provides a one-click system for analyzing
      loop-resolution Hi-C experiments. Cell Systems. (.hic format, hicstraw)
  Abdennur and Mirny 2020. Cooler: scalable storage for Hi-C data and other
      genomically labeled arrays. Bioinformatics.
  ENCODE Project Consortium et al. 2020. Expanded encyclopaedias of DNA
      elements in the human and mouse genomes. Nature. (intact Hi-C data)
  Frankish et al. 2023. GENCODE: reference annotation for the human and mouse
      genomes in 2023. Nucleic Acids Research.
  Phipson and Smyth 2010. Permutation P-values should never be zero.
      Statistical Applications in Genetics and Molecular Biology. (relate)
  Morales et al. 2022. A joint NCBI and EMBL-EBI transcript set for clinical
      genomics and research. Nature. (MANE Select)
  GTEx Consortium 2020. The GTEx Consortium atlas of genetic regulatory
      effects across human tissues. Science.
  Geeven et al. 2018. peakC: a flexible, non-parametric peak calling package
      for 4C and Capture-C data. Nucleic Acids Research.
  Love et al. 2014. Moderated estimation of fold change and dispersion for
      RNA-seq data with DESeq2. Genome Biology. (input of `relate`)
  Sen 1968. Estimates of the regression coefficient based on Kendall's tau.
      Journal of the American Statistical Association. (Theil-Sen line)
"""

from __future__ import annotations

import argparse
import gzip
import math
import sys
import time
from pathlib import Path

import numpy as np

WINDOW = 1_000_000
RESOLUTION = 10_000
ANCHOR_HALF = 10_000
GROUP_SPAN = 2_000_000
MIN_RAW = 50
LOCATION_COLUMNS = "s3_url,download_url,path"


# ---------------------------------------------------------------------------
# small helpers

def read_tsv(path) -> list[dict]:
    rows = []
    with open(path) as fh:
        head = fh.readline().rstrip("\r\n").split("\t")
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            f = line.split("\t")
            f += [""] * (len(head) - len(f))
            rows.append(dict(zip(head, f)))
    return rows


def load_tss(gtf: Path, wanted: set) -> list[dict]:
    """MANE Select TSS where available, otherwise the gene record."""
    best = {}
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] not in ("gene", "transcript"):
                continue
            a = f[8]
            i = a.find('gene_name "')
            if i < 0:
                continue
            name = a[i + 11:a.find('"', i + 11)]
            if name not in wanted:
                continue
            mane = "MANE_Select" in a
            rank = 2 if (f[2] == "transcript" and mane) else (1 if f[2] == "gene" else 0)
            if rank == 0:
                continue
            tss = int(f[3]) - 1 if f[6] == "+" else int(f[4])
            prev = best.get(name)
            if prev is None or rank > prev["rank"]:
                best[name] = {"gene": name, "chrom": f[0], "tss": tss,
                              "strand": f[6], "rank": rank}
    return sorted(best.values(), key=lambda g: (g["chrom"], g["tss"]))


def group_anchors(items, span: int):
    """Consecutive anchors on one chromosome whose span stays under `span`.

    Neighbouring viewpoints share the same blocks of the file, so one query
    per group instead of one per viewpoint saves most of the reading time.
    """
    groups, cur = [], []
    for it in items:
        if cur and it[1] - cur[0][1] > span:
            groups.append(cur)
            cur = []
        cur.append(it)
    if cur:
        groups.append(cur)
    return groups


def n_anchor_bins(pos, res, chrom_len, anchor_half):
    """Number of anchor bins inside the chromosome, as in profiles_for_group."""
    lo, hi = pos - anchor_half - res, pos + anchor_half + res
    return sum(1 for b in range(int(math.ceil(max(0, lo) / res)),
                                int(min(chrom_len, hi) // res) + 1)
               if lo < b * res < hi)


def rankdata(a) -> np.ndarray:
    """Ranks with ties averaged: mtDNA copy number is identical for every
    library of one tissue, so ties are the rule."""
    a = np.asarray(a, dtype=float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), dtype=float)
    sa = a[order]
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and sa[j + 1] == sa[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks


def spearman(a, b) -> float:
    ra, rb = rankdata(a), rankdata(b)
    ra -= ra.mean()
    rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d else 0.0


def partial_spearman(a, b, c) -> float:
    """Spearman of a and b with c taken out of both."""
    def z(v):
        r = rankdata(v) - rankdata(v).mean()
        return r / np.sqrt((r ** 2).sum())
    za, zb, zc = z(a), z(b), z(c)
    ra, rb = za - (za @ zc) * zc, zb - (zb @ zc) * zc
    d = np.sqrt((ra @ ra) * (rb @ rb))
    return float((ra @ rb) / d) if d else 0.0


def theil_sen(x, y):
    """Median of the pairwise slopes; returns (slope, intercept)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    sl = []
    for i in range(len(x)):
        dx = x[i + 1:] - x[i]
        ok = dx != 0
        if ok.any():
            sl.append((y[i + 1:][ok] - y[i]) / dx[ok])
    if not sl:
        return 0.0, float(np.median(y))
    sl = np.concatenate(sl)
    b = float(np.median(sl))
    return b, float(np.median(y - b * x))


# ---------------------------------------------------------------------------
# readers: .hic through hicstraw, .cool / .mcool through cooler
#
# records(chrom, g0, g1, w0, w1) returns (x, y, counts), x and y being the
# genomic start of the two bins of a contact, for every stored contact with
# one bin start in [g0, g1] and the other in [w0, w1] (bp, both ends
# included), each contact once. This is what hicstraw getRecords returns for
# an intra-chromosomal query.

def is_remote(loc: str) -> bool:
    return loc.startswith("http")


def guess_format(loc: str, declared: str = "") -> str:
    declared = (declared or "").strip().lower()
    if declared:
        return declared
    base = loc.split("?")[0].lower()
    for ext in ("mcool", "cool", "hic"):
        if base.endswith("." + ext):
            return ext
    return "hic"


class HicReader:
    def __init__(self, loc, res, retries):
        import hicstraw
        last = None
        for attempt in range(retries):
            try:
                self.h = hicstraw.HiCFile(str(loc))
                last = None
                break
            except Exception as exc:                      # noqa: BLE001
                last = exc
                print(f"    open attempt {attempt+1}: "
                      f"{type(exc).__name__}: {exc}", flush=True)
                if attempt + 1 < retries:
                    time.sleep(5 * 2 ** attempt)
        if last is not None:
            raise last
        self.res = res
        self.resolutions = set(self.h.getResolutions())
        self.chroms = {c.name: c.length for c in self.h.getChromosomes()}
        self._mzd = {}

    def prepare(self, chrom):
        self._mzd[chrom] = self.h.getMatrixZoomData(
            chrom, chrom, "observed", "NONE", "BP", self.res)

    def records(self, chrom, g0, g1, w0, w1):
        recs = self._mzd[chrom].getRecords(int(g0), int(g1), int(w0), int(w1))
        if not recs:
            return None
        x = np.array([r.binX for r in recs], dtype=np.float64)
        y = np.array([r.binY for r in recs], dtype=np.float64)
        c = np.array([r.counts for r in recs], dtype=np.float64)
        return x, y, c


class CoolReader:
    def __init__(self, loc, res, is_mcool):
        import cooler
        if is_mcool:
            have = {int(p.rsplit("/", 1)[-1])
                    for p in cooler.fileops.list_coolers(str(loc))
                    if p.rsplit("/", 1)[-1].isdigit()}
            self.resolutions = have
            if res not in have:
                self.chroms = {}
                return
            self.c = cooler.Cooler(f"{loc}::/resolutions/{int(res)}")
        else:
            self.c = cooler.Cooler(str(loc))
            self.resolutions = {int(self.c.binsize)}
        self.res = int(self.c.binsize)
        self.chroms = {k: int(v) for k, v in self.c.chromsizes.items()}
        self._sel = self.c.matrix(balance=False, as_pixels=True, join=True)

    def prepare(self, chrom):
        pass

    def records(self, chrom, g0, g1, w0, w1):
        # The whole window against itself: a cooler stores the upper triangle
        # only, so an anchor x window query would return the contacts on one
        # side of the anchor and miss those on the other.
        r = self.res
        start = int(math.ceil(w0 / r)) * r
        end = min(self.chroms[chrom], (int(w1) // r + 1) * r)
        if end <= start:
            return None
        reg = f"{chrom}:{start}-{end}"
        px = self._sel.fetch(reg, reg)
        if px is None or len(px) == 0:
            return None
        x = px["start1"].to_numpy(dtype=np.float64)
        y = px["start2"].to_numpy(dtype=np.float64)
        c = px["count"].to_numpy(dtype=np.float64)
        # the selection of hicstraw getRecords: bin starts inside the
        # requested intervals (bp, both ends included), in either orientation
        keep = (((x >= g0) & (x <= g1) & (y >= w0) & (y <= w1))
                | ((y >= g0) & (y <= g1) & (x >= w0) & (x <= w1)))
        if not keep.any():
            return None
        return x[keep], y[keep], c[keep]


def open_library(row, args):
    """First readable location among --location-columns, in order.

    Returns (reader, column name) or (None, None)."""
    for key in args.location_columns.split(","):
        key = key.strip()
        loc = (row.get(key) or "").strip()
        if not loc:
            continue
        remote = is_remote(loc)
        if not remote and not Path(loc).exists():
            print(f"    {key}: {loc} not found", flush=True)
            continue
        fmt = guess_format(loc, row.get(args.format_column, ""))
        try:
            if fmt == "hic":
                return HicReader(loc, args.resolution,
                                 args.retries if remote else 1), key
            if fmt in ("cool", "mcool"):
                if remote:
                    print(f"    {key}: cooler reads local files only",
                          flush=True)
                    continue
                return CoolReader(loc, args.resolution, fmt == "mcool"), key
            print(f"    {key}: unknown format {fmt!r} (hic, cool or mcool)",
                  flush=True)
        except Exception as exc:                          # noqa: BLE001
            print(f"    {key}: {type(exc).__name__}: {exc}", flush=True)
    return None, None


# ---------------------------------------------------------------------------
# extract

def profiles_for_group(reader, chrom, anchors, chrom_len: int, res: int,
                       window: int, anchor_half: int):
    """One query for the whole group, profiles sliced out locally.

    Returns {viewpoint index: vector indexed by offset in bins}, the counts
    divided by the number of anchor bins.
    """
    n = 2 * (window // res) + 1
    g0 = max(0, anchors[0][1] - anchor_half)
    g1 = min(chrom_len, anchors[-1][1] + anchor_half)
    w0 = max(0, anchors[0][1] - window)
    w1 = min(chrom_len, anchors[-1][1] + window)
    if w1 - w0 < 10 * res:
        return {}
    try:
        got = reader.records(chrom, g0, g1, w0, w1)
    except Exception:                                     # noqa: BLE001
        return {}
    if got is None:
        return {}
    bx, by, cc = got
    finite = np.isfinite(cc)
    bx, by, cc = bx[finite], by[finite], cc[finite]

    out = {}
    for idx, tss in anchors:
        lo, hi = tss - anchor_half - res, tss + anchor_half + res
        centre = (tss // res) * res
        acc = np.zeros(n)
        # divide by the number of anchor bins inside the chromosome, not by
        # the number of records seen: sparse formats omit the zeros
        n_abins = n_anchor_bins(tss, res, chrom_len, anchor_half)
        for pa, pb in ((bx, by), (by, bx)):
            m = (pa > lo) & (pa < hi) & ~((pb > lo) & (pb < hi))
            if not m.any():
                continue
            k = ((pb[m] - centre) // res + window // res).astype(int)
            ok = (k >= 0) & (k < n)
            np.add.at(acc, k[ok], cc[m][ok])
        vec = np.full(n, np.nan)
        a0 = max(0, tss - window)
        a1 = min(chrom_len, tss + window)
        klo = int((a0 - centre) // res + window // res)
        khi = int((a1 - centre) // res + window // res)
        inside = np.zeros(n, dtype=bool)
        inside[max(0, klo):min(n, khi + 1)] = True
        # the anchor band itself carries no usable value
        band = np.zeros(n, dtype=bool)
        for b in range(n):
            if lo < centre + (b - window // res) * res < hi:
                band[b] = True
        if n_abins:
            vec[inside & ~band] = acc[inside & ~band] / n_abins
        out[idx] = vec
    return out


def load_anchors(args):
    if args.anchors:
        rows = read_tsv(args.anchors)
        return [{"gene": r["name"], "chrom": r["chrom"], "tss": int(r["pos"]),
                 "set": r.get("set") or "partner"} for r in rows]
    wanted_m = set(Path(args.genes).read_text().split())
    wanted_c = set(Path(args.controls).read_text().split()) if args.controls else set()
    tss = load_tss(Path(args.gtf), wanted_m | wanted_c)
    for g in tss:
        g["set"] = "mito" if g["gene"] in wanted_m else "control"
    return tss


def chrom_in_file(name, chroms):
    for c in (name, name.replace("chr", ""), "chr" + name):
        if c in chroms:
            return c
    return None


def cmd_extract(args):
    rows = read_tsv(args.tissues)
    if args.index is not None:
        if not 1 <= args.index <= len(rows):
            sys.exit(f"--index {args.index} outside 1..{len(rows)}")
        rows = [rows[args.index - 1]]
    anchors = load_anchors(args)
    if not anchors:
        sys.exit("no anchor; check --gtf with --genes, or --anchors")
    half = args.window // args.resolution
    n_bins = 2 * half + 1
    offsets = np.arange(-half, half + 1) * args.resolution
    print(f"{len(anchors)} anchors, {len(rows)} files, raw counts")
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    for row in rows:
        acc = row[args.accession_column]
        dest = Path(args.outdir) / f"{acc}.npz"
        if dest.exists() and not args.force:
            print(f"[skip] {acc}")
            continue
        t0 = time.time()
        reader, via = open_library(row, args)
        if reader is None:
            print(f"[{acc}] no readable location in "
                  f"{args.location_columns}, skipped")
            continue
        if args.resolution not in reader.resolutions:
            print(f"[{acc}] {args.resolution} bp absent "
                  f"(has {sorted(reader.resolutions)[:6]}), skipped")
            continue
        chroms = reader.chroms
        mat = np.full((len(anchors), n_bins), np.nan, dtype=np.float32)
        by_chrom = {}
        for i, g in enumerate(anchors):
            c = chrom_in_file(g["chrom"], chroms)
            if c is not None:
                by_chrom.setdefault(c, []).append((i, g["tss"]))
        for chrom, items in by_chrom.items():
            items.sort(key=lambda x: x[1])
            try:
                reader.prepare(chrom)
            except Exception as exc:                      # noqa: BLE001
                print(f"    {chrom}: {type(exc).__name__}", flush=True)
                continue
            for grp in group_anchors(items, args.group_span):
                prof = profiles_for_group(reader, chrom, grp, chroms[chrom],
                                          args.resolution, args.window,
                                          args.anchor_half)
                for idx, vec in prof.items():
                    k = n_anchor_bins(anchors[idx]["tss"], args.resolution,
                                      chroms[chrom], args.anchor_half)
                    mat[idx] = np.round(vec * k).astype(np.float32)
        raw = np.nansum(mat, axis=1)
        mat[raw < args.min_raw] = np.nan
        np.savez_compressed(
            dest, counts=mat, raw_total=raw,
            names=np.array([g["gene"] for g in anchors]),
            chrom=np.array([g["chrom"] for g in anchors]),
            pos=np.array([g["tss"] for g in anchors], dtype=np.int64),
            set=np.array([g["set"] for g in anchors]),
            offsets=offsets, accession=acc,
            term=row.get(args.term_column, ""),
            gtex=row.get(args.gtex_column, ""),
            mtcn=float(row.get(args.mtcn_column, "nan") or "nan"),
            resolution=args.resolution, anchor_half=args.anchor_half, via=via)
        kept = int(np.isfinite(mat).any(1).sum())
        print(f"[{acc}] {kept}/{len(anchors)} anchors, via {via}, "
              f"{time.time()-t0:.0f}s", flush=True)


# ---------------------------------------------------------------------------
# select

def cmd_select(args):
    info = {r[args.accession_column]: r for r in read_tsv(args.tissues)}
    rows = []
    for f in sorted(Path(args.indir).glob("*.npz")):
        d = np.load(f, allow_pickle=True)
        acc = str(d["accession"])
        sel = d["set"] == args.depth_set
        depth = float(np.nansum(d["raw_total"][sel]))
        r = info.get(acc, {})
        rows.append({"accession": acc,
                     "term": r.get(args.term_column) or str(d["term"]),
                     "gtex": r.get(args.gtex_column) or str(d["gtex"]),
                     "mtcn": float(d["mtcn"]), "depth": depth})
    rows = [r for r in rows if np.isfinite(r["mtcn"]) and r["depth"] > 0]
    if not rows:
        sys.exit("no library with mtCN and contacts")
    lg = np.log10([r["depth"] for r in rows])
    med = float(np.median(lg))
    sd = float(np.std(lg, ddof=1)) if len(lg) > 1 else 0.0
    cut = med - args.k_sd * sd
    keep = lg >= cut
    rule = f"log10 depth >= median - {args.k_sd:g} SD = {cut:.2f}"
    if keep.sum() < args.min_keep:
        order = np.argsort(-lg)
        keep = np.zeros(len(rows), bool)
        keep[order[:args.min_keep]] = True
        rule = (f"median - {args.k_sd:g} SD would keep {int((lg >= cut).sum())}; "
                f"the {args.min_keep} deepest kept instead "
                f"(log10 depth >= {lg[keep].min():.2f})")
    print(f"{len(rows)} libraries, log10 depth median {med:.2f}, SD {sd:.2f}")
    print(f"rule: {rule}")
    print(f"kept {int(keep.sum())}, dropped {int((~keep).sum())}")
    for r, k, l in sorted(zip(rows, keep, lg), key=lambda x: x[2]):
        if not k:
            print(f"  dropped {r['accession']}  {l:.2f}  {r['term']}")
    terms = {}
    for r, k in zip(rows, keep):
        if k:
            terms.setdefault(r["term"], 0)
            terms[r["term"]] += 1
    lost = sorted({r["term"] for r in rows} - set(terms))
    print(f"{len(terms)} terms remain, {sum(v >= 2 for v in terms.values())} "
          f"with replicates; terms lost entirely: "
          f"{', '.join(lost) if lost else 'none'}")
    with open(args.out, "w") as fh:
        fh.write("accession\tterm\tgtex\tmtcn\tdepth\tlog10_depth\tkept\n")
        for r, k, l in zip(rows, keep, lg):
            fh.write(f"{r['accession']}\t{r['term']}\t{r['gtex']}\t"
                     f"{r['mtcn']:.1f}\t{r['depth']:.0f}\t{l:.4f}\t{int(k)}\n")


# ---------------------------------------------------------------------------
# export

def cmd_export(args):
    sets = [s.strip() for s in args.set.split(",") if s.strip()]
    bad = set(sets) - {"mito", "partner", "control"}
    if bad or not sets:
        sys.exit(f"--set takes mito, partner, control or a comma-separated "
                 f"list of them, got {args.set!r}")
    kept = {r["accession"] for r in read_tsv(args.libraries) if r["kept"] == "1"}
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    decay_rows = []
    n = 0
    for f in sorted(Path(args.indir).glob("*.npz")):
        d = np.load(f, allow_pickle=True)
        acc = str(d["accession"])
        if acc not in kept:
            continue
        sel = np.isin(d["set"], sets)
        C = d["counts"][sel]
        offs = d["offsets"]
        with gzip.open(out / f"{acc}.tsv.gz", "wt") as fh:
            fh.write("name\tchrom\tpos\t" + "\t".join(str(int(o)) for o in offs)
                     + "\n")
            for name, ch, p, row in zip(d["names"][sel], d["chrom"][sel],
                                        d["pos"][sel], C):
                fh.write(f"{name}\t{ch}\t{int(p)}\t" + "\t".join(
                    "NA" if not np.isfinite(v) else str(int(v)) for v in row)
                    + "\n")
        # the expected count at a distance is taken on the control viewpoints
        if "mito" in sets or "control" in sets:
            ctrl = d["counts"][d["set"] == "control"]
            decay_rows.append((acc, np.nansum(ctrl, axis=0), offs))
        n += 1
    if decay_rows:
        with open(out / "decay.tsv", "w") as fh:
            fh.write("accession\t" + "\t".join(str(int(o)) for o in
                                               decay_rows[0][2]) + "\n")
            for acc, v, _ in decay_rows:
                fh.write(acc + "\t" + "\t".join(f"{x:.0f}" for x in v) + "\n")
    print(f"{n} libraries exported to {out}/ (set {','.join(sets)})")
    missing = kept - {p.name[:-len(".tsv.gz")] for p in out.glob("*.tsv.gz")}
    if missing:
        print(f"  kept but not found in {args.indir}: {', '.join(sorted(missing))}")


# ---------------------------------------------------------------------------
# relate

def perm_test(stat_fn, group_of_term, mt_group, perms, rng):
    """Two-sided; mtCN moved between GTEx tissues, the terms of one tissue
    moving together."""
    g = np.array(sorted(set(group_of_term)))
    pos = {t: i for i, t in enumerate(g)}
    gi = np.array([pos[t] for t in group_of_term])
    obs = stat_fn(mt_group[gi])
    null = np.array([stat_fn(mt_group[rng.permutation(len(g))][gi])
                     for _ in range(perms)])
    return obs, (1.0 + (np.abs(null) >= abs(obs)).sum()) / (perms + 1.0)


def cmd_relate(args):
    libs = [r for r in read_tsv(args.libraries) if r["kept"] == "1"]
    terms = sorted({r["term"] for r in libs})
    ti = {t: i for i, t in enumerate(terms)}
    T = len(terms)
    mt = np.array([np.median([float(r["mtcn"]) for r in libs if r["term"] == t])
                   for t in terms])
    gtex = [next(r["gtex"] for r in libs if r["term"] == t) for t in terms]
    depth = np.array([sum(float(r["depth"]) for r in libs if r["term"] == t)
                      for t in terms])

    allowed = None
    if args.min_terms > 1:
        rows = read_tsv(Path(args.loops) / "loops.tsv")
        allowed = {r["loop_id"] for r in rows
                   if int(r["n_terms"]) >= args.min_terms}
        print(f"{len(allowed)} of {len(rows)} loops called in at least "
              f"{args.min_terms} terms")

    pres_rows = read_tsv(Path(args.loops) / "presence.tsv")
    P = np.array([[int(r.get(t, "0")) for t in terms] for r in pres_rows]).T
    shared = P.astype(int) @ P.T.astype(int)

    size = np.zeros((T, T)); signed = np.zeros((T, T)); events = np.zeros((T, T))
    types = {k: np.zeros((T, T)) for k in ("gain", "loss", "stronger", "weaker")}
    min_l2 = math.log2(args.min_fold)
    with gzip.open(Path(args.loops) / "deseq_pairs.tsv.gz", "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        c = {h: i for i, h in enumerate(head)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if allowed is not None and f[c["loop_id"]] not in allowed:
                continue
            a, b = ti[f[c["term_a"]]], ti[f[c["term_b"]]]
            lfc = float(f[c["log2fc"]])
            if not np.isfinite(lfc):
                continue
            size[a, b] += abs(lfc); size[b, a] += abs(lfc)
            signed[a, b] += lfc; signed[b, a] -= lfc
            padj = f[c["padj"]]
            if padj not in ("NA", "") and float(padj) < args.fdr \
                    and abs(lfc) >= min_l2:
                events[a, b] += 1; events[b, a] += 1
                ca, cb = f[c["called_a"]] == "1", f[c["called_b"]] == "1"
                hi_is_a = mt[a] >= mt[b]
                up_in_hi = (lfc > 0) == hi_is_a
                called_hi, called_lo = (ca, cb) if hi_is_a else (cb, ca)
                if up_in_hi:
                    kind = "gain" if called_hi and not called_lo else "stronger"
                else:
                    kind = "loss" if called_lo and not called_hi else "weaker"
                types[kind][a, b] += 1; types[kind][b, a] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        norm = np.where(shared > 0, shared, np.nan)
        size /= norm; signed /= norm; events /= norm
        for k in types:
            types[k] /= norm

    iu = np.triu_indices(T, 1)
    ok = np.isfinite(size[iu])
    g_sorted = sorted(set(gtex))
    mt_group = np.array([np.median([mt[t] for t in range(T) if gtex[t] == g])
                         for g in g_sorted])
    power = np.log(2.0 / (1.0 / depth[:, None] + 1.0 / depth[None, :]))

    def sym(D):
        d = D[iu][ok]
        return lambda m: spearman(d, np.abs(m[:, None] - m[None, :])[iu][ok])

    def direction(m):
        a, b = iu[0][ok], iu[1][ok]
        hi = m[a] > m[b]
        keep = m[a] != m[b]
        net = np.where(hi, signed[a, b], signed[b, a])[keep]
        return spearman(net, np.abs(m[a] - m[b])[keep])

    res = []
    for name, fn, D in (("size", sym(size), size), ("events", sym(events), events),
                        ("direction", direction, None)):
        # the same permutations for every statistic
        r, p = perm_test(fn, gtex, mt_group, args.perms,
                         np.random.default_rng(args.seed))
        if D is not None:
            M = np.abs(mt[:, None] - mt[None, :])
            rp = partial_spearman(D[iu][ok], M[iu][ok], power[iu][ok])
            rw = spearman(D[iu][ok], power[iu][ok])
        else:
            rp = rw = np.nan
        res.append((name, r, p, rp, rw))

    print(f"{T} terms, {int(ok.sum())} pairs, median {np.median(shared[iu]):.0f} "
          f"shared viewpoints, {len(g_sorted)} mtCN values")
    print(f"size per shared viewpoint: median {np.nanmedian(size[iu]):.4f}")
    print(f"significant events per shared viewpoint: median "
          f"{np.nanmedian(events[iu]):.4f}")
    print(f"\n  {'':<11}{'rho':>7}{'p':>9}{'power-partial':>15}{'vs power':>10}")
    for name, r, p, rp, rw in res:
        print(f"  {name:<11}{r:>+7.3f}{p:>9.4f}{rp:>+15.3f}{rw:>+10.3f}")
    print("  size: sum of |shrunken log2 fold change| over the loops, per shared "
          "viewpoint, against |delta mtCN|. events: significant changes "
          f"(FDR {args.fdr}, fold >= {args.min_fold}). direction: signed sum "
          "oriented towards the higher-mtCN term, against delta mtCN.")
    M = np.abs(mt[:, None] - mt[None, :])
    print("  types against |delta mtCN| (descriptive): " + ", ".join(
        f"{k} {spearman(types[k][iu][ok], M[iu][ok]):+.2f}" for k in types))

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    with open(out / "relate_stats.tsv", "w") as fh:
        fh.write("statistic\trho\tp_perm\trho_power_partial\trho_vs_power\n")
        for name, r, p, rp, rw in res:
            fh.write(f"{name}\t{r:.4f}\t{p:.4f}\t{rp:.4f}\t{rw:.4f}\n")
    with open(out / "relate_pairs.tsv", "w") as fh:
        fh.write("term_a\tterm_b\tmtcn_a\tmtcn_b\tshared_viewpoints\t"
                 "size_per_vp\tsigned_a_over_b\tevents_per_vp\t"
                 + "\t".join(f"{k}_per_vp" for k in types) + "\n")
        for a, b in zip(*iu):
            fh.write(f"{terms[a]}\t{terms[b]}\t{mt[a]:.0f}\t{mt[b]:.0f}\t"
                     f"{shared[a, b]}\t{size[a, b]:.5f}\t{signed[a, b]:.5f}\t"
                     f"{events[a, b]:.5f}\t"
                     + "\t".join(f"{types[k][a, b]:.5f}" for k in types) + "\n")
    draw(out / "relate_size.svg", M[iu][ok], size[iu][ok], res[0], T)
    print(f"\nWritten to {out}/")


def draw(path, dmt, d, res, T):
    W, H, pad = 720, 520, 80
    pw, ph = W - 2 * pad, H - 2 * pad
    xr = float(dmt.max()) or 1.0
    lo, hi = float(d.min()), float(d.max())
    my = 0.06 * (hi - lo + 1e-12)
    sx = lambda v: pad + v / xr * pw
    sy = lambda v: pad + ph - (v - lo + my) / (hi - lo + 2 * my) * ph
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}" font-family="Helvetica,Arial" font-size="11">',
         f'<rect width="{W}" height="{H}" fill="white"/>',
         f'<clipPath id="plot"><rect x="{pad}" y="{pad}" width="{pw}" '
         f'height="{ph}"/></clipPath>',
         f'<rect x="{pad}" y="{pad}" width="{pw}" height="{ph}" fill="none" '
         f'stroke="#999"/>']
    for i in range(len(d)):
        s.append(f'<circle cx="{sx(dmt[i]):.1f}" cy="{sy(d[i]):.1f}" r="2.6" '
                 f'fill="#c44" fill-opacity="0.45"/>')
    b, a = theil_sen(dmt, d)
    s.append(f'<line x1="{sx(0):.1f}" y1="{sy(a):.1f}" x2="{sx(xr):.1f}" '
             f'y2="{sy(a + b * xr):.1f}" stroke="#333" stroke-dasharray="6 4" '
             f'stroke-width="1.5" clip-path="url(#plot)"/>')
    for t in range(5):
        v = xr * t / 4
        s.append(f'<text x="{sx(v):.1f}" y="{pad+ph+16}" '
                 f'text-anchor="middle">{v:,.0f}</text>')
        v = lo + (hi - lo) * t / 4
        s.append(f'<text x="{pad-6}" y="{sy(v)+4:.1f}" '
                 f'text-anchor="end">{v:.3f}</text>')
    _, r, p, _, _ = res
    ptxt = f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.1e}"
    s.append(f'<text x="{pad}" y="{pad-26}" font-size="12">rho {r:+.2f}, '
             f'{ptxt}, {T} terms, one point per pair</text>')
    s.append(f'<text x="{pad}" y="{pad-10}" fill="#666">summed |shrunken log2 '
             f'fold change| of the loops, per shared viewpoint</text>')
    s.append(f'<text x="{pad+pw/2:.0f}" y="{H-24}" text-anchor="middle">'
             f'difference in mtDNA copy number between the two terms</text>')
    s.append(f'<text x="24" y="{pad+ph/2:.0f}" text-anchor="middle" '
             f'transform="rotate(-90 24 {pad+ph/2:.0f})">change in looping '
             f'between the two terms</text>')
    s.append("</svg>")
    path.write_text("\n".join(s))


# ---------------------------------------------------------------------------

def add_columns(p, which):
    if "accession" in which:
        p.add_argument("--accession-column", default="accession",
                       help="accession column of --tissues")
    if "term" in which:
        p.add_argument("--term-column", default="term_name",
                       help="ENCODE term column of --tissues")
    if "gtex" in which:
        p.add_argument("--gtex-column", default="gtex_tissue",
                       help="GTEx tissue column of --tissues")
    if "mtcn" in which:
        p.add_argument("--mtcn-column", default="mtcn",
                       help="mtDNA copy number column of --tissues")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    D = argparse.ArgumentDefaultsHelpFormatter

    p = sub.add_parser("extract", formatter_class=D,
                       help="virtual-4C profiles, one npz per library")
    p.add_argument("--tissues", required=True, help="library table (TSV)")
    p.add_argument("--gtf", help="GENCODE GTF, gzipped or plain")
    p.add_argument("--genes", help="mitochondrial gene symbols, one per line")
    p.add_argument("--controls", help="control gene symbols, one per line")
    p.add_argument("--anchors",
                   help="TSV name/chrom/pos[/set], instead of the GTF")
    p.add_argument("--outdir", required=True)
    p.add_argument("--index", type=int,
                   help="1-based row of --tissues to run alone (job arrays)")
    p.add_argument("--resolution", type=int, default=RESOLUTION,
                   help="bin size in bp")
    p.add_argument("--window", type=int, default=WINDOW,
                   help="profile half-width in bp")
    p.add_argument("--anchor-half", type=int, default=ANCHOR_HALF,
                   help="anchor half-width around the TSS in bp")
    p.add_argument("--group-span", type=int, default=GROUP_SPAN,
                   help="viewpoints within this span (bp) share one query")
    p.add_argument("--min-raw", type=float, default=MIN_RAW,
                   help="profiles with fewer raw contacts are set to NaN")
    p.add_argument("--location-columns", default=LOCATION_COLUMNS,
                   help="columns of --tissues holding a URL or a local path, "
                        "tried in this order")
    p.add_argument("--format-column", default="format",
                   help="optional column of --tissues giving hic, cool or "
                        "mcool; otherwise the extension decides")
    p.add_argument("--retries", type=int, default=4,
                   help="attempts to open a remote .hic before the next "
                        "location column")
    add_columns(p, ("accession", "term", "gtex", "mtcn"))
    p.add_argument("--force", action="store_true",
                   help="overwrite existing npz files")

    p = sub.add_parser("select", formatter_class=D,
                       help="library depth and depth filter")
    p.add_argument("--tissues", required=True, help="library table (TSV)")
    p.add_argument("--indir", required=True, help="npz files of extract")
    p.add_argument("--out", default="libraries.tsv")
    p.add_argument("--k-sd", type=float, default=1.0,
                   help="drop libraries below median - k SD of log10 depth")
    p.add_argument("--min-keep", type=int, default=55,
                   help="keep at least this many libraries (the deepest); "
                        "with fewer libraries than this, all are kept, so "
                        "set it to suit the size of the library table")
    p.add_argument("--depth-set", default="control",
                   choices=["control", "mito", "partner"],
                   help="viewpoint set whose contacts define library depth")
    add_columns(p, ("accession", "term", "gtex"))

    p = sub.add_parser("export", formatter_class=D,
                       help="one table per kept library for R")
    p.add_argument("--indir", required=True, help="npz files of extract")
    p.add_argument("--libraries", required=True, help="libraries.tsv of select")
    p.add_argument("--set", required=True,
                   help="mito, partner, control, or a comma-separated list")
    p.add_argument("--outdir", required=True)

    p = sub.add_parser("relate", formatter_class=D,
                       help="pairwise route from the DESeq2 results")
    p.add_argument("--libraries", required=True, help="libraries.tsv of select")
    p.add_argument("--loops", required=True,
                   help="output folder of v4c_loops.R (loops.tsv, "
                        "presence.tsv, deseq_pairs.tsv.gz)")
    p.add_argument("--outdir", default="relate")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--min-fold", type=float, default=1.5,
                   help="minimum fold change of a significant event")
    p.add_argument("--min-terms", type=int, default=1,
                   help="keep only the loops called in at least this many "
                        "terms; 1 keeps every loop")
    p.add_argument("--perms", type=int, default=5000)
    p.add_argument("--seed", type=int, default=0)

    args = ap.parse_args()
    if args.cmd == "extract" and not args.anchors and not (args.gtf and args.genes):
        sys.exit("extract needs --anchors, or --gtf with --genes")
    {"extract": cmd_extract, "select": cmd_select, "export": cmd_export,
     "relate": cmd_relate}[args.cmd](args)


if __name__ == "__main__":
    main()
