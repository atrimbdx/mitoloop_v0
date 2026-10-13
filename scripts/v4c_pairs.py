#!/usr/bin/env python3
"""
v4c_pairs.py - virtual-4C profiles at mitochondrial promoters in paired Hi-C
libraries, two metabolic states of one cell type per experiment

Purpose
-------
Each experiment of the manifest holds two states of the same cell type that
differ in mitochondrial content:

    mito+   the state with the higher mitochondrial content
            (myotube, cardiomyocyte, activated lymphocyte, adipocyte)
    mito-   the other state of the same experiment
            (myoblast, iPSC, resting lymphocyte, preadipocyte)

The question is whether the loops called at mitochondrial promoters are the
same in the two states, and whether they differ more than the loops of control
promoters do. It is the within-experiment counterpart of the across-tissue
analysis of v4c_pipeline.py: same cells, same laboratory, same protocol. This
script extracts and exports the profiles; the comparison is v4c_pairs.R.

Subcommands
-----------
  extract  One virtual-4C profile per viewpoint and per library: observed
           counts, no normalisation, summed over the anchor bins
           (TSS +/- --anchor-half) against every bin within +/- --window.
           Reads .hic (local file or URL, through hicstraw) and .cool /
           .mcool (through cooler), so that datasets distributed in either
           format go through the same code. One npz per library; --index
           runs one library, for array jobs.

  export   One table per library for R, a viewpoint per row and a bin per
           column, plus libraries_pairs.tsv with the depth of each library.

Inputs
------
  --pairs     manifest, tab separated, one line per library; lines starting
              with '#' are ignored. Columns:
                experiment  name of the experiment (pairs the two states)
                state       --plus-state or --minus-state (default plus, minus)
                label       unique; becomes the npz and tsv.gz file name
                path        local file, or an http(s) URL for .hic
                format      hic | cool | mcool  (cool/mcool need cooler)
                genome      key of a --gtf entry; may be empty if --genome
                            is given
                resolution  bin size in bp, must exist in the file; empty
                            means --resolution
              Example: config/pairs.example.tsv
                experiment  state  label        path                         format  genome  resolution
                myogenesis  minus  myoblast_r1  data/pairs/myoblast_r1.mcool mcool   hg19    10000
                myogenesis  plus   myotube_r1   data/pairs/myotube_r1.mcool  mcool   hg19    10000
  --gtf       GENOME=PATH, repeated once per genome build; GTF (optionally
              gzipped) with gene_name attributes. The MANE Select transcript
              start is the TSS where tagged, otherwise the gene start.
  --genes     mitochondrial gene symbols, one per line
  --controls  control gene symbols, one per line (optional); a symbol in both
              lists is kept as mito

Gene symbols are matched without case, so human symbols find their mouse
counterparts (NDUFA1 -> Ndufa1). Orthology is not checked; a mouse dataset is
a supporting observation, not a replicate of a human one.

Outputs
-------
  extract  <outdir>/<label>.npz with keys
             counts       viewpoints x offsets, raw counts (NaN outside the
                          chromosome, in the anchor band, or for viewpoints
                          below --min-raw contacts)
             raw_total    contacts per viewpoint before the --min-raw filter
             names, symbol, chrom, pos, set ('mito' or 'control')
             offsets      bin offsets in bp from the TSS bin
             label, experiment, state, genome, resolution, anchor_half
  export   <outdir>/<label>.tsv.gz with columns name, set, chrom, pos, then
           one column per offset in bp (NA where missing);
           <outdir>/libraries_pairs.tsv with columns label, experiment,
           state, genome, resolution, depth_control, depth_mito,
           n_viewpoints (depths are contacts summed over the viewpoints of
           each set).

Method
------
For each viewpoint, contacts between the anchor bins (TSS +/- --anchor-half,
plus one bin) and every bin of the window are summed; contacts within the
anchor band are excluded. Neighbouring viewpoints (within --group-span) are
read in one query and sliced locally. Raw counts are kept because peakC
normalises each profile itself and the binomial thinning and splitting in
v4c_pairs.R need integer counts. A .cool file stores the upper triangle only,
so the whole window is read against itself and the anchor is selected
afterwards; hicstraw returns both orientations of an intra-chromosomal query.

Usage
-----
    python3 scripts/v4c_pairs.py extract --pairs config/pairs.tsv \\
        --gtf hg38=data/gencode.v44.gtf.gz --gtf hg19=data/gencode.v19.gtf.gz \\
        --gtf mm10=data/gencode.vM25.gtf.gz \\
        --genes genes.txt --controls controls.txt --outdir results/raw_pairs
    python3 scripts/v4c_pairs.py export --pairs config/pairs.tsv \\
        --indir results/raw_pairs --outdir results/export_pairs
    Rscript scripts/v4c_pairs.R compare \\
        --libraries results/export_pairs/libraries_pairs.tsv \\
        --export results/export_pairs --wsize 3 --out results/pairs_out

References
----------
  Durand et al. 2016. Juicer provides a one-click system for analyzing
    loop-resolution Hi-C experiments. Cell Systems. (hicstraw, .hic format)
  Abdennur and Mirny 2020. Cooler: scalable storage for Hi-C data and other
    genomically labeled arrays. Bioinformatics. (.cool / .mcool format)
  Geeven et al. 2018. peakC: a flexible, non-parametric peak calling package
    for 4C and Capture-C data. Nucleic Acids Research. (downstream caller)
  Frankish et al. 2023. GENCODE: reference annotation for the human and
    mouse genomes in 2023. Nucleic Acids Research. doi:10.1093/nar/gkac1071
    (GTF annotation)
  Morales et al. 2022. A joint NCBI and EMBL-EBI transcript set for clinical
    genomics and research. Nature. (MANE Select)
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
GROUP_SPAN = 400_000
MIN_RAW = 100
RETRIES = 4
PLUS_STATE = "plus"
MINUS_STATE = "minus"


# ---------------------------------------------------------------------------
# small helpers

def read_tsv(path):
    """Rows of a tab-separated table as dicts; '#' lines are comments."""
    rows = []
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            f += [""] * (len(head) - len(f))
            rows.append(dict(zip(head, f)))
    return rows


def load_tss(gtf: Path, wanted_upper: set):
    """MANE Select TSS where available, otherwise the gene record.

    Symbols are compared in upper case, so human symbols find mouse genes.
    """
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
            key = name.upper()
            if key not in wanted_upper:
                continue
            mane = "MANE_Select" in a
            rank = 2 if (f[2] == "transcript" and mane) else (1 if f[2] == "gene" else 0)
            if rank == 0:
                continue
            tss = int(f[3]) - 1 if f[6] == "+" else int(f[4])
            prev = best.get(key)
            if prev is None or rank > prev["rank"]:
                best[key] = {"gene": key, "symbol": name, "chrom": f[0],
                             "tss": tss, "strand": f[6], "rank": rank}
    return sorted(best.values(), key=lambda g: (g["chrom"], g["tss"]))


def group_anchors(items, span: int):
    """Split (index, tss) items sorted by tss into groups read in one query."""
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
    """Number of anchor bins inside the chromosome, as profiles_for_group."""
    lo, hi = pos - anchor_half - res, pos + anchor_half + res
    return sum(1 for b in range(int(math.ceil(max(0, lo) / res)),
                                int(min(chrom_len, hi) // res) + 1)
               if lo < b * res < hi)


# ---------------------------------------------------------------------------
# two readers behind one interface
#
# fetch(chrom, g0, g1, w0, w1) -> (x, y, counts), x and y being the genomic
# start of the two bins of a contact, in base pairs.

class HicReader:
    kind = "hic"

    def __init__(self, path, res, retries=RETRIES):
        import hicstraw
        # remote .hic files are opened over http; transient failures are retried
        last = None
        for attempt in range(retries):
            try:
                self.h = hicstraw.HiCFile(str(path))
                last = None
                break
            except Exception as exc:                      # noqa: BLE001
                last = exc
                print(f"    open attempt {attempt+1}: "
                      f"{type(exc).__name__}: {exc}", flush=True)
                time.sleep(5 * 2 ** attempt)
        if last is not None:
            raise last
        self.res = res
        self.resolutions = set(self.h.getResolutions())
        self.chroms = {c.name: c.length for c in self.h.getChromosomes()}
        self._zoom = {}

    def _mzd(self, chrom):
        if chrom not in self._zoom:
            self._zoom[chrom] = self.h.getMatrixZoomData(
                chrom, chrom, "observed", "NONE", "BP", self.res)
        return self._zoom[chrom]

    def fetch(self, chrom, g0, g1, w0, w1):
        recs = self._mzd(chrom).getRecords(int(g0), int(g1), int(w0), int(w1))
        if not recs:
            return None
        x = np.array([r.binX for r in recs], dtype=np.float64)
        y = np.array([r.binY for r in recs], dtype=np.float64)
        c = np.array([r.counts for r in recs], dtype=np.float64)
        return x, y, c


class CoolReader:
    kind = "cool"

    def __init__(self, path, res, is_mcool):
        import cooler
        uri = f"{path}::/resolutions/{int(res)}" if is_mcool else str(path)
        self.c = cooler.Cooler(uri)
        self.res = int(self.c.binsize)
        self.resolutions = {self.res}
        self.chroms = dict(self.c.chromsizes.astype(int))
        self._sel = self.c.matrix(balance=False, as_pixels=True, join=True)

    def fetch(self, chrom, g0, g1, w0, w1):
        # The whole window against itself, not the anchor against the window:
        # a cooler stores the upper triangle only, so a rectangular anchor x
        # window query returns the contacts on one side of the anchor and
        # drops those on the other. The caller selects the anchor itself.
        # Cooler regions are half-open; one more base includes the bin that
        # starts at w1, the last offset of the profile, as hicstraw does.
        end = min(int(w1) + 1, int(self.chroms[chrom]))
        r2 = f"{chrom}:{int(w0)}-{end}"
        try:
            px = self._sel.fetch(r2, r2)
        except Exception:                                 # noqa: BLE001
            return None
        if px is None or len(px) == 0:
            return None
        x = px["start1"].to_numpy(dtype=np.float64)
        y = px["start2"].to_numpy(dtype=np.float64)
        c = px["count"].to_numpy(dtype=np.float64)
        return x, y, c


def open_matrix(row, res, retries=RETRIES):
    fmt = (row.get("format") or "hic").strip().lower()
    path = row["path"].strip()
    if fmt == "hic":
        return HicReader(path, res, retries=retries)
    if fmt in ("cool", "mcool"):
        return CoolReader(path, res, is_mcool=(fmt == "mcool"))
    raise ValueError(f"unknown format {fmt!r} (hic, cool or mcool)")


# ---------------------------------------------------------------------------
# profiles

def profiles_for_group(reader, chrom, anchors, chrom_len, res, window,
                       anchor_half):
    """One fetch for a group of neighbouring viewpoints, sliced locally.

    Returns {viewpoint index: vector indexed by offset in bins}. Counts are
    divided by the number of anchor bins inside the chromosome, not by the
    number of records seen: a sparse format omits the zeros, and dividing by
    what was seen would inflate every value.
    """
    n = 2 * (window // res) + 1
    g0 = max(0, anchors[0][1] - anchor_half)
    g1 = min(chrom_len, anchors[-1][1] + anchor_half)
    w0 = max(0, anchors[0][1] - window)
    w1 = min(chrom_len, anchors[-1][1] + window)
    if w1 - w0 < 10 * res:
        return {}
    got = reader.fetch(chrom, g0, g1, w0, w1)
    if got is None:
        return {}
    bx, by, cc = got
    finite = np.isfinite(cc)
    bx, by, cc = bx[finite], by[finite], cc[finite]
    if bx.size == 0:
        return {}

    out = {}
    for idx, tss in anchors:
        lo, hi = tss - anchor_half - res, tss + anchor_half + res
        centre = (tss // res) * res
        acc = np.zeros(n)
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
        band = np.zeros(n, dtype=bool)
        for b in range(n):
            if lo < centre + (b - window // res) * res < hi:
                band[b] = True
        if n_abins:
            vec[inside & ~band] = acc[inside & ~band] / n_abins
        out[idx] = vec
    return out


# ---------------------------------------------------------------------------
# extract

def check_manifest(rows, plus=PLUS_STATE, minus=MINUS_STATE):
    bad = [r for r in rows if r.get("state") not in (plus, minus)]
    if bad:
        sys.exit(f"state must be {plus!r} or {minus!r}; offending labels: "
                 + ", ".join(r.get("label", "?") for r in bad[:5]))
    seen = {}
    for r in rows:
        seen.setdefault(r["experiment"], set()).add(r["state"])
    lone = [e for e, s in seen.items() if s != {plus, minus}]
    if lone:
        sys.exit("these experiments do not have both states: "
                 + ", ".join(sorted(lone)))
    labels = [r["label"] for r in rows]
    if len(set(labels)) != len(labels):
        sys.exit("labels must be unique across the manifest")


def cmd_extract(args):
    rows = read_tsv(args.pairs)
    check_manifest(rows, args.plus_state, args.minus_state)
    if args.index is not None:
        if not 1 <= args.index <= len(rows):
            sys.exit(f"--index {args.index} outside 1..{len(rows)}")
        print(f"array task: library {args.index} of {len(rows)}")
        rows = [rows[args.index - 1]]

    gtfs = {}
    for spec in args.gtf:
        if "=" not in spec:
            sys.exit(f"--gtf wants genome=path, got {spec!r}")
        g, p = spec.split("=", 1)
        gtfs[g.strip()] = Path(p.strip())

    wanted_m = {g.upper() for g in Path(args.genes).read_text().split()}
    wanted_c = ({g.upper() for g in Path(args.controls).read_text().split()}
                if args.controls else set())
    if wanted_m & wanted_c:
        print(f"note: {len(wanted_m & wanted_c)} symbols are in both lists, "
              f"kept as mito")
        wanted_c -= wanted_m

    anchors_by_genome = {}

    def anchors_for(genome):
        if genome in anchors_by_genome:
            return anchors_by_genome[genome]
        if genome not in gtfs:
            sys.exit(f"no --gtf given for genome {genome!r}")
        tss = load_tss(gtfs[genome], wanted_m | wanted_c)
        for g in tss:
            g["set"] = "mito" if g["gene"] in wanted_m else "control"
        n_m = sum(1 for g in tss if g["set"] == "mito")
        print(f"[{genome}] {len(tss)} viewpoints "
              f"({n_m} mito, {len(tss)-n_m} control) from {gtfs[genome].name}")
        if n_m == 0:
            sys.exit(f"no mito gene found in {gtfs[genome]}; check the build")
        anchors_by_genome[genome] = tss
        return tss

    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    for row in rows:
        label = row["label"]
        dest = Path(args.outdir) / f"{label}.npz"
        if dest.exists() and not args.force:
            print(f"[skip] {label}")
            continue
        res = int(row.get("resolution") or args.resolution)
        if res != args.resolution:
            print(f"[{label}] note: the manifest gives {res} bp, "
                  f"--resolution is {args.resolution} bp; {res} bp is used")
        genome = (row.get("genome") or "").strip() or args.genome
        if not genome:
            sys.exit(f"[{label}] no genome in the manifest and no --genome")
        anchors = anchors_for(genome)
        # the matrix width follows the resolution of this library, so that a
        # manifest resolution different from --resolution still fits
        n_bins = 2 * (args.window // res) + 1
        offsets = np.arange(-(args.window // res),
                            args.window // res + 1) * res
        t0 = time.time()
        try:
            reader = open_matrix(row, res, retries=args.retries)
        except Exception as exc:                          # noqa: BLE001
            print(f"[{label}] cannot open: {type(exc).__name__}: {exc}")
            continue
        if res not in reader.resolutions:
            print(f"[{label}] {res} bp absent from the file "
                  f"(has {sorted(reader.resolutions)[:6]}...), skipped")
            continue

        mat = np.full((len(anchors), n_bins), np.nan, dtype=np.float32)
        by_chrom = {}
        for i, g in enumerate(anchors):
            c = g["chrom"]
            if c not in reader.chroms:
                c = c.replace("chr", "")
            if c not in reader.chroms:
                c = "chr" + g["chrom"]
            if c in reader.chroms:
                by_chrom.setdefault(c, []).append((i, g["tss"]))
        if not by_chrom:
            print(f"[{label}] no chromosome name in common, skipped "
                  f"(file has {list(reader.chroms)[:4]})")
            continue
        for chrom, items in by_chrom.items():
            items.sort(key=lambda x: x[1])
            for grp in group_anchors(items, args.group_span):
                try:
                    prof = profiles_for_group(reader, chrom, grp,
                                              reader.chroms[chrom], res,
                                              args.window, args.anchor_half)
                except Exception as exc:                  # noqa: BLE001
                    print(f"    {chrom}: {type(exc).__name__}", flush=True)
                    continue
                for idx, vec in prof.items():
                    k = n_anchor_bins(anchors[idx]["tss"], res,
                                      reader.chroms[chrom], args.anchor_half)
                    mat[idx] = np.round(vec * k).astype(np.float32)
        raw = np.nansum(mat, axis=1)
        mat[raw < args.min_raw] = np.nan
        np.savez_compressed(
            dest, counts=mat, raw_total=raw,
            names=np.array([g["gene"] for g in anchors]),
            symbol=np.array([g["symbol"] for g in anchors]),
            chrom=np.array([g["chrom"] for g in anchors]),
            pos=np.array([g["tss"] for g in anchors], dtype=np.int64),
            set=np.array([g["set"] for g in anchors]),
            offsets=offsets, label=label,
            experiment=row["experiment"], state=row["state"],
            genome=genome, resolution=res, anchor_half=args.anchor_half)
        kept = int(np.isfinite(mat).any(1).sum())
        print(f"[{label}] {kept}/{len(anchors)} viewpoints, "
              f"{time.time()-t0:.0f}s", flush=True)


# ---------------------------------------------------------------------------
# export

def cmd_export(args):
    man = {r["label"]: r for r in read_tsv(args.pairs)}
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    lib_rows = []
    for f in sorted(Path(args.indir).glob("*.npz")):
        d = np.load(f, allow_pickle=True)
        label = str(d["label"])
        if label not in man:
            print(f"  {label} not in the manifest, skipped")
            continue
        counts = d["counts"]
        sets = d["set"]
        ctrl_depth = float(np.nansum(d["raw_total"][sets == "control"]))
        mito_depth = float(np.nansum(d["raw_total"][sets == "mito"]))
        offs = d["offsets"]
        with gzip.open(out / f"{label}.tsv.gz", "wt") as fh:
            fh.write("name\tset\tchrom\tpos\t"
                     + "\t".join(str(int(o)) for o in offs) + "\n")
            for name, st, ch, p, row in zip(d["names"], sets, d["chrom"],
                                            d["pos"], counts):
                fh.write(f"{name}\t{st}\t{ch}\t{int(p)}\t" + "\t".join(
                    "NA" if not np.isfinite(v) else str(int(v)) for v in row)
                    + "\n")
        lib_rows.append({
            "label": label, "experiment": str(d["experiment"]),
            "state": str(d["state"]), "genome": str(d["genome"]),
            "resolution": int(d["resolution"]),
            "depth_control": ctrl_depth, "depth_mito": mito_depth,
            "n_viewpoints": int(np.isfinite(counts).any(1).sum())})
    if not lib_rows:
        sys.exit(f"nothing exported; check {args.indir}")
    with open(out / "libraries_pairs.tsv", "w") as fh:
        cols = ["label", "experiment", "state", "genome", "resolution",
                "depth_control", "depth_mito", "n_viewpoints"]
        fh.write("\t".join(cols) + "\n")
        for r in sorted(lib_rows, key=lambda x: (x["experiment"], x["state"],
                                                 x["label"])):
            fh.write("\t".join(f"{r[c]:.0f}" if isinstance(r[c], float)
                               else str(r[c]) for c in cols) + "\n")
    print(f"{len(lib_rows)} libraries exported to {out}/")
    by_exp = {}
    for r in lib_rows:
        by_exp.setdefault(r["experiment"], {}).setdefault(r["state"], 0)
        by_exp[r["experiment"]][r["state"]] += 1
    for e, s in sorted(by_exp.items()):
        ok = args.plus_state in s and args.minus_state in s
        print(f"  {e}: " + ", ".join(f"{k} x{v}" for k, v in sorted(s.items()))
              + ("" if ok else "   <- incomplete, will be skipped"))


# ---------------------------------------------------------------------------

def add_state_options(parser):
    parser.add_argument("--plus-state", default=PLUS_STATE,
                        help="value of the manifest 'state' column for the "
                             "state with the higher mitochondrial content "
                             "(default %(default)s)")
    parser.add_argument("--minus-state", default=MINUS_STATE,
                        help="value of the manifest 'state' column for the "
                             "other state (default %(default)s)")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("extract", help="virtual 4C profiles, one npz per library")
    e.add_argument("--pairs", required=True, help="library manifest (TSV)")
    e.add_argument("--gtf", action="append", default=[],
                   metavar="GENOME=PATH",
                   help="GTF of one genome build; repeat once per build")
    e.add_argument("--genome", default=None,
                   help="genome key used for manifest rows whose 'genome' "
                        "is empty (default: none, the column is required)")
    e.add_argument("--genes", required=True,
                   help="mitochondrial gene symbols, one per line")
    e.add_argument("--controls", default=None,
                   help="control gene symbols, one per line")
    e.add_argument("--outdir", required=True, help="directory for the npz files")
    e.add_argument("--index", type=int,
                   help="1-based library (data row, comments excluded) of the "
                        "manifest to run alone, for array jobs")
    e.add_argument("--window", type=int, default=WINDOW,
                   help="half-width of the profile around the TSS, bp "
                        "(default %(default)s)")
    e.add_argument("--resolution", type=int, default=RESOLUTION,
                   help="bin size in bp for rows with no resolution in the "
                        "manifest (default %(default)s)")
    e.add_argument("--anchor-half", type=int, default=ANCHOR_HALF,
                   help="half-width of the viewpoint anchor around the TSS, "
                        "bp (default %(default)s)")
    e.add_argument("--group-span", type=int, default=GROUP_SPAN,
                   help="viewpoints closer than this are read in one query, "
                        "bp (default %(default)s)")
    e.add_argument("--min-raw", type=float, default=MIN_RAW,
                   help="viewpoints with fewer contacts in the window are "
                        "set to NA (default %(default)s)")
    e.add_argument("--retries", type=int, default=RETRIES,
                   help="attempts to open a .hic file (default %(default)s)")
    add_state_options(e)
    e.add_argument("--force", action="store_true",
                   help="overwrite existing npz files")
    e.set_defaults(func=cmd_extract)

    x = sub.add_parser("export", help="one table per library for R")
    x.add_argument("--pairs", required=True, help="library manifest (TSV)")
    x.add_argument("--indir", required=True, help="npz files of `extract`")
    x.add_argument("--outdir", required=True, help="directory for the tables")
    add_state_options(x)
    x.set_defaults(func=cmd_export)

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
