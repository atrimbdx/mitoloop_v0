#!/usr/bin/env python3
"""
make_synthetic_v4c.py - synthetic Hi-C libraries with planted promoter loops
for the test of v4c_pipeline.py and v4c_loops.R.

Purpose
-------
Build a small library set that goes through the whole loop pipeline in a few
minutes: Hi-C libraries in ENCODE-like terms grouped into GTEx-like tissues of
known mtDNA copy number, mitochondrial and control viewpoints on two
chromosomes, one planted loop at every mitochondrial viewpoint, and for half
of them ("trend" genes) a loop strength that rises with mtDNA copy number.

Inputs
------
None; every value is an option.

Outputs (in --outdir)
---------------------
  libs/<accession>.cool     one cooler per library (raw counts)
  libs/<accession>.pairs    with --hic-pairs: the contacts of that library in
                            the Juicer short format with score, for
                            `juicer_tools pre`
  libs/<accession>.mcool    with --mcool: that library as a multi-resolution
                            cooler
  chrom.sizes               chromosome lengths
  gencode.gtf.gz            gene and MANE Select transcript records
  mito.txt, controls.txt    viewpoint gene symbols, one per line
  tissues.tsv               accession, term_name, gtex_tissue, mtcn, s3_url,
                            download_url, path (relative to --outdir)
  planted.tsv               gene, set, chrom, tss, partner (bin centre of the
                            planted partner), fold_min, fold_max
  sets.tsv                  source, set, gene (trend and flat genes)

Method
------
Expected contacts between bins i and j of one chromosome follow
--amplitude * depth / (1 + |i - j|), with a per-library depth factor, gamma
noise per pixel (shape --shape) and Poisson sampling, up to --band bins from
the diagonal. A loop multiplies a 3 x 3 dot centred on (TSS bin, partner bin):
by the fold at its centre, by 1 + (fold - 1) / 2 at the four side bins and by
1 + (fold - 1) / 4 at the corners. The fold is --fold for flat genes and
--fold * (mtcn / --mtcn-ref) ^ --beta for trend genes.

Usage
-----
    python3 make_synthetic_v4c.py --outdir synth --hic-pairs SYN03 \\
        --mcool SYN05

References
----------
  Abdennur and Mirny 2020. Cooler: scalable storage for Hi-C data and other
      genomically labeled arrays. Bioinformatics.
  Durand et al. 2016. Juicer provides a one-click system for analyzing
      loop-resolution Hi-C experiments. Cell Systems.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--res", type=int, default=10_000)
    ap.add_argument("--chrom-len", type=int, default=12_800_000)
    ap.add_argument("--n-chroms", type=int, default=2)
    ap.add_argument("--first-vp", type=int, default=1_100_000)
    ap.add_argument("--vp-spacing", type=int, default=800_000)
    ap.add_argument("--vp-per-chrom", type=int, default=14)
    ap.add_argument("--band", type=int, default=110,
                    help="bins from the diagonal that carry contacts")
    ap.add_argument("--amplitude", type=float, default=120.0)
    ap.add_argument("--shape", type=float, default=30.0,
                    help="gamma shape of the per-pixel noise")
    ap.add_argument("--fold", type=float, default=6.0)
    ap.add_argument("--beta", type=float, default=0.4)
    ap.add_argument("--mtcn-ref", type=float, default=2000.0)
    ap.add_argument("--libs-per-term", type=int, default=2)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--hic-pairs", default="",
                    help="accession also written as Juicer pairs")
    ap.add_argument("--mcool", default="",
                    help="accession also written as .mcool")
    a = ap.parse_args()

    import cooler

    rng = np.random.default_rng(a.seed)
    out = Path(a.outdir)
    (out / "libs").mkdir(parents=True, exist_ok=True)
    res = a.res
    nbin = a.chrom_len // res
    chroms = [f"chr{i + 1}" for i in range(a.n_chroms)]

    # GTEx-like tissues with their mtDNA copy number, ENCODE-like terms
    tissues = [("tissue_A", 500.0), ("tissue_B", 1000.0), ("tissue_C", 2000.0),
               ("tissue_D", 4000.0), ("tissue_E", 8000.0)]
    term_tissue = [("term 01", 0), ("term 02", 0), ("term 03", 1),
                   ("term 04", 2), ("term 05", 2), ("term 06", 3),
                   ("term 07", 4), ("term 08", 4)]

    # viewpoints: every third one a control, the others mitochondrial
    offsets = [15, -20, 25, -12, 30, -18, 22, -27, 14, -24]
    genes = []
    n_m = n_c = 0
    for c in chroms:
        for k in range(a.vp_per_chrom):
            tss = a.first_vp + k * a.vp_spacing + int(rng.integers(0, res))
            strand = "+" if rng.random() < 0.5 else "-"
            if k % 3 == 2:
                n_c += 1
                genes.append(dict(gene=f"CTRL{n_c:02d}", set="control",
                                  chrom=c, tss=tss, strand=strand))
            else:
                n_m += 1
                d = offsets[(n_m - 1) % len(offsets)]
                genes.append(dict(gene=f"MITO{n_m:02d}",
                                  set="trend" if n_m % 2 else "flat",
                                  chrom=c, tss=tss, strand=strand,
                                  partner_bin=tss // res + d))

    # annotation
    with gzip.open(out / "gencode.gtf.gz", "wt") as fh:
        fh.write("##description: synthetic annotation\n")
        for g in genes:
            if g["strand"] == "+":
                s, e = g["tss"] + 1, g["tss"] + 20_000
            else:
                s, e = g["tss"] - 20_000, g["tss"]
            attr = (f'gene_id "ID_{g["gene"]}"; gene_type "protein_coding"; '
                    f'gene_name "{g["gene"]}";')
            fh.write(f'{g["chrom"]}\tSYN\tgene\t{s}\t{e}\t.\t{g["strand"]}\t.\t'
                     f'{attr}\n')
            fh.write(f'{g["chrom"]}\tSYN\ttranscript\t{s}\t{e}\t.\t'
                     f'{g["strand"]}\t.\t{attr} transcript_id "T_{g["gene"]}"; '
                     f'tag "MANE_Select";\n')
    (out / "mito.txt").write_text(
        "".join(g["gene"] + "\n" for g in genes if g["set"] != "control"))
    (out / "controls.txt").write_text(
        "".join(g["gene"] + "\n" for g in genes if g["set"] == "control"))
    (out / "chrom.sizes").write_text(
        "".join(f"{c}\t{a.chrom_len}\n" for c in chroms))
    with open(out / "sets.tsv", "w") as fh:
        fh.write("source\tset\tgene\n")
        for g in genes:
            if g["set"] != "control":
                fh.write(f"synthetic\t{g['set']}\t{g['gene']}\n")

    def fold(g, mtcn):
        if g["set"] == "trend":
            return a.fold * (mtcn / a.mtcn_ref) ** a.beta
        return a.fold

    with open(out / "planted.tsv", "w") as fh:
        fh.write("gene\tset\tchrom\ttss\tpartner\tfold_min\tfold_max\n")
        for g in genes:
            if g["set"] == "control":
                continue
            fs = [fold(g, m) for _, m in tissues]
            fh.write(f"{g['gene']}\t{g['set']}\t{g['chrom']}\t{g['tss']}\t"
                     f"{g['partner_bin'] * res + res // 2}\t{min(fs):.3f}\t"
                     f"{max(fs):.3f}\n")

    bins = pd.DataFrame({
        "chrom": np.repeat(chroms, nbin),
        "start": np.tile(np.arange(nbin) * res, len(chroms)),
    })
    bins["end"] = bins["start"] + res

    # pixel layout shared by every library: (i, i + d), d = 0 .. band
    ii, dd = [], []
    for d in range(a.band + 1):
        ii.append(np.arange(nbin - d))
        dd.append(np.full(nbin - d, d))
    ii = np.concatenate(ii); dd = np.concatenate(dd)
    jj = ii + dd
    order = np.lexsort((jj, ii))
    ii, jj, dd = ii[order], jj[order], dd[order]
    base = a.amplitude / (1.0 + dd)
    base[dd == 0] *= 2.0
    pix_index = {(int(i), int(j)): k for k, (i, j) in
                 enumerate(zip(ii, jj)) if abs(int(i) - int(j)) >= 5}

    rows = []
    n = 0
    for term, ti in term_tissue:
        gtex, mtcn = tissues[ti]
        for _ in range(a.libs_per_term):
            n += 1
            acc = f"SYN{n:02d}"
            depth = rng.uniform(0.7, 1.4)
            counts_all, b1_all, b2_all = [], [], []
            for ci, c in enumerate(chroms):
                lam = base * depth
                for g in genes:
                    if g["chrom"] != c or g["set"] == "control":
                        continue
                    f = fold(g, mtcn)
                    t, p = g["tss"] // res, g["partner_bin"]
                    for di in (-1, 0, 1):
                        for dj in (-1, 0, 1):
                            m = 1 + (f - 1) / (2 ** (abs(di) + abs(dj)))
                            x, y = sorted((t + di, p + dj))
                            k = pix_index.get((x, y))
                            if k is not None:
                                lam[k] *= m
                lam = lam * rng.gamma(a.shape, 1.0 / a.shape, size=lam.size)
                cnt = rng.poisson(lam)
                nz = cnt > 0
                b1_all.append(ii[nz] + ci * nbin)
                b2_all.append(jj[nz] + ci * nbin)
                counts_all.append(cnt[nz])
            pixels = pd.DataFrame({"bin1_id": np.concatenate(b1_all),
                                   "bin2_id": np.concatenate(b2_all),
                                   "count": np.concatenate(counts_all)})
            path = out / "libs" / f"{acc}.cool"
            cooler.create_cooler(str(path), bins, pixels, ordered=True,
                                 dtypes={"count": np.int32})
            if acc == a.hic_pairs:
                write_pairs(out / "libs" / f"{acc}.pairs", bins, pixels, res)
            if acc == a.mcool:
                cooler.zoomify_cooler(str(path), str(out / "libs" / f"{acc}.mcool"),
                                      resolutions=[res, 2 * res],
                                      chunksize=10_000_000)
            rows.append((acc, term, gtex, mtcn, f"libs/{acc}.cool"))

    with open(out / "tissues.tsv", "w") as fh:
        fh.write("accession\tterm_name\tgtex_tissue\tmtcn\ts3_url\t"
                 "download_url\tpath\n")
        for acc, term, gtex, mtcn, path in rows:
            fh.write(f"{acc}\t{term}\t{gtex}\t{mtcn:.0f}\t\t\t{path}\n")
    print(f"{len(rows)} libraries, {len(genes)} viewpoints "
          f"({n_m} mito, {n_c} control) written to {out}/")


def write_pairs(path, bins, pixels, res):
    """Juicer short format with score, sorted by chromosome and position."""
    b1 = bins.iloc[pixels["bin1_id"].to_numpy()]
    b2 = bins.iloc[pixels["bin2_id"].to_numpy()]
    df = pd.DataFrame({"c1": b1["chrom"].to_numpy(),
                       "p1": b1["start"].to_numpy() + res // 2,
                       "c2": b2["chrom"].to_numpy(),
                       "p2": b2["start"].to_numpy() + res // 2,
                       "n": pixels["count"].to_numpy()})
    df = df.sort_values(["c1", "p1", "p2"], kind="mergesort")
    with open(path, "w") as fh:
        for c1, p1, c2, p2, k in df.itertuples(index=False):
            fh.write(f"0 {c1} {p1} 0 16 {c2} {p2} 1 {k}\n")


if __name__ == "__main__":
    main()
