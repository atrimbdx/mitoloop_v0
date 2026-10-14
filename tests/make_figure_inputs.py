#!/usr/bin/env python3
"""
make_figure_inputs.py - synthetic inputs for the figure scripts

Purpose
    Builds a small contact map with planted loops and every table the figure
    scripts read, so that locus_figure.py, folding_figure.py and
    svg_to_pptx.py can be run end to end without any download.

Inputs
    none; the layout is set by the constants below and the noise is seeded

Outputs (in OUTDIR)
    synthetic.cool           one chromosome, 10 kb bins, distance decay plus
                             planted loops, Poisson counts, no balancing
    pairs.txt, chrom.sizes   the same contacts as Juicer short-with-score
                             lines, for `juicer_tools pre`
    genes.gtf                gene and transcript records (gene_name,
                             gene_type, a MANE_Select tag on transcripts)
    loops.tsv                gene, chrom, vp_pos, summit, start, end,
                             n_terms, called, loop_id, for three viewpoints
    re2g_elements.tsv        gene, set, chrom, start, end, distance,
                             prevalence
    extract.npz              keys written by `v4c_pipeline.py extract`
    profiles.npz             keys profiles, genes, gene_set, offsets,
                             accession, tissue, resolution (no position)
    planted.tsv              gene, summit: the planted loops, for assertions

Method
    Expected contacts between bins i and j are 200 / (1 + |i - j|) ** 0.9;
    each planted loop multiplies a 3 x 3 block around (viewpoint, summit) by
    ENRICHMENT. Counts are drawn from a Poisson law with a set
    seed. The virtual 4C profiles of the npz files are the viewpoint rows of
    that matrix, summed over the anchor bins, over +/- WINDOW.

Usage
    python3 tests/make_figure_inputs.py tmp/figures

References
    Abdennur N, Mirny LA. 2020. Cooler: scalable storage for Hi-C data and
        other genomically labeled arrays. Bioinformatics.
    Durand NC et al. 2016. Juicer provides a one-click system for analyzing
        loop-resolution Hi-C experiments. Cell Systems.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import cooler

RES = 10_000
CHROM = "chr1"
LENGTH = 3_000_000
WINDOW = 500_000
ANCHOR_HALF = 10_000
MAX_SEP = 120                   # bins of separation stored
ENRICHMENT = 6.0

# name, TSS (0-based), strand, gene_type
GENES = [
    ("GENEA", 1_000_000, "+", "protein_coding"),
    ("GENEB", 1_400_000, "-", "protein_coding"),
    ("GENEC", 1_800_000, "+", "protein_coding"),
    ("GENED", 1_150_000, "+", "protein_coding"),     # carries no loop
    ("LNCE", 1_600_000, "+", "lncRNA"),
]
# viewpoint gene, summit, n_terms
LOOPS = [
    ("GENEA", 1_250_000, 12),
    ("GENEA", 750_000, 5),
    ("GENEB", 1_650_000, 9),
    ("GENEC", 1_400_000, 4),
]
TERMS = ["psoas muscle", "gastrocnemius medialis", "heart left ventricle",
         "liver", "spleen", "lung", "thyroid gland", "stomach", "adrenal gland",
         "pancreas", "testis", "ovary"]


def tss_of(name):
    return next(g[1] for g in GENES if g[0] == name)


def build_matrix(rng):
    n = LENGTH // RES
    i, j = np.triu_indices(n)
    keep = (j - i) <= MAX_SEP
    i, j = i[keep], j[keep]
    lam = 200.0 / (1.0 + (j - i)) ** 0.9
    lam = lam.astype(float)
    boost = np.ones(len(i))
    for gene, summit, _ in LOOPS:
        a, b = sorted((tss_of(gene) // RES, summit // RES))
        boost[(np.abs(i - a) <= 1) & (np.abs(j - b) <= 1)] = ENRICHMENT
    counts = rng.poisson(lam * boost)
    return n, i, j, counts


def write_cool(path, n, i, j, counts):
    bins = pd.DataFrame({"chrom": CHROM, "start": np.arange(n) * RES,
                         "end": np.minimum((np.arange(n) + 1) * RES, LENGTH)})
    nz = counts > 0
    pixels = pd.DataFrame({"bin1_id": i[nz], "bin2_id": j[nz],
                           "count": counts[nz].astype(np.int32)})
    cooler.create_cooler(str(path), bins, pixels, ordered=True)


def write_pairs(out, i, j, counts):
    nz = counts > 0
    with open(out / "pairs.txt", "w") as fh:
        for a, b, c in zip(i[nz], j[nz], counts[nz]):
            fh.write(f"0 {CHROM} {a * RES + RES // 2} 0 "
                     f"16 {CHROM} {b * RES + RES // 2} 1 {c}\n")
    (out / "chrom.sizes").write_text(f"{CHROM}\t{LENGTH}\n")


def write_gtf(path):
    with open(path, "w") as fh:
        fh.write("##description: synthetic annotation\n")
        for name, tss, strand, kind in GENES:
            s, e = (tss + 1, tss + 30_000) if strand == "+" else \
                (tss - 30_000 + 1, tss)
            attr = (f'gene_id "{name}.1"; gene_type "{kind}"; '
                    f'gene_name "{name}";')
            fh.write(f"{CHROM}\tSYN\tgene\t{s}\t{e}\t.\t{strand}\t.\t{attr}\n")
            fh.write(f"{CHROM}\tSYN\ttranscript\t{s}\t{e}\t.\t{strand}\t.\t"
                     f'{attr} transcript_id "{name}-201"; tag "MANE_Select";\n')


def write_loops(path):
    with open(path, "w") as fh:
        fh.write("gene\tchrom\tvp_pos\tsummit\tstart\tend\tn_terms\tcalled\t"
                 "loop_id\n")
        for gene, summit, nt in LOOPS:
            fh.write(f"{gene}\t{CHROM}\t{tss_of(gene)}\t{summit}\t"
                     f"{summit - RES}\t{summit + RES}\t{nt}\t"
                     f"{'|'.join(TERMS[:nt])}\t{gene}:{summit}\n")


def write_elements(path):
    with open(path, "w") as fh:
        fh.write("gene\tset\tchrom\tstart\tend\tdistance\tprevalence\n")
        # replicable element on the first GENEA anchor, rare one elsewhere
        fh.write(f"GENEA\tmito\t{CHROM}\t1248000\t1252000\t250000\t80\n")
        fh.write(f"GENEB\tmito\t{CHROM}\t1648000\t1651000\t250000\t10\n")


def profiles(n, i, j, counts):
    dense = np.zeros((n, n))
    dense[i, j] = counts
    dense[j, i] = counts
    half = WINDOW // RES
    offsets = np.arange(-half, half + 1) * RES
    names = [g for g, *_ in GENES[:3]]
    mat = np.full((len(names), len(offsets)), np.nan, dtype=np.float32)
    for k, g in enumerate(names):
        vp = tss_of(g)
        lo = (vp - ANCHOR_HALF) // RES
        hi = (vp + ANCHOR_HALF) // RES + 1
        row = dense[lo:hi].sum(0)
        for m, off in enumerate(offsets):
            b = (vp + off) // RES
            if 0 <= b < n:
                mat[k, m] = row[b]
    return names, offsets, mat


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "figure_inputs")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(11)
    n, i, j, counts = build_matrix(rng)
    write_cool(out / "synthetic.cool", n, i, j, counts)
    write_pairs(out, i, j, counts)
    write_gtf(out / "genes.gtf")
    write_loops(out / "loops.tsv")
    write_elements(out / "re2g_elements.tsv")

    names, offsets, mat = profiles(n, i, j, counts)
    np.savez_compressed(
        out / "extract.npz", counts=mat, raw_total=np.nansum(mat, axis=1),
        names=np.array(names), chrom=np.array([CHROM] * len(names)),
        pos=np.array([tss_of(g) for g in names], dtype=np.int64),
        set=np.array(["mito"] * len(names)), offsets=offsets,
        accession="SYN0001", term="psoas muscle", gtex="Muscle - Skeletal",
        mtcn=1500.0, resolution=RES, anchor_half=ANCHOR_HALF, via="local")
    np.savez_compressed(
        out / "profiles.npz", profiles=mat, genes=np.array(names),
        gene_set=np.array(["mito"] * len(names)), offsets=offsets,
        accession="SYN0001", tissue="psoas muscle", resolution=RES)

    with open(out / "planted.tsv", "w") as fh:
        fh.write("gene\tsummit\n")
        for gene, summit, _ in LOOPS:
            fh.write(f"{gene}\t{summit}\n")
    print(f"synthetic inputs written to {out}")


if __name__ == "__main__":
    main()
