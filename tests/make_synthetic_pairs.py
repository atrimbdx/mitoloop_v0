#!/usr/bin/env python3
"""
make_synthetic_pairs.py - synthetic paired Hi-C libraries for testing
v4c_pairs.py and v4c_pairs.R

Purpose
-------
Builds two experiments of two states each, on one synthetic chromosome:

  planted   a fraction (--frac) of the mito promoters gain a loop in the
            plus state
  null      the two states differ only by depth and by Poisson noise

A correct pipeline shows a loss of concordance, or an excess of changing
loops, for mito against control promoters in `planted`, and neither in
`null`.

Inputs
------
None; everything is generated from the options.

Outputs (in --outdir)
---------------------
  <experiment>_<state>.cool   one library each, raw integer counts
  pairs.tsv                   manifest for v4c_pairs.py (format cool), columns
                              experiment, state, label, path, format, genome,
                              resolution
  synth.gtf.gz                one 'gene' record per viewpoint
  mito.txt, controls.txt      the two gene lists, one symbol per line
  with --hic-pairs, also:
  <experiment>_<state>.pairs.txt  Juicer Tools 'short with score' input
                              (str chr pos frag str chr pos frag count), one
                              line per non-zero pixel, sorted by position
  chrom.sizes                 for `juicer_tools pre`
  pairs_hic.tsv               the same manifest pointing at
                              <experiment>_<state>.hic (format hic), to be
                              built with `juicer_tools pre`

Method
------
Expected contacts between two bins at distance d (in bins) follow
--decay-scale / (1 + d)^--decay-exponent, up to --band bins. Every promoter
has two shared loops (offsets --shared, in bins) in both states; a loop is a
3 x 3 block of pixels multiplied by --boost. In the planted experiment the
first round(--frac * --n-genes) mito promoters receive one extra loop in the
plus state, at offsets taken in turn from --planted-offsets. The plus state is
sampled --depth-ratio times deeper than the minus state, so that the depth
matching of v4c_pairs.R is exercised. Counts are Poisson draws.

Usage
-----
    python3 tests/make_synthetic_pairs.py --outdir results/synthetic_pairs
    python3 tests/make_synthetic_pairs.py --outdir results/synthetic_pairs \\
        --n-genes 12 --frac 0.25 --hic-pairs

References
----------
  Abdennur and Mirny 2020. Cooler: scalable storage for Hi-C data and other
    genomically labeled arrays. Bioinformatics.
  Durand et al. 2016. Juicer provides a one-click system for analyzing
    loop-resolution Hi-C experiments. Cell Systems.
"""

import argparse
import gzip
import os

import numpy as np
import pandas as pd
import cooler


def int_list(text):
    return [int(x) for x in text.split(",") if x.strip()]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="synthetic_pairs",
                   help="output directory (default %(default)s)")
    p.add_argument("--n-genes", type=int, default=6,
                   help="viewpoints per set, mito and control (default %(default)s)")
    p.add_argument("--frac", type=float, default=0.5,
                   help="fraction of mito viewpoints gaining a loop in the "
                        "plus state of `planted` (default %(default)s)")
    p.add_argument("--resolution", type=int, default=10_000,
                   help="bin size, bp (default %(default)s)")
    p.add_argument("--spacing", type=int, default=1_500_000,
                   help="distance between consecutive viewpoints, bp "
                        "(default %(default)s)")
    p.add_argument("--band", type=int, default=130,
                   help="largest simulated distance, in bins; a little over "
                        "the 1 Mb profile window at 10 kb (default %(default)s)")
    p.add_argument("--chrom", default="chr1",
                   help="chromosome name (default %(default)s)")
    p.add_argument("--genome", default="synth",
                   help="genome key written to the manifest (default %(default)s)")
    p.add_argument("--decay-scale", type=float, default=260.0,
                   help="expected contacts at distance 0 (default %(default)s)")
    p.add_argument("--decay-exponent", type=float, default=0.9,
                   help="exponent of the distance decay (default %(default)s)")
    p.add_argument("--shared", type=int_list, default=[20, -35],
                   help="offsets, in bins, of the loops present in both states "
                        "(default 20,-35)")
    p.add_argument("--planted-offsets", type=int_list,
                   default=[50, -60, 45, -70, 55, -45, 65, -55],
                   help="offsets, in bins, of the planted loops, used in turn "
                        "(default 50,-60,45,-70,55,-45,65,-55)")
    p.add_argument("--boost", type=float, default=4.0,
                   help="multiplier of the expected contacts at a loop "
                        "(default %(default)s)")
    p.add_argument("--depth-ratio", type=float, default=1.6,
                   help="depth of the plus state relative to the minus state "
                        "(default %(default)s)")
    p.add_argument("--seed", type=int, default=7,
                   help="seed of the Poisson draws (default %(default)s)")
    p.add_argument("--hic-pairs", action="store_true",
                   help="also write Juicer Tools input, chrom.sizes and "
                        "pairs_hic.tsv")
    return p.parse_args()


def main():
    args = parse_args()
    out = args.outdir
    os.makedirs(out, exist_ok=True)
    res = args.resolution
    n_genes = args.n_genes
    chrom_len = 2 * args.spacing + 2 * n_genes * args.spacing
    n_bin = chrom_len // res
    band = args.band
    rng = np.random.default_rng(args.seed)

    genes = []
    for i in range(n_genes):
        genes.append((f"MITOG{i+1}", args.spacing * (i + 1), "mito"))
    for i in range(n_genes):
        genes.append((f"CTRLG{i+1}", args.spacing * (n_genes + i + 1), "control"))

    offs = args.planted_offsets
    planted = {f"MITOG{i+1}": offs[i % len(offs)]
               for i in range(max(1, int(round(args.frac * n_genes))))}

    def lam_band(extra):
        """lam[d][i] = expected contacts between bin i and bin i + d."""
        lam = {}
        for d in range(1, band + 1):
            lam[d] = np.full(n_bin - d,
                             args.decay_scale / (1.0 + d) ** args.decay_exponent)
        for i, j, mult in extra:
            d = abs(j - i)
            a = min(i, j)
            if 1 <= d <= band and 0 <= a < n_bin - d:
                lam[d][a] *= mult
        return lam

    def boosts(with_planted):
        cells = []
        for name, tss, _s in genes:
            g = tss // res
            loop_offsets = list(args.shared)
            if with_planted and name in planted:
                loop_offsets.append(planted[name])
            for o in loop_offsets:
                p = g + o
                for a in range(g - 1, g + 2):
                    for b in range(p - 1, p + 2):
                        cells.append((a, b, args.boost))
        return cells

    def write_cool(path, lam, scale):
        bins = pd.DataFrame({
            "chrom": [args.chrom] * n_bin,
            "start": np.arange(n_bin) * res,
            "end": (np.arange(n_bin) + 1) * res})
        b1, b2, cnt = [], [], []
        for d, v in lam.items():
            c = rng.poisson(v * scale)
            nz = np.nonzero(c)[0]
            if nz.size == 0:
                continue
            b1.append(nz)
            b2.append(nz + d)
            cnt.append(c[nz])
        px = pd.DataFrame({"bin1_id": np.concatenate(b1),
                           "bin2_id": np.concatenate(b2),
                           "count": np.concatenate(cnt).astype(int)})
        cooler.create_cooler(path, bins, px, dtypes={"count": int},
                             ordered=False, triucheck=False)
        return px

    def write_juicer_pairs(path, px):
        px = px.sort_values(["bin1_id", "bin2_id"])
        p1 = px["bin1_id"].to_numpy() * res + res // 2
        p2 = px["bin2_id"].to_numpy() * res + res // 2
        with open(path, "w") as fh:
            for a, b, c in zip(p1, p2, px["count"].to_numpy()):
                fh.write(f"0 {args.chrom} {a} 0 0 {args.chrom} {b} 1 {c}\n")

    rows = []
    for exp in ("planted", "null"):
        for state in ("minus", "plus"):
            lam = lam_band(boosts(exp == "planted" and state == "plus"))
            scale = 1.0 if state == "minus" else args.depth_ratio
            label = f"{exp}_{state}"
            path = os.path.join(out, f"{label}.cool")
            px = write_cool(path, lam, scale)
            print(f"{exp:8s} {state:5s} {int(px['count'].sum()):>12,d} "
                  f"contacts -> {path}")
            rows.append((exp, state, label, path, "cool", args.genome, res))
            if args.hic_pairs:
                write_juicer_pairs(os.path.join(out, f"{label}.pairs.txt"), px)

    head = "experiment\tstate\tlabel\tpath\tformat\tgenome\tresolution\n"
    with open(os.path.join(out, "pairs.tsv"), "w") as fh:
        fh.write(head)
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    if args.hic_pairs:
        with open(os.path.join(out, "pairs_hic.tsv"), "w") as fh:
            fh.write(head)
            for exp, state, label, _p, _f, genome, r in rows:
                fh.write("\t".join(str(x) for x in (
                    exp, state, label, os.path.join(out, f"{label}.hic"),
                    "hic", genome, r)) + "\n")
        with open(os.path.join(out, "chrom.sizes"), "w") as fh:
            fh.write(f"{args.chrom}\t{chrom_len}\n")

    with gzip.open(os.path.join(out, "synth.gtf.gz"), "wt") as fh:
        for name, tss, _s in genes:
            fh.write(f"{args.chrom}\tsynth\tgene\t{tss+1}\t{tss+2000}\t.\t+\t.\t"
                     f'gene_id "{name}"; gene_name "{name}";\n')

    with open(os.path.join(out, "mito.txt"), "w") as fh:
        fh.write("\n".join(n for n, _t, s in genes if s == "mito") + "\n")
    with open(os.path.join(out, "controls.txt"), "w") as fh:
        fh.write("\n".join(n for n, _t, s in genes if s == "control") + "\n")

    print(f"\nmanifest, GTF and gene lists in {out}/")
    print(f"planted loops (plus state only, offsets in bins): {planted}")


if __name__ == "__main__":
    main()
