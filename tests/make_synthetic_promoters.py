#!/usr/bin/env python3
"""
make_synthetic_promoters.py - synthetic inputs for promoter_features.py

Purpose
    Build a small genome with known promoters: target ("mito") promoters are
    CpG-rich and carry an NRF1 site; a third of the control promoters are
    CpG-rich without it, the others CpG-poor. Some target genes have a
    bidirectional partner, some loop onto another promoter, and one target
    promoter is methylated everywhere except in muscle.

Inputs
    none; every quantity is an option.

Outputs (in --out)
    genome.fa                 one chromosome, chr1, 60 bp lines
    gencode.gtf.gz            gene and MANE Select transcript per gene, plus
                              bidirectional partners (lncRNA) and hub genes
    genes.txt, controls.txt   gene lists
    cpg_islands.txt.gz        UCSC cpgIslandExt layout (leading bin column)
    loops.tsv                 loops.tsv format
    re2g_elements.tsv         re2g_panel.py format, one element per target
                              gene
    meth/                     muscle_1.bedMethyl.gz, muscle_2.bedMethyl.gz
                              (ENCODE bedMethyl), liver.cov.gz (Bismark),
                              brain.bedGraph.gz (fraction, no coverage),
                              heart.bedGraph.gz (percent, coverage)
    methylation.tsv           --meth-table
    tpm.gct                   GTEx median TPM layout
    truth.tsv                 gene, set, cpg_rich, nrf1, bidirectional,
                              promoter_loop, muscle_specific

Method
    Background sequence has 40% GC with most CpGs mutated to TpG; CpG-rich
    promoters (TSS - 1 kb to TSS + 1 kb) have 65% GC without depletion. The
    NRF1 site GCGCATGCGC is written 200 bp upstream of each target TSS, after
    every other NRF1 match of the genome has been destroyed. Methylation is
    written at every CpG of every promoter and element: 0.05 at target
    promoters, 0.1 at CpG-rich control promoters, 0.8 at CpG-poor ones, 0.3
    at elements; MITOG3 is at 0.8 outside muscle.

Usage
    python3 make_synthetic_promoters.py --out synth
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np

SPACING = 200_000
NRF1 = "GCGCATGCGC"


def rand_seq(rng, n, gc, deplete):
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    s = rng.choice(list("ACGT"), size=n, p=p)
    if deplete:
        cg = np.flatnonzero((s[:-1] == "C") & (s[1:] == "G"))
        s[cg[rng.random(len(cg)) < 0.8] + 1] = "A"
    return s


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True)
    p.add_argument("--n-mito", type=int, default=30)
    p.add_argument("--n-control", type=int, default=30)
    p.add_argument("--seed", type=int, default=5)
    a = p.parse_args()
    rng = np.random.default_rng(a.seed)
    out = Path(a.out)
    (out / "meth").mkdir(parents=True, exist_ok=True)

    genes = [(f"MITOG{i + 1}", "mito") for i in range(a.n_mito)] + \
        [(f"CTRLG{i + 1}", "control") for i in range(a.n_control)]
    L = SPACING * (len(genes) + 2)
    seq = rand_seq(rng, L, 0.40, True)

    info = []
    for i, (g, s) in enumerate(genes):
        tss = SPACING * (i + 1)
        strand = "+" if i % 2 == 0 else "-"
        k = int(g[5:])
        rich = s == "mito" or k <= a.n_control // 3
        info.append({"gene": g, "set": s, "tss": tss, "strand": strand,
                     "cpg_rich": int(rich), "nrf1": int(s == "mito"),
                     "bidirectional": int(s == "mito" and k <= 10),
                     "promoter_loop": int(s == "mito" and k <= 5),
                     "muscle_specific": int(g == "MITOG3")})
        if rich:
            seq[tss - 1000:tss + 1000] = rand_seq(rng, 2000, 0.65, False)
        if s == "mito":
            # a CpG-rich element 50 kb away, so that it has CpGs to measure
            seq[tss + 49_600:tss + 50_400] = rand_seq(rng, 800, 0.65, False)

    # no NRF1 site anywhere but where it is planted
    text = "".join(seq)
    text = re.sub(NRF1, "GCGCAAGCGC", text)
    text = re.sub("GCGCACGCGC", "GCGCAAGCGC", text)
    seq = np.array(list(text))
    for r in info:
        if r["nrf1"]:
            pos = r["tss"] - 200 if r["strand"] == "+" else r["tss"] + 190
            seq[pos:pos + len(NRF1)] = list(NRF1)
    text = "".join(seq)
    with open(out / "genome.fa", "w") as fh:
        fh.write(">chr1\n")
        for j in range(0, len(text), 60):
            fh.write(text[j:j + 60] + "\n")

    with gzip.open(out / "gencode.gtf.gz", "wt") as fh:
        n = 0

        def gene(name, tss, strand, gtype):
            nonlocal n
            n += 1
            attr = (f'gene_id "ENSG{n:011d}.1"; gene_type "{gtype}"; '
                    f'gene_name "{name}";')
            if strand == "+":
                gs, ge = tss + 1, tss + 20000
            else:
                gs, ge = tss - 20000, tss
            # a gene record that starts away from the TSS: the MANE Select
            # transcript must win
            fh.write(f"chr1\tSYN\tgene\t{gs - 3000 if strand == '+' else gs}"
                     f"\t{ge if strand == '+' else ge + 3000}\t.\t{strand}\t."
                     f"\t{attr}\n")
            fh.write(f"chr1\tSYN\ttranscript\t{gs}\t{ge}\t.\t{strand}\t.\t"
                     f'{attr} transcript_id "ENST{n:011d}"; '
                     'tag "MANE_Select";\n')
        for r in info:
            gene(r["gene"], r["tss"], r["strand"], "protein_coding")
            if r["bidirectional"]:
                other = "-" if r["strand"] == "+" else "+"
                d = -300 if r["strand"] == "+" else 300
                gene(f"BIDIR_{r['gene']}", r["tss"] + d, other, "lncRNA")
            if r["promoter_loop"]:
                gene(f"HUB_{r['gene']}", r["tss"] + 60_000, "+",
                     "protein_coding")

    Path(out / "genes.txt").write_text(
        "".join(f"{r['gene']}\n" for r in info if r["set"] == "mito"))
    Path(out / "controls.txt").write_text(
        "".join(f"{r['gene']}\n" for r in info if r["set"] == "control"))

    with gzip.open(out / "cpg_islands.txt.gz", "wt") as fh:
        for r in info:
            if r["cpg_rich"]:
                fh.write(f"585\tchr1\t{r['tss'] - 900}\t{r['tss'] + 900}\t"
                         f"CpG: 150\n")

    with open(out / "loops.tsv", "w") as fh:
        fh.write("gene\tchrom\tvp_pos\tsummit\tstart\tend\tn_terms\tcalled\t"
                 "loop_id\n")
        for r in info:
            if r["set"] != "mito":
                continue
            far = r["tss"] - 90_000
            fh.write(f"{r['gene']}\tchr1\t{r['tss']}\t{far}\t{far - 5000}\t"
                     f"{far + 5000}\t3\tpsoas muscle\t{r['gene']}:{far}\n")
            if r["promoter_loop"]:
                sm = r["tss"] + 62_000
                fh.write(f"{r['gene']}\tchr1\t{r['tss']}\t{sm}\t{sm - 5000}\t"
                         f"{sm + 5000}\t3\tpsoas muscle\t{r['gene']}:{sm}\n")

    with open(out / "re2g_elements.tsv", "w") as fh:
        fh.write("gene\tset\tchrom\tstart\tend\tdistance\tprevalence\n")
        for r in info:
            if r["set"] == "mito":
                m = r["tss"] + 50_000
                fh.write(f"{r['gene']}\tmito\tchr1\t{m - 400}\t{m + 400}\t"
                         f"50000\t6\n")

    # methylation at every CpG of the promoters and elements
    def level(r, group):
        if r["muscle_specific"]:
            return 0.05 if group == "muscle" else 0.8
        if r["set"] == "mito":
            return 0.05
        return 0.1 if r["cpg_rich"] else 0.8

    groups = {"muscle_1": "muscle", "muscle_2": "muscle", "liver": "liver",
              "brain": "brain", "heart": "heart"}
    calls = {s: [] for s in groups}
    for r in info:
        for lo, hi, el in ((r["tss"] - 1500, r["tss"] + 1500, False),
                           (r["tss"] + 49_600, r["tss"] + 50_400, True)):
            if el and r["set"] != "mito":
                continue
            for c in [m.start() for m in re.finditer("CG", text[lo:hi])]:
                for smp, grp in groups.items():
                    b = 0.3 if el else level(r, grp)
                    calls[smp].append((lo + c, b))
    for smp in ("muscle_1", "muscle_2"):
        with gzip.open(out / "meth" / f"{smp}.bedMethyl.gz", "wt") as fh:
            fh.write('track name="synthetic" itemRgb="On"\n')
            for pos, b in calls[smp]:
                for q, strand in ((pos, "+"), (pos + 1, "-")):
                    fh.write(f"chr1\t{q}\t{q + 1}\t.\t10\t{strand}\t{q}\t"
                             f"{q + 1}\t0,255,0\t10\t{round(100 * b)}\n")
            # a CpG under the coverage threshold, fully methylated
            pos = info[0]["tss"]
            fh.write(f"chr1\t{pos}\t{pos + 1}\t.\t1\t+\t{pos}\t{pos + 1}\t"
                     f"0,255,0\t1\t100\n")
    with gzip.open(out / "meth" / "liver.cov.gz", "wt") as fh:
        for pos, b in calls["liver"]:
            m = round(20 * b)
            fh.write(f"chr1\t{pos + 1}\t{pos + 1}\t{100 * m / 20:.2f}\t{m}\t"
                     f"{20 - m}\n")
    with gzip.open(out / "meth" / "brain.bedGraph.gz", "wt") as fh:
        for pos, b in calls["brain"]:
            fh.write(f"chr1\t{pos}\t{pos + 2}\t{b:.3f}\n")
    with gzip.open(out / "meth" / "heart.bedGraph.gz", "wt") as fh:
        for pos, b in calls["heart"]:
            fh.write(f"chr1\t{pos}\t{pos + 2}\t{100 * b:.1f}\t12\n")
    with open(out / "methylation.tsv", "w") as fh:
        fh.write("sample\tgroup\tpath\tformat\tgtex_tissue\n")
        fh.write("muscle_1\tskeletal muscle\tmeth/muscle_1.bedMethyl.gz\t\t"
                 "Muscle - Skeletal\n")
        fh.write("muscle_2\tskeletal muscle\tmeth/muscle_2.bedMethyl.gz\t"
                 "bedmethyl\tMuscle - Skeletal\n")
        fh.write("liver\tliver\tmeth/liver.cov.gz\t\tLiver\n")
        fh.write("brain\tbrain\tmeth/brain.bedGraph.gz\t\tBrain - Cortex\n")
        fh.write("heart\theart\tmeth/heart.bedGraph.gz\t\t"
                 "Heart - Left Ventricle\n")

    tissues = ["Brain - Cortex", "Heart - Left Ventricle", "Liver",
               "Muscle - Skeletal"]
    with open(out / "tpm.gct", "w") as fh:
        fh.write("#1.2\n")
        fh.write(f"{len(info)}\t{len(tissues)}\n")
        fh.write("Name\tDescription\t" + "\t".join(tissues) + "\n")
        for n, r in enumerate(info):
            if r["muscle_specific"]:
                v = [1, 3, 2, 200]
            else:
                v = list(rng.uniform(10, 50, len(tissues)).round(2))
            # Ensembl ids follow the order of gene() calls; match by symbol
            fh.write(f"ENSGX{n:011d}.1\t{r['gene']}\t"
                     + "\t".join(str(x) for x in v) + "\n")

    with open(out / "truth.tsv", "w") as fh:
        cols = ["gene", "set", "tss", "cpg_rich", "nrf1", "bidirectional",
                "promoter_loop", "muscle_specific"]
        fh.write("\t".join(cols) + "\n")
        for r in info:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
