#!/usr/bin/env python3
"""
make_synthetic_l1.py - synthetic inputs for l1_score.py

Purpose
    Build a tissue atlas with RNA counts and protein abundances measured in
    the same samples, where the mode of control of every gene is known, a
    second atlas without replicates, and a patient fibroblast cohort with
    planted RNA outliers, plus external per-gene tables.

Inputs
    none

Outputs (in --out)
    gencode.gtf.gz, genes.txt, complexes.tsv
    atlas_rna.gct, atlas_protein.tsv, atlas_samples.tsv
                                 16 tissues x 4 donors; ischemic time column
                                 SMTSISCH
    hpa_rna.tsv, hpa_protein.tsv, hpa_samples.tsv
                                 the same tissues, one sample each, TPM
    patients_rna.tsv, patients_protein.tsv, patients_samples.tsv
                                 60 individuals
    buffering.tsv, expected_high.txt, expected_low.txt
    truth.tsv                    gene, class

Method
    Classes of the 80 target genes (MITOG1-80):
      transcriptional  1-15, 73-80: protein follows the mRNA profile
      buffered        16-30: mRNA varies across tissues, protein does not
      post            31-40: protein varies, mRNA does not
      flat            41-50: neither varies
      noisy           51-60: transcriptional, protein measured with large
                      noise (the correction for attenuation must recover it)
      subunit         61-72: members of one complex; mRNA of each subunit
                      varies on its own around the complex profile, protein
                      follows the complex profile
    Every target gene also moves with a tissue mitochondrial content
    (+/- 2 log2) in RNA and protein, which the content adjustment must
    remove. MITOG73-80 proteins depend on the ischemic time of the sample.
    Two mtDNA genes take up to half of the reads in tissues rich in
    mitochondria. Patients: one individual per gene of MITOG1-15 and
    MITOG16-30 has its mRNA divided by 8; the protein follows for the first
    group only (missing for MITOG2). Patient proteins carry a batch shift
    (three TMT plexes) that must be removed.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=7)
    a = p.parse_args()
    rng = np.random.default_rng(a.seed)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    mito = [f"MITOG{i}" for i in range(1, 81)]
    other = [f"OTHERG{i}" for i in range(1, 201)]
    mt = ["MT-ATP8", "MT-CO1"]
    genes = mito + other
    cls = {}
    for i, g in enumerate(mito, 1):
        cls[g] = ("transcriptional" if i <= 15 or i >= 73 else
                  "buffered" if i <= 30 else "post" if i <= 40 else
                  "flat" if i <= 50 else "noisy" if i <= 60 else "subunit")
    for i, g in enumerate(other):
        cls[g] = "transcriptional" if i % 2 == 0 else "buffered"
    ens = {g: f"ENSG{n:011d}" for n, g in enumerate(genes + mt)}
    with gzip.open(out / "gencode.gtf.gz", "wt") as fh:
        for g in genes + mt:
            c = "chrM" if g in mt else "chr1"
            fh.write(f"{c}\tSYN\tgene\t1\t1000\t.\t+\t.\tgene_id "
                     f'"{ens[g]}.2"; gene_type "protein_coding"; '
                     f'gene_name "{g}";\n')
    (out / "genes.txt").write_text("\n".join(mito) + "\n")
    with open(out / "complexes.tsv", "w") as fh:
        fh.write("complex\tgene\n")
        for g in mito[60:72]:
            fh.write(f"CI subunits\t{g}\n")
    with open(out / "truth.tsv", "w") as fh:
        fh.write("gene\tclass\n")
        for g in genes:
            fh.write(f"{g}\t{cls[g]}\n")
    (out / "expected_high.txt").write_text("\n".join(mito[:15]) + "\n")
    (out / "expected_low.txt").write_text("\n".join(mito[15:30]) + "\n")
    with open(out / "buffering.tsv", "w") as fh:
        fh.write("gene\tbuffering\n")
        for g in genes:
            b = {"transcriptional": 0.1, "buffered": 0.8}.get(cls[g], 0.5)
            fh.write(f"{g}\t{b + rng.normal(0, 0.1):.3f}\n")

    T = 16
    tissues = [f"tissue_{t:02d}" for t in range(T)]
    content = rng.uniform(-2, 2, T)
    complex_prof = rng.normal(0, 1.2, T)
    base = {g: rng.uniform(6, 12) for g in genes}
    prof_r, prof_p = {}, {}
    for g in genes:
        c = cls[g]
        x = rng.normal(0, 1.2, T)
        if c in ("transcriptional", "noisy"):
            prof_r[g], prof_p[g] = x, x
        elif c == "buffered":
            prof_r[g], prof_p[g] = x, rng.normal(0, 0.05, T)
        elif c == "post":
            prof_r[g], prof_p[g] = rng.normal(0, 0.05, T), x
        elif c == "flat":
            prof_r[g], prof_p[g] = (rng.normal(0, 0.05, T),
                                    rng.normal(0, 0.05, T))
        else:
            prof_r[g] = complex_prof + rng.normal(0, 0.8, T)
            prof_p[g] = complex_prof + rng.normal(0, 0.1, T)
        if g in mito:
            prof_r[g] = prof_r[g] + content
            prof_p[g] = prof_p[g] + content

    def atlas(donors, prefix, units):
        rows_r = {g: [] for g in genes + mt}
        rows_p = {g: [] for g in genes}
        samp = []
        for t in range(T):
            isch_mean = rng.uniform(0, 3)
            for d in range(donors):
                sid = f"{prefix}{t:02d}_{d}"
                isch = isch_mean + rng.normal(0, 1.5)
                depth = rng.uniform(0.5, 2)
                samp.append((tissues[t], f"R{sid}", f"P{sid}", f"D{d}",
                             isch))
                for g in genes:
                    lr = base[g] + prof_r[g][t] + rng.normal(0, 0.2)
                    noise = 1.5 if cls[g] == "noisy" else 0.25
                    lp = base[g] + prof_p[g][t] + rng.normal(0, noise)
                    if g in mito[72:]:
                        lp += -0.6 * isch
                    rows_r[g].append(2 ** lr * depth)
                    rows_p[g].append(lp)
                tot = sum(rows_r[g][-1] for g in genes)
                for g in mt:
                    rows_r[g].append(tot * 0.5 * 2 ** (content[t] - 2))
        return rows_r, rows_p, samp

    rr, pp, samp = atlas(4, "A", "counts")
    with open(out / "atlas_rna.gct", "w") as fh:
        fh.write(f"#1.2\n{len(rr)}\t{len(samp)}\n")
        fh.write("Name\tDescription\t" + "\t".join(s[1] for s in samp) + "\n")
        for g, v in rr.items():
            fh.write(f"{ens[g]}.2\t{g}\t"
                     + "\t".join(str(int(round(x))) for x in v) + "\n")
    with open(out / "atlas_protein.tsv", "w") as fh:
        fh.write("gene\t" + "\t".join(s[2] for s in samp) + "\n")
        for g, v in pp.items():
            vals = [f"{x:.4f}" for x in v]
            if g == "OTHERG7":
                vals[0] = "NA"
            fh.write(f"{g}\t" + "\t".join(vals) + "\n")
    with open(out / "atlas_samples.tsv", "w") as fh:
        fh.write("tissue\trna_sample\tprotein_sample\tdonor\tSMTSISCH\n")
        for s in samp:
            fh.write(f"{s[0]}\t{s[1]}\t{s[2]}\t{s[3]}\t{s[4]:.3f}\n")

    rr, pp, samp = atlas(1, "H", "tpm")
    with open(out / "hpa_rna.tsv", "w") as fh:
        fh.write("gene\t" + "\t".join(s[1] for s in samp) + "\n")
        cols = list(zip(*rr.values()))
        tot = [sum(c) for c in cols]
        for g, v in rr.items():
            fh.write(f"{g}\t" + "\t".join(f"{x / tot[k] * 1e6:.4f}"
                                          for k, x in enumerate(v)) + "\n")
    with open(out / "hpa_protein.tsv", "w") as fh:
        fh.write("gene\t" + "\t".join(s[2] for s in samp) + "\n")
        for g, v in pp.items():
            fh.write(f"{g}\t" + "\t".join(f"{2 ** x:.2f}" for x in v) + "\n")
    with open(out / "hpa_samples.tsv", "w") as fh:
        fh.write("tissue\trna_sample\tprotein_sample\n")
        for s in samp:
            fh.write(f"{s[0]}\t{s[1]}\t{s[2]}\n")

    # patients
    n = 60
    drop = {}
    for k, g in enumerate(mito[:30]):
        drop[g] = k * 2 % n
    with open(out / "patients_rna.tsv", "w") as fr, \
            open(out / "patients_protein.tsv", "w") as fp, \
            open(out / "patients_samples.tsv", "w") as fs:
        fr.write("gene_id\t" + "\t".join(f"FR{k}" for k in range(n)) + "\n")
        fp.write("gene\t" + "\t".join(f"FP{k}" for k in range(n)) + "\n")
        fs.write("rna_sample\tprotein_sample\tbatch\n")
        for k in range(n):
            fs.write(f"FR{k}\tFP{k}\tTMT{k % 3}\n")
        shift = {g: rng.normal(0, 1.5, 3) for g in genes}
        depth = rng.uniform(0.7, 1.4, n)
        for g in genes + mt:
            lr = base.get(g, 14) + rng.normal(0, 0.25, n)
            lp = base.get(g, 14) + 0.8 * (lr - base.get(g, 14)) \
                + rng.normal(0, 0.2, n)
            if g in shift:
                lp += shift[g][np.arange(n) % 3]
            if g in drop:
                k = drop[g]
                lr[k] -= 3
                if cls[g] == "transcriptional":
                    lp[k] -= 3
            fr.write(f"{ens[g]}\t" + "\t".join(
                str(int(round(2 ** x * d))) for x, d in zip(lr, depth))
                + "\n")
            if g in mt:
                continue
            vals = [f"{2 ** x:.1f}" for x in lp]
            if g == "MITOG2":
                vals[drop[g]] = "NA"
            fp.write(f"{g}\t" + "\t".join(vals) + "\n")


if __name__ == "__main__":
    main()
