#!/usr/bin/env python3
"""
make_synthetic_gtex.py - small GTEx-like inputs for gtex_eqtl.py and
make_controls.py.

Purpose
-------
Writes a GENCODE-like GTF, a GTEx-like median TPM GCT, a GTEx-like tar of
conditionally independent cis-eQTL tables, the same eQTL table in a
per-tissue file with different column names, and a gene list, with planted
signal that a correct analysis must recover.

Inputs
------
  outdir   positional; created if absent.

Outputs (in outdir)
-------------------
  gencode.gtf.gz              "gene" lines with gene_id, gene_name, gene_type
  gtex_median_tpm.gct.gz      #1.2 header, dimensions line, Name, Description,
                              one column per tissue
  GTEx_eQTL_independent.tar   GTEx_Analysis_v8_eQTL_independent/
                              <tissue>.v8.independent_eqtls.txt.gz, columns
                              gene_id, variant_id, tss_distance, maf,
                              pval_nominal, slope, slope_se, rank
  muscle_eqtl_renamed.txt.gz  the Muscle_Skeletal table with columns
                              phenotype_id, start_distance, slope
  genes.txt                   gene-set symbols, plus one symbol absent from
                              the annotation

Method
------
Four kinds of protein-coding gene are placed on two chromosomes, every TSS at
least 300 kb from any gene-set TSS except one decoy:
  MITO*   gene set: expressed near 100 TPM in every tissue
  BROAD*  expressed near 100 TPM in every tissue (correct matches)
  SPEC*   alternating 1000 and 10 TPM across tissues: similar median, high
          cross-tissue variability (wrong matches for make_controls.py)
  LOW*    near 1 TPM everywhere (wrong matches for both scripts)
  NEAR1   features identical to MITO1, TSS 50 kb from it (removed by
          --exclusion)
  PAD*    absent from the TPM table; they flank each chromosome so that
          every expressed gene has the same local gene density
Non-coding genes and a "_PAR_Y" copy carrying a gene-set symbol are added and
must be ignored. In Muscle_Skeletal, MITO* signals have |slope| in
[0.6, 1.0] and |distance| in [50, 500] kb, every other gene |slope| in
[0.05, 0.3] and |distance| below 19 kb. In Heart_Left_Ventricle all genes
share one distribution.

Usage
-----
    python3 make_synthetic_gtex.py synth

References
----------
  GTEx Consortium 2020. The GTEx Consortium atlas of genetic regulatory
    effects across human tissues. Science. (file formats imitated)
"""

from __future__ import annotations

import gzip
import io
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

N_MITO = 20
N_BROAD = 60
N_SPEC = 60
N_LOW = 60
SIGNALS_PER_GENE = 4
SPACING = 400_000
# genes without expression data flanking each chromosome, so that every
# expressed gene has the same local density
N_PAD = 3
TISSUES = ["Adipose - Subcutaneous", "Muscle - Skeletal",
           "Heart - Left Ventricle", "Liver", "Lung", "Brain - Cortex",
           "Skin - Sun Exposed (Lower leg)", "Testis", "Thyroid",
           "Whole Blood"]
EQTL_TISSUES = {"Muscle_Skeletal": "planted", "Heart_Left_Ventricle": "null"}


def main() -> None:
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(11)

    genes = []
    kinds = ([("MITO", N_MITO), ("BROAD", N_BROAD), ("SPEC", N_SPEC),
              ("LOW", N_LOW)])
    names = [f"{k}{i + 1}" for k, n in kinds for i in range(n)]
    rng.shuffle(names)
    # irregular spacing gives every kind of gene the same spread of density
    next_tss = {"chr1": SPACING * (N_PAD + 1), "chr2": SPACING * (N_PAD + 1)}
    for idx, name in enumerate(names):
        chrom = "chr1" if idx % 2 == 0 else "chr2"
        tss = next_tss[chrom]
        next_tss[chrom] += int(rng.integers(SPACING // 2, 3 * SPACING // 2))
        length = int(rng.integers(28_000, 32_000))
        strand = "+" if rng.random() < 0.5 else "-"
        genes.append({"name": name, "chrom": chrom, "tss": tss,
                      "length": length, "strand": strand})
    mito1 = next(g for g in genes if g["name"] == "MITO1")
    genes.append({"name": "NEAR1", "chrom": mito1["chrom"],
                  "tss": mito1["tss"] + 50_000, "length": mito1["length"],
                  "strand": mito1["strand"]})
    pads = []
    for chrom in ("chr1", "chr2"):
        last = max(g["tss"] for g in genes if g["chrom"] == chrom)
        for k in range(1, N_PAD + 1):
            for tss in (SPACING * k, last + SPACING * k):
                pads.append({"name": f"PAD{len(pads) + 1}", "chrom": chrom,
                             "tss": tss, "length": 30_000, "strand": "+"})
    for i, g in enumerate(genes + pads):
        g["id"] = f"ENSG{i + 1:011d}"

    # GTF
    lines = ["##description: synthetic annotation"]
    for g in genes + pads:
        if g["strand"] == "+":
            start, end = g["tss"], g["tss"] + g["length"]
        else:
            start, end = g["tss"] - g["length"], g["tss"]
        lines.append(gtf_line(g["chrom"], start + 1, end, g["strand"],
                              g["id"] + ".7", g["name"], "protein_coding"))
        lines.append(gtf_line(g["chrom"], start + 1, start + 500, g["strand"],
                              g["id"] + ".7", g["name"], "protein_coding",
                              feature="exon"))
    for i in range(5):
        lines.append(gtf_line("chr1", 10_000 + i * 1000, 10_500 + i * 1000,
                              "+", f"ENSG9{i:010d}.1", f"LNC{i + 1}", "lncRNA"))
    lines.append(gtf_line("chrY_PAR_Y", 1000, 30_000, "+",
                          f"ENSG{1:011d}.7_PAR_Y", "MITO1", "protein_coding"))
    with gzip.open(out / "gencode.gtf.gz", "wt") as fh:
        fh.write("\n".join(lines) + "\n")

    # TPM
    rows = []
    for g in genes:
        kind = g["name"].rstrip("0123456789")
        if kind in ("MITO", "BROAD"):
            v = 100 * np.exp(rng.normal(0, 0.1, len(TISSUES)))
        elif kind == "SPEC":
            v = np.where(np.arange(len(TISSUES)) % 2 == 1, 1000.0, 10.0)
            v = v * np.exp(rng.normal(0, 0.1, len(TISSUES)))
        elif kind == "LOW":
            v = 1 * np.exp(rng.normal(0, 0.1, len(TISSUES)))
        else:
            v = None
        g["tpm"] = v
    near = next(g for g in genes if g["name"] == "NEAR1")
    near["tpm"] = mito1["tpm"].copy()
    for g in genes:
        rows.append([g["id"] + ".3", g["name"]] + [f"{x:.3f}" for x in g["tpm"]])
    rows.append(["ENSG90000000000.1", "LNC1"] + ["5.000"] * len(TISSUES))
    with gzip.open(out / "gtex_median_tpm.gct.gz", "wt") as fh:
        fh.write("#1.2\n")
        fh.write(f"{len(rows)}\t{len(TISSUES)}\n")
        fh.write("\t".join(["Name", "Description"] + TISSUES) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    # eQTL
    tar_path = out / "GTEx_eQTL_independent.tar"
    with tarfile.open(tar_path, "w") as tf:
        for tissue, mode in EQTL_TISSUES.items():
            table = eqtl_table(genes, mode, rng)
            if tissue == "Muscle_Skeletal":
                renamed = table.rename(columns={"gene_id": "phenotype_id",
                                                "tss_distance": "start_distance"})
                renamed[["phenotype_id", "variant_id", "start_distance",
                         "slope"]].to_csv(out / "muscle_eqtl_renamed.txt.gz",
                                          sep="\t", index=False)
            buf = io.BytesIO()
            with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
                gz.write(table.to_csv(sep="\t", index=False).encode())
            info = tarfile.TarInfo(
                f"GTEx_Analysis_v8_eQTL_independent/{tissue}"
                ".v8.independent_eqtls.txt.gz")
            info.size = buf.tell()
            buf.seek(0)
            tf.addfile(info, buf)

    gene_set = [f"MITO{i + 1}" for i in range(N_MITO)] + ["NOTAGENE"]
    (out / "genes.txt").write_text("\n".join(gene_set) + "\n")
    print(f"synthetic GTEx inputs written to {out}")


def gtf_line(chrom, start, end, strand, gid, name, gtype, feature="gene"):
    attr = (f'gene_id "{gid}"; gene_type "{gtype}"; gene_name "{name}"; '
            f'level 2;')
    return "\t".join([chrom, "SYNTH", feature, str(start), str(end), ".",
                      strand, ".", attr])


def eqtl_table(genes, mode, rng) -> pd.DataFrame:
    rows = []
    for g in genes:
        is_mito = g["name"].startswith("MITO")
        for rank in range(1, SIGNALS_PER_GENE + 1):
            if mode == "planted" and is_mito:
                slope = rng.uniform(0.6, 1.0)
                dist = int(rng.integers(50_000, 500_000))
            elif mode == "planted":
                slope = rng.uniform(0.05, 0.3)
                dist = int(rng.integers(0, 19_000))
            else:
                slope = rng.uniform(0.05, 1.0)
                dist = int(rng.integers(0, 500_000))
            sign_s = 1 if rng.random() < 0.5 else -1
            sign_d = 1 if rng.random() < 0.5 else -1
            pos = g["tss"] + sign_d * dist
            rows.append({"gene_id": g["id"] + ".3",
                         "variant_id": f"{g['chrom']}_{pos}_A_G_b38",
                         "tss_distance": sign_d * dist,
                         "maf": round(rng.uniform(0.05, 0.5), 3),
                         "pval_nominal": 1e-8, "slope": sign_s * slope,
                         "slope_se": 0.05, "rank": rank})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
