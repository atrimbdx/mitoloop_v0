#!/usr/bin/env python3
"""
make_synthetic_regulome.py - synthetic inputs for regulome.py

Purpose
    Build MitoCarta-like sheets, a small Gene Ontology and its annotations,
    a GTF, GTEx-like expression, protein, eQTL and Vg tables, and per-gene
    feature tables, with a known pathway per gene and a known degree of
    transcriptional control.

Inputs
    none; every quantity is an option.

Outputs (in --out)
    mitocarta_A.tsv, mitocarta_C.tsv   sheets A and C of MitoCarta3.0: 50 of
                                       the 60 target genes and 20 MitoCarta
                                       genes outside the set; pathways given
                                       as full hierarchies or as leaf names
    go.obo, goa.gaf.gz                 ontology and annotations
    genes.txt                          the 60 target genes (MITOG3 under its
                                       former symbol OLDMITOG3)
    gencode.gtf.gz                     gene records of 400 genes
    tpm.gct, protein.tsv, tissue_map.tsv, eqtl.txt.gz, vg.tsv
    promoter_features.tsv, methylation_by_group.tsv, re2g_per_gene.tsv,
    loops.tsv, assoc_loops.tsv         feature tables of the target genes
    truth.tsv                          gene, section, z (transcriptional
                                       control, 0 to 1)

Method
    Target gene i belongs to pathway i mod 4 (mtDNA replication, CI
    subunits, TIM23 presequence pathway, Fusion); MITOG7 also to mtDNA
    replication. Every gene with a pathway carries the GO term of that
    pathway, the term 'mitochondrion' and two random terms; MITOG50 has no
    MitoPathway in sheet A and MITOG51-60 are absent from MitoCarta, so
    their pathway must come from GO. A NOT annotation and an ND annotation
    to a wrong pathway must be ignored.

    Transcriptional control z is 1 for MITOG1-20, 0.5 for MITOG21-40 and 0
    for MITOG41-60 (uniform for other genes): log2 TPM across tissues varies
    with sd proportional to z, protein follows mRNA with weight z, eQTL
    signals and Vg grow with z. Mean expression is independent of z. The
    number of rE2G elements and of loops grows with z; promoter and
    methylation features do not depend on it.

Usage
    python3 make_synthetic_regulome.py --out synth
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np

PATHS = [
    ("Mitochondrial central dogma > mtDNA maintenance", "mtDNA replication"),
    ("OXPHOS > Complex I", "CI subunits"),
    ("Protein import, sorting and homeostasis > Protein import and sorting",
     "TIM23 presequence pathway"),
    ("Mitochondrial dynamics and surveillance", "Fusion"),
]
SECTIONS = [p[0].split(" > ")[0] for p in PATHS]


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=3)
    a = p.parse_args()
    rng = np.random.default_rng(a.seed)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    mito = [f"MITOG{i}" for i in range(1, 61)]
    mconly = [f"MCONLY{i}" for i in range(1, 21)]
    others = [f"OTHERG{i}" for i in range(1, 321)]
    universe = mito + mconly + others
    pw = {g: [int(g[5:]) % 4] for g in mito}
    pw["MITOG7"].append(0)
    for i, g in enumerate(mconly):
        pw[g] = [i % 4]

    # ------------------------------------------------------------ MitoCarta
    with open(out / "mitocarta_A.tsv", "w") as fh:
        fh.write("Symbol\tSynonyms\tMitoCarta3.0_SubMitoLocalization\t"
                 "MitoCarta3.0_MitoPathways\n")
        for g in mito[:50] + mconly:
            k = int(g.lstrip("MITOGCNLY"))
            if g == "MITOG50":
                paths = "0"
            elif k % 2:
                paths = " | ".join(f"{PATHS[j][0]} > {PATHS[j][1]}"
                                   for j in pw[g])
            else:
                paths = " | ".join(PATHS[j][1] for j in pw[g])
            syn = "OLDMITOG3|M3" if g == "MITOG3" else "-"
            fh.write(f"{g}\t{syn}\tMatrix\t{paths}\n")
    with open(out / "mitocarta_C.tsv", "w") as fh:
        fh.write("MitoPathways Hierarchy\tMitoPathway\tGenes\n")
        for j, (h, leaf) in enumerate(PATHS):
            genes = [g for g in mito[:49] + mconly if j in pw[g]]
            fh.write(f"{h}\t{leaf}\t{', '.join(genes)}\n")
    (out / "genes.txt").write_text(
        "".join(f"{'OLDMITOG3' if g == 'MITOG3' else g}\n" for g in mito)
        + "FXN_GAA\n")

    # ------------------------------------------------------------------- GO
    sec_t = [f"GO:{100000 + j:07d}" for j in range(4)]
    leaf_t = [f"GO:{110000 + j:07d}" for j in range(4)]
    noise = [f"GO:{200000 + j:07d}" for j in range(60)]
    with open(out / "go.obo", "w") as fh:
        fh.write("format-version: 1.2\n\n")

        def term(t, ns, parents=(), part=(), obs=False):
            fh.write(f"[Term]\nid: {t}\nname: {t}\nnamespace: {ns}\n")
            for x in parents:
                fh.write(f"is_a: {x} ! parent\n")
            for x in part:
                fh.write(f"relationship: part_of {x} ! whole\n")
            if obs:
                fh.write("is_obsolete: true\n")
            fh.write("\n")
        term("GO:0008150", "biological_process")
        term("GO:0005575", "cellular_component")
        term("GO:0005739", "cellular_component", ["GO:0005575"])
        for j in range(4):
            term(sec_t[j], "biological_process", ["GO:0008150"])
            term(leaf_t[j], "biological_process", [], [sec_t[j]])
        for t in noise:
            term(t, "biological_process", ["GO:0008150"])
        term("GO:9999999", "biological_process", ["GO:0008150"], obs=True)
        fh.write("[Typedef]\nid: part_of\nname: part of\n")

    def gaf(fh, g, t, ev="IDA", qual="", aspect="P"):
        fh.write(f"UniProtKB\tP{abs(hash(g)) % 99999:05d}\t{g}\t{qual}\t{t}\t"
                 f"PMID:1\t{ev}\t\t{aspect}\t{g} protein\t{g}|syn_{g}\t"
                 "protein\ttaxon:9606\t20240101\tSYN\t\t\n")
    with gzip.open(out / "goa.gaf.gz", "wt") as fh:
        fh.write("!gaf-version: 2.2\n")
        for g in universe:
            for t in rng.choice(noise, 2 if g in pw else 3, replace=False):
                gaf(fh, g, t)
            if g in pw:
                for j in pw[g]:
                    gaf(fh, g, leaf_t[j])
                gaf(fh, g, "GO:0005739", aspect="C")
        gaf(fh, "MITOG52", leaf_t[(52 + 1) % 4], qual="NOT")
        gaf(fh, "MITOG53", leaf_t[(53 + 1) % 4], ev="ND")
        gaf(fh, "MITOG54", "GO:9999999")

    # ------------------------------------------------------------ expression
    z = {}
    for i, g in enumerate(mito):
        z[g] = 1.0 if i < 20 else (0.5 if i < 40 else 0.0)
    for g in mconly + others:
        z[g] = float(rng.uniform())
    base = {g: float(rng.uniform(2, 8)) for g in universe}
    tissues = [f"Tissue {k:02d}" for k in range(1, 19)] + \
        ["Muscle - Skeletal", "Cells - Cultured fibroblasts"]
    ens = {g: f"ENSG{n:011d}" for n, g in enumerate(universe)}
    with gzip.open(out / "gencode.gtf.gz", "wt") as fh:
        for n, g in enumerate(universe):
            s = 1_000_000 + n * 100_000
            fh.write(f"chr1\tSYN\tgene\t{s}\t{s + 10000}\t.\t+\t.\t"
                     f'gene_id "{ens[g]}.3"; gene_type "protein_coding"; '
                     f'gene_name "{g}";\n')
        fh.write('chr1\tSYN\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG99"; '
                 'gene_type "lncRNA"; gene_name "LNC1";\n')
    logs = {}
    with open(out / "tpm.gct", "w") as fh:
        fh.write(f"#1.2\n{len(universe)}\t{len(tissues)}\n")
        fh.write("Name\tDescription\t" + "\t".join(tissues) + "\n")
        for g in universe:
            prof = rng.normal(0, 1, len(tissues))
            lg = base[g] + 2.5 * z[g] * prof + rng.normal(0, 0.05,
                                                          len(tissues))
            logs[g] = lg
            tpm = np.maximum(2 ** lg - 1, 0)
            fh.write(f"{ens[g]}.3\t{g}\t"
                     + "\t".join(f"{x:.3f}" for x in tpm) + "\n")
    ptiss = tissues[:15]
    with open(out / "tissue_map.tsv", "w") as fh:
        fh.write("protein_tissue\tgtex_tissue\n")
        for k, t in enumerate(ptiss):
            fh.write(f"prot_{k}\t{t}\n")
    with open(out / "protein.tsv", "w") as fh:
        fh.write("gene\t" + "\t".join(f"prot_{k}" for k in range(15)) + "\n")
        for g in universe:
            lg = logs[g][:15]
            lg = (lg - lg.mean()) / (lg.std() + 1e-9)
            pr = z[g] * lg + (1 - z[g]) * rng.normal(0, 1, 15) \
                + rng.normal(0, 0.1, 15)
            key = ens[g] if g.startswith("OTHERG") else g
            fh.write(f"{key}\t" + "\t".join(f"{x:.3f}" for x in pr) + "\n")
    with gzip.open(out / "eqtl.txt.gz", "wt") as fh:
        fh.write("variant_id\tgene_id\ttss_distance\tslope\trank\n")
        for g in universe:
            k = int(round(3 * z[g] + rng.normal(0, 0.3)))
            for r in range(max(0, k)):
                fh.write(f"chr1_{r}_A_G_b38\t{ens[g]}.3\t{1000 * r}\t"
                         f"{(0.2 + 0.6 * z[g]) * (-1) ** r:.3f}\t{r + 1}\n")
    with open(out / "vg.tsv", "w") as fh:
        fh.write("gene\tVg\n")
        for g in universe:
            fh.write(f"{ens[g]}.3\t{10 ** (-2 + 2 * z[g] + rng.normal(0, 0.2)):.5f}\n")

    # -------------------------------------------------------------- features
    with open(out / "re2g_per_gene.tsv", "w") as fh:
        fh.write("gene\tset\tchrom\ttss\tn_ccre_window\tn_union\t"
                 "n_intersection\tmax_prevalence\n")
        for g in mito:
            nu = max(0, int(round(3 + 8 * z[g] + rng.normal(0, 1))))
            fh.write(f"{g}\tmito\tchr1\t1\t{rng.integers(1500, 2500)}\t{nu}\t"
                     f"{rng.integers(0, 3)}\t{rng.integers(1, 64)}\n")
    with open(out / "loops.tsv", "w") as fh, \
            open(out / "assoc_loops.tsv", "w") as fa:
        fh.write("gene\tchrom\tvp_pos\tsummit\tstart\tend\tn_terms\tcalled\t"
                 "loop_id\n")
        fa.write("gene\tloop_id\tlog2fc_per_doubling\tt\n")
        for g in mito:
            for r in range(int(round(4 * z[g] + rng.uniform(0, 0.9)))):
                fh.write(f"{g}\tchr1\t1\t{r}\t{r}\t{r}\t"
                         f"{rng.integers(3, 7)}\tpsoas muscle\t{g}:{r}\n")
                fa.write(f"{g}\t{g}:{r}\t{rng.normal(0.02, 0.01):.4f}\t1\n")
    with open(out / "promoter_features.tsv", "w") as fh:
        fh.write("gene\tgene_id\tset\tchrom\ttss\tstrand\tprom_start\t"
                 "prom_end\tgc\tcpg_oe\tpromoter_class\tmotif_NRF1\t"
                 "bidirectional\n")
        for g in mito:
            fh.write(f"{g}\t{ens[g]}\tmito\tchr1\t1\t+\t0\t1\t"
                     f"{rng.uniform(0.4, 0.7):.3f}\t{rng.uniform(0.2, 1.1):.3f}"
                     f"\t{rng.choice(['HCP', 'LCP'])}\t{rng.integers(0, 2)}\t"
                     f"{rng.integers(0, 2)}\n")
    with open(out / "methylation_by_group.tsv", "w") as fh:
        fh.write("gene\tset\tmuscle\tliver\n")
        for g in mito:
            fh.write(f"{g}\tmito\t{rng.uniform(0, 0.2):.3f}\t"
                     f"{rng.uniform(0, 0.2):.3f}\n")
    with open(out / "truth.tsv", "w") as fh:
        fh.write("gene\tsection\tz\n")
        for g in mito:
            fh.write(f"{g}\t{SECTIONS[pw[g][0]]}\t{z[g]}\n")


if __name__ == "__main__":
    main()
