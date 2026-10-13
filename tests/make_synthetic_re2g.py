#!/usr/bin/env python3
"""
make_synthetic_re2g.py - synthetic inputs for re2g_panel.py and
peaks_vs_elements.py

Purpose
    Build a small, fully known data set: a GTF, an ENCODE-like cCRE
    registry, ENCODE-rE2G link files in the column format re2g_panel.py
    parses, target and control gene lists, an alias table and two loop
    tables. Target ("mito") genes carry --n-core elements linked in every
    muscle biosample; every gene carries one element linked in a single
    biosample. Rows that the pipeline must discard are planted as well.

Inputs
    none; every quantity is an option.

Outputs (in --out)
    links/<accession>.bed.gz   one link file per biosample; header
                               chr, start, end, name, class, TargetGene,
                               TargetGeneTSS, TargetGeneEnsemblID, CellType,
                               isSelfPromoter, eight '<name>.Feature'
                               columns of the dhs_megamap model, Score (last)
    ccre.bed.gz                registry: 300 bp elements every 5 kb within
                               +/- 1 Mb of every TSS
    gencode.gtf.gz             gene, MANE Select transcript and one other
                               transcript per gene, plus one chrM gene
    genes.txt, controls.txt    gene lists; genes.txt holds an outdated
                               symbol resolved by the built-in alias table
                               (TAZ), one resolved only by aliases.tsv
                               (OLDG2), one non-gene entry and one chrM gene
    aliases.tsv                extra alias table (old, new) with a header
    loops_planted.tsv          loops.tsv format; per target gene, loops on
                               two of its core elements, one loop elsewhere
                               (called in a muscle term) and one decoy
                               called only in liver
    loops_null.tsv             loops.tsv format; per target gene, three
                               loops at random distances, independent of the
                               elements
    truth.tsv                  gene, set, chrom, tss, n_core, core offsets

Method
    Genes sit 3 Mb apart on chr1, alternating strands. The registry element
    at offset k * 5 kb from the TSS has midpoint TSS + k * 5 kb, so element
    distances are exact multiples of 5 kb. Link-file elements are the
    registry elements with a few bp of jitter per biosample, which
    re2g_panel.py must undo by mapping onto the registry. Discarded rows:
    a self-promoter link, a link 5 kb from the TSS (not distal), a link
    under the score threshold, and a link beyond the 1 Mb window.

Usage
    python3 make_synthetic_re2g.py --out synth
    python3 make_synthetic_re2g.py --out synth_wide --n-mito 4 \\
        --n-control 0 --n-muscle 70 --n-other 0 --n-core 1 --core-fraction 0.5

References
    Gschwind AR et al. 2023. An encyclopedia of enhancer-gene regulatory
        interactions in the human genome. bioRxiv.
    ENCODE Project Consortium et al. 2020. Expanded encyclopaedias of DNA
        elements in the human and mouse genomes (Registry of candidate
        cis-regulatory elements). Nature.
"""

import argparse
import gzip
from pathlib import Path

import numpy as np

STEP = 5_000            # spacing of registry elements, bp
HALF = 150              # half-width of a registry element, bp
SPAN = 200              # registry elements on each side of a TSS
SPACING = 3_000_000     # distance between genes
MUSCLE_TYPES = ("psoas muscle", "gastrocnemius medialis",
                "skeletal muscle myoblast")
FEATURES = ("ABC.Score", "3DContact", "distanceToTSS", "normalizedDNase_prom",
            "numCandidateEnhGene", "numTSSEnhGene", "sumNearbyEnhancers",
            "ubiquitousExpressedGene")
LOOP_COLS = ("gene", "chrom", "vp_pos", "summit", "start", "end", "n_terms",
             "called", "loop_id")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True)
    p.add_argument("--n-mito", type=int, default=20)
    p.add_argument("--n-control", type=int, default=10)
    p.add_argument("--n-control-empty", type=int, default=2,
                   help="control genes that receive no link at all")
    p.add_argument("--n-muscle", type=int, default=6)
    p.add_argument("--n-other", type=int, default=2,
                   help="non-muscle biosamples, without the core links")
    p.add_argument("--n-core", type=int, default=3)
    p.add_argument("--core-fraction", type=float, default=1.0,
                   help="share of target genes carrying core elements")
    p.add_argument("--seed", type=int, default=11)
    a = p.parse_args()

    rng = np.random.default_rng(a.seed)
    out = Path(a.out)
    (out / "links").mkdir(parents=True, exist_ok=True)

    genes = []                       # (name, set, tss, strand, ensembl id)
    for i in range(a.n_mito):
        name = "TAFAZZIN" if i == 3 else f"MITOG{i + 1}"
        genes.append((name, "mito"))
    for i in range(a.n_control):
        genes.append((f"CTRLG{i + 1}", "control"))
    genes = [(g, s, SPACING * (i + 1), "+" if i % 2 == 0 else "-",
              f"ENSG{i + 1:011d}")
             for i, (g, s) in enumerate(genes)]

    n_core_genes = int(round(a.core_fraction * a.n_mito))
    core = {}
    for i, (g, s, _t, _st, _e) in enumerate(genes):
        if s == "mito" and i < n_core_genes:
            core[g] = sorted(rng.choice(np.arange(80, SPAN), a.n_core,
                                        replace=False).tolist())
        else:
            core[g] = []
    empty = {g for g, s, *_ in genes if s == "control"}
    empty = set(sorted(empty, key=lambda x: int(x[5:]))[-a.n_control_empty:]
                if a.n_control_empty else [])

    # GTF: a gene record that does not start at the TSS, so the MANE Select
    # transcript must be preferred
    with gzip.open(out / "gencode.gtf.gz", "wt") as fh:
        fh.write("##description: synthetic\n")
        for g, _s, tss, strand, ens in genes:
            attr = (f'gene_id "{ens}.1"; gene_type "protein_coding"; '
                    f'gene_name "{g}";')
            if strand == "+":
                gs, ge, ts, te = tss - 4999, tss + 20000, tss + 1, tss + 20000
                os_, oe = tss + 2001, tss + 20000
            else:
                gs, ge, ts, te = tss - 20000, tss + 5000, tss - 20000, tss
                os_, oe = tss - 20000, tss - 2000
            fh.write(f"chr1\tSYN\tgene\t{gs}\t{ge}\t.\t{strand}\t.\t{attr}\n")
            fh.write(f"chr1\tSYN\ttranscript\t{os_}\t{oe}\t.\t{strand}\t.\t"
                     f'{attr} transcript_id "{ens}T2"; tag "basic";\n')
            fh.write(f"chr1\tSYN\ttranscript\t{ts}\t{te}\t.\t{strand}\t.\t"
                     f'{attr} transcript_id "{ens}T1"; tag "basic"; '
                     'tag "MANE_Select";\n')
        fh.write('chrM\tSYN\tgene\t3307\t4262\t.\t+\t.\tgene_id '
                 '"ENSG99999999999.1"; gene_type "protein_coding"; '
                 'gene_name "MT-ND1";\n')

    with gzip.open(out / "ccre.bed.gz", "wt") as fh:
        for _g, _s, tss, _st, _e in genes:
            for k in range(-SPAN, SPAN + 1):
                m = tss + k * STEP
                fh.write(f"chr1\t{m - HALF}\t{m + HALF}\tEH{tss}_{k}\n")

    head = (["chr", "start", "end", "name", "class", "TargetGene",
             "TargetGeneTSS", "TargetGeneEnsemblID", "CellType",
             "isSelfPromoter"] + [f"{f}.Feature" for f in FEATURES]
            + ["Score"])
    n_bio = a.n_muscle + a.n_other
    for b in range(n_bio):
        muscle = b < a.n_muscle
        ctype = MUSCLE_TYPES[b % len(MUSCLE_TYPES)] if muscle else "liver"
        rows = []
        for g, _s, tss, _st, ens in genes:
            if g in empty:
                continue
            # the link files may carry an older symbol; the Ensembl id finds it
            label = "OLDMITOG5" if g == "MITOG5" else g
            links = [(0, "TRUE", 0.9),           # self promoter
                     (1, "FALSE", 0.8),          # 5 kb, not distal
                     (-2, "FALSE", 0.15),        # under the threshold
                     (240, "FALSE", 0.9),        # beyond the window
                     (-(3 + b), "FALSE", 0.4)]   # this biosample only
            if muscle:
                links += [(k, "FALSE", 0.5) for k in core[g]]
            for k, selfp, score in links:
                m = tss + k * STEP
                j1, j2 = rng.integers(0, 60, 2)
                s, e = m - HALF - int(j1), m + HALF + int(j2)
                feats = [f"{rng.random():.4f}" for _ in FEATURES]
                feats[2] = f"{abs((s + e) / 2 - tss):.1f}"
                rows.append(["chr1", s, e, f"chr1:{s}-{e}", "intergenic",
                             label, tss, ens, ctype, selfp] + feats
                            + [f"{score:.6f}"])
        with gzip.open(out / "links" / f"ENCFF{b:04d}SYN.bed.gz", "wt") as fh:
            fh.write("\t".join(head) + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")

    def gene_file(path, names):
        Path(path).write_text("\n".join(names) + "\n")

    listed = []
    for g, s, *_ in genes:
        if s != "mito":
            continue
        listed.append({"TAFAZZIN": "TAZ", "MITOG2": "OLDG2"}.get(g, g))
    gene_file(out / "genes.txt", listed + ["FXN_GAA", "MT-ND1"])
    gene_file(out / "controls.txt", [g for g, s, *_ in genes if s == "control"])
    (out / "aliases.tsv").write_text("old\tnew\nOLDG2\tMITOG2\n")

    def loop(g, tss, summit, called):
        return [g, "chr1", tss, summit, summit - 5000, summit + 5000,
                called.count("|") + 1, called, f"{g}:{summit}"]

    planted, null = [], []
    for g, s, tss, _st, _e in genes:
        if s != "mito":
            continue
        for k in core[g][:2]:
            planted.append(loop(g, tss, tss + k * STEP + 3000,
                                "psoas muscle|liver"))
        planted.append(loop(g, tss, tss + int(rng.uniform(400e3, 1e6)),
                            "gastrocnemius medialis"))
        planted.append(loop(g, tss, tss + int(rng.uniform(400e3, 1e6)),
                            "liver"))
        for _ in range(3):
            null.append(loop(g, tss, tss + int(rng.uniform(400e3, 1e6)),
                             "psoas muscle"))
    for name, rows in (("loops_planted.tsv", planted),
                       ("loops_null.tsv", null)):
        with open(out / name, "w") as fh:
            fh.write("\t".join(LOOP_COLS) + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")

    with open(out / "truth.tsv", "w") as fh:
        fh.write("gene\tset\tchrom\ttss\tn_core\tcore_offsets\n")
        for g, s, tss, _st, _e in genes:
            fh.write(f"{g}\t{s}\tchr1\t{tss}\t{len(core[g])}\t"
                     + ",".join(str(k * STEP) for k in core[g]) + "\n")


if __name__ == "__main__":
    main()
