#!/usr/bin/env python3
"""
make_synthetic_mitoshow.py - synthetic inputs for mitoshow_data.py

Purpose
    Build MitoCarta-like sheets with the sub-mitochondrial localisations in
    the proportions of MitoCarta3.0, PanelApp genes outside MitoCarta, a
    Human Protein Atlas subcellular table, GO annotations, STRING links and
    a curated location table, at the size of the real data so that the
    viewer can also be tried on them.

Inputs
    none

Outputs (in --out)
    mitocarta_A.tsv, mitocarta_C.tsv, genes.txt, hpa.tsv, goa.gaf,
    extra_locations.tsv, string_links.txt.gz (detailed channels),
    string_info.txt.gz, structures.tsv, complex.cif (a three-chain mmCIF),
    hgnc.tsv (symbol, uniprot_ids), truth.tsv (gene, expected locations)

Method
    MitoCarta genes MC1..MC1136: Matrix 525, MIM 340, IMS 50, MOM 110,
    Membrane 41, unknown 70; matrix genes of the mtDNA pathways go to the
    nucleoid. MC5 is also nuclear and MC6 cytosolic in HPA (Supported);
    MC7 is nuclear in HPA with Uncertain reliability (ignored). PanelApp
    genes PA1..PA110: PA1-60 placed by HPA, PA61-80 by GO only, PA81-110
    nowhere. MC8 is given the ER by the curated table, MC9 a nested OXPHOS
    place. MT-TL1 and MT-TK (PanelApp) have no location. STRING links join
    genes of the same pathway; one link per pair is under --min-string. The
    mmCIF places MC530, MC531 and MC532 in "Complex I": the first two in
    contact, the third 60 A away.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np

PATHWAYS = [
    ("Mitochondrial central dogma > mtDNA maintenance", "mtDNA replication"),
    ("Mitochondrial central dogma > Translation", "Mitochondrial ribosome"),
    ("Mitochondrial central dogma > mtRNA metabolism", "Transcription"),
    ("OXPHOS > Complex I", "CI subunits"),
    ("OXPHOS > Complex IV", "CIV assembly factors"),
    ("Metabolism > Carbohydrate metabolism", "TCA cycle"),
    ("Metabolism > Lipid metabolism", "Fatty acid oxidation"),
    ("Protein import, sorting and homeostasis > Protein import and sorting",
     "TIM23 presequence pathway"),
    ("Mitochondrial dynamics and surveillance", "Fusion"),
    ("Small molecule transport", "SLC25A family"),
    ("Signaling", "Apoptosis"),
]
LOCS = [("Matrix", 525), ("MIM", 340), ("IMS", 50), ("MOM", 110),
        ("Membrane", 41), ("unknown", 70)]


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=2)
    a = p.parse_args()
    rng = np.random.default_rng(a.seed)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    locs = [loc for loc, n in LOCS for _ in range(n)]
    mc = [f"MC{i + 1}" for i in range(len(locs))]
    pw = {}
    for i, g in enumerate(mc):
        if locs[i] == "Matrix":
            pw[g] = int(rng.choice([0, 1, 2, 5, 6]))
        elif locs[i] == "MIM":
            pw[g] = int(rng.choice([3, 4, 7, 9]))
        elif locs[i] == "MOM":
            pw[g] = int(rng.choice([7, 8, 10]))
        else:
            pw[g] = int(rng.integers(0, len(PATHWAYS)))
    pw["MC1"] = 0                       # matrix, mtDNA: nucleoid
    with open(out / "mitocarta_A.tsv", "w") as fh:
        fh.write("Symbol\tSynonyms\tMitoCarta3.0_SubMitoLocalization\t"
                 "MitoCarta3.0_MitoPathways\n")
        for i, g in enumerate(mc):
            h, leaf = PATHWAYS[pw[g]]
            path = leaf if i % 2 else f"{h} > {leaf}"
            syn = "OLDMC2|X2" if g == "MC2" else "-"
            fh.write(f"{g}\t{syn}\t{locs[i]}\t{path}\n")
    with open(out / "mitocarta_C.tsv", "w") as fh:
        fh.write("MitoPathways Hierarchy\tMitoPathway\tGenes\n")
        for j, (h, leaf) in enumerate(PATHWAYS):
            fh.write(f"{h}\t{leaf}\t"
                     + ", ".join(g for g in mc if pw[g] == j) + "\n")
    pa = [f"PA{i}" for i in range(1, 111)]
    (out / "genes.txt").write_text("\n".join(pa + ["OLDMC2", "FXN_GAA",
                                                    "MT-TL1", "MT-TK"])
                                   + "\n")
    hpa_locs = ["Nucleoplasm", "Cytosol", "Plasma membrane",
                "Endoplasmic reticulum", "Mitochondria", "Golgi apparatus"]
    truth = {}
    with open(out / "hpa.tsv", "w") as fh:
        fh.write("Gene\tGene name\tReliability\tMain location\t"
                 "Additional location\n")
        fh.write("ENSG5\tMC5\tSupported\tMitochondria\tNucleoplasm\n")
        fh.write("ENSG6\tMC6\tApproved\tCytosol;Mitochondria\t\n")
        fh.write("ENSG7\tMC7\tUncertain\tNucleoplasm\t\n")
        for k, g in enumerate(pa[:60]):
            main = hpa_locs[k % len(hpa_locs)]
            fh.write(f"ENSGPA{k}\t{g}\tEnhanced\t{main}\t\n")
            truth[g] = {"Nucleoplasm": "nucleus", "Cytosol": "cytosol",
                        "Plasma membrane": "plasma_membrane",
                        "Endoplasmic reticulum": "er",
                        "Mitochondria": "mito_unknown",
                        "Golgi apparatus": "other"}[main]
    with open(out / "goa.gaf", "w") as fh:
        fh.write("!gaf-version: 2.2\n")
        for k, g in enumerate(pa[60:80]):
            term = ["GO:0005743", "GO:0005634", "GO:0005829",
                    "GO:0005777"][k % 4]
            truth[g] = {"GO:0005743": "imm", "GO:0005634": "nucleus",
                        "GO:0005829": "cytosol",
                        "GO:0005777": "peroxisome"}[term]
            fh.write(f"UniProtKB\tQ{k}\t{g}\tlocated_in\t{term}\tPMID:1\tIDA"
                     f"\t\tC\t{g}\t\tprotein\ttaxon:9606\t20240101\tGO\t\t\n")
        fh.write("UniProtKB\tQ99\tPA81\tNOT|located_in\tGO:0005634\tPMID:1\t"
                 "IDA\t\tC\tPA81\t\tprotein\ttaxon:9606\t20240101\tGO\t\t\n")
    for g in pa[80:]:
        truth[g] = "unlocalized"
    (out / "extra_locations.tsv").write_text(
        "gene\tcompartment\tpath\tsource\nMC8\ter\t\tcurated\n"
        "MC9\toxphos\tSupercomplexes > Respirasome assembly\tpaper\n")
    # a three-chain structure, entities mapped through UniProt accessions
    with open(out / "hgnc.tsv", "w") as fh:
        fh.write("hgnc_id\tsymbol\tuniprot_ids\n")
        for k, g in enumerate(("MC530", "MC531", "MC532")):
            fh.write(f"HGNC:{k}\t{g}\tP0000{k}\n")
    with open(out / "complex.cif", "w") as fh:
        fh.write("data_SYN\n#\nloop_\n_struct_ref.id\n_struct_ref.db_name\n"
                 "_struct_ref.db_code\n_struct_ref.pdbx_db_accession\n"
                 "_struct_ref.entity_id\n")
        for k in range(3):
            fh.write(f"{k + 1} UNP X{k}_HUMAN P0000{k} {k + 1}\n")
        fh.write("#\n_entity.id 1\n_entity.pdbx_description 'a ; quoted'\n#\n")
        fh.write("loop_\n_atom_site.group_PDB\n_atom_site.id\n"
                 "_atom_site.label_atom_id\n_atom_site.label_entity_id\n"
                 "_atom_site.auth_asym_id\n_atom_site.Cartn_x\n"
                 "_atom_site.Cartn_y\n_atom_site.Cartn_z\n")
        n = 0
        for k, (ch, x0) in enumerate((("A", 0.0), ("B", 6.0), ("C", 60.0))):
            for r in range(20):
                for atom in ("N", "CA"):
                    n += 1
                    fh.write(f"ATOM {n} {atom} {k + 1} {ch} {x0 + r * 0.3:.2f} "
                             f"{r * 1.5:.2f} {k * 2.0:.2f}\n")
    (out / "structures.tsv").write_text(
        "complex\tfile\tprojection\nComplex I\t" + str(out / "complex.cif")
        + "\tpca\n")

    names = mc + pa
    with gzip.open(out / "string_info.txt.gz", "wt") as fh:
        fh.write("#string_protein_id\tpreferred_name\tprotein_size\t"
                 "annotation\n")
        for k, g in enumerate(names):
            fh.write(f"9606.ENSP{k:011d}\t{g}\t300\t-\n")
    with gzip.open(out / "string_links.txt.gz", "wt") as fh:
        fh.write("protein1 protein2 neighborhood fusion cooccurence "
                 "coexpression experimental database textmining "
                 "combined_score\n")
        by = {}
        for g in mc:
            by.setdefault(pw[g], []).append(g)
        idx = {g: k for k, g in enumerate(names)}
        for j, gs in by.items():
            for _ in range(len(gs) * 3):
                x, y = rng.choice(gs, 2, replace=False)
                sc = int(rng.integers(450, 999))
                exp = int(sc * rng.uniform(0.2, 1))
                for u, v in ((x, y), (y, x)):
                    fh.write(f"9606.ENSP{idx[u]:011d} 9606.ENSP{idx[v]:011d} "
                             f"0 0 0 {int(sc * 0.3)} {exp} {int(sc * 0.5)} "
                             f"{int(sc * 0.2)} {sc}\n")
        fh.write(f"9606.ENSP{idx['MC1']:011d} 9606.ENSP{idx['MC3']:011d} "
                 "0 0 0 0 0 0 150 150\n")
        fh.write(f"9606.ENSP{idx['MC1']:011d} 9606.ENSP{idx['PA1']:011d} "
                 "0 0 0 0 900 0 0 900\n")
    with open(out / "truth.tsv", "w") as fh:
        fh.write("gene\tlocation\n")
        for g, c in truth.items():
            fh.write(f"{g}\t{c}\n")


if __name__ == "__main__":
    main()
