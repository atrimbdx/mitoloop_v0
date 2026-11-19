#!/usr/bin/env python3
"""JEU DE TEST CONSTRUIT À LA MAIN (séquences et gènes fictifs, aucune donnée réelle ni patient).

Sert uniquement à vérifier la logique de annoter_chemin1.py : détection des sites 8mer / 7mer-m8 /
7mer-A1, gains et pertes, hypergéométrique, lncRNA antisens, partenaires Pol III, REF discordante.
"""
import csv
import gzip
import os
import sys

d = sys.argv[1]
os.makedirs(d, exist_ok=True)


def tsv(nom, header, rows):
    op = gzip.open if nom.endswith(".gz") else open
    with op(os.path.join(d, nom), "wt") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# miARN fictif : mature REF = U GAGCUAC CAGGUAGCUUAGCC (graine nt 2-8 = GAGCUAC)
mature = "TGAGCTACCAGGTAGCTTAGCC"
mature_alt = mature[:3] + "T" + mature[4:]          # graine pos 4 G>T : GAGCUAC -> GATCUAC
s8 = rc(mature[1:8]) + "A"                           # 8mer REF
s7a1 = rc(mature[1:7]) + "A"                         # 7mer-A1 REF
salt = rc(mature_alt[1:8]) + "A"                     # 8mer ALT
fond = "CCCCCCCCCCCCCCCCCCCC"
tsv("utr3.tsv.gz", ["gene", "transcrit", "chrom", "brin", "longueur", "blocs_0based", "sequence"], [
    ["MITO1", "T1", "chrT", "+", 48, ".", fond + s8 + fond],
    ["MITO2", "T2", "chrT", "+", 48, ".", fond + s7a1 + fond],
    ["NONM1", "T3", "chrT", "+", 40, ".", fond + fond],
    ["NONM2", "T4", "chrT", "+", 48, ".", fond + salt + fond],
])
pre = "GGGGGGGGGGGGGGGGGGGG" + mature + "AAAAAAAAAA" + rc(mature)[:20] + "GGGG"
# précurseur sur chrT:1001-..., brin +, bras mature à l'offset 20 (0-based)
tsv("mirna.tsv.gz", ["gene", "transcrit", "chrom", "brin", "debut0", "fin", "miRBase_hairpin",
                     "matures_nom:offset0:seq", "sequence_transcrit_flancs"],
    [["MIRTEST", "TX", "chrT", "+", 1000, 1000 + len(pre), "hsa-mir-test", "hsa-miR-test-5p:20:" + mature, pre]])
lnc = "ACGT" * 25
rmrp = "GATC" * 25
tsv("elements_ncrna.tsv.gz", ["gene", "classe", "sous_classe", "chrom", "brin", "debut0", "fin", "transcrits",
                              "antisens_chevauchant_codant", "sequence_brin_plus"], [
    ["LNCTEST", "lncRNA", "lncRNA", "chrT", "-", 5000, 5100, "TXL", "MITO1,NONM1", lnc],
    ["RMRP", "ARN_misc", "PolIII_RMRP", "chrT", "-", 8000, 8100, "TXR", ".", rmrp],
])
tsv("genes_mito.tsv", ["gene", "sources"], [["MITO1", "test"], ["MITO2", "test"], ["PNPT1", "test"]])
tsv("partenaires_pol3_arnt_arnr.tsv", ["arn", "partenaire", "ensemble_mito", "agit_sur_arn_nucleaire", "lien", "source"],
    [["RMRP", "PNPT1", "test", "oui (test)", "import (test)", "test"]])
# variants de test (1-based) : graine pos 4 (offset 20+3 -> pos 1024), REF discordante, lncRNA, RMRP, hors élément
tsv("variants_test.tsv", ["chr", "pos", "ref", "alt", "attendu"], [
    ["chrT", 1024, "G", "T", "miARN graine : perdues MITO1 MITO2, gagnée NONM2"],
    ["chrT", 1024, "A", "T", "REF discordante"],
    ["chrT", 5050, "G", "A", "lncRNA antisens de MITO1 -> faible"],
    ["chrT", 8050, "T", "C", "RMRP -> PNPT1 hypothese"],
    ["chrT", 20000, "A", "G", "hors élément -> aucun"],
])
print("mini paquet de test écrit dans", d)
