#!/usr/bin/env python3
"""Construit le TSV des témoins (chr/pos/ref/alt hg38) à partir de la nomenclature n. du catalogue,
en la convertissant sur le transcrit GENCODE v27 de l'élément (paquet : elements_ncrna.tsv.gz).

La numérotation n. publiée suit le transcrit RefSeq (NR_) ou l'épingle miRBase ; on suppose qu'elle
coïncide avec le début du transcrit GENCODE et on CONTRÔLE la base de référence attendue : toute
discordance est écrite dans la colonne controle_ref et le témoin est à vérifier, pas à interpréter.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from commun import lire_variants, ouvrir, revcomp

# gene, n (1-based sur le transcrit), ref (brin du transcrit), alt, étiquette, source
TEMOINS = [
    ("MIR204", 37, "C", "T", "temoin_publie", "catalogue : n.37C>T (graine), Conte 2015 doi:10.1073/pnas.1401464112"),
    ("MIR184", 57, "C", "T", "temoin_publie", "catalogue : c.57C>T (graine)"),
    ("RNU12", 84, "C", "T", "temoin_publie", "catalogue : n.84C>T homozygote (SCAR33)"),
    ("RNU4ATAC", 51, "G", "A", "temoin_publie", "catalogue : n.51G>A (boucle 5', NR_023343.1)"),
    ("RNU6ATAC", 43, "G", "A", "temoin_publie", "catalogue : n.43G>A"),
    ("RNU2-2", 4, "G", "A", "temoin_publie", "catalogue : n.4G>A (dominant)"),
    ("RNU2-2", 35, "A", "G", "temoin_publie", "catalogue : n.35A>G (dominant)"),
    ("RNU1-1", 3, "A", "G", "temoin_somatique", "catalogue : r.3A>G somatique (médulloblastome)"),
    ("RMRP", 71, "A", "G", "temoin_publie", "catalogue : n.71A>G fondateur (numérotation à vérifier : n.72 récent)"),
    ("RMRP", 72, "A", "G", "temoin_publie", "catalogue : n.72A>G (nomenclature récente du même fondateur)"),
]


def main():
    paquet, out = sys.argv[1], sys.argv[2]
    el = {}
    with ouvrir(os.path.join(paquet, "elements_ncrna.tsv.gz")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            el.setdefault(r["gene"], r)
    mir = {}
    with ouvrir(os.path.join(paquet, "mirna.tsv.gz")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mir.setdefault(r["gene"], r)
    rows = []
    for g, n, ref, alt, etiq, src in TEMOINS:
        e = el.get(g)
        if e is None:
            rows.append(["NA", "NA", ref, alt, g, etiq, src, "element_absent"])
            continue
        s, f, brin = int(e["debut0"]), int(e["fin"]), e["brin"]
        seq = e["sequence_brin_plus"]
        if brin == "+":
            pos = s + n
            gref, galt = ref, alt
        else:
            pos = f - n + 1
            gref, galt = revcomp(ref), revcomp(alt)
        i = pos - 1 - s
        obs = seq[i] if 0 <= i < len(seq) else "hors_element(longueur=%d)" % len(seq)
        rows.append([e["chrom"], pos, gref, galt, g, etiq, src, "ok" if obs == gref else "REF_DISCORDANTE(obs=%s)" % obs])
    # témoin publié avec coordonnée génomique donnée par le catalogue
    rows.append(["chr12", 120291839, "T", "TA", "RNU4-2", "temoin_publie",
                 "catalogue : n.64_65insT = chr12:120291839 T>TA (ReNU)", "coordonnee_catalogue"])
    # MIR96 : le résumé de Mencía 2009 (doi:10.1038/ng.355) ne donne pas les positions ; on construit
    # TOUTES les substitutions de la graine (pos 2-8) de miR-96-5p : JEU CONSTRUIT, pas les variants publiés.
    m = mir["MIR96"]
    arms = {x.split(":")[0]: (int(x.split(":")[1]), x.split(":")[2]) for x in m["matures_nom:offset0:seq"].split(";")}
    o, mat = arms["hsa-miR-96-5p"]
    s0, f0, brin = int(m["debut0"]), int(m["fin"]), m["brin"]
    for p in range(2, 9):
        off = o + p - 1
        b = mat[p - 1]
        for a in "ACGT":
            if a == b:
                continue
            pos = s0 + off + 1 if brin == "+" else f0 - off
            gr, ga = (b, a) if brin == "+" else (revcomp(b), revcomp(a))
            rows.append([m["chrom"], pos, gr, ga, "MIR96", "construit_graine",
                         "miR-96-5p graine pos %d %s>%s (construit)" % (p, b, a), "construit"])
    with open(out, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["chr", "pos", "ref", "alt", "gene", "etiquette", "source", "controle_ref"])
        w.writerows(rows)
    print("%d témoins écrits dans %s" % (len(rows), out))


if __name__ == "__main__":
    main()
