#!/usr/bin/env bash
# locus_figure.py --ag2-dir : pistes de la campagne AlphaGenome 2 (délétion entière des éléments
# et des ancres, perte de contact), sans carte Hi-C (profil 4C virtuel d'un npz), données synthétiques
set -euo pipefail
cd "$(dirname "$0")/.."
python3 -c "import matplotlib" 2>/dev/null || { echo "SKIP test_locus_ag2 (matplotlib absent)"; exit 0; }
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
python3 - "$T" <<'PY'
import numpy as np, os, sys
T = sys.argv[1]; os.makedirs(f"{T}/ag")
off = np.arange(-500000, 500001, 10000)
prof = np.exp(-np.abs(off) / 80000) * 50
np.savez(f"{T}/v4c.npz", names=np.array(["GX"]), chrom=np.array(["chr1"]), pos=np.array([1000000]),
         profiles=np.array([prof]), offsets=off, resolution=10000, anchor_half=10000, accession="X", term="m")
w = lambda f, rows: open(f"{T}/ag/{f}", "w").write("".join("\t".join(map(str, r)) + "\n" for r in rows))
w("choix_conditions.tsv", [("type", "condition", "nom", "expression", "acces"), ("muscle", "X", "m", "rna_seq", "")])
w("convention_ag2.tsv", [("cle", "valeur"), ("alt_realigne", "non"), ("unite_arn", "ln")])
w("cibles.tsv", [("cible", "gene", "chrom", "debut", "fin", "tss", "brin", "classe", "type", "exon1"),
                 ("E1", "GX", "chr1", 1030000, 1030400, 1000000, "+", "element_hors_ancre", "muscle", "")])
w("multi_muscle.tsv", [("cible", "echelle", "gene", "type", "chrom", "debut", "fin", "taille", "variant", "sortie",
                        "gene_arn", "piste", "imputation", "raw_score", "quantile_score"),
                       ("E1", "element", "GX", "muscle", "chr1", 1030000, 1030400, 400, ".", "rna_seq", "GX", "p", "o", -0.6931, "."),
                       ("A1", "ancre_10kb", "GX", "muscle", "chr1", 1195000, 1205000, 10000, ".", "rna_seq", "GX", "p", "o", 0.0, ".")])
w("lecture_multi_muscle.tsv", [("cible", "echelle", "gene", "type", "taille", "mesure", "ref", "alt", "log2", "delta"),
                               ("A1", "ancre_10kb", "GX", "muscle", 10000, "contact_promoteur_1200kb_flanc_d", 1, 0.5, ".", -0.5),
                               ("A1", "ancre_10kb", "GX", "muscle", 10000, "contact_promoteur_1200kb_3x3", 1, "nan", ".", "nan")])
sys.path.insert(0, "scripts")
from locus_figure import ag2_regions
r = sorted(ag2_regions(f"{T}/ag", "GX", "muscle", "chr1"))
assert [(x[2], round(x[3], 3) if x[3] is not None else None, x[4]) for x in r] == \
    [("element_hors_ancre", -1.0, None), ("ancre", 0.0, -0.5)], r
PY
python3 scripts/locus_figure.py --npz "$T/v4c.npz" --gene GX --ag2-dir "$T/ag" --ag2-type muscle --out "$T/f" --dpi 60 > "$T/log" 2>&1 \
    || { cat "$T/log"; echo "FAIL test_locus_ag2"; exit 1; }
[ -s "$T/f.png" ] && grep -q "1 éléments, 1 ancres" "$T/log" || { cat "$T/log"; echo "FAIL test_locus_ag2 content"; exit 1; }
echo "PASS test_locus_ag2"
