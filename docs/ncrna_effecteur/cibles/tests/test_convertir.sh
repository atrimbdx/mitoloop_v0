#!/bin/sh
# Test de convertir_cibles.py sur des extraits SYNTHETIQUES étiquetés (bases inaccessibles depuis la session).
# Usage : sh test_convertir.sh HAIRPIN_FA HGNC_TXT
set -e
D=$(dirname "$0"); P="$D/../convertir_cibles.py"; O=$(mktemp -d)
python3 -I "$P" mirtarbase "$D/mirtarbase_synthetique.csv" --hairpin "$1" --hgnc "$2" --version synthetique > "$O/mtb.tsv"
grep -q "MIR210	.*ISCU	arnm	repression_graine	.	experimentale_humaine" "$O/mtb.tsv"
test "$(grep -c '^mirna_canonique	MIR181A[12]	.*COX11' "$O/mtb.tsv")" -eq 2
grep -q "MIRLET7A1	.*ATP5MC1	.*base_de_donnees" "$O/mtb.tsv"
! grep -q "Cox1" "$O/mtb.tsv"
! grep -q "TFAM" "$O/mtb.tsv"
grep -v '^#' "$D/generique_synthetique.tsv" > "$O/gen_in.tsv"
python3 -I "$P" generique "$O/gen_in.tsv" --type-id tdr_trf_tirna --col-element tdr_id --col-cible gene_cible --col-position site --type-cible arnm --nature repression_graine --preuve prediction --source "SYNTHETIQUE" --reference "synthetique" > "$O/gen.tsv"
test "$(wc -l < "$O/gen.tsv")" -eq 3
python3 -I "$P" fusion --hgnc "$2" --out "$O/f.tsv" "$O/mtb.tsv" "$O/gen.tsv"
test "$(head -1 "$O/f.tsv")" = "$(printf 'type_id\telement\telement_id\tcible\ttype_cible\tnature\tposition_cible\tpreuve\tsource\treference\tregion_guide')"
echo "OK : $(($(wc -l < "$O/f.tsv")-1)) lignes fusionnées"; rm -r "$O"
