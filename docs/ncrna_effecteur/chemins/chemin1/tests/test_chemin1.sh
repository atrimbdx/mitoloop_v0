#!/usr/bin/env bash
# Tests du chemin 1.
# 1) Logique, sur un JEU DE TEST CONSTRUIT À LA MAIN (tests/mini_paquet.py) : toujours exécuté.
# 2) Témoins du catalogue, sur un vrai paquet si PAQUET=<dossier> est défini (mesure, pas de seuil ajusté).
set -euo pipefail
ICI="$(cd "$(dirname "$0")" && pwd)"
TMP="$(mktemp -d)"
python3 "$ICI/mini_paquet.py" "$TMP/mini"
python3 "$ICI/../vm/annoter_chemin1.py" --variants "$TMP/mini/variants_test.tsv" --paquet "$TMP/mini" --out "$TMP/mini_out.tsv"
echec=0
verifie() { if grep -q -- "$2" <(sed -n "$1p" "$TMP/mini_out.tsv"); then echo "OK   ligne $1 : $2"; else echo "ECHEC ligne $1 : $2"; sed -n "$1p" "$TMP/mini_out.tsv"; echec=1; fi; }
verifie 2 "perdues 2 dont 2 mito (P=0.17"
verifie 2 "gagnées 1 dont 0 mito"
verifie 2 "|hypothese|"
verifie 3 "allèle REF discordant"
verifie 4 "|faible|LNCTEST (lncRNA) : antisens chevauchant le(s) gène(s) mito MITO1"
verifie 5 "|hypothese|RMRP : partenaires protéiques de l'ensemble mito : PNPT1"
verifie 6 "0|aucun|variant hors des éléments du chemin 1"
python3 - <<PY
import sys; sys.path.insert(0, "$ICI/..")
from commun import hypergeom_sf
assert abs(hypergeom_sf(2, 4, 2, 2) - 1/6) < 1e-12
assert abs(hypergeom_sf(1, 10, 3, 2) - (1 - 21/45)) < 1e-12
print("OK   hypergéométrique")
PY
if [ -n "${PAQUET:-}" ]; then
  python3 "$ICI/construire_temoins.py" "$PAQUET" "$TMP/temoins.tsv"
  python3 "$ICI/../vm/annoter_chemin1.py" --variants "$TMP/temoins.tsv" --paquet "$PAQUET" --out "${SORTIE_TEMOINS:-$TMP/temoins_annotes.tsv}"
  echo "Témoins annotés : ${SORTIE_TEMOINS:-$TMP/temoins_annotes.tsv}"
fi
exit $echec
