#!/usr/bin/env bash
# Test du chemin 2 sur des temoins CONSTRUITS A LA MAIN (variants fictifs, etiquetes) et un temoin du catalogue.
# Usage : bash tests/test_chemin2.sh   (depuis le dossier chemin2)
set -euo pipefail
ICI="$(cd "$(dirname "$0")/.." && pwd)"
SORTIE="${TMPDIR:-/tmp}/chemin2_test_$$.tsv"
python3 "$ICI/scripts/annoter_chemin2.py" --variants "$ICI/tests/temoins_chemin2.tsv" --out "$SORTIE"
python3 - "$SORTIE" <<'PY'
import sys
ok = 0; ko = 0
with open(sys.argv[1]) as fh:
    tete = fh.readline().rstrip("\n").split("\t")
    for l in fh:
        r = dict(zip(tete, l.rstrip("\n").split("\t")))
        score, niveau, just = r["NC_chemin2"].split("|", 2)
        bon = niveau == r["attendu_niveau"]
        assert niveau in ("aucun", "hypothese", "faible"), niveau
        assert "ACMG" in just or niveau == "aucun"
        if r["id"] == "T10":
            bon = bon and "REF discordant" in just
        print(("OK " if bon else "ECHEC ") + r["id"], niveau, score, "|", just[:160])
        ok += bon; ko += (not bon)
print(f"{ok} OK, {ko} ECHEC")
sys.exit(1 if ko else 0)
PY
cp "$SORTIE" "$ICI/tests/sortie_temoins_chemin2.tsv"
rm -f "$SORTIE"
