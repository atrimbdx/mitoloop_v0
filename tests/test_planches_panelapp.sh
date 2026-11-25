#!/usr/bin/env bash
# planches_panelapp.py : liste des gènes PanelApp (niveaux, renommage), fenêtre par gène, librairie
# la plus profonde, assemblage de la planche, index et reprise ; sans réseau ni .hic (le tracé
# du locus échoue ici et doit être consigné sans arrêter la boucle), données synthétiques
set -euo pipefail
cd "$(dirname "$0")/.."
python3 -c "import matplotlib, PIL" 2>/dev/null || { echo "SKIP test_planches_panelapp (matplotlib ou Pillow absent)"; exit 0; }
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
fail() { echo "FAIL test_planches_panelapp: $1"; exit 1; }
mkdir -p "$T/L/muscle"
printf '# PanelApp Mitochondrial disorders id 112 version 1\ngene\tconfidence_level\tentity_type\nGA\t3\tgene\nOLDB\t3\tgene\nGC\t2\tgene\nGD\t3\tstr\nGA\t3\tgene\n' > "$T/panel.tsv"
printf 'old\tnew\nOLDB\tGB\n' > "$T/ren.tsv"
printf 'accession\tterm\tdepth\tkept\nLIB1\tmuscle\t100\t1\nLIB2\tmuscle\t900\t1\nLIB3\tmuscle\t5000\t0\n' > "$T/L/muscle/libraries.tsv"
printf 'gene\tchrom\tvp_pos\tsummit\nGA\tchr1\t1000000\t1400000\nGA\tchr1\t1000000\t1100000\n' > "$T/L/muscle/loops.tsv"
S=scripts/planches_panelapp.py
python3 $S --panelapp "$T/panel.tsv" --renommes "$T/ren.tsv" --tissu muscle --dir "$T/L" --out "$T/o" --liste > "$T/liste" 2> "$T/l.log" \
    || { cat "$T/l.log"; fail liste; }
[ "$(cut -f1,3 "$T/liste" | tr '\n' ';')" = "GA	oui;GB	non;" ] || { cat "$T/liste"; fail "niveaux, renommage, doublons"; }
grep -q "2 gènes PanelApp (niveaux 3) ; 1 avec boucle" "$T/l.log" || { cat "$T/l.log"; fail "décompte"; }
python3 $S --panelapp "$T/panel.tsv" --renommes "$T/ren.tsv" --niveaux 3,2 --tissu muscle --dir "$T/L" --out "$T/o" --liste > "$T/liste2" 2>/dev/null
[ "$(wc -l < "$T/liste2")" = "3" ] || fail "--niveaux 3,2"
python3 - "$T" <<'PY'
import sys
sys.path.insert(0, "scripts")
from pathlib import Path
import planches_panelapp as pp
T = Path(sys.argv[1])
assert pp.librairie(T / "L", "muscle") == "LIB2"                       # la plus profonde des gardées
f = pp.fenetres(T / "L", "muscle", 1_000_000)
assert f == {"GA": 550_000}, f                                         # 400 kb + 150 kb
assert pp.fenetres(T / "L", "muscle", 500_000) == {"GA": 500_000}      # plafond
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
for n, (w, h) in {"a": (3, 6), "b": (8, 3)}.items():
    fig = plt.figure(figsize=(w, h)); fig.savefig(T / f"{n}.png", dpi=50); plt.close(fig)
pp.planche(T / "a.png", T / "b.png", T / "ab.png")
from PIL import Image
a, b, ab = (Image.open(T / x) for x in ("a.png", "b.png", "ab.png"))
assert ab.size == (a.width + b.width + 20, a.height), (a.size, b.size, ab.size)
pp.planche(T / "a.png", None, T / "a2.png")
assert Image.open(T / "a2.png").size == a.size
PY
# sans .hic : chaque gène en échec, consigné, la boucle continue ; table des URL absente ou trouvée seule
! (cd "$T" && python3 "$OLDPWD/$S" --panelapp panel.tsv --tissu muscle --dir L --out o 2> /dev/null) || fail "table introuvable"
! python3 $S --panelapp "$T/panel.tsv" --tissu muscle --dir "$T/L" --bibliotheque TABLE --out "$T/o" 2> /dev/null || fail "--bibliotheque absente"
printf 'accession\ts3_url\nLIB2\t\n' > "$T/bib.tsv"
mkdir -p "$T/data" && cp "$T/bib.tsv" "$T/data/libs.tsv"
python3 - "$T" <<'PY'
import sys; sys.path.insert(0, "scripts")
from pathlib import Path
import planches_panelapp as pp
T = Path(sys.argv[1])
assert pp.trouver_bibliotheque(T, "LIB2").endswith("data/libs.tsv")
assert pp.trouver_bibliotheque(T, "LIBX") == ""
PY
python3 $S --panelapp "$T/panel.tsv" --renommes "$T/ren.tsv" --tissu muscle --dir "$T/L" --bibliotheque "$T/bib.tsv" \
    --gtf "$T/absent.gtf" --out "$T/o" 2> "$T/r.log" || { cat "$T/r.log"; fail "boucle"; }
[ "$(awk -F'\t' 'NR>1{print $1, $4}' "$T/o/index.tsv" | tr '\n' ';')" = "GA echec;GB echec;" ] || { cat "$T/o/index.tsv" "$T/r.log"; fail index; }
grep -q "0 complètes, 0 sans zoom, 2 en échec" "$T/r.log" || { cat "$T/r.log"; fail bilan; }
[ -z "$(ls "$T/o" | grep -v index.tsv)" ] || fail "fichiers intermédiaires laissés"
# tranches (tableau SLURM) puis fusion des index
for k in 0 1; do
    python3 $S --panelapp "$T/panel.tsv" --renommes "$T/ren.tsv" --tissu muscle --dir "$T/L" --bibliotheque "$T/bib.tsv" \
        --gtf "$T/absent.gtf" --tranche $k/2 --out "$T/t" 2> "$T/t$k.log" || { cat "$T/t$k.log"; fail "tranche $k"; }
done
[ "$(awk -F'\t' 'NR>1{print $1}' "$T/t/index_0.tsv")" = "GA" ] && [ "$(awk -F'\t' 'NR>1{print $1}' "$T/t/index_1.tsv")" = "GB" ] \
    || { cat "$T"/t/index_*.tsv; fail "partage des tranches"; }
python3 $S --tissu muscle --out "$T/t" --fusionner 2> "$T/f.log" || { cat "$T/f.log"; fail fusion; }
[ "$(awk -F'\t' 'NR>1{print $1, $4}' "$T/t/index.tsv" | tr '\n' ';')" = "GA echec;GB echec;" ] || { cat "$T/t/index.tsv"; fail "index fusionné"; }
echo "PASS test_planches_panelapp"
