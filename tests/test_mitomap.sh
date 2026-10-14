#!/usr/bin/env bash
# test_mitomap.sh - end-to-end test of mitomap_data.py and app/mitomap.html
#
# Purpose      Build MitoCarta-like sheets, PanelApp genes, HPA and GO
#              locations, STRING links and a curated table at the size of
#              the real data (tests/make_synthetic_mitomap.py), run
#              mitomap_data.py and check the locations, pathways and links
#              against the planted truth. When Playwright and a Chromium are
#              available, load the result and a score table in the viewer and
#              check that it draws without error. Offline; prints
#              "PASS test_mitomap" on success.
# Inputs       none
# Outputs      none kept; everything is written to a temporary directory
# Method       1. MitoCarta sub-localisations (nucleoid for mtDNA matrix
#                 genes), HPA locations outside mitochondria added (Uncertain
#                 ignored), curated table, PanelApp genes placed by HPA, then
#                 GO, else unlocalized; old symbol found through MitoCarta
#                 synonyms; non-gene panel entries dropped.
#              2. STRING: links above --min-string between genes of the set,
#                 with their evidence channels, mapped through protein.info.
#              3. viewer (optional): data and score loaded, every label inside
#                 the drawing, scale chosen, unmatched identifiers reported.
# Usage        bash tests/test_mitomap.sh
# References   see the header of scripts/mitomap_data.py

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/../scripts"
APP="$(cd "$HERE/../app" && pwd)/mitomap.html"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

fail() { echo "FAIL test_mitomap: $*" >&2; exit 1; }

python3 "$HERE/make_synthetic_mitomap.py" --out "$WORK/s"
cd "$WORK/s"
python3 "$SCRIPTS/mitomap_data.py" --mitocarta mitocarta_A.tsv \
    --mitopathways mitocarta_C.tsv --genes genes.txt --hpa hpa.tsv \
    --gaf goa.gaf --extra-locations extra_locations.tsv \
    --string-links string_links.txt.gz --string-info string_info.txt.gz \
    --out mitomap.json > build.log
grep -q "MitoCarta: 1136 genes" build.log || fail "MitoCarta genes"
grep -q "110 genes outside MitoCarta added" build.log || fail "PanelApp genes"

python3 - <<'PY' || fail "mitomap.json"
import json
import pandas as pd
d = json.load(open("mitomap.json"))
g = {x["s"]: x for x in d["genes"]}
assert len(g) == 1246 and "OLDMC2" not in g and "FXN_GAA" not in g
assert g["MC1"]["loc"] == ["nucleoid"]
assert g["MC1"]["pw"] == ["Mitochondrial central dogma > mtDNA maintenance > mtDNA replication"]
assert g["MC2"]["pw"][0].startswith("Mitochondrial central dogma > Translation")
assert g["MC5"]["loc"] == ["matrix", "nucleus"] and g["MC5"]["ls"] == ["MitoCarta3.0", "HPA"]
assert g["MC6"]["loc"] == ["matrix", "cytosol"]
assert g["MC7"]["loc"] == ["matrix"]                    # Uncertain HPA ignored
assert g["MC8"]["loc"][1:] == ["er"] and g["MC8"]["ls"][1] == "curated"
assert g["MC600"]["loc"] == ["imm"] and g["MC950"]["loc"] == ["omm"]
assert g["MC870"]["loc"] == ["ims"] and g["MC1050"]["loc"] == ["mito_membrane"]
assert g["MC1100"]["loc"] == ["mito_unknown"]
t = pd.read_csv("truth.tsv", sep="\t")
for _, r in t.iterrows():
    assert g[r.gene]["loc"] == [r.location], (r.gene, g[r.gene]["loc"])
    assert g[r.gene]["src"] == "panelapp"
assert g["PA81"]["loc"] == ["unlocalized"]               # NOT annotation
assert d["meta"]["sources"] == ["MitoCarta3.0", "HPA", "GO", "STRING"]
idx = {x["s"]: i for i, x in enumerate(d["genes"])}
pairs = {(a, b) for a, b, *_ in d["links"]}
assert (idx["MC1"], idx["MC3"]) not in pairs             # under --min-string
assert (min(idx["MC1"], idx["PA1"]), max(idx["MC1"], idx["PA1"])) in pairs
assert all(len(l) == 7 and 400 <= l[2] <= 999 and l[0] < l[1] for l in d["links"])
assert len(pairs) == len(d["links"])                     # one row per pair
exp = [l for l in d["links"] if (l[0], l[1]) == (min(idx["MC1"], idx["PA1"]), max(idx["MC1"], idx["PA1"]))][0]
assert exp[3] == 900                                     # experimental channel
PY

# --------------------------------------------------- viewer (optional)
CHROME="${CHROME:-}"
if [ -z "$CHROME" ] && [ -x /opt/pw-browsers/chromium-1194/chrome-linux/chrome ]; then
    CHROME=/opt/pw-browsers/chromium-1194/chrome-linux/chrome
fi
if python3 -c "import playwright" 2>/dev/null; then
    printf 'gene\tvalue\nMC5\t2.5\nMC6\t-1\nPA1\t0.5\nOLDMC2\t1\nUNKNOWN1\t3\n' > scores.tsv
    python3 "$HERE/mitomap_smoke.py" "$APP" "$WORK/s/mitomap.json" \
        "$WORK/s/scores.tsv" $CHROME > smoke.json || { cat smoke.json; fail "viewer errors"; }
    python3 - <<'PY' || fail "viewer"
import json
r = json.load(open("smoke.json"))
assert r["genes"] == 1246 and r["instances"] == 1249 and r["multi"] == 3
assert r["scored"] == 4, r                # OLDMC2 found through its synonym
assert r["scale"] == "div" and r["outside"] == 0
assert r["unmatched"].startswith("1 ")
assert r["selected"] == "MC5" and r["info"].startswith("MC5")
assert r["labels"] > 100 and r["demo_genes"] > 100
PY
else
    echo "viewer check skipped (Playwright not installed)"
fi

echo "PASS test_mitomap"
