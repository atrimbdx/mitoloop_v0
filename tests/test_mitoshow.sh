#!/usr/bin/env bash
# test_mitoshow.sh - end-to-end test of mitoshow_data.py and app/mitoshow.html
#
# Purpose      Build MitoCarta-like sheets, PanelApp genes, HPA and GO
#              locations, STRING links and a curated table at the size of
#              the real data (tests/make_synthetic_mitoshow.py), run
#              mitoshow_data.py and check the locations, pathways and links
#              against the planted truth. When Playwright and a Chromium are
#              available, load the result and a score table in the viewer and
#              check that it draws without error. Offline; prints
#              "PASS test_mitoshow" on success.
# Inputs       none
# Outputs      none kept; everything is written to a temporary directory
# Method       1. MitoCarta sub-localisations (nucleoid for mtDNA matrix
#                 genes), HPA locations outside mitochondria added (Uncertain
#                 ignored), curated places with a nested path, mtDNA genes
#                 without location, PanelApp genes placed by HPA, then
#                 GO, else unlocalized; old symbol found through MitoCarta
#                 synonyms; non-gene panel entries dropped.
#              2. STRING: links above --min-string between genes of the set,
#                 with their evidence channels, mapped through protein.info;
#                 a structure (mmCIF) gives subunit positions and contacts.
#              3. viewer (optional): data and score loaded, every label inside
#                 the drawing, scale chosen, unmatched identifiers reported.
# Usage        bash tests/test_mitoshow.sh
# References   see the header of scripts/mitoshow_data.py

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/../scripts"
APP="$(cd "$HERE/../app" && pwd)/mitoshow.html"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

fail() { echo "FAIL test_mitoshow: $*" >&2; exit 1; }

python3 "$HERE/make_synthetic_mitoshow.py" --out "$WORK/s"
cd "$WORK/s"
python3 "$SCRIPTS/mitoshow_data.py" --mitocarta mitocarta_A.tsv \
    --mitopathways mitocarta_C.tsv --genes genes.txt --hpa hpa.tsv \
    --gaf goa.gaf --extra-locations extra_locations.tsv \
    --structures structures.tsv --id-map hgnc.tsv \
    --string-links string_links.txt.gz --string-info string_info.txt.gz \
    --out mitoshow.json > build.log
grep -q "MitoCarta: 1136 genes" build.log || fail "MitoCarta genes"
grep -q "112 genes outside MitoCarta added" build.log || fail "PanelApp genes"

python3 - <<'PY' || fail "mitoshow.json"
import json
import pandas as pd
d = json.load(open("mitoshow.json"))
g = {x["s"]: x for x in d["genes"]}
assert len(g) == 1248 and "OLDMC2" not in g and "FXN_GAA" not in g
assert g["MC1"]["loc"] == ["matrix"]          # the viewer derives the nucleoid
assert g["MT-TL1"]["loc"] == [] and g["MT-TK"]["src"] == "panelapp"
assert g["MC1"]["pw"] == ["Mitochondrial central dogma > mtDNA maintenance > mtDNA replication"]
assert g["MC2"]["pw"][0].startswith("Mitochondrial central dogma > Translation")
assert g["MC5"]["loc"] == ["matrix", "nucleus"] and g["MC5"]["ls"] == ["MitoCarta3.0", "HPA"]
assert g["MC6"]["loc"] == ["matrix", "cytosol"]
assert g["MC7"]["loc"] == ["matrix"]                    # Uncertain HPA ignored
assert g["MC8"]["extra"] == [{"c": "er", "p": [], "src": "curated"}]
assert g["MC9"]["extra"] == [{"c": "oxphos", "p": ["Supercomplexes", "Respirasome assembly"], "src": "paper"}]
xy = {s: g[s]["xy"]["Complex I"] for s in ("MC530", "MC531", "MC532")}
assert all(0.04 <= v <= 0.96 for p in xy.values() for v in p), xy
assert abs(xy["MC532"][0] - xy["MC530"][0]) > 0.5 or abs(xy["MC532"][1] - xy["MC530"][1]) > 0.5
assert g["MC600"]["loc"] == ["imm"] and g["MC950"]["loc"] == ["omm"]
assert g["MC870"]["loc"] == ["ims"] and g["MC1050"]["loc"] == ["mito_membrane"]
assert g["MC1100"]["loc"] == ["mito_unknown"]
t = pd.read_csv("truth.tsv", sep="\t")
for _, r in t.iterrows():
    assert g[r.gene]["loc"] == [r.location], (r.gene, g[r.gene]["loc"])
    assert g[r.gene]["src"] == "panelapp"
assert g["PA81"]["loc"] == ["unlocalized"]               # NOT annotation
assert d["meta"]["sources"] == ["MitoCarta3.0", "HPA", "GO", "STRING"]
i0, i1, i2 = (next(i for i, x in enumerate(d["genes"]) if x["s"] == s) for s in ("MC530", "MC531", "MC532"))
st = {(l[0], l[1]): l[7] for l in d["links"] if l[7]}
assert (min(i0, i1), max(i0, i1)) in st and (min(i0, i2), max(i0, i2)) not in st
idx = {x["s"]: i for i, x in enumerate(d["genes"])}
pairs = {(a, b) for a, b, *_ in d["links"]}
assert (idx["MC1"], idx["MC3"]) not in pairs             # under --min-string
assert (min(idx["MC1"], idx["PA1"]), max(idx["MC1"], idx["PA1"])) in pairs
assert all(len(l) == 8 and 400 <= l[2] <= 999 and l[0] < l[1] for l in d["links"])
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
    python3 "$HERE/mitoshow_smoke.py" "$APP" "$WORK/s/mitoshow.json" \
        "$WORK/s/scores.tsv" $CHROME > smoke.json || { cat smoke.json; fail "viewer errors"; }
    python3 - <<'PY' || fail "viewer"
import json
r = json.load(open("smoke.json"))
assert r["genes"] == 1248 and r["multi"] == 4, r
assert r["scored"] == 4, r                # OLDMC2 found through its synonym
assert r["scale"] == "div" and r["outside"] == 0
assert r["unmatched"].startswith("1 ")
assert r["selected"] == "MC5" and r["info"].startswith("MC5")
assert r["labels"] >= 4 and r["demo_genes"] > 300
PY
else
    echo "viewer check skipped (Playwright not installed)"
fi

echo "PASS test_mitoshow"
