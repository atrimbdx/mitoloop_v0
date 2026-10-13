#!/usr/bin/env bash
# test_figures.sh - end-to-end test of the figure scripts
#
# Purpose      Runs locus_figure.py, folding_figure.py and svg_to_pptx.py on
#              synthetic inputs with planted loops and checks their outputs.
# Inputs       none; tests/make_figure_inputs.py writes them in a temporary
#              directory. JUICER_TOOLS (path to juicer_tools.jar) is optional:
#              when it is set, java is present and hicstraw is importable, a
#              .hic is built and read as well; otherwise that part is skipped.
# Outputs      none kept; prints PASS test_figures on success, exits non-zero
#              on the first failed assertion.
# Method       Builds a .cool (and a .hic), npz files in both key schemes,
#              loops.tsv with several viewpoints, a GTF and re2g_elements.tsv;
#              draws the locus figure with and without the map, with arcs and
#              gene track, the folding figure with --gtf and --elements, and
#              converts SVGs to .pptx. Checks that files exist and are not
#              empty, that the virtual 4C peaks at the planted summits, that
#              the profiles from the map and from the npz agree, and that
#              each .pptx opens with python-pptx and holds shapes.
# Usage        bash tests/test_figures.sh
#              JUICER_TOOLS=tools/juicer_tools.jar bash tests/test_figures.sh
# References   Durand NC et al. 2016. Juicer provides a one-click system for
#                  analyzing loop-resolution Hi-C experiments. Cell Systems.
#              Abdennur N, Mirny LA. 2020. Cooler: scalable storage for Hi-C
#                  data and other genomically labeled arrays. Bioinformatics.
#              python-pptx (software). https://github.com/scanny/python-pptx

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/../scripts"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

fail() { echo "FAIL test_figures: $*" >&2; exit 1; }
nonempty() { for f in "$@"; do [[ -s "$f" ]] || fail "missing or empty: $f"; done; }
run() { "$@" > "$TMP/last.log" 2>&1 || { cat "$TMP/last.log" >&2; fail "command failed: $*"; }; }

IN="$TMP/in"
OUT="$TMP/out"
mkdir -p "$OUT"
run python3 "$HERE/make_figure_inputs.py" "$IN"
nonempty "$IN/synthetic.cool" "$IN/loops.tsv" "$IN/genes.gtf" \
         "$IN/re2g_elements.tsv" "$IN/extract.npz" "$IN/profiles.npz"

# --- locus figure: map from .cool, arcs, gene track ------------------------
run python3 "$SCRIPTS/locus_figure.py" --hic "$IN/synthetic.cool" \
    --gene GENEA --gtf "$IN/genes.gtf" --loops "$IN/loops.tsv" \
    --norm NONE --arcs --genes-track --out "$OUT/locus_cool"
nonempty "$OUT/locus_cool".{svg,png,tsv}
grep -q "2 called loops inside the window" "$TMP/last.log" \
    || fail "expected 2 called loops at GENEA"

# raw contacts, square map, options exposed for the drawing
run python3 "$SCRIPTS/locus_figure.py" --hic "$IN/synthetic.cool" \
    --gene GENEB --gtf "$IN/genes.gtf" --loops "$IN/loops.tsv" --raw \
    --norm NONE --font-scale 1.3 --marker-size 7 --vmin 1 --dpi 100 \
    --width 5 --out "$OUT/locus_raw"
nonempty "$OUT/locus_raw".{svg,png,tsv}

# --- locus figure without the map: npz in both key schemes ----------------
run python3 "$SCRIPTS/locus_figure.py" --npz "$IN/extract.npz" --gene GENEA \
    --loops "$IN/loops.tsv" --out "$OUT/v4c_extract"
nonempty "$OUT/v4c_extract".{svg,png,tsv}
run python3 "$SCRIPTS/locus_figure.py" --npz "$IN/profiles.npz" --gene GENEB \
    --loops "$IN/loops.tsv" --out "$OUT/v4c_profiles"
nonempty "$OUT/v4c_profiles".{svg,png,tsv}
grep -q "position taken from loops.tsv" "$TMP/last.log" \
    || fail "position of GENEB should come from loops.tsv"
# no position in the npz and none in loops.tsv for this call: from the GTF
run python3 "$SCRIPTS/locus_figure.py" --npz "$IN/profiles.npz" --gene GENEC \
    --gtf "$IN/genes.gtf" --hic "$IN/synthetic.cool" --loops "$IN/loops.tsv" \
    --arcs --genes-track --out "$OUT/v4c_map_npz"
nonempty "$OUT/v4c_map_npz".{svg,png,tsv}
grep -q "the rest at GENEB" "$TMP/last.log" \
    || fail "arc panel should include the loop of the neighbouring viewpoint"

# a library table without a URL must be refused with a message
printf "accession\tterm_name\tgtex_tissue\tmtcn\tdownload_url\nSYN0001\tpsoas muscle\tMuscle\t1500\tnone\n" \
    > "$IN/tissues.tsv"
if python3 "$SCRIPTS/locus_figure.py" --tissues "$IN/tissues.tsv" \
        --term psoas --gene GENEA --gtf "$IN/genes.gtf" \
        --out "$OUT/no_url" > "$TMP/last.log" 2>&1; then
    fail "a term without URL should fail"
fi
grep -q "has no URL or existing file" "$TMP/last.log" \
    || fail "missing location message"

# a library table pointing to a local file through its `path` column
printf "accession\tterm_name\tgtex_tissue\tmtcn\tpath\nSYN0002\tpsoas muscle\tMuscle\t1500\t%s\n" \
    "$IN/synthetic.cool" > "$IN/tissues_local.tsv"
run python3 "$SCRIPTS/locus_figure.py" --tissues "$IN/tissues_local.tsv" \
    --term psoas --gene GENEA --gtf "$IN/genes.gtf" --norm NONE \
    --out "$OUT/local_path"
nonempty "$OUT/local_path.svg" "$OUT/local_path.png"

# planted signal and agreement between the map and the npz profile
python3 - "$OUT" "$IN/planted.tsv" <<'EOF' || fail "profile checks"
import sys
import pandas as pd
out, planted = sys.argv[1], pd.read_csv(sys.argv[2], sep="\t")
cool = pd.read_csv(f"{out}/locus_cool.tsv", sep="\t")
npz = pd.read_csv(f"{out}/v4c_extract.tsv", sep="\t")
assert list(cool.columns) == ["chrom", "start", "end", "distance",
                              "virtual_4c"], cool.columns
both = cool.merge(npz, on=["chrom", "start", "end"], suffixes=("_m", "_n"))
both = both.dropna()
assert len(both) > 50
assert (abs(both.virtual_4c_m - both.virtual_4c_n) < 1e-3).all(), \
    "map and npz profiles differ"
for s in planted.loc[planted.gene == "GENEA", "summit"]:
    at = cool.loc[(cool.start <= s) & (cool.end > s), "virtual_4c"].iloc[0]
    far = cool.loc[(abs(cool.start - s) > 60_000) &
                   (abs(cool.distance) > 100_000), "virtual_4c"].median()
    assert at > 3 * far, (s, at, far)
EOF

# --- .hic, when Juicer Tools and hicstraw are available -------------------
if [[ -n "${JUICER_TOOLS:-}" && -f "${JUICER_TOOLS}" ]] && command -v java >/dev/null \
        && python3 -c "import hicstraw" 2>/dev/null; then
    sort -k3,3n -k7,7n "$IN/pairs.txt" > "$IN/pairs.sorted.txt"
    run java -Xmx2g -jar "$JUICER_TOOLS" pre -r 10000 -n \
        "$IN/pairs.sorted.txt" "$IN/synthetic.hic" "$IN/chrom.sizes"
    run python3 "$SCRIPTS/locus_figure.py" --hic "$IN/synthetic.hic" \
        --gene GENEA --gtf "$IN/genes.gtf" --loops "$IN/loops.tsv" \
        --norm NONE --arcs --out "$OUT/locus_hic"
    nonempty "$OUT/locus_hic".{svg,png,tsv}
    cmp -s <(cut -f1-5 "$OUT/locus_hic.tsv") <(cut -f1-5 "$OUT/locus_cool.tsv") \
        || fail ".hic and .cool profiles differ"
else
    echo "  skip: .hic part (JUICER_TOOLS, java or hicstraw unavailable)"
fi

# --- folding figure --------------------------------------------------------
run python3 "$SCRIPTS/folding_figure.py" --loops "$IN/loops.tsv" \
    --gene GENEA --gtf "$IN/genes.gtf" --elements "$IN/re2g_elements.tsv" \
    --out "$OUT/folding"
nonempty "$OUT/folding".{svg,png,tsv}
grep -q "1 of them carrying a replicable element" "$TMP/last.log" \
    || fail "the +250 kb anchor should carry the replicable element"
head -1 "$OUT/folding.tsv" | grep -q "^petal	anchor_a	anchor_b" \
    || fail "folding.tsv header"
awk -F'\t' 'NR>1 && $3=="+250 kb" && $8=="enh"{ok=1} END{exit !ok}' \
    "$OUT/folding.tsv" || fail "enhancer anchor missing from the rosette"
run python3 "$SCRIPTS/folding_figure.py" --loops "$IN/loops.tsv" \
    --gene GENEB --lang fr --seed 5 --min-terms 5 --font-size 9 \
    --out "$OUT/folding_fr"
nonempty "$OUT/folding_fr".{svg,png,tsv}
grep -q "3 loops" "$TMP/last.log" || fail "--min-terms 5 should keep 3 loops"

# --- SVG to PowerPoint -----------------------------------------------------
run python3 "$SCRIPTS/svg_to_pptx.py" "$OUT/locus_cool.svg" "$OUT/locus.pptx"
run python3 "$SCRIPTS/svg_to_pptx.py" "$OUT/folding.svg" "$OUT/folding.pptx"
nonempty "$OUT/locus.pptx" "$OUT/folding.pptx"
python3 - "$OUT" <<'EOF' || fail "pptx checks"
import sys
from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE_TYPE
out = sys.argv[1]
for name, want_picture in (("locus", True), ("folding", False)):
    prs = Presentation(f"{out}/{name}.pptx")
    shapes = list(prs.slides[0].shapes)
    texts = [s.text_frame.text for s in shapes if s.has_text_frame
             and s.text_frame.text]
    pics = [s for s in shapes if s.shape_type == MSO_SHAPE_TYPE.PICTURE]
    assert len(shapes) > 20, (name, len(shapes))
    assert texts, f"{name}: no text box"
    assert "GENEA" in " ".join(texts), f"{name}: gene label missing"
    if want_picture:
        assert pics, f"{name}: contact map picture missing"
EOF

echo "PASS test_figures"
