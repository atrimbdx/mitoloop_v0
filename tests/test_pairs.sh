#!/usr/bin/env bash
# test_pairs.sh - end-to-end test of v4c_pairs.py and v4c_pairs.R
#
# Purpose
#   Builds a planted and a null paired experiment with make_synthetic_pairs.py,
#   runs `v4c_pairs.py extract` and `export`, then `v4c_pairs.R compare`, and
#   checks that the planted experiment shows a loss of concordance or an
#   excess of changing loops for mito against control promoters, and that the
#   null experiment shows neither.
#
# Inputs
#   None. Optional environment: JUICER_TOOLS (path to juicer_tools.jar) to
#   also build .hic libraries and check that the .hic reader (hicstraw)
#   returns the same profiles as the .cool reader. Without it, or without
#   java or hicstraw, that part is skipped with a message.
#
# Outputs
#   Prints PASS pairs on success; exits non-zero on the first failed check.
#   Everything is written to a temporary directory removed on exit.
#
# Method
#   Checks file presence, npz keys, table columns, the planted loops in
#   pairs_peaks.tsv (called in the plus state only, contacts higher in plus),
#   the mito-minus-control statistics of pairs_meta.tsv in both experiments,
#   and that relabelling the states (--plus-state / --minus-state) leaves the
#   results unchanged.
#
# Usage
#   bash tests/test_pairs.sh
#   JUICER_TOOLS=tools/juicer_tools.jar bash tests/test_pairs.sh
#
# References
#   Geeven et al. 2018. peakC: a flexible, non-parametric peak calling package
#     for 4C and Capture-C data. Nucleic Acids Research.
#   Durand et al. 2016. Juicer provides a one-click system for analyzing
#     loop-resolution Hi-C experiments. Cell Systems.
#   Abdennur and Mirny 2020. Cooler: scalable storage for Hi-C data and other
#     genomically labeled arrays. Bioinformatics.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="$ROOT/scripts/v4c_pairs.py"
RS="$ROOT/scripts/v4c_pairs.R"
GEN="$HERE/make_synthetic_pairs.py"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cd "$TMP"

fail() { echo "FAIL pairs: $*" >&2; exit 1; }

command -v Rscript >/dev/null || fail "Rscript not found"
Rscript -e 'quit(status = !requireNamespace("peakC", quietly = TRUE))' \
    || fail "R package peakC not installed"

python3 -m py_compile "$PY" "$GEN"
Rscript -e "invisible(parse('$RS'))"

# --- synthetic data ----------------------------------------------------------
HIC=0
if [[ -n "${JUICER_TOOLS:-}" && -f "${JUICER_TOOLS}" ]] && command -v java >/dev/null \
        && python3 -c "import hicstraw" 2>/dev/null; then
    HIC=1
fi
GEN_OPTS=(--outdir synth --n-genes 6 --frac 0.5 --resolution 10000 --seed 7)
[[ $HIC == 1 ]] && GEN_OPTS+=(--hic-pairs)
python3 "$GEN" "${GEN_OPTS[@]}" > gen.log
for f in pairs.tsv synth.gtf.gz mito.txt controls.txt \
         planted_minus.cool planted_plus.cool null_minus.cool null_plus.cool; do
    [[ -s synth/$f ]] || fail "generator did not write $f"
done

# --- extract and export (.cool) ------------------------------------------------
EXTRACT_OPTS=(--gtf synth=synth/synth.gtf.gz --genes synth/mito.txt
              --controls synth/controls.txt --window 1000000
              --resolution 10000 --anchor-half 10000 --min-raw 100)
python3 "$PY" extract --pairs synth/pairs.tsv "${EXTRACT_OPTS[@]}" \
    --outdir raw > extract.log
python3 "$PY" export --pairs synth/pairs.tsv --indir raw --outdir export \
    > export.log

python3 - <<'EOF' || fail "extract/export outputs"
import gzip
import numpy as np
keys = {"counts", "raw_total", "names", "symbol", "chrom", "pos", "set",
        "offsets", "label", "experiment", "state", "genome", "resolution",
        "anchor_half"}
for lab in ("planted_minus", "planted_plus", "null_minus", "null_plus"):
    d = np.load(f"raw/{lab}.npz")
    assert set(d.files) == keys, (lab, sorted(d.files))
    assert d["counts"].shape == (12, 201), d["counts"].shape
    assert set(d["set"]) == {"mito", "control"}
    assert np.isfinite(d["counts"]).any(1).all(), "a viewpoint is empty"
    with gzip.open(f"export/{lab}.tsv.gz", "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        assert head[:4] == ["name", "set", "chrom", "pos"], head[:4]
        assert head[4] == "-1000000" and head[-1] == "1000000", head[4::196]
        assert sum(1 for _ in fh) == 12
lib = open("export/libraries_pairs.tsv").read().splitlines()
assert lib[0].split("\t") == ["label", "experiment", "state", "genome",
    "resolution", "depth_control", "depth_mito", "n_viewpoints"], lib[0]
assert len(lib) == 5
EOF

# --- .hic route --------------------------------------------------------------
if [[ $HIC == 1 ]]; then
    for lab in planted_minus planted_plus null_minus null_plus; do
        java -Xmx2g -jar "$JUICER_TOOLS" pre -r 10000 -n \
            "synth/$lab.pairs.txt" "synth/$lab.hic" synth/chrom.sizes \
            > "juicer_$lab.log" 2>&1 || fail "juicer pre failed on $lab"
    done
    python3 "$PY" extract --pairs synth/pairs_hic.tsv "${EXTRACT_OPTS[@]}" \
        --outdir raw_hic > extract_hic.log
    python3 - <<'EOF' || fail ".hic and .cool profiles differ"
import numpy as np
for lab in ("planted_minus", "planted_plus", "null_minus", "null_plus"):
    a = np.load(f"raw/{lab}.npz")["counts"]
    b = np.load(f"raw_hic/{lab}.npz")["counts"]
    assert np.array_equal(a, b, equal_nan=True), lab
EOF
    echo ".hic route: profiles identical to .cool"
else
    echo "skipping the .hic part (needs JUICER_TOOLS, java and hicstraw)"
fi

# --- comparison ----------------------------------------------------------------
COMPARE_OPTS=(--wsize 3 --tol 20000 --excess-quantile 0.95 --thin-on control
              --seed 1 --cores 1)
Rscript "$RS" compare --libraries export/libraries_pairs.tsv --export export \
    --out out "${COMPARE_OPTS[@]}" > compare.log 2>&1 \
    || { cat compare.log >&2; fail "v4c_pairs.R compare"; }
for f in pairs_viewpoints.tsv pairs_peaks.tsv pairs_summary.tsv pairs_meta.tsv \
         pairs.svg; do
    [[ -s out/$f ]] || fail "compare did not write $f"
done

python3 - <<'EOF' || fail "planted signal or null behaviour"
import pandas as pd
res, tol = 10000, 20000
# "null" is an experiment name here, not a missing value
na = dict(keep_default_na=False, na_values=["", "NA"])
planted = {"MITOG1": 50, "MITOG2": -60, "MITOG3": 45}

pk = pd.read_csv("out/pairs_peaks.tsv", sep="\t", **na)
for c in ("name", "summit", "distance", "called_minus", "called_plus",
          "count_minus", "count_plus", "log2_ratio", "experiment", "set"):
    assert c in pk.columns, c
for g, off in planted.items():
    want = off * res + res // 2
    hit = pk[(pk.experiment == "planted") & (pk.name == g)
             & ((pk.distance - want).abs() <= tol)]
    assert len(hit) == 1, (g, "planted loop not found")
    h = hit.iloc[0]
    assert bool(h.called_plus) and not bool(h.called_minus), (g, h.to_dict())
    assert h.log2_ratio > 0.5, (g, h.log2_ratio)
    gone = pk[(pk.experiment == "null") & (pk.name == g)
              & ((pk.distance - want).abs() <= tol)]
    assert gone.empty, (g, "loop called in the null experiment")

m = pd.read_csv("out/pairs_meta.tsv", sep="\t", **na).set_index("experiment")
for c in ("d_concordance_loss", "d_excess", "d_direction"):
    assert c in m.columns, c
p, n = m.loc["planted"], m.loc["null"]
print(f"planted: d_concordance_loss {p.d_concordance_loss:+.3f}, "
      f"d_excess {p.d_excess:+.3f}")
print(f"null:    d_concordance_loss {n.d_concordance_loss:+.3f}, "
      f"d_excess {n.d_excess:+.3f}")
assert p.d_concordance_loss > 0.05 or p.d_excess > n.d_excess + 0.05, \
    "planted experiment shows no change at mito promoters"
assert p.d_concordance_loss > n.d_concordance_loss
assert p.d_excess > n.d_excess
assert abs(n.d_concordance_loss) < 0.05, "null shows a loss of concordance"
assert n.d_excess < 0.1, "null shows an excess"
EOF

# --- state labels are options ------------------------------------------------
awk 'BEGIN{FS=OFS="\t"} NR>1{$2 = ($2 == "plus") ? "high" : "low"} {print}' \
    synth/pairs.tsv > synth/pairs_relabel.tsv
python3 "$PY" extract --pairs synth/pairs_relabel.tsv "${EXTRACT_OPTS[@]}" \
    --plus-state high --minus-state low --outdir raw_relabel > /dev/null
python3 "$PY" export --pairs synth/pairs_relabel.tsv --indir raw_relabel \
    --outdir export_relabel --plus-state high --minus-state low > /dev/null
Rscript "$RS" compare --libraries export_relabel/libraries_pairs.tsv \
    --export export_relabel --out out_relabel "${COMPARE_OPTS[@]}" \
    --plus-state high --minus-state low > compare_relabel.log 2>&1 \
    || { cat compare_relabel.log >&2; fail "compare with relabelled states"; }
cmp -s out/pairs_summary.tsv out_relabel/pairs_summary.tsv \
    || fail "relabelled states change the results"

echo "PASS pairs"
