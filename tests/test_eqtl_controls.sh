#!/usr/bin/env bash
# test_eqtl_controls.sh - end-to-end test of make_controls.py and gtex_eqtl.py
#
# Purpose      Builds small synthetic GTEx-like inputs with planted signal
#              (tests/make_synthetic_gtex.py), runs both scripts offline and
#              checks their outputs.
# Inputs       none; everything is generated in a temporary directory.
# Outputs      none kept; prints "PASS eqtl_controls" on success and exits
#              non-zero on the first failed assertion.
# Method       make_controls.py must pick broadly expressed controls, reject
#              tissue-specific and weakly expressed decoys, and drop a perfect
#              decoy lying within --exclusion of a gene-set TSS.
#              gtex_eqtl.py must match on expression, recover larger effect
#              sizes and distances for the gene set in the planted tissue,
#              report no shift in the null tissue, and give identical results
#              from the release tar and from a per-tissue table with renamed
#              columns.
# Usage        bash tests/test_eqtl_controls.sh
# References   see the headers of scripts/make_controls.py and
#              scripts/gtex_eqtl.py

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/../scripts"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

fail() { echo "FAIL eqtl_controls: $*" >&2; exit 1; }

python3 "$HERE/make_synthetic_gtex.py" "$TMP/data" > /dev/null
D="$TMP/data"

# ---- make_controls.py ------------------------------------------------------
cd "$TMP"
python3 "$SCRIPTS/make_controls.py" --gtf "$D/gencode.gtf.gz" \
    --tpm "$D/gtex_median_tpm.gct.gz" --genes "$D/genes.txt" \
    --out controls.txt --pairs control_pairs.csv \
    --report controls_check.tsv > mc.log || fail "make_controls.py exited non-zero"

for f in controls.txt control_pairs.csv controls_check.tsv; do
    [[ -s $f ]] || fail "make_controls.py did not write $f"
done
[[ $(wc -l < controls.txt) -eq 40 ]] || fail "expected 40 controls (20 x 2)"
grep -qvE '^BROAD[0-9]+$' controls.txt && fail "non-BROAD gene among controls: $(grep -vE '^BROAD' controls.txt | tr '\n' ' ')"
head -1 control_pairs.csv | grep -qx 'mito_gene,control_gene,distance' \
    || fail "control_pairs.csv header"
head -1 controls_check.tsv | grep -qx $'set\tn\tmedian_log_tpm\tmedian_tpm_sd\tmedian_length_kb\tmedian_density' \
    || fail "controls_check.tsv header"

# without exclusion the perfect decoy next to MITO1 must be taken
python3 "$SCRIPTS/make_controls.py" --gtf "$D/gencode.gtf.gz" \
    --tpm "$D/gtex_median_tpm.gct.gz" --genes "$D/genes.txt" \
    --exclusion 0 --out controls_noexcl.txt --pairs pairs_noexcl.csv \
    --report check_noexcl.tsv > mc0.log || fail "make_controls.py --exclusion 0"
grep -qx NEAR1 controls_noexcl.txt || fail "NEAR1 missing with --exclusion 0"

# a pool too small leaves gene-set genes unmatched, flagged by an empty field
python3 "$SCRIPTS/make_controls.py" --gtf "$D/gencode.gtf.gz" \
    --tpm "$D/gtex_median_tpm.gct.gz" --genes "$D/genes.txt" \
    --controls 15 --out c15.txt --pairs p15.csv --report r15.tsv > mc15.log \
    || fail "make_controls.py --controls 15"
grep -q "short of 15 controls" mc15.log || fail "pool exhaustion not reported"

python3 - <<'PY' || fail "make_controls.py value checks"
import pandas as pd
rep = pd.read_csv("controls_check.tsv", sep="\t").set_index("set")
assert rep.loc["mito", "n"] == 20, rep
assert rep.loc["control", "n"] == 40, rep
assert abs(rep.loc["mito", "median_log_tpm"] - rep.loc["control", "median_log_tpm"]) < 0.1
assert abs(rep.loc["mito", "median_tpm_sd"] - rep.loc["control", "median_tpm_sd"]) < 0.1
p = pd.read_csv("control_pairs.csv")
assert p.groupby("mito_gene").size().eq(2).all()
assert not p["mito_gene"].str.contains("NOTAGENE").any()
p15 = pd.read_csv("p15.csv", keep_default_na=False)
assert (p15["control_gene"] == "").sum() > 0
PY

# ---- gtex_eqtl.py ----------------------------------------------------------
python3 "$SCRIPTS/gtex_eqtl.py" run --gtf "$D/gencode.gtf.gz" \
    --tpm "$D/gtex_median_tpm.gct.gz" --eqtl "$D/GTEx_eQTL_independent.tar" \
    --genes "$D/genes.txt" --outdir eqtl_tar > ge.log \
    || fail "gtex_eqtl.py run (tar) exited non-zero"
for f in gtex_matched_pairs.csv gtex_matching_check.csv gtex_eqtl_signals.csv \
         gtex_summary.csv gtex_per_gene.csv; do
    [[ -s eqtl_tar/$f ]] || fail "gtex_eqtl.py did not write $f"
done

# same inputs found through --datadir
python3 "$SCRIPTS/gtex_eqtl.py" run --datadir "$D" --genes "$D/genes.txt" \
    --outdir eqtl_datadir > ge_dd.log || fail "gtex_eqtl.py run --datadir"

# per-tissue table with tensorQTL-style column names
python3 "$SCRIPTS/gtex_eqtl.py" run --gtf "$D/gencode.gtf.gz" \
    --tpm "$D/gtex_median_tpm.gct.gz" --eqtl "$D/muscle_eqtl_renamed.txt.gz" \
    --eqtl-gene-column phenotype_id --eqtl-distance-column start_distance \
    --genes "$D/genes.txt" --outdir eqtl_renamed > ge_rn.log \
    || fail "gtex_eqtl.py run (renamed table)"

# null tissue
python3 "$SCRIPTS/gtex_eqtl.py" run --datadir "$D" --genes "$D/genes.txt" \
    --tissue Heart_Left_Ventricle --tpm-column "Heart - Left Ventricle" \
    --overlap-bins 20 \
    --outdir eqtl_null > ge_null.log || fail "gtex_eqtl.py run (null tissue)"

# fetch leaves files already present untouched and needs no network then
python3 "$SCRIPTS/gtex_eqtl.py" fetch --datadir "$D" > fetch.log \
    || fail "gtex_eqtl.py fetch with all files present"
[[ $(grep -c '\[skip\]' fetch.log) -eq 3 ]] || fail "fetch did not skip present files"

# a missing input must stop the offline step with a message
if python3 "$SCRIPTS/gtex_eqtl.py" run --genes "$D/genes.txt" \
        --outdir eqtl_none > ge_none.log 2>&1; then
    fail "gtex_eqtl.py run without inputs should fail"
fi
grep -q -- "--gtf is required" ge_none.log || fail "missing-input message"

python3 - <<'PY' || fail "gtex_eqtl.py value checks"
import pandas as pd
s = pd.read_csv("eqtl_tar/gtex_summary.csv").set_index("quantity")
assert list(s.columns) == ["n_mito", "n_control", "median_mito",
                           "median_control", "cohens_d",
                           "prob_superiority", "overlap"], s.columns
assert s.loc["abs_slope", "n_mito"] == 80 and s.loc["abs_slope", "n_control"] == 240
for q in ("abs_slope", "abs_dist"):
    assert s.loc[q, "prob_superiority"] == 1.0, s
    assert s.loc[q, "overlap"] < 0.05, s
    assert s.loc[q, "cohens_d"] > 2, s

pairs = pd.read_csv("eqtl_tar/gtex_matched_pairs.csv")
assert len(pairs) == 60 and pairs["control_gene"].nunique() == 60
chk = pd.read_csv("eqtl_tar/gtex_matching_check.csv").set_index("set")
assert abs(chk.loc["mito", "median_tpm"] / chk.loc["control", "median_tpm"] - 1) < 0.2, chk

sig = pd.read_csv("eqtl_tar/gtex_eqtl_signals.csv")
assert list(sig.columns) == ["gene_id", "abs_slope", "abs_dist",
                             "tss_distance", "slope", "set"]
assert not sig["gene_id"].str.contains(r"\.").any()

pg = pd.read_csv("eqtl_tar/gtex_per_gene.csv")
assert list(pg.columns) == ["set", "gene_id", "n_signals", "max_dist", "frac_distal"]
assert (pg.loc[pg.set == "mito", "frac_distal"] == 1).all()
assert (pg.loc[pg.set == "control", "frac_distal"] == 0).all()

for other in ("eqtl_datadir", "eqtl_renamed"):
    o = pd.read_csv(f"{other}/gtex_summary.csv").set_index("quantity")
    pd.testing.assert_frame_equal(s, o)

n = pd.read_csv("eqtl_null/gtex_summary.csv").set_index("quantity")
for q in ("abs_slope", "abs_dist"):
    assert 0.35 < n.loc[q, "prob_superiority"] < 0.65, n
    assert n.loc[q, "overlap"] > 0.6, n
PY

echo "PASS eqtl_controls"
