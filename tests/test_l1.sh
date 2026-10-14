#!/usr/bin/env bash
# test_l1.sh - end-to-end test of l1_score.py
#
# Purpose      Build a tissue atlas with matched RNA counts and protein, a
#              second atlas without replicates, a patient fibroblast cohort
#              and external tables (tests/make_synthetic_l1.py), run `score`
#              and `validate`, and check them against the planted modes of
#              control. Offline; prints "PASS test_l1" on success and exits
#              non-zero on the first failed assertion.
# Inputs       none
# Outputs      none kept; everything is written to a temporary directory
# Method       1. score: transcriptional genes score near 1, buffered and
#                 post-transcriptional genes near 0, flat genes are
#                 'flat' and unscored, noisy genes are recovered by the
#                 correction for attenuation; complex subunits follow their
#                 complex at the protein level more than at the mRNA level.
#              2. without the mitochondrial content adjustment, buffered
#                 genes look transcriptional; without the ischemic-time
#                 covariate, the genes it affects score lower.
#              3. an atlas without replicates gives uncorrected scores.
#              4. validate: replicate, external table with its sign, expected
#                 sets, patient RNA outliers reaching the protein for
#                 high-L1 genes.
#              5. regulome.py compare reads l1.tsv.
# Usage        bash tests/test_l1.sh
# References   see the header of scripts/l1_score.py

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/../scripts"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

fail() { echo "FAIL test_l1: $*" >&2; exit 1; }

python3 "$HERE/make_synthetic_l1.py" --out "$WORK/s"
cd "$WORK/s"

run_score() {
    python3 "$SCRIPTS/l1_score.py" score --rna atlas_rna.gct \
        --protein atlas_protein.tsv --samples atlas_samples.tsv \
        --gtf gencode.gtf.gz --genes genes.txt --complexes complexes.tsv \
        --boot 100 "$@"
}

# ---------------------------------------------------------------- 1. score
run_score --covariates SMTSISCH --out l1 > score.log
run_score --out l1_nocov > /dev/null
run_score --covariates SMTSISCH --no-content-adjustment --out l1_raw > /dev/null
grep -q "64 matched samples" score.log || fail "sample matching"
grep -q "(counts)" score.log || fail "RNA units"
python3 - <<'PY' || fail "score"
import pandas as pd
t = pd.read_csv("truth.tsv", sep="\t")
def load(d):
    m = pd.read_csv(f"{d}/l1.tsv", sep="\t").merge(t, on="gene")
    return m[m.set == "mito"].set_index("gene")
m = load("l1")
med = m.groupby("class").score.median()
assert med["transcriptional"] > 0.9 and med["noisy"] > 0.9, med
assert abs(med["buffered"]) < 0.3 and abs(med["post"]) < 0.3, med
assert (m[m["class"] == "flat"].status == "flat").all()
assert m[m["class"] == "flat"].score.isna().all()
assert (m[m["class"] != "flat"].status == "estimable").all()
noisy = m[m["class"] == "noisy"]
assert (noisy.r_corrected - noisy.r).median() > 0.05      # attenuation
assert noisy.rel_protein.median() < m[m["class"] == "transcriptional"].rel_protein.median()
e = m[m.status == "estimable"]
assert (e.ci_low <= e.ci_high).all() and e.ci_low.between(-1, 1).all()
assert ((e.ci_low <= e.score + 0.05) & (e.ci_high >= e.score - 0.05)).mean() > 0.9
sub = m[m["class"] == "subunit"]
assert (sub.complex == "CI subunits").all()
assert (sub.protein_coupling > sub.rna_coupling).mean() > 0.9
assert sub.complex_score.iloc[0] > sub.score.median() + 0.1
assert "MT-ATP8" not in m.index
raw = load("l1_raw")
assert raw.groupby("class").score.median()["buffered"] > 0.4
nocov = load("l1_nocov")
aff = [f"MITOG{i}" for i in range(73, 81)]
assert m.loc[aff, "score"].mean() > nocov.loc[aff, "score"].mean()
c = pd.read_csv("l1/content_index.tsv", sep="\t")
assert c.groupby("tissue").content_rna.mean().std() > 0.5
s = pd.read_csv("l1/l1_complexes.tsv", sep="\t")
assert s.loc[0, "complex"] == "CI subunits" and s.loc[0, "n_subunits"] == 12
PY

# ------------------------------------------------- 3. atlas without replicates
python3 "$SCRIPTS/l1_score.py" score --rna hpa_rna.tsv \
    --protein hpa_protein.tsv --samples hpa_samples.tsv --gtf gencode.gtf.gz \
    --genes genes.txt --out l1_hpa > hpa.log
grep -q "(tpm)" hpa.log || fail "TPM detected"
grep -q "(linear)" hpa.log || fail "linear protein detected"
python3 - <<'PY' || fail "uncorrected atlas"
import pandas as pd
m = pd.read_csv("l1_hpa/l1.tsv", sep="\t")
assert set(m.status) <= {"estimable_uncorrected", "flat"}
assert m.rel_rna.isna().all()
PY

# ------------------------------------------------------------- 4. validate
python3 "$SCRIPTS/l1_score.py" validate --l1 l1/l1.tsv \
    --replicate l1_hpa/l1.tsv --external buffering=buffering.tsv:buffering:- \
    --expected-high expected_high.txt --expected-low expected_low.txt \
    --patients-rna patients_rna.tsv --patients-protein patients_protein.tsv \
    --patients-samples patients_samples.tsv --gtf gencode.gtf.gz \
    --out val > validate.log
grep -q "protein centred within 3 batches" validate.log || fail "TMT batches"
python3 - <<'PY' || fail "validate"
import pandas as pd
v = pd.read_csv("val/validation.tsv", sep="\t").set_index("test")
for t in ("replicate", "replicate_mito", "buffering", "buffering_mito",
          "expected_sets", "outlier_protein_z",
          "fibroblast_rna_protein_r_mito"):
    assert v.loc[t, "agrees"] == 1, t
    assert v.loc[t, "p"] < 0.01, (t, v.loc[t, "p"])
assert v.loc["buffering", "value"] < 0
e = pd.read_csv("val/patient_events.tsv", sep="\t")
planted = {f"MITOG{i}" for i in range(1, 31)}
assert planted <= set(e.gene)
hi = e[e.gene.isin([f"MITOG{i}" for i in range(1, 16)])]
lo = e[e.gene.isin([f"MITOG{i}" for i in range(16, 31)]) & (e.rna_z < -6)]
assert hi.propagated.mean() > 0.9 and lo.propagated.mean() < 0.2
assert e[(e.gene == "MITOG2") & (e.rna_z < -6)].protein_missing.iloc[0] == 1
PY

# ------------------------------------------------- 5. regulome.py compare
python3 - <<'PY'
import numpy as np, pandas as pd
m = pd.read_csv("l1/l1.tsv", sep="\t")
rng = np.random.default_rng(1)
g = m[m.set == "mito"].gene
pd.DataFrame({"gene": g, "set": "mito", "n_union": np.round(
    3 + 8 * m.set_index("gene").loc[g, "score"].fillna(0).to_numpy()
    + rng.normal(0, 1, len(g)))}).to_csv("re2g_per_gene.tsv", sep="\t",
                                         index=False)
PY
python3 "$SCRIPTS/regulome.py" compare --scores l1/l1.tsv \
    --re2g re2g_per_gene.tsv --perms 0 --out cmp > cmp.log
python3 - <<'PY' || fail "compare on l1.tsv"
import pandas as pd
c = pd.read_csv("cmp/compare.tsv", sep="\t").set_index("feature")
assert c.loc["n_union", "rho"] > 0.7, c.loc["n_union", "rho"]
PY

echo "PASS test_l1"
