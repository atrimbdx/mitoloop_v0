#!/usr/bin/env bash
# test_promoters.sh - end-to-end test of promoter_features.py
#
# Purpose      Build a synthetic genome, GTF, CpG islands, loops, elements,
#              methylation files in four formats and a TPM table with
#              tests/make_synthetic_promoters.py, run promoter_features.py and
#              check its outputs against the planted truth. Offline; prints
#              "PASS test_promoters" on success and exits non-zero on the
#              first failed assertion.
# Inputs       none
# Outputs      none kept; everything is written to a temporary directory
# Method       1. sequence: target promoters CpG-rich (HCP) with an NRF1 site,
#                 a third of the controls CpG-rich without it; NRF1 enriched
#                 also within promoter classes, CpG island overlap not.
#              2. architecture: 10 bidirectional target promoters, 5 with a
#                 promoter-promoter loop.
#              3. methylation: ENCODE bedMethyl (a CpG under the coverage
#                 threshold), Bismark coverage, bedGraph in fraction and in
#                 percent; MITOG3 unmethylated only in muscle, found as the
#                 one variable target promoter, anti-correlated with its
#                 expression; elements at 0.3.
# Usage        bash tests/test_promoters.sh
# References   see the header of scripts/promoter_features.py

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/../scripts"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

fail() { echo "FAIL test_promoters: $*" >&2; exit 1; }

# value of column $2 in the row of $1 whose first column is $3
cell() {
    python3 - "$1" "$2" "$3" <<'PY'
import sys
import pandas as pd
t = pd.read_csv(sys.argv[1], sep="\t", dtype=str, keep_default_na=False)
r = t[t.iloc[:, 0] == sys.argv[3]]
if r.empty:
    sys.exit(f"no row {sys.argv[3]} in {sys.argv[1]}")
print(r[sys.argv[2]].iloc[0])
PY
}

check() {
    python3 -c "import sys; a=float(sys.argv[1]); b=float(sys.argv[2]); sys.exit(0 if ($3) else 1)" "$1" "$2" || fail "$4 (a=$1 b=$2)"
}

python3 "$HERE/make_synthetic_promoters.py" --out "$WORK/s"
cd "$WORK/s"
python3 "$SCRIPTS/promoter_features.py" --gtf gencode.gtf.gz \
    --genes genes.txt --controls controls.txt --fasta genome.fa \
    --cpg-islands cpg_islands.txt.gz --loops loops.tsv \
    --meth-table methylation.tsv --elements re2g_elements.tsv \
    --min-prevalence 6 --tpm tpm.gct --min-tissues 3 --focus muscle \
    --perms 200 --seed 1 --out pf > run.log

for f in promoter_features.tsv promoter_feature_tests.tsv methylation_long.tsv \
         methylation_by_group.tsv elements_methylation_by_group.tsv \
         methylation_summary.tsv methylation_variable.tsv \
         methylation_expression.tsv; do
    [ -s "pf/$f" ] || fail "missing pf/$f"
done

# ------------------------------------------------------------ 1. sequence
F=pf/promoter_features.tsv
[ "$(cell $F promoter_class MITOG7)" = HCP ] || fail "target promoter class"
[ "$(cell $F promoter_class CTRLG30)" = LCP ] || fail "CpG-poor control class"
[ "$(cell $F promoter_class CTRLG2)" = HCP ] || fail "CpG-rich control class"
[ "$(cell $F motif_NRF1 MITOG7)" = 1 ] || fail "planted NRF1 site (+ strand)"
[ "$(cell $F motif_NRF1 MITOG8)" = 1 ] || fail "planted NRF1 site (- strand)"
[ "$(cell $F motif_NRF1 CTRLG2)" = 0 ] || fail "NRF1 in a control"
[ "$(cell $F tss MITOG1)" = 200000 ] || fail "MANE TSS (+)"
[ "$(cell $F tss MITOG2)" = 400000 ] || fail "MANE TSS (-)"
[ "$(cell $F cpg_island MITOG1)" = 1 ] || fail "CpG island"
[ "$(cell $F cpg_island CTRLG30)" = 0 ] || fail "no CpG island"
T=pf/promoter_feature_tests.tsv
check "$(cell $T p motif_NRF1)" 1e-6 "a<b" "NRF1 enrichment"
check "$(cell $T p_within_class motif_NRF1)" 0.05 "a<b" "NRF1 within classes"
check "$(cell $T p cpg_island)" 1e-4 "a<b" "CpG island enrichment"
check "$(cell $T p_within_class cpg_island)" 0.5 "a>b" "CpG islands explained by class"
check "$(cell $T target cpg_oe)" 0.75 "a>b" "target CpG o/e"

# -------------------------------------------------------- 2. architecture
[ "$(cell $F bidirectional MITOG4)" = 1 ] || fail "bidirectional pair"
[ "$(cell $F bidirectional_partner MITOG4)" = BIDIR_MITOG4 ] || fail "partner"
[ "$(cell $F bidirectional MITOG14)" = 0 ] || fail "no bidirectional pair"
check "$(cell $T target bidirectional)" 0.3333 "abs(a-b)<0.001" "bidirectional share"
[ "$(cell $F n_promoter_loops MITOG2)" = 1 ] || fail "promoter-promoter loop"
[ "$(cell $F loop_partners MITOG2)" = HUB_MITOG2 ] || fail "loop partner"
[ "$(cell $F n_promoter_loops MITOG9)" = 0 ] || fail "loop without promoter"
[ "$(cell $F n_loops MITOG9)" = 1 ] || fail "loop count"

# --------------------------------------------------------- 3. methylation
# experimental annotations: ChIP peaks (one long peak nesting a short one),
# CAGE peak widths, ChromHMM states at the TSS, SCREEN classes in column 6
python3 - <<'PY'
import gzip
import pandas as pd
f = pd.read_csv("pf/promoter_features.tsv", sep="\t")
tss = dict(zip(f.gene, f.tss))
with open("chip.bed", "w") as fh:
    for g in f[f.set == "mito"].gene:
        fh.write(f"chr1\t{tss[g] - 100}\t{tss[g] + 100}\tpeak\n")
    fh.write(f"chr1\t{tss['CTRLG2'] - 50000}\t{tss['CTRLG2'] - 20000}\tlong\n")
    fh.write(f"chr1\t{tss['CTRLG2'] - 30000}\t{tss['CTRLG2'] - 29000}\tshort\n")
with open("cage.bed", "w") as fh:
    fh.write(f"chr1\t{tss['MITOG1'] - 10}\t{tss['MITOG1'] + 10}\tp1\n")
    fh.write(f"chr1\t{tss['MITOG2'] - 150}\t{tss['MITOG2'] + 150}\tp2\n")
    fh.write(f"chr1\t{tss['MITOG2'] + 400}\t{tss['MITOG2'] + 405}\tp3\n")
with gzip.open("states.bed.gz", "wt") as fh:
    rows = sorted([(tss["MITOG1"] - 600, tss["MITOG1"] + 400, "1_TssA"),
                   (tss["CTRLG30"] - 500_000, tss["CTRLG30"] + 500_000, "15_Quies")])
    for s, e, n in rows:
        fh.write(f"chr1\t{s}\t{e}\t{n}\n")
with open("ccre.bed", "w") as fh:
    fh.write(f"chr1\t{tss['MITOG7'] - 150}\t{tss['MITOG7'] + 150}\tEH38D1\tEH38E1\tPLS,CTCF-bound\n")
    fh.write(f"chr1\t{tss['MITOG8'] - 900}\t{tss['MITOG8'] - 600}\tEH38D2\tEH38E2\tpELS\n")
PY
python3 "$SCRIPTS/promoter_features.py" --gtf gencode.gtf.gz \
    --genes genes.txt --controls controls.txt --perms 50 --out pa \
    --annotation chip=chip.bed --annotation cage=cage.bed:width:200 \
    --annotation tssA=states.bed.gz:label=TssA:0 \
    --annotation pls=ccre.bed:label=^PLS --annotation absent=nofile.bed > pa.log
A=pa/promoter_features.tsv
[ "$(cell $A exp_chip MITOG3)" = 1 ] || fail "ChIP peak on a target"
[ "$(cell $A exp_chip CTRLG2)" = 0 ] || fail "ChIP peak away from the promoter"
[ "$(cell $A expw_cage MITOG1)" = 20 ] || fail "sharp CAGE peak"
[ "$(cell $A expw_cage MITOG2)" = 300 ] || fail "broad CAGE peak (nearest of two)"
[ "$(cell $A expw_cage MITOG3)" = NA ] || fail "no CAGE peak"
[ "$(cell $A exp_tssA MITOG1)" = 1 ] || fail "TssA state"
[ "$(cell $A exp_tssA CTRLG30)" = 0 ] || fail "quiescent state"
[ "$(cell $A exp_pls MITOG7)" = 1 ] || fail "PLS cCRE"
[ "$(cell $A exp_pls MITOG8)" = 0 ] || fail "pELS is not PLS"
grep -q "absent: nofile.bed absent, skipped" pa.log || fail "absent annotation file"
check "$(cell pa/promoter_feature_tests.tsv p exp_chip)" 1e-6 "a<b" "ChIP enrichment test"

G=pf/methylation_by_group.tsv
for grp in "skeletal muscle" liver brain heart; do
    check "$(cell $G "$grp" MITOG1)" 0.05 "abs(a-b)<0.01" "MITOG1 $grp"
    check "$(cell $G "$grp" CTRLG30)" 0.8 "abs(a-b)<0.03" "CTRLG30 $grp"
done
check "$(cell $G liver MITOG3)" 0.8 "abs(a-b)<0.03" "MITOG3 outside muscle"
check "$(cell $G "skeletal muscle" MITOG3)" 0.05 "abs(a-b)<0.01" "MITOG3 in muscle"
check "$(cell pf/elements_methylation_by_group.tsv heart chr1:249600-250400)" \
    0.3 "abs(a-b)<0.01" "element methylation"
[ "$(tail -n +2 pf/methylation_variable.tsv | awk -F'\t' '$2=="mito"' | wc -l)" = 1 ] \
    || fail "one variable target promoter"
[ "$(cell pf/methylation_variable.tsv focus_minus_other MITOG3)" = -0.75 ] \
    || fail "focus difference"
check "$(cell pf/methylation_expression.tsv rho MITOG3)" 0 "a<-0.5" "methylation-expression"
[ "$(tail -n +2 pf/methylation_expression.tsv | wc -l)" = 1 ] \
    || fail "correlations only where methylation varies"
python3 - <<'PY' || fail "methylation summary"
import pandas as pd
s = pd.read_csv("pf/methylation_summary.tsv", sep="\t")
p = s[(s.kind == "promoter") & (s.group == "liver")].set_index("set")
assert abs(p.loc["mito", "share_unmethylated"] - 29 / 30) < 1e-3
assert abs(p.loc["control", "share_methylated"] - 20 / 30) < 1e-3
assert p.loc["mito", "p"] < 1e-6
PY

# --controls-all: every other protein-coding gene (30 controls and 5 hubs)
python3 "$SCRIPTS/promoter_features.py" --gtf gencode.gtf.gz \
    --genes genes.txt --controls-all --fasta genome.fa --perms 0 \
    --out pf_all > run_all.log
[ "$(cell pf_all/promoter_feature_tests.tsv n_control motif_NRF1)" = 35 ] \
    || fail "--controls-all"

# the .fai reader used without pysam
python3 - "$SCRIPTS" <<'PY' || fail "FASTA reader without pysam"
import sys
sys.path.insert(0, sys.argv[1])
sys.modules["pysam"] = None
import promoter_features as pf
fa = pf.Fasta("genome.fa")
assert fa.pysam is None
ref = "".join(l.strip() for l in open("genome.fa") if not l.startswith(">"))
for s, e in ((0, 10), (59, 61), (199_000, 200_500), (len(ref) - 5, len(ref))):
    assert fa.fetch("chr1", s, e) == ref[s:e], (s, e)
assert fa.fetch("1", 100, 130) == ref[100:130]
PY

echo "PASS test_promoters"
