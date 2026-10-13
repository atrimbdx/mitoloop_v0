#!/usr/bin/env bash
# test_re2g.sh - end-to-end test of re2g_panel.py and peaks_vs_elements.py
#
# Purpose      Build synthetic ENCODE-rE2G link files, a cCRE registry, a GTF,
#              gene lists and loop tables with tests/make_synthetic_re2g.py,
#              run both scripts and check their outputs against the planted
#              truth. Offline; prints "PASS test_re2g" on success and exits
#              non-zero on the first failed assertion.
# Inputs       none
# Outputs      none kept; everything is written to a temporary directory
# Method       1. 20 target genes with 3 elements linked in each of 6 muscle
#                 biosamples, 10 control genes with single-biosample links,
#                 2 liver biosamples excluded by --keep; rows the script must
#                 discard (self promoter, proximal, under threshold, beyond
#                 the window); symbols resolved by the built-in and the extra
#                 alias tables and by Ensembl id.
#              2. --no-builtin-aliases leaves the outdated symbol unresolved.
#              3. 70 biosamples (more than one 64-bit word per element) for
#                 the core curve, with --model.
#              4. peaks_vs_elements.py on loops planted on the replicable
#                 elements (small distance-matched p) and on loops placed
#                 independently of them (no small p).
# Usage        bash tests/test_re2g.sh
# References   see the headers of scripts/re2g_panel.py and
#              scripts/peaks_vs_elements.py

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/../scripts"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

fail() { echo "FAIL test_re2g: $*" >&2; exit 1; }

# value of column $2 in the row of $1 whose first column is $3
cell() {
    python3 - "$1" "$2" "$3" <<'EOF'
import sys
import pandas as pd
t = pd.read_csv(sys.argv[1], sep="\t", dtype=str, keep_default_na=False)
r = t[t.iloc[:, 0] == sys.argv[3]]
if r.empty:
    sys.exit(f"no row {sys.argv[3]} in {sys.argv[1]}")
print(r[sys.argv[2]].iloc[0])
EOF
}

# float comparison: check "<expr>" with a and b bound
check() {
    python3 -c "import sys; a=float(sys.argv[1]); b=float(sys.argv[2]); sys.exit(0 if ($3) else 1)" "$1" "$2" || fail "$4 (a=$1 b=$2)"
}

KEEP="psoas muscle|gastrocnemius|skeletal muscle"

# ---------------------------------------------------------------- 1. main run
python3 "$HERE/make_synthetic_re2g.py" --out "$WORK/s" >/dev/null
cd "$WORK/s"
python3 "$SCRIPTS/re2g_panel.py" --links-dir links --ccre ccre.bed.gz \
    --gtf gencode.gtf.gz --genes genes.txt --controls controls.txt \
    --aliases aliases.tsv --keep "$KEEP" --out re2g --cores 2 \
    --curve-draws 20 --seed 1 > run1.log

for f in re2g_elements.tsv re2g_per_gene.tsv re2g_summary.tsv \
         re2g_sentence.txt re2g_core_curve.tsv re2g_core_curve.svg \
         re2g_elements_per_gene.svg genes_without_tss.txt; do
    [ -s "re2g/$f" ] || fail "missing re2g/$f"
done
[ "$(head -1 re2g/re2g_elements.tsv)" = \
  "$(printf 'gene\tset\tchrom\tstart\tend\tdistance\tprevalence')" ] \
    || fail "re2g_elements.tsv header"
grep -q "threshold 0.201" run1.log || fail "dhs_megamap threshold not detected"
grep -q "6 files" run1.log || fail "--keep did not select the 6 muscle files"
grep -q "Across 6 biosamples of 3 types" re2g/re2g_sentence.txt \
    || fail "sentence content"

# gene set: 20 target genes (TAZ and OLDG2 renamed, FXN_GAA and MT-ND1
# dropped) and 10 controls
[ "$(cell re2g/re2g_summary.tsv n_genes mito)" = 20 ] || fail "mito genes"
[ "$(cell re2g/re2g_summary.tsv n_genes control)" = 10 ] || fail "control genes"
grep -qx "MT-ND1" re2g/genes_without_tss.txt || fail "chrM gene not left out"
for g in TAFAZZIN MITOG2 MITOG5; do
    [ "$(cell re2g/re2g_per_gene.tsv n_intersection $g)" = 3 ] \
        || fail "$g: alias, extra alias or Ensembl lookup"
done

# planted counts: 3 core + 6 single-biosample elements per target gene;
# self-promoter, proximal, sub-threshold and out-of-window rows discarded
[ "$(cell re2g/re2g_per_gene.tsv n_union MITOG1)" = 9 ] || fail "union"
[ "$(cell re2g/re2g_per_gene.tsv links_per_biosample_median MITOG1)" = 4.0000 ] \
    || fail "links per biosample"
[ "$(cell re2g/re2g_per_gene.tsv max_prevalence MITOG1)" = 6 ] || fail "prevalence"
[ "$(cell re2g/re2g_per_gene.tsv n_ccre_window MITOG1)" = 401 ] || fail "registry count"
[ "$(cell re2g/re2g_per_gene.tsv tss MITOG1)" = 3000000 ] || fail "MANE TSS (+)"
cut -f1 re2g/re2g_per_gene.tsv | grep -qxE "OLDG2|TAZ|OLDMITOG5" \
    && fail "old symbol in output"
[ "$(cell re2g/re2g_per_gene.tsv tss MITOG2)" = 6000000 ] || fail "MANE TSS (-)"
[ "$(cell re2g/re2g_per_gene.tsv n_intersection CTRLG1)" = 0 ] || fail "control intersection"
[ "$(cell re2g/re2g_summary.tsv pct_genes_without_any_link control)" = 20.0000 ] \
    || fail "control genes without link"

# planted core elements: right distance (mapped onto the registry, so an
# exact multiple of 5 kb) and prevalence 6
python3 - <<'EOF' || fail "core element distances"
import pandas as pd
el = pd.read_csv("re2g/re2g_elements.tsv", sep="\t")
truth = pd.read_csv("truth.tsv", sep="\t", keep_default_na=False)
for _, t in truth[truth.n_core > 0].iterrows():
    want = sorted(int(x) for x in t.core_offsets.split(","))
    got = sorted(el[(el.gene == t.gene) & (el.prevalence == 6)].distance)
    assert got == want, (t.gene, got, want)
assert (el.distance.abs() >= 10_000).all()
assert (el.distance.abs() <= 1_000_000).all()
EOF
check "$(cell re2g/re2g_core_curve.tsv subsampled_mito 6)" 1 "abs(a-b)<0.0001" "core curve mito"
check "$(cell re2g/re2g_core_curve.tsv subsampled_control 2)" 0 "abs(a-b)<0.0001" "core curve control"

# ------------------------------------------------------- 2. no built-in table
python3 "$SCRIPTS/re2g_panel.py" --links-dir links --ccre ccre.bed.gz \
    --gtf gencode.gtf.gz --genes genes.txt --aliases aliases.tsv \
    --no-builtin-aliases --keep "$KEEP" --out re2g_noalias --cores 1 \
    --curve-draws -1 > run2.log
grep -qx "TAZ" re2g_noalias/genes_without_tss.txt || fail "--no-builtin-aliases"
[ "$(cell re2g_noalias/re2g_summary.tsv n_genes mito)" = 19 ] || fail "genes without built-in aliases"

# ------------------------------------------- 3. more than 64 biosamples, --model
python3 "$HERE/make_synthetic_re2g.py" --out "$WORK/w" --n-mito 4 \
    --n-control 0 --n-muscle 70 --n-other 0 --n-core 1 --core-fraction 0.5 \
    >/dev/null
cd "$WORK/w"
python3 "$SCRIPTS/re2g_panel.py" --links-dir links --links-glob "*.bed.gz" \
    --ccre ccre.bed.gz --gtf gencode.gtf.gz --genes genes.txt \
    --aliases aliases.tsv --model dhs_intact_hic --keep "$KEEP" \
    --out re2g --cores 2 \
    --curve-draws 5 > run3.log
grep -q "threshold 0.243" run3.log || fail "--model threshold"
check "$(cell re2g/re2g_core_curve.tsv subsampled_mito 70)" 0.5 "abs(a-b)<0.0001" "70-biosample curve, k=70"
check "$(cell re2g/re2g_core_curve.tsv subsampled_mito 1)" 1 "abs(a-b)<0.0001" "70-biosample curve, k=1"
check "$(cell re2g/re2g_core_curve.tsv observed_mito 70)" 0.5 "abs(a-b)<0.0001" "70-biosample observed curve"

# ------------------------------------------------------ 4. peaks vs elements
cd "$WORK/s"
python3 "$SCRIPTS/peaks_vs_elements.py" --elements re2g/re2g_elements.tsv \
    --loops loops_planted.tsv --list-terms > terms.log
grep -q "liver" terms.log || fail "--list-terms"
for L in planted null; do
    python3 "$SCRIPTS/peaks_vs_elements.py" --elements re2g/re2g_elements.tsv \
        --loops "loops_$L.tsv" --min-prevalence 6 --min-distance 50000 \
        --called-in "psoas muscle|gastrocnemius" --perms 200 --seed 3 \
        --out "pv_$L" > "pv_$L.log"
    [ -s "pv_$L/peaks_vs_elements.tsv" ] || fail "missing pv_$L output"
done
grep -q "60 of 80 loops kept" pv_planted.log || fail "--called-in filter"
T=pv_planted/peaks_vs_elements.tsv
[ "$(cell $T replicable_elements mito)" = 60 ] || fail "replicable elements"
check "$(cell $T elements_on_a_peak mito)" 40 "a>=b" "planted elements on a peak"
check "$(cell $T observed_elements_on_a_peak_by_distance mito)" \
      "$(cell $T elements_on_a_peak mito)" "abs(a-b)<=3" "distance vs coordinate counts"
check "$(cell $T p_distance_elements mito)" 0.01 "a<=b" "planted, distance-matched p (elements)"
check "$(cell $T p_distance_peaks mito)" 0.01 "a<=b" "planted, distance-matched p (peaks)"
check "$(cell $T p_uniform_elements mito)" 0.01 "a<=b" "planted, uniform p"
T=pv_null/peaks_vs_elements.tsv
check "$(cell $T p_distance_elements mito)" 0.1 "a>b" "unplanted, distance-matched p (elements)"
check "$(cell $T p_distance_peaks mito)" 0.1 "a>b" "unplanted, distance-matched p (peaks)"

echo "PASS test_re2g"
