#!/usr/bin/env bash
# test_choisir_locus.sh - choice of the locus to illustrate (choisir_locus.py)
#                         and the rE2G element track of locus_figure.py
#
# Purpose      a loop ending on an element of muscle biosamples is ranked;
#              a loop on an element of heart biosamples only is not; the
#              locus figure draws the element track, one element on a called
#              loop anchor. Prints "PASS test_choisir_locus".
# Usage        bash tests/test_choisir_locus.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
S="$HERE/../scripts"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
fail() { echo "FAIL test_choisir_locus: $*" >&2; exit 1; }
cd "$WORK"
python3 - <<'PY'
import numpy as np, os
os.makedirs("T/muscle",exist_ok=True)
open("T/loops_muscle.tsv","w").write("gene\tchrom\tvp_pos\tsummit\tstart\tend\tn_fichiers\tn_presents\tfichiers\ttermes\tloop_id\n"
"GA\tchr1\t1000000\t1200000\t1190000\t1210000\t3\t3\tM1|M2|M3\tpsoas\tGA:1200000\n"
"GA\tchr1\t1000000\t900000\t890000\t910000\t3\t1\tM1\tpsoas\tGA:900000\n"
"GB\tchr2\t500000\t530000\t520000\t540000\t3\t2\tM2|M3\tpsoas\tGB:530000\n"
"GC\tchr3\t500000\t800000\t790000\t810000\t3\t3\tM1|M2|M3\tpsoas\tGC:800000\n")
open("T/muscle/peaks_mito.tsv","w").write("term\tgene\tchrom\tvp_pos\tn_libs\tstart\tend\tsummit\tsummit_ratio\treciprocal\n"
"muscle\tGA\tchr1\t1000000\t3\t1190000\t1210000\t1200000\t3.5\t1\nmuscle\tGB\tchr2\t500000\t3\t520000\t540000\t530000\t5\t1\n")
open("T/muscle/loops.tsv","w").write("gene\tchrom\tvp_pos\tsummit\tstart\tend\tn_terms\tcalled\tloop_id\n"
"GA\tchr1\t1000000\t1200000\t1190000\t1210000\t1\tmuscle\tGA:1200000\nGA\tchr1\t1000000\t900000\t890000\t910000\t1\tmuscle\tGA:900000\n")
open("el.tsv","w").write("gene\tset\tchrom\tstart\tend\tdistance\tprevalence\tcelltypes\tbiosamples\n"
"GA\tmito\tchr1\t1199500\t1200300\t0\t3\tx\ta=skeletal muscle myoblast;b=psoas muscle;c=heart left ventricle\n"
"GA\tmito\tchr1\t1100000\t1100500\t0\t1\tx\ta=skeletal muscle myoblast\n"
"GB\tmito\tchr2\t530000\t530500\t0\t2\tx\ta=skeletal muscle myoblast;b=psoas muscle\n"
"GC\tmito\tchr3\t800000\t800500\t0\t2\tx\tc=heart left ventricle\n")
off=np.arange(-500000,500001,10000); prof=np.exp(-np.abs(off)/1e5)*100; prof[off==200000]+=60
np.savez("M1.npz",names=np.array(["GA"]),chrom=np.array(["chr1"]),pos=np.array([1000000]),counts=prof[None,:],offsets=off,accession="M1",term="psoas",resolution=10000,anchor_half=10000)
PY
python3 "$S/choisir_locus.py" --dir T --tissu muscle --elements el.tsv 2> c.log || { cat c.log; fail choisir; }
[ "$(cut -f1 T/locus_muscle.tsv | tail -n +2 | tr '\n' ' ')" = "GA GB " ] || { cat T/locus_muscle.tsv; fail "ranking"; }
[ "$(awk -F'\t' 'NR==2{print $6, $11, $12}' T/locus_muscle.tsv)" = "1200000 2 3" ] || { cat T/locus_muscle.tsv; fail "best loop of GA"; }
grep -q -- "--elements-tissu muscle" T/locus_muscle.tsv || fail "command"
# --l1: GB transmitted (95% interval above 0.5), GA uncertain -> GB alone; tous keeps both
printf 'gene\tscore\tci_low\tci_high\tstatus\nGA\t0.5\t0.2\t0.8\testimable\nGB\t0.8\t0.6\t0.9\testimable\n' > l1.tsv
python3 "$S/choisir_locus.py" --dir T --tissu muscle --elements el.tsv --l1 l1.tsv --out l1.out 2> c2.log || { cat c2.log; fail "--l1"; }
[ "$(cut -f1-3 l1.out | tail -n +2)" = "$(printf 'GB\ttransmitted\t0.80')" ] || { cat l1.out; fail "L1 filter"; }
python3 "$S/choisir_locus.py" --dir T --tissu muscle --elements el.tsv --l1 l1.tsv --l1-classes tous --out l1b.out 2> c3.log
[ "$(cut -f1,2 l1b.out | tail -n +2 | tr '\t\n' ': ')" = "GA:uncertain GB:transmitted " ] || { cat l1b.out; fail "L1 tous"; }
python3 "$S/locus_figure.py" --npz M1.npz --gene GA --loops T/muscle/loops.tsv --arcs \
    --elements el.tsv --elements-tissu muscle --window 500000 --out fig/GA > f.log 2>&1 || { cat f.log; fail figure; }
grep -q "2 rE2G elements of GA in the window (muscle), 1 on a called loop anchor" f.log || { cat f.log; fail "element track"; }
[ -s fig/GA.svg ] || fail svg
grep -q "biosamples rE2G par tissu (seuil) : muscle 2 (2)" c.log || { cat c.log; fail "frac-tissu threshold"; }
# figure with the 80 % rule: muscle has 2 biosamples, an element needs 2 -> the
# 1,100,000 element (1 muscle biosample) is left out
python3 "$S/locus_figure.py" --npz M1.npz --gene GA --loops T/muscle/loops.tsv --arcs \
    --elements el.tsv --elements-tissu muscle --elements-frac 0.8 --window 500000 --out fig/GA8 > f8.log 2>&1 || { cat f8.log; fail "figure frac"; }
grep -q "an element needs 2" f8.log && grep -q "1 rE2G elements of GA in the window (muscle), 1 on a called loop anchor" f8.log \
    || { cat f8.log; fail "element track frac"; }
# anchor rule shared with the reannot AG2 campaign: the peakC region (start-end) widened by one
# bin; an element 6 kb past the summit is on the anchor, one 15 kb past the region end is not
python3 - "$S" <<'PY' || fail "anchor rule"
import sys
sys.path.insert(0, sys.argv[1])
import choisir_locus as c
el = [("chr1", 36000, 36100, {}), ("chr1", 45000, 45100, {})]
assert c.sur_ancre(el, "chr1", 30200, 10000, 30000, 30400) == [el[0]]
assert c.sur_ancre(el, "chr1", 30200, 10000, *c.region({})) == [el[0]]       # summit +/- one bin
assert c.region({"start": "", "end": ""}) == (None, None)
PY
echo "PASS test_choisir_locus"
