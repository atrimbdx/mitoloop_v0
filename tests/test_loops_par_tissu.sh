#!/usr/bin/env bash
# test_loops_par_tissu.sh - per-tissue loops (loops_par_tissu.py)
#
# Purpose      librairies: kept libraries grouped by tissue, one folder per
#              tissue, adrenal gland not kidney, smooth muscle not muscle;
#              partenaires: partner anchors of all tissues missing from an
#              existing export; fusionner: two export folders merged library
#              by library; separer: loops_<tissue>.tsv from each tissue's
#              loops.tsv, with the per-library control count (call_lib).
#              The R steps (call, reciprocal) are not run here. Prints
#              "PASS test_loops_par_tissu".
# Usage        bash tests/test_loops_par_tissu.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
S="$HERE/../scripts"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
fail() { echo "FAIL test_loops_par_tissu: $*" >&2; exit 1; }
cd "$WORK"
cat > libs.tsv <<'T'
accession	term	gtex	mtcn	depth	log10_depth	kept
M1	gastrocnemius medialis	x	1	1	1	1
M2	psoas muscle	x	1	1	1	1
H1	heart left ventricle	x	1	1	1	1
H2	heart left ventricle	x	1	1	1	1
H3	right cardiac atrium	x	1	1	1	1
K1	kidney	x	1	1	1	1
A1	adrenal gland	x	1	1	1	1
S1	smooth muscle of the stomach	x	1	1	1	1
X1	gastrocnemius medialis	x	1	1	1	0
T
python3 "$S/loops_par_tissu.py" librairies --libraries libs.tsv --out T 2> l.log || { cat l.log; fail librairies; }
[ "$(cut -f1 T/muscle/libraries.tsv | tail -n +2 | tr '\n' ' ')" = "M1 M2 " ] || fail "muscle libraries"
[ "$(cut -f1 T/rein/libraries.tsv | tail -n +2 | tr '\n' ' ')" = "K1 " ] || fail "kidney: adrenal excluded"
[ "$(cut -f2 T/cardiomyocyte/libraries.tsv | tail -n +2 | sort -u)" = "cardiomyocyte" ] || fail "term = tissue"
head -1 T/muscle/libraries.tsv | grep -q "term_encode" || fail "term_encode column"
[ ! -e T/neurone ] || fail "no folder without library"
grep -q "rein           1 librairies (kidney) -> appel sur une seule librairie" l.log || { cat l.log; fail "single library flag"; }

# partner anchors: P_chr1_300000 already exported, P_chr2_900000 new
printf 'name\tchrom\tpos\tset\nP_chr1_300000\tchr1\t300000\tpartner\n' > T/muscle/partners.tsv
printf 'name\tchrom\tpos\tset\nP_chr1_300000\tchr1\t300000\tpartner\nP_chr2_900000\tchr2\t900000\tpartner\n' > T/rein/partners.tsv
mkdir -p exp_old exp_new
for a in M1 K1; do
    printf 'name\tchrom\tpos\t-10000\t0\t10000\nP_chr1_300000\tchr1\t300000\t1\t2\t3\n' | gzip > exp_old/$a.tsv.gz
    printf 'name\tchrom\tpos\t-10000\t0\t10000\nP_chr2_900000\tchr2\t900000\t4\t5\t6\nP_chr1_300000\tchr1\t300000\t9\t9\t9\n' | gzip > exp_new/$a.tsv.gz
done
python3 "$S/loops_par_tissu.py" partenaires --dir T --export exp_old --out manq.tsv 2> p.log || { cat p.log; fail partenaires; }
[ "$(cut -f1 manq.tsv | tail -n +2)" = "P_chr2_900000" ] || { cat manq.tsv; fail "missing partners"; }
# --genes: only the partners of the peaks of the wanted genes
printf 'term\tgene\tchrom\tvp_pos\tn_libs\tstart\tend\tsummit\tsummit_ratio\nrein\tGB\tchr2\t70000\t1\t890000\t910000\t900000\t2\nrein\tCTRL\tchr1\t10000\t1\t290000\t310000\t300000\t2\n' > T/rein/peaks_mito.tsv
printf 'GB\n' > mito.txt
python3 "$S/loops_par_tissu.py" partenaires --dir T --export exp_old --genes mito.txt --out manq2.tsv 2> p2.log || { cat p2.log; fail "partenaires --genes"; }
grep -q "^2 partenaires\|^1 partenaires" p2.log && [ "$(cut -f1 manq2.tsv | tail -n +2)" = "P_chr2_900000" ] || { cat p2.log manq2.tsv; fail "partenaires --genes content"; }
python3 "$S/loops_par_tissu.py" fusionner --export exp_old exp_new --out exp_all 2> f.log || { cat f.log; fail fusionner; }
[ "$(zcat exp_all/K1.tsv.gz | cut -f1,4 | tr '\t\n' ' ;')" = "name -10000;P_chr1_300000 1;P_chr2_900000 4;" ] \
    || { zcat exp_all/K1.tsv.gz; fail "merged export"; }

# separer: each tissue's loops.tsv, control count from peaks_lib.tsv
printf 'gene\tchrom\tvp_pos\tsummit\tstart\tend\tn_terms\tcalled\tloop_id\nGA\tchr1\t50000\t100000\t90000\t110000\t1\tmuscle\tGA:100000\nGZ\tchr1\t50000\t100000\t90000\t110000\t1\tmuscle\tGZ:100000\n' > T/muscle/loops.tsv
printf 'gene\tchrom\tvp_pos\tsummit\tstart\tend\tn_terms\tcalled\tloop_id\nGB\tchr2\t70000\t500000\t490000\t510000\t1\trein\tGB:500000\n' > T/rein/loops.tsv
printf 'accession\tterm\tgene\tchrom\tvp_pos\tstart\tend\tsummit\tsummit_ratio\nM1\tgastrocnemius medialis\tGA\tchr1\t50000\t95000\t105000\t100000\t3\nM2\tpsoas muscle\tGA\tchr1\t50000\t300000\t310000\t305000\t2\n' > peaks_lib.tsv
printf 'GA\nGB\n' > genes.txt
python3 "$S/loops_par_tissu.py" separer --dir T --genes genes.txt --peaks-lib peaks_lib.tsv 2> s.log \
    || { cat s.log; fail separer; }
[ "$(cut -f1,4,7,8,9,10 T/loops_muscle.tsv | tail -n +2)" = "$(printf 'GA\t100000\t2\t1\tM1\tgastrocnemius medialis|psoas muscle')" ] \
    || { cat T/loops_muscle.tsv; fail "muscle loops"; }
[ "$(cut -f1,7 T/loops_rein.tsv | tail -n +2)" = "$(printf 'GB\t1')" ] || { cat T/loops_rein.tsv; fail "kidney loops"; }
[ -f T/loops_neurone.tsv ] && [ "$(tail -n +2 T/loops_neurone.tsv | wc -l)" = "0" ] || fail "neurone empty"
grep -q "^rein	1	1	1" T/resume_tissus.tsv || { cat T/resume_tissus.tsv; fail resume; }
echo "PASS test_loops_par_tissu"
