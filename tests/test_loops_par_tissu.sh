#!/usr/bin/env bash
# test_loops_par_tissu.sh - per-tissue replication of loops (loops_par_tissu.py)
#
# Purpose      Two muscle libraries, four heart libraries, one kidney library
#              and one smooth-muscle library (excluded from muscle). A loop is
#              kept in 2/2 muscle libraries, in >= 3 of the 4 heart ones, and
#              from the single kidney library; summits 10 kb apart are one
#              loop. Prints "PASS test_loops_par_tissu".
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
H3	right atrium auricular region	x	1	1	1	1
H4	heart right ventricle	x	1	1	1	1
K1	kidney	x	1	1	1	1
S1	smooth muscle of the stomach	x	1	1	1	1
X1	gastrocnemius medialis	x	1	1	1	0
T
{
printf 'accession\tterm\tgene\tchrom\tvp_pos\tstart\tend\tsummit\tsummit_ratio\n'
for l in M1 M2; do printf '%s\tt\tGA\tchr1\t50000\t95000\t105000\t100000\t3\n' $l; done
printf 'M2\tt\tGA\tchr1\t50000\t300000\t310000\t305000\t2\n'          # 1/2 : dropped
printf 'S1\tt\tGA\tchr1\t50000\t300000\t310000\t305000\t2\n'          # smooth muscle: not muscle
printf 'X1\tt\tGA\tchr1\t50000\t300000\t310000\t305000\t2\n'          # library not kept
for l in H1 H2; do printf '%s\tt\tGA\tchr1\t50000\t195000\t205000\t200000\t3\n' $l; done
printf 'H3\tt\tGA\tchr1\t50000\t205000\t215000\t210000\t4\n'          # 10 kb away: same loop
printf 'H1\tt\tGA\tchr1\t50000\t395000\t405000\t400000\t2\n'
printf 'H4\tt\tGA\tchr1\t50000\t395000\t405000\t400000\t2\n'          # 2/4 : dropped
printf 'K1\tt\tGB\tchr2\t70000\t495000\t505000\t500000\t2\n'
printf 'K1\tt\tGZ\tchr2\t70000\t495000\t505000\t500000\t2\n'          # not in --genes
} > peaks_lib.tsv
printf 'GA\nGB\n' > genes.txt
python3 "$S/loops_par_tissu.py" --peaks peaks_lib.tsv --libraries libs.tsv --genes genes.txt --out o 2> log \
    || { cat log; fail run; }
[ "$(cut -f1,4,7,8 o/loops_muscle.tsv | tail -n +2)" = "$(printf 'GA\t100000\t2\t2')" ] || { cat o/loops_muscle.tsv; fail muscle; }
[ "$(cut -f1,4,7,8,9 o/loops_cardiomyocyte.tsv | tail -n +2)" = "$(printf 'GA\t200000\t4\t3\tH1|H2|H3')" ] \
    || { cat o/loops_cardiomyocyte.tsv; fail heart; }
[ "$(cut -f1,4,7,8 o/loops_rein.tsv | tail -n +2)" = "$(printf 'GB\t500000\t1\t1')" ] || { cat o/loops_rein.tsv; fail kidney; }
[ -f o/loops_neurone.tsv ] && [ "$(tail -n +2 o/loops_neurone.tsv | wc -l)" = "0" ] || fail neurone
grep -q "^neurone	0	0	0	0" o/resume_tissus.tsv || { cat o/resume_tissus.tsv; fail resume; }
echo "PASS test_loops_par_tissu"
