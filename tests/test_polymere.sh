#!/usr/bin/env bash
# test_polymere.sh - preparation of a MitoLoop locus for simpol path A and the
#                    analysis of the deletions (polymere/), without OpenMM:
#
# Purpose      preparer: window, anchors and hub from cibles.tsv and the
#              tissue's loops; oriented CTCF sites found under synthetic peaks;
#              SANSELEMENT removes the site inside the element, one
#              SANSANCRE<kb> per distal loop summit removes the site there, and a deletion that removes no
#              site is flagged null by construction; plan: the deletion plan
#              of the chosen parameter set; contact_deletion: contact
#              frequencies from synthetic conformations (element bead close
#              to the promoter in SAUVAGE only). Prints "PASS test_polymere".
# Usage        bash tests/test_polymere.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
P="$HERE/../polymere"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
fail() { echo "FAIL test_polymere: $*" >&2; exit 1; }
cd "$WORK"
python3 - <<'PY'
import random
random.seed(1)
# chr1, 1.2 Mb of random sequence with three CTCF motifs (MA0139 consensus)
seq = [random.choice("ACGT") for _ in range(1_200_000)]
motif = "CCGCGAGGGGGCGCTGTGG"                 # reverse complement of a strong MA0139 match
fwd = "CCACCAGGGGGCGCTA"
def put(pos, s):
    seq[pos:pos + len(s)] = list(s)
put(500_100, "TGGCCACCAGGGGGCGCTA")            # promoter
put(620_100, "TGGCCACCAGGGGGCGCTA")            # element (on the anchor)
put(700_100, "TAGCGCCCCCTGGTGGCCA")            # distal summit, minus strand
open("g.fa", "w").write(">chr1\n" + "\n".join("".join(seq[i:i + 60]) for i in range(0, len(seq), 60)) + "\n")
with open("pics.bed", "w") as fo:
    for s in (500_000, 620_000, 700_000):
        fo.write(f"chr1\t{s}\t{s + 300}\tp\t100\t.\t50\t-1\t-1\t150\n")
with open("cibles.tsv", "w") as fo:
    fo.write("cible\tgene\tchrom\tdebut\tfin\ttss\tbrin\tclasse\ttype\n")
    fo.write("promoteur:GA\tGA\tchr1\t499750\t500100\t500000\t+\tpromoteur\ttous\n")
    fo.write("element:neurone:GA\tGA\tchr1\t619900\t620400\tGA\t+\telement_sur_ancre\tneurone\n".replace("\tGA\t+", "\t500000\t+"))
    fo.write("element:rein:GA\tGA\tchr1\t640000\t640300\t500000\t+\telement_sur_ancre\trein\n")
with open("loops.tsv", "w") as fo:
    fo.write("gene\tchrom\tvp_pos\tsummit\tstart\tend\tn_terms\tcalled\tloop_id\n")
    fo.write("GA\tchr1\t500000\t620000\t610000\t630000\t1\tneurone\tGA:620000\n")
    fo.write("GA\tchr1\t500000\t700000\t690000\t710000\t1\tneurone\tGA:700000\n")
PY
python3 "$P/preparer_locus.py" preparer --gene GA --type neurone --cibles cibles.tsv --loops loops.tsv \
    --ctcf-pics pics.bed --genome g.fa --marge 100000 --flanc 50000 --replicats-ajustement 2 --out L > p.log 2>&1 \
    || { cat p.log; fail preparer; }
grep -q "locus ga : chr1:450,000-750,000" p.log && grep -q "hub : promoteur_GA, element_1" p.log \
    || { cat p.log; fail "window or hub"; }
[ "$(wc -l < L/ctcf_SAUVAGE.bed)" = "3" ] || { cat p.log L/ctcf_SAUVAGE.bed; fail "CTCF sites"; }
grep -q "SANSELEMENT      1 site(s) CTCF retiré(s) sur 3 : 620,109(+)" p.log \
    && grep -q "SANSANCRE620     1 site(s) CTCF retiré(s) sur 3 : 620,109(+)" p.log \
    && grep -q "SANSANCRE700     1 site(s) CTCF retiré(s) sur 3 : 700,109(-)" p.log \
    || { cat p.log; fail "deletions"; }
[ "$(wc -l < L/plan_ajustement.tsv)" = "37" ] || fail "fitting plan"
python3 "$P/preparer_locus.py" plan --dir L --jeu 5 --replicats 4 > q.log || { cat q.log; fail plan; }
[ "$(cut -f1,3 L/plan_deletion.tsv | sort -u | tail -n +2 | wc -l)" = "4" ] && grep -q "16 trajectoires" q.log \
    || { cat q.log L/plan_deletion.tsv; fail "deletion plan"; }
# synthetic conformations: 60 beads (450,000 - 100,000 ... 5 kb), promoter bead
# and element bead in contact in SAUVAGE, apart in SANSELEMENT
SIMPOL_LOCI=L/loci.json P="$P" python3 - <<'PY'
import os, sys, json
import numpy as np
sys.path.insert(0, os.environ["P"] + "/simpol")
import config as cfg
n, d0, f0 = cfg.geometrie("ga")
idx = cfg.indices_ancres("ga")
os.makedirs("L/resultats", exist_ok=True)
rng = np.random.default_rng(0)
for cond, proche in (("SAUVAGE", True), ("SANSELEMENT", False), ("SANSANCRE700", True)):
    for r in range(4):
        conf = np.cumsum(rng.normal(0, 1, (200, n, 3)), axis=1) * 3
        if proche:
            conf[:, idx["element_1"]] = conf[:, idx["promoteur_GA"]] + 0.5
        np.save(f"L/resultats/ga_{cond}_A_jeu005_rep{r:02d}.npy", conf.astype(np.float32))
PY
SIMPOL_LOCI=L/loci.json python3 "$P/contact_deletion.py" --dir L --locus ga --tirages 200 > c.log 2>&1 \
    || { cat c.log; fail contact_deletion; }
python3 - <<'PY' || { cat c.log L/contacts_deletion.tsv; fail "contact_deletion content"; }
import csv
r = {(x["condition"], x["ancre"]): x for x in csv.DictReader(open("L/contacts_deletion.tsv"), delimiter="\t")}
e = r[("SANSELEMENT", "element_1")]
assert float(e["sauvage"]) == 1.0 and float(e["deletion"]) < 0.2 and float(e["ci_haut"]) < -2, e
assert abs(float(r[("SANSANCRE700", "element_1")]["log2_ratio"])) < 0.01
PY
[ -s L/contacts_deletion.png ] || fail png
# witnesses: equal arrest forces, CTCF sites at random positions
python3 "$P/preparer_locus.py" temoins --dir L --jeu 5 --jeux 5 --tirages 2 --replicats 2 > t.log || { cat t.log; fail temoins; }
grep -q "plan_temoins.tsv : 6 trajectoires" t.log && [ "$(cut -f7 L/ctcf_UNIFORME.bed | sort -u | wc -l)" = "1" ] \
    && [ "$(cut -f6 L/ctcf_ALEATOIRE1.bed | sort | tr '\n' ' ')" = "$(cut -f6 L/ctcf_SAUVAGE.bed | sort | tr '\n' ' ')" ] \
    || { cat t.log L/ctcf_UNIFORME.bed L/ctcf_ALEATOIRE1.bed; fail "temoins content"; }
# observed / expected agreement: measured matrix with a promoter-element peak
SIMPOL_LOCI=L/loci.json P="$P" python3 - <<'PY'
import os, sys
sys.path.insert(0, os.environ["P"] + "/simpol")
import config as cfg
n, d0, f0 = cfg.geometrie("ga")
with open("L/matrice.tsv", "w") as fo:
    fo.write("startA\tendA\tstartB\tendB\tN\n")
    for i in range(450_000, 750_000, 10_000):
        for j in range(i, 750_000, 10_000):
            v = 100 / (1 + (j - i) / 10_000)
            if i == 500_000 and j == 620_000:
                v *= 6
            fo.write(f"{i}\t{i + 10000}\t{j}\t{j + 10000}\t{v}\n")
PY
SIMPOL_LOCI=L/loci.json python3 "$P/evaluer_4c.py" --dir L --locus ga > e.log 2>&1 || { cat e.log; fail evaluer_4c; }
python3 - <<'PY' || { cat e.log; fail "evaluer_4c content"; }
import csv
r = {(x["condition"], x["jeu"]): x for x in csv.DictReader(open("L/evaluation_4c.tsv"), delimiter="\t")}
assert float(r[("MESURE", ".")]["oe_promoteur_GA__element_1"]) > 3
assert float(r[("SAUVAGE", "5")]["oe_promoteur_GA__element_1"]) > float(r[("SANSELEMENT", "5")]["oe_promoteur_GA__element_1"])
PY
# deletions in the UNIFORME setting, and the reading against UNIFORME
python3 "$P/preparer_locus.py" uniforme --dir L --jeu 5 --replicats 2 > un.log || { cat un.log; fail uniforme; }
grep -q "plan_uniforme.tsv : 8 trajectoires" un.log && [ "$(cut -f7 L/ctcf_UNIF_SANSANCRE700.bed | sort -u | wc -l)" = "1" ] \
    || { cat un.log L/plan_uniforme.tsv; fail "uniforme content"; }
SIMPOL_LOCI=L/loci.json P="$P" python3 - <<'PY'
import os, sys
import numpy as np
sys.path.insert(0, os.environ["P"] + "/simpol")
import config as cfg
n, d0, f0 = cfg.geometrie("ga")
idx = cfg.indices_ancres("ga")
rng = np.random.default_rng(1)
for cond, proche in (("UNIFORME", True), ("UNIF_SANSELEMENT", False)):
    for r in range(3):
        conf = np.cumsum(rng.normal(0, 1, (200, n, 3)), axis=1) * 3
        if proche:
            conf[:, idx["element_1"]] = conf[:, idx["promoteur_GA"]] + 0.5
        np.save(f"L/resultats/ga_{cond}_A_jeu005_rep{r:02d}.npy", conf.astype(np.float32))
PY
SIMPOL_LOCI=L/loci.json python3 "$P/contact_deletion.py" --dir L --locus ga --reference UNIFORME --conditions conditions_uniforme.tsv --tirages 100 > cu.log 2>&1 \
    || { cat cu.log; fail "contact_deletion UNIFORME"; }
grep -q "UNIF_SANSELEMENT promoteur_GA ↔ element_1" cu.log || { cat cu.log; fail "uniforme reading"; }
echo "PASS test_polymere"
