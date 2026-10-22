#!/usr/bin/env bash
# test_promoter_elements.sh - core promoter elements and promoter types in
#                             the MitoLoop package (mitoloop_paquet.py,
#                             ajout_mitoloop.py)
#
# Purpose      A small genome with three genes: GENET (+, TATA box at -30,
#              initiator at +1, G-quadruplex at +30, expressed in one
#              tissue: type tissu), GENEH (+, CpG-rich, expressed everywhere:
#              menage), GENED (+, CpG-rich, Polycomb-repressed in every
#              ChromHMM tissue: developpement); an anchor of GENEH that holds
#              the promoter of another gene. Checks the type of each TSS, the
#              arguments TATA_detruite, Inr_detruit, G4_detruit and
#              zone_pause, MITO_contact_promoteur, MITO_type_promoteur, and
#              the priority rule for distal elements by promoter type.
#              Offline; prints "PASS test_promoter_elements".
# Usage        bash tests/test_promoter_elements.sh

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
R="$HERE/../reannot"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail() { echo "FAIL test_promoter_elements: $*" >&2; exit 1; }
cd "$WORK"

python3 - <<'PY'
import gzip, random
random.seed(11)
N = 30000
seq = [random.choice("AT") if random.random() < 0.6 else random.choice("ACGT") for _ in range(N)]
def put(p, s):
    seq[p:p + len(s)] = list(s)
# GENET, + strand, TSS 5000 (0-based): TATA at -30, Inr at -2..+5, G4 at +30
put(5000 - 30, "TATAAAAG")
put(5000 - 2, "CCAGTCC")
put(5000 + 30, "GGGAGGGAGGGAGGG")
# GENEH and GENED: CpG-rich promoters (-1000..+500)
for t in (12000, 20000):
    for p in range(t - 1000, t + 500):
        seq[p] = random.choice("CG") if random.random() < 0.8 else random.choice("AT")
seq = "".join(seq)
open("genome.fa", "w").write(">chr1\n" + "\n".join(seq[i:i + 60] for i in range(0, N, 60)) + "\n")
open("genes.txt", "w").write("GENET\nGENEH\nGENED\n")
open("l1.tsv", "w").write("gene\tgene_id\tset\tstatus\tscore\tci_low\tci_high\tcomplex\t"
                          "protein_coupling\trna_coupling\tl1_class\n")
with open("regions.bed", "w") as fh:
    fh.write("#chrom\tstart\tend\tgene\tclass\tsource\ttss\tstrand\tdetail\n")
    fh.write("chr1\t4500\t5600\tGENET\tpromoter\tgencode\t5000\t+\t\n")
    fh.write("chr1\t8000\t8400\tGENET\telement\tre2g\t5000\t+\t\n")
    fh.write("chr1\t11000\t12600\tGENEH\tpromoter\tgencode\t12000\t+\t\n")
    fh.write("chr1\t15000\t15400\tGENEH\tanchor\tloops\t12000\t+\t\n")
    fh.write("chr1\t19000\t20600\tGENED\tpromoter\tgencode\t20000\t+\t\n")
# EPDnew: -49..+10, spaces
with open("epd.bed", "w") as fh:
    for g, t in (("GENET", 5000), ("GENEH", 12000), ("GENED", 20000)):
        fh.write(f"chr1 {t - 49} {t + 11} {g}_1 900 +\n")
with gzip.open("cage.bed.gz", "wt") as fh:
    fh.write("chr1\t4998\t5004\tp1\t1\t+\n")          # sharp
    fh.write("chr1\t11960\t12060\tp2\t1\t+\n")        # broad
# ChromHMM: GENED TSS Polycomb-repressed in both tissues
for t in ("muscle", "heart"):
    with gzip.open(f"chm_{t}.bed.gz", "wt") as fh:
        fh.write("chr1\t0\t19900\t15_Quies\nchr1\t19900\t20100\t13_ReprPC\nchr1\t20100\t30000\t15_Quies\n")
with gzip.open("genes.gtf.gz", "wt") as fh:
    for g, s, e, st, ty in (("GENET", 5001, 7000, "+", "protein_coding"),
                            ("GENEH", 12001, 14000, "+", "protein_coding"),
                            ("GENED", 20001, 22000, "+", "protein_coding"),
                            ("OTHERP", 15201, 17000, "+", "protein_coding")):
        fh.write(f'chr1\tx\tgene\t{s}\t{e}\t.\t{st}\t.\tgene_id "ENSG_{g}.1"; gene_type "{ty}"; gene_name "{g}";\n')
with open("tpm.gct", "w") as fh:
    fh.write("#1.2\n3\t5\nName\tDescription\tT1\tT2\tT3\tT4\tT5\n")
    fh.write("ENSG_GENET\tGENET\t200\t0\t0\t0\t0\n")
    fh.write("ENSG_GENEH\tGENEH\t50\t60\t55\t40\t45\n")
    fh.write("ENSG_GENED\tGENED\t5\t6\t4\t5\t6\n")
comp = {"A": "C", "C": "A", "G": "T", "T": "G"}
with open("run.final.tsv", "w") as fh:
    fh.write("transcrit\tvariant_entree\tchrom\tpos\tref\talt\tgene\tstatut\n")
    def row(name, p0, alt=None):
        b = seq[p0]
        fh.write(f"T\t{name}\tchr1\t{p0 + 1}\t{b}\t{alt or comp[b]}\t.\tok\n")
    row("tata", 5000 - 28)                 # second T of TATAAAAG
    row("inr", 5000, "G")                  # A at +1
    row("g4", 5000 + 31)                   # inside the first G run
    row("anchor", 15100)                   # in the anchor holding OTHERP's promoter
PY

python3 "$R/mitoloop_paquet.py" --genes genes.txt --l1 l1.tsv --regions regions.bed \
    --epd epd.bed --cage cage.bed.gz --chromhmm muscle=chm_muscle.bed.gz \
    --chromhmm heart=chm_heart.bed.gz --gtf genes.gtf.gz --fasta genome.fa --tpm tpm.gct \
    --out pkg > paquet.log 2>&1 || { cat paquet.log; fail "paquet"; }
grep -q "type de promoteur (TSS EPDnew) : developpement 1, menage 1, tissu 1 ; 3 gènes" paquet.log \
    || { cat paquet.log; fail "promoter types"; }
grep -q "ancres portant le promoteur d'un autre gène : 1" paquet.log || { cat paquet.log; fail "anchor"; }
python3 - <<'PY' || fail "genes.tsv"
import csv
g = {r["gene"]: r for r in csv.DictReader(open("pkg/genes.tsv"), delimiter="\t")}
assert g["GENET"]["type_promoteur"] == "tissu", g["GENET"]
assert g["GENEH"]["type_promoteur"] == "menage", g["GENEH"]
assert g["GENED"]["type_promoteur"] == "developpement", g["GENED"]
PY

python3 "$R/ajout_mitoloop.py" --data pkg --fasta genome.fa --samtools no-samtools-here \
    run.final.tsv 2> ajout.log || { cat ajout.log; fail "ajout"; }
python3 - "$R" <<'PY' || fail "columns"
import csv, sys
from types import SimpleNamespace
v = {r["variant_entree"]: r for r in csv.DictReader(open("run.final.mitoloop.tsv"), delimiter="\t")}
t = v["tata"]
assert "TATA_detruite" in t["MITO_arguments"] and "TATA_detruite[fort]" in t["MITO_evenements"], t
assert t["MITO_priorite"] == "1_fort", t
assert t["MITO_type_promoteur"].startswith("tissu(CpG=") and "TATA=oui" in t["MITO_type_promoteur"] \
    and "CAGE=6pb" in t["MITO_type_promoteur"], t["MITO_type_promoteur"]
i = v["inr"]
assert "Inr_detruit" in i["MITO_arguments"] and "Inr_detruit[fort]" in i["MITO_evenements"], i   # promoteur minimal
g4 = v["g4"]
assert "G4_detruit" in g4["MITO_arguments"] and "zone_pause" in g4["MITO_arguments"], g4
assert g4["MITO_evenements"] == "G4_detruit[modere]" and g4["MITO_priorite"] == "2_modere", g4
a = v["anchor"]
assert a["MITO_contact_promoteur"] == "GENEH<-promoteur:OTHERP", a["MITO_contact_promoteur"]
assert a["MITO_type_promoteur"].startswith("menage"), a["MITO_type_promoteur"]
# priority of a moderate event in a distal element, by promoter type
sys.path.insert(0, sys.argv[1])
from ajout_mitoloop import priority
opt = SimpleNamespace(spliceai_seuil=0.2, spliceai2_seuil=0.1, alphagenome_seuil=20.0)
base = {"MITO_arguments": "site_NRF1", "MITO_L1_classe": "not_estimable", "_classes": {"element"}}
ev = [("modere", "site_NRF1(IC1.5,d-1.2)")]
assert priority({}, dict(base, _type="tissu"), ev, opt)[1] == "1_fort"
assert priority({}, dict(base, _type="developpement"), ev, opt)[1] == "1_fort"
assert priority({}, dict(base, _type="menage"), ev, opt)[1] == "2_modere"
assert priority({}, dict(base, _type="tissu", _classes={"element", "promoter"}), ev, opt)[1] == "2_modere"
PY
echo "PASS test_promoter_elements"
