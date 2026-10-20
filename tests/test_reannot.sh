#!/usr/bin/env bash
# test_reannot.sh - the MitoLoop module of the reannot pipeline
#
# Purpose      Build a small genome with two genes (one per strand), an NRF1
#              site planted 30 bp upstream of the + gene, EPDnew, CAGE,
#              ReMap and ChromHMM files, a JASPAR matrix and an l1.tsv; make
#              the package with mitoloop_paquet.py, then annotate a reannot
#              table with ajout_mitoloop.py and check every column.
#              Offline; prints "PASS test_reannot".
# Inputs       none
# Outputs      none kept; everything is written to a temporary directory
# Method       variants: a SNV at an informative position of the NRF1 site
#              (site_NRF1, negative delta), a SNV at the TSS of the - gene
#              (promoteur_minimal, distance 0, oriented), a SNV outside every
#              region (hors_regions_mito), an indel (no motif), a SNV whose
#              REF does not match the genome (REF_discordante), a row
#              without position (NA).
# Usage        bash tests/test_reannot.sh

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
R="$HERE/../reannot"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail() { echo "FAIL test_reannot: $*" >&2; exit 1; }
cd "$WORK"

python3 - <<'PY'
import gzip, random
random.seed(5)
seq = [random.choice("ACGT") for _ in range(20000)]
site = "GCGCATGCGC"                         # NRF1, palindrome
s0 = 5000 - 30                              # 0-based start, + gene TSS at 5000 (0-based)
seq[s0:s0 + 10] = list(site)
seq = "".join(seq)
open("genome.fa", "w").write(">chr1\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
open("genes.txt", "w").write("GENEP\nGENEM\nGENEX\n")
with open("regions.bed", "w") as fh:
    fh.write("#chrom\tstart\tend\tgene\tclass\tsource\ttss\tstrand\tdetail\n")
    fh.write("chr1\t4000\t5500\tGENEP\tpromoter\tgencode\t5000\t+\t\n")
    fh.write("chr1\t5500\t9000\tGENEP\tgene_body\tgencode\t5000\t+\t\n")
    fh.write("chr1\t12500\t14000\tGENEM\tpromoter\tgencode\t13000\t-\t\n")
    fh.write("chr1\t13600\t13800\tGENEM\telement\tre2g\t13000\t-\t\n")
with open("l1.tsv", "w") as fh:
    fh.write("gene\tgene_id\tset\tstatus\tscore\tci_low\tci_high\tcomplex\tprotein_coupling\trna_coupling\tl1_class\n")
    fh.write("GENEP\tENSG00000000001\tmito\testimable\t0.9\t0.7\t0.99\tCI subunits\t0.95\t0.5\ttransmitted\n")
    fh.write("GENEM\tENSG00000000002\tmito\tflat\tNA\tNA\tNA\t\t\t\tnot_estimable\n")
open("patient_genes.tsv", "w").write("gene\tfibroblast_r\nGENEP\t0.41\n")
# EPDnew: spaces, -49..+10 around the TSS
open("epd.bed", "w").write("chr1 4951 5011 GENEP_1 900 +\nchr1 12990 13050 GENEM_1 900 -\n")
with gzip.open("cage.bed.gz", "wt") as fh:
    fh.write("chr1\t4990\t5010\tp1\t1\t+\n")
with gzip.open("remap_NRF1.bed.gz", "wt") as fh:
    fh.write("chr1\t4900\t5100\tNRF1:K-562,HepG2,GM12878,HEK293,MCF-7,A549,IMR-90\t0\t.\n")
    fh.write("chr1\t12900\t13100\tNRF1:K-562\t0\t.\n")
    fh.write("chr1\t30000\t30100\tNRF1:K-562\t0\t.\n")          # outside: dropped
with gzip.open("e108.bed.gz", "wt") as fh:
    fh.write("chr1\t0\t4500\t15_Quies\nchr1\t4500\t5600\t1_TssA\nchr1\t5600\t20000\t15_Quies\n")
with open("jaspar.txt", "w") as fh:
    m = {"A": [0,0,0,0,20,0,0,0,0,0], "C": [0,20,0,20,0,0,0,20,0,20],
         "G": [20,0,20,0,0,0,20,0,20,0], "T": [0,0,0,0,0,20,0,0,0,0]}
    fh.write(">MA0506.3\tNRF1\n" + "".join(f"{b}  [ {' '.join(map(str, m[b]))} ]\n" for b in "ACGT"))
    fh.write(">MA0079.5\tSP1\nA  [ 1 1 ]\nC  [ 1 1 ]\nG  [ 9 9 ]\nT  [ 1 1 ]\n")
# gnomAD v4.1 constraint (MANE Select row first), GISMO (20 genes, CSV)
with open("gnomad.tsv", "w") as fh:
    fh.write("gene\tgene_id\ttranscript\tcanonical\tmane_select\tlof.obs\tlof.exp\tlof.oe\t"
             "lof.pLI\tlof.oe_ci.upper\tmis.oe\tmis.z_score\n")
    fh.write("GENEP\t1234\tNM_1\ttrue\ttrue\t1\t20\t0.05\t0.99\t0.12\t0.5\t3.1\n")
    fh.write("GENEP\tENSG00000000001\tENST9\tfalse\tfalse\t9\t10\t0.9\t0\t1.5\t1\t0\n")
    fh.write("GENEP\tENSG00000000001\tENST1\ttrue\ttrue\t2\t25\t0.08\t0.97\t0.21\t0.61\t2.8\n")
    fh.write("GENEM\tENSG00000000002\tENST2\ttrue\ttrue\t30\t28\t1.07\t0\t1.4\t1.02\t-0.2\n")
with open("gismo.csv", "w") as fh:
    fh.write("ensembl_gene_id,hgnc_symbol,GISMO,GISMO-mis\n")
    fh.write("ENSG00000000001,GENEP,0.001,0.05\n")
    for i in range(3, 21):
        fh.write(f"ENSG{i:011d},G{i},{0.01 * i:.3f},{0.1 * i:.3f}\n")
    fh.write("ENSG00000000002,GENEM,0.9,0.95\n")
# phyloP: 3.5 over the NRF1 site and its flanks, -1 elsewhere
import pyBigWig
bw = pyBigWig.open("phylop.bw", "w")
bw.addHeader([("chr1", 20000)])
bw.addEntries(["chr1", "chr1", "chr1"], [0, s0 - 10, s0 + 20], ends=[s0 - 10, s0 + 20, 20000],
              values=[-1.0, 3.5, -1.0])
bw.close()
# reannot table: site position 2 of the motif (C at 0-based s0+1) -> 1-based s0+2
ref_tss_m = seq[13000]
with open("run.final.tsv", "w") as fh:
    fh.write("transcrit\tvariant_entree\tchrom\tpos\tref\talt\tgene\tstatut\n")
    fh.write(f"ENST1\tv1\tchr1\t{s0 + 2}\tC\tA\tGENEP\tok\n")
    fh.write(f"ENST2\tv2\tchr1\t13001\t{ref_tss_m}\t{'A' if ref_tss_m != 'A' else 'C'}\tGENEM\tok\n")
    fh.write(f"ENST3\tv3\tchr1\t11000\t{seq[10999]}\tA\tGENEX\tok\n")
    fh.write(f"ENST4\tv4\tchr1\t{s0 + 2}\tC\tCT\tGENEP\tok\n")
    bad = "A" if seq[s0 + 6] != "A" else "C"
    fh.write(f"ENST5\tv5\tchr1\t{s0 + 7}\t{bad}\tG\tGENEP\tok\n")
    fh.write("ENST6\tv6\tNA\tNA\tNA\tNA\tNA\tnon_converti\n")
PY

python3 "$R/mitoloop_paquet.py" --genes genes.txt --l1 l1.tsv \
    --fibroblast patient_genes.tsv --regions regions.bed --epd epd.bed \
    --cage cage.bed.gz --remap NRF1=remap_NRF1.bed.gz --chromhmm muscle=e108.bed.gz \
    --jaspar jaspar.txt --tfs NRF1,YY1 --gnomad gnomad.tsv --gismo gismo.csv \
    --phylop phylop.bw --out pkg > paquet.log || { cat paquet.log; fail "paquet"; }
grep -q "gnomAD : 2 gènes ; GISMO : 2 gènes" paquet.log || { cat paquet.log; fail "gene constraints"; }
grep -q "GISMO : gismo.csv, 20 gènes, colonnes 'GISMO', 'GISMO-mis', clé 'ensembl_gene_id'" paquet.log \
    || { cat paquet.log; fail "GISMO table"; }
[ -s pkg/annotations/phylop.bedgraph.gz.tbi ] || { cat paquet.log; fail "phyloP index"; }
grep -q "régions : 3 (element, promoter)" paquet.log || { cat paquet.log; fail "regions"; }
grep -q "remap_NRF1 : 2 intervalles" paquet.log || { cat paquet.log; fail "ReMap cut to the regions"; }
grep -q "motifs JASPAR : NRF1 (absents : YY1)" paquet.log || { cat paquet.log; fail "JASPAR subset"; }
[ "$(zcat pkg/annotations/epd.bed.gz | head -1 | cut -f2,4)" = "$(printf '5000\t+')" ] || fail "EPD TSS (+)"
[ "$(zcat pkg/annotations/epd.bed.gz | tail -1 | cut -f2,4)" = "$(printf '13000\t-')" ] || fail "EPD TSS (-)"

python3 "$R/ajout_mitoloop.py" --data pkg --fasta genome.fa --samtools no-samtools-here \
    run.final.tsv 2> ajout.log || { cat ajout.log; fail "ajout"; }
[ -s run.final.mitoloop.tsv ] || fail "output"
python3 - <<'PY' || fail "columns"
import csv
rows = list(csv.DictReader(open("run.final.mitoloop.tsv"), delimiter="\t"))
v = {r["variant_entree"]: r for r in rows}
a = v["v1"]
assert a["MITO_statut"] == "annote" and a["MITO_gene"] == "GENEP:promoter", a
assert a["MITO_L1_classe"] == "transmitted" and a["MITO_L1_IC"] == "[0.7;0.99]", a
assert "sous-unité tenue par CI subunits" in a["MITO_L1_note"] and "r_fibroblastes=0.41" in a["MITO_L1_note"], a
assert a["MITO_dist_TSS"] == "-29", a["MITO_dist_TSS"]
assert a["MITO_fixation"] == "NRF1(7)", a
assert a["MITO_motifs"].startswith("NRF1:p"), a["MITO_motifs"]
assert "/d-" in a["MITO_motifs"], a["MITO_motifs"]
args = a["MITO_arguments"].split(",")
for x in ("L1_transmis", "promoteur_minimal", "fixation_NRF1", "site_NRF1", "chromatine_active_muscle"):
    assert x in args, (x, args)
assert "pic_CAGE" not in args
assert a["MITO_n_arguments"] == str(len(args))
assert a["MITO_phyloP"] == "3.50" and "conservation_mammiferes" in args, a
assert a["MITO_gnomAD"] == "LoF=2/25;oe_LoF=0.08;LOEUF=0.21;pLI=0.97;oe_mis=0.61;z_mis=2.8", a["MITO_gnomAD"]
assert a["MITO_GISMO"] == "GISMO=0.001(d1);GISMO-mis=0.05(d1)", a["MITO_GISMO"]
assert v["v2"]["MITO_phyloP"] == "-1.00" and "conservation" not in v["v2"]["MITO_arguments"]
assert v["v2"]["MITO_GISMO"] == "GISMO=0.9(d10);GISMO-mis=0.95(d5)", v["v2"]["MITO_GISMO"]
b = v["v2"]
assert b["MITO_gene"] == "GENEM:promoter" and b["MITO_dist_TSS"] == "0", b
assert b["MITO_L1_classe"] == "not_estimable" and b["MITO_fixation"] == "NRF1(1)", b
assert "fixation_NRF1" not in b["MITO_arguments"] and "promoteur_minimal" in b["MITO_arguments"], b
assert b["MITO_chromatine"] == "muscle:15_Quies", b
assert v["v3"]["MITO_statut"] == "hors_regions_mito"
assert v["v4"]["MITO_motifs"] == "." and v["v4"]["MITO_statut"] == "annote"
assert v["v5"]["MITO_motifs"] == "REF_discordante", v["v5"]
assert v["v6"]["MITO_statut"] == "NA"
assert list(rows[0])[:8] == ["transcrit", "variant_entree", "chrom", "pos", "ref", "alt", "gene", "statut"]
PY
# validation on fine-mapped eQTL: in 30 signals of GENEP, the causal variant
# breaks the NRF1 site, its companions lie elsewhere in the promoter
python3 - <<'PY'
import gzip, random
random.seed(2)
seq = "".join(l.strip() for l in open("genome.fa") if not l.startswith(">"))
comp = {"A": "G", "G": "A", "C": "T", "T": "C"}
with gzip.open("dapg.txt.gz", "wt") as fh:
    fh.write("tissue_id\tgene_id\tcluster_id\tcluster_pip\tvariant_id\tvariant_pip\n")
    for k in range(30):
        p = 4970 + 1 + random.choice([0, 2, 6, 8])          # informative NRF1 positions
        fh.write(f"T{k}\tENSG00000000001.5\t1\t1\tchr1_{p}_{seq[p-1]}_{comp[seq[p-1]]}_b38\t0.9\n")
        for _ in range(3):
            q = random.randint(4100, 4800)
            fh.write(f"T{k}\tENSG00000000001.5\t1\t1\tchr1_{q}_{seq[q-1]}_{comp[seq[q-1]]}_b38\t0.01\n")
        fh.write(f"T{k}\tENSG00000000009.1\t1\t1\tchr1_4975_A_G_b38\t0.9\n")   # another gene
PY
python3 "$R/validation_eqtl.py" --data pkg --dapg dapg.txt.gz --fasta genome.fa \
    --samtools no-samtools-here --out val.tsv > val.log || { cat val.log; fail "validation"; }
grep -q "30 avec au moins un variant causal" val.log || { cat val.log; fail "signals"; }
python3 - <<'PY' || { cat val.log; fail "validation statistics"; }
import csv
v = {r["argument"]: r for r in csv.DictReader(open("val.tsv"), delimiter="\t")}
assert float(v["site_NRF1"]["causal"]) == 1 and float(v["site_NRF1"]["background"]) == 0
assert float(v["site_NRF1"]["p"]) < 1e-6 and float(v["site_NRF1"]["z"]) > 0
assert float(v["n_arguments"]["difference"]) > 0
PY

echo "PASS test_reannot"
