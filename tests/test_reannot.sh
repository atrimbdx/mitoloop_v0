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
grep -q "GISMO : gismo_mis <- gismo.csv, 20 lignes, colonne 'GISMO-mis', décile 'calculé', clé 'ensembl_gene_id'" paquet.log \
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

# input preparation: Excel with HGVS (NM_ prefix), VCF with indels and the
# MANE Select transcript taken from the GTF; one-command chain with stubs
python3 - <<'PY'
import gzip, openpyxl
wb = openpyxl.Workbook(); ws = wb.active
ws.append(["Patient", "Transcript", "Variant", "Note"])
ws.append(["P1", "ENST00000269305.9", "NM_000546.6:c.-29+2T>G", "x"])
ws.append(["P2", "ENST00000269305", "chr17:g.7676154G>C", ""])
ws.append(["P3", "", "", ""])
wb.save("variants.xlsx")
with gzip.open("ens.gtf.gz", "wt") as fh:
    fh.write('1\tens\ttranscript\t5001\t9000\t.\t+\t.\tgene_id "G1"; transcript_id "ENST00000000011"; tag "MANE_Select";\n')
    fh.write('1\tens\ttranscript\t10000\t13000\t.\t-\t.\tgene_id "G2"; transcript_id "ENST00000000022"; tag "MANE_Select";\n')
    fh.write('1\tens\ttranscript\t5001\t9500\t.\t+\t.\tgene_id "G1"; transcript_id "ENST00000000033";\n')
with open("in.vcf", "w") as fh:
    fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    fh.write("chr1\t4972\t.\tC\tA,T\t.\t.\t.\n")
    fh.write("chr1\t13100\t.\tATG\tA\t.\t.\t.\n")
    fh.write("chr1\t13200\t.\tA\tACC\t.\t.\t.\n")
    fh.write("chr1\t13300\t.\tA\t<DEL>\t.\t.\t.\n")
PY
python3 "$R/preparer_entree.py" variants.xlsx e1.tsv > p1.log || { cat p1.log; fail "xlsx"; }
[ "$(cat e1.tsv)" = "$(printf 'ENST00000269305.9\tc.-29+2T>G\nENST00000269305\tchr17:g.7676154G>C')" ] \
    || { cat e1.tsv p1.log; fail "xlsx content"; }
python3 "$R/preparer_entree.py" in.vcf e2.tsv --gtf ens.gtf.gz > p2.log || { cat p2.log; fail "vcf"; }
[ "$(cat e2.tsv)" = "$(printf 'ENST00000000011\tchr1:g.4972C>A\nENST00000000011\tchr1:g.4972C>T\nENST00000000022\tchr1:g.13101_13102del\nENST00000000022\tchr1:g.13200_13201insCC')" ] \
    || { cat e2.tsv p2.log; fail "vcf content"; }
# the chain, with stubs of the reannot scripts (the real ones need VEP)
mkdir -p base/scripts base/data/ref base/data/spliceai2
echo 115 > base/VEP_RELEASE
cp ens.gtf.gz base/data/ref/Homo_sapiens.GRCh38.115.gtf.gz
cp -r pkg base/data/mitoloop
cp "$R/preparer_entree.py" "$R/ajout_mitoloop.py" "$R/reannot_complet.sh" base/scripts/
cat > base/scripts/reannot_alphagenome.sh <<'SH'
#!/usr/bin/env bash
# stub: one final.tsv row per input line, the position of a g. notation
awk -F'\t' 'BEGIN { print "transcrit\tvariant_entree\tchrom\tpos\tref\talt\tgene\tstatut" }
  { split($2, a, /[:.>]/); g = substr($2, index($2, "g.") + 2); match(g, /[0-9]+/); p = substr(g, RSTART, RLENGTH)
    r = substr($2, length($2) - 2, 1); x = substr($2, length($2))
    print $1 "\t" $2 "\t" a[1] "\t" p "\t" r "\t" x "\t.\tok" }' "$1" > "$2.final.tsv"
SH
cat > base/scripts/ajout_spliceai2.sh <<'SH'
#!/usr/bin/env bash
awk -F'\t' -v OFS='\t' 'NR == 1 { print $0, "SpliceAI2_statut"; next } { print $0, "non_couvert" }' "$1" > "${1%.tsv}.spliceai2.tsv"
SH
touch base/data/spliceai2/spliceai2.tsv.gz
printf '#chrom\tpos\tref\talt\nchr1\t4972\tC\tA\n' > direct.tsv
REANNOT_BASE="$PWD/base" FASTA="$PWD/genome.fa" SAMTOOLS=no-samtools-here \
    bash base/scripts/reannot_complet.sh in.vcf run2 > chain.log 2>&1 || { cat chain.log; fail "chain"; }
[ -s run2.final.spliceai2.mitoloop.tsv ] || { cat chain.log; fail "chain output"; }
head -1 run2.final.spliceai2.mitoloop.tsv | grep -q "SpliceAI2_statut.*MITO_gene" || fail "chain columns"
grep -q "ENST00000000011.*GENEP:promoter" run2.final.spliceai2.mitoloop.tsv || { cat run2.final.spliceai2.mitoloop.tsv; fail "chain MitoLoop"; }
python3 "$R/ajout_mitoloop.py" --data pkg --fasta genome.fa --samtools no-samtools-here in.vcf 2> v.log \
    || { cat v.log; fail "VCF in ajout_mitoloop"; }
[ "$(awk -F'\t' 'NR > 1 { print $NF }' in.mitoloop.tsv | sort | uniq -c | tr -s ' ' | tr '\n' ';')" = " 4 annote;" ] \
    || { cat in.mitoloop.tsv; fail "VCF statuses"; }

# VM install: three updates on a reannot tree; the downloaded data must not
# move, the previous versions are archived, two archives kept
V=vm/reannot; mkdir -p $V/scripts $V/envs/vep/bin
for d in alphagenome annovar morfee ref spliceai spliceai2 vep_cache vep_plugins; do
    mkdir -p $V/data/$d; echo "$d" > $V/data/$d/donnees.txt; done
echo "pipeline" > $V/scripts/reannot_alphagenome.sh
( cd $V/data && find . -type f -exec sha256sum {} + | sort ) > avant.sha
for k in 1 2 3; do
    rm -rf transfert; cp -r pkg transfert
    cp "$R/ajout_mitoloop.py" "$R/preparer_entree.py" "$R/reannot_complet.sh" "$R/07_maj_mitoloop_vm.sh" \
       "$R/reannot_alphagenome.sh" "$R/run_morfee.R" "$R/verifier_morfee.sh" transfert/
    echo "mitoloop_v0 test$k" > transfert/VERSION
    ( cd transfert && sha256sum $(find . -type f ! -name SHA256SUMS | sed 's|^\./||' | sort) > SHA256SUMS )
    REANNOT_BASE="$PWD/$V" bash transfert/07_maj_mitoloop_vm.sh > install$k.log 2>&1 \
        || { cat install$k.log; fail "install $k"; }
done
( cd $V/data && find . -path ./mitoloop -prune -o -type f -exec sha256sum {} + | sort ) > apres.sha
diff avant.sha apres.sha > /dev/null || fail "the pipeline data changed"
[ "$(cat $V/scripts/reannot_alphagenome.sh)" = pipeline ] || fail "a pipeline script changed"
[ "$(cat $V/data/mitoloop/VERSION)" = "mitoloop_v0 test3" ] || fail "installed version"
[ "$(ls $V/archives | tr '\n' ' ')" = "mitoloop_mitoloop_v0_test1.tar.gz mitoloop_mitoloop_v0_test2.tar.gz " ] \
    || { ls $V/archives; fail "archives"; }
tar -tzf $V/archives/mitoloop_mitoloop_v0_test2.tar.gz | grep -q "scripts/ajout_mitoloop.py" || fail "scripts archived"
grep -q "Données du pipeline non modifiées : alphagenome annovar" install3.log || { cat install3.log; fail "install report"; }
# a locally modified reannot_alphagenome.sh is left alone
grep -q "version locale modifiée, non remplacée" install3.log || { cat install3.log; fail "AG script guard"; }
[ ! -e $V/archives/reannot_alphagenome.sh.origine ] || fail "AG original archived without replacement"
# the accelerated version already in place is recognised
cp "$R/reannot_alphagenome.sh" $V/scripts/reannot_alphagenome.sh
REANNOT_BASE="$PWD/$V" bash transfert/07_maj_mitoloop_vm.sh > install4.log 2>&1 || { cat install4.log; fail "install 4"; }
grep -q "déjà la version accélérée" install4.log || { cat install4.log; fail "AG script recognised"; }

echo "PASS test_reannot"
