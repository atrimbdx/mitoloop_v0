#!/usr/bin/env bash
# test_sv.sh - structural variants: preparer_entree.py --sv, the gnomAD-SV
#              file of the package (mitoloop_paquet.py) and ajout_mitoloop_sv.py
#
# Purpose      Two genes with exons, promoter, distal element and loop anchor;
#              an Excel table type / chromosome / start / end without header
#              (the format of the WGS SV exports); a mock gnomAD-SV BED and a
#              SNV table of the same patient. Checks the VCF written, each
#              effect (deletion of exons and TSS, of an element, internal and
#              whole-gene duplication, inversion through the gene or
#              separating it from its element), the gnomAD-SV downgrade, the
#              second allele and the statuses. Offline; prints "PASS test_sv".
# Usage        bash tests/test_sv.sh

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
R="$HERE/../reannot"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail() { echo "FAIL test_sv: $*" >&2; exit 1; }
cd "$WORK"

python3 - <<'PY'
import gzip, openpyxl
open("genes.txt", "w").write("GENEA\nGENEB\n")
open("l1.tsv", "w").write("gene\tgene_id\tset\tstatus\tscore\tci_low\tci_high\tcomplex\t"
                          "protein_coupling\trna_coupling\tl1_class\n")
with open("regions.bed", "w") as fh:
    fh.write("#chrom\tstart\tend\tgene\tclass\tsource\ttss\tstrand\tdetail\n")
    fh.write("chr1\t9000\t10500\tGENEA\tpromoter\tgencode\t10000\t+\t\n")
    fh.write("chr1\t30000\t30500\tGENEA\telement\tre2g\t10000\t+\t\n")
    fh.write("chr1\t40000\t41000\tGENEA\tanchor\tloops\t10000\t+\t\n")
    fh.write("chr1\t69500\t71000\tGENEB\tpromoter\tgencode\t69999\t-\t\n")
with gzip.open("ens.gtf.gz", "wt") as fh:
    def tx(g, tid, s, e, st, exons):
        a = f'gene_id "ENSG_{g}"; transcript_id "{tid}"; gene_name "{g}"; gene_type "protein_coding";'
        fh.write(f"1\tens\tgene\t{s}\t{e}\t.\t{st}\t.\t{a}\n")
        fh.write(f'1\tens\ttranscript\t{s}\t{e}\t.\t{st}\t.\t{a} tag "MANE_Select";\n')
        for es, ee in exons:
            fh.write(f"1\tens\texon\t{es}\t{ee}\t.\t{st}\t.\t{a}\n")
    tx("GENEA", "ENST00000000001", 10001, 20000, "+", [(10001, 10500), (14001, 14200), (19001, 20000)])
    tx("GENEB", "ENST00000000002", 60001, 70000, "-", [(60001, 61000), (65001, 65500), (69001, 70000)])
with gzip.open("gnomad_sv.bed.gz", "wt") as fh:
    fh.write("#chrom\tstart\tend\tname\tsvtype\tFILTER\tAF\n")
    fh.write("chr1\t64990\t66000\tsv1\tDEL\tPASS\t0.05\n")      # common deletion of GENEB exon 2
    fh.write("chr1\t9499\t12000\tsv2\tDEL\tUNRESOLVED\t0.3\n")  # not PASS: dropped
wb = openpyxl.Workbook(); ws = wb.active
for r in (["DEL", 1, 9500, 12000],       # exon 1, TSS and promoter of GENEA
          ["DEL", 1, 30101, 30400],      # GENEA element
          ["DUP", 1, 13001, 15000],      # inside GENEA, exon 2 copied
          ["DUP", 1, 5001, 25000],       # whole GENEA
          ["INV", 1, 15001, 35000],      # breakpoint inside GENEA
          ["inversion", 1, 25001, 35000],  # separates GENEA from its element
          ["DEL", 1, 500001, 501000],    # far from every gene
          ["DEL", 1, 65001, 66000]):     # GENEB exon 2, common in gnomAD-SV
    ws.append(r)
wb.save("sv.xlsx")
with open("snv.tsv", "w") as fh:
    fh.write("transcrit\tvariant_entree\tchrom\tpos\tref\talt\tgene\tMITO_gene\tMITO_priorite\n")
    fh.write("ENST00000000001\tc.200A>G\tchr1\t10200\tA\tG\tGENEA\t.\t.\n")
    fh.write("ENST00000000001\tc.-30C>T\tchr1\t9970\tC\tT\t.\tGENEA:promoter\t1_fort\n")
PY

python3 "$R/preparer_entree.py" sv.xlsx sv.vcf --sv > prep.log || { cat prep.log; fail "preparer --sv"; }
grep -q "colonnes type 1, chromosome 2, début 3, fin 4 ; 8 SV" prep.log || { cat prep.log; fail "SV columns"; }
grep -q $'chr1\t9500\tligne1_DEL\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=12000;SVLEN=-2501' sv.vcf || { cat sv.vcf; fail "VCF line"; }
grep -q "SVTYPE=INV;END=35000" sv.vcf || fail "inversion synonym"

python3 "$R/mitoloop_paquet.py" --genes genes.txt --l1 l1.tsv --regions regions.bed \
    --gtf ens.gtf.gz --gnomad-sv gnomad_sv.bed.gz --out pkg > paquet.log 2>&1 \
    || { cat paquet.log; fail "paquet"; }
grep -q "gnomAD-SV : 1 SV PASS" paquet.log || { cat paquet.log; fail "gnomAD-SV PASS"; }

python3 "$R/ajout_mitoloop_sv.py" --data pkg --gtf ens.gtf.gz --snv snv.tsv sv.vcf 2> ajout.log \
    || { cat ajout.log; fail "ajout_mitoloop_sv"; }
python3 - <<'PY' || { cat sv.mitoloop.tsv; fail "columns"; }
import csv
v = {(r["pos"], r["svtype"]): r for r in csv.DictReader(open("sv.mitoloop.tsv"), delimiter="\t")}
d = v[("9500", "DEL")]
assert d["MITO_SV_genes"] == "GENEA:exons 1/3,TSS,cassure_dans_le_gene,promoter", d["MITO_SV_genes"]
assert d["MITO_evenements"] == "GENEA:DEL:exons 1/3,TSS[fort]" and d["MITO_priorite"] == "1_fort", d
assert d["MITO_second_allele"] == "GENEA:chr1:10200A>G,chr1:9970C>T[1_fort]", d["MITO_second_allele"]
e = v[("30101", "DEL")]
assert e["MITO_evenements"] == "GENEA:DEL:element[modere]" and e["MITO_priorite"] == "2_modere", e
i = v[("13001", "DUP")]
assert i["MITO_evenements"] == "GENEA:DUP:interne exons 1/3[fort]", i["MITO_evenements"]
w = v[("5001", "DUP")]
assert w["MITO_SV_genes"].startswith("GENEA:gene_entier") and w["MITO_evenements"] == "." \
    and w["MITO_priorite"] == "3_contexte", w
c = v[("15001", "INV")]
assert c["MITO_evenements"] == "GENEA:INV:coupe_le_gene[fort]", c["MITO_evenements"]
s = v[("25001", "INV")]
assert s["MITO_evenements"] == "GENEA:INV:separe_element[modere]", s["MITO_evenements"]
assert v[("500001", "DEL")]["MITO_statut"] == "hors_genes_mito"
b = v[("65001", "DEL")]
assert b["MITO_SV_gnomAD"].startswith("AF=0.05(DEL,RO="), b["MITO_SV_gnomAD"]
assert "GENEB:DEL:exons 1/3[fort]" in b["MITO_evenements"] and "gnomAD-SV_frequent(AF=0.05)" in b["MITO_evenements"]
assert b["MITO_priorite"] == "2_modere", b                    # one level down
assert b["MITO_L1_classe"] == "GENEB:not_estimable" and b["MITO_second_allele"] == ".", b
PY
# one command on a reannot tree, without AnnotSV
mkdir -p base/scripts base/data/ref
echo 116 > base/VEP_RELEASE
cp ens.gtf.gz base/data/ref/Homo_sapiens.GRCh38.116.gtf.gz
cp -r pkg base/data/mitoloop
cp "$R/preparer_entree.py" "$R/ajout_mitoloop.py" "$R/ajout_mitoloop_sv.py" "$R/reannot_sv.sh" base/scripts/
REANNOT_BASE="$PWD/base" bash base/scripts/reannot_sv.sh sv.xlsx run1 --snv snv.tsv > chain.log 2>&1 \
    || { cat chain.log; fail "reannot_sv.sh"; }
grep -q "AnnotSV non installé" chain.log || { cat chain.log; fail "AnnotSV message"; }
[ "$(awk -F'\t' 'NR > 1 { print $NF }' run1.sv.mitoloop.tsv | sort | uniq -c | tr -s ' ' | tr '\n' ';')" \
    = " 7 annote; 1 hors_genes_mito;" ] || { cat run1.sv.mitoloop.tsv; fail "chain statuses"; }
echo "PASS test_sv"
