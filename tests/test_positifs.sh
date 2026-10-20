#!/usr/bin/env bash
# test_positifs.sh - test of positifs_mitoloop.py, offline
#
# Purpose      Mock MANE summary, ClinVar VCF, Ensembl answers (cache) and a
#              supplementary table; check which variants are kept, the merge
#              of sources, the unconverted list, and that preparer_entree.py
#              reads the output as chrom/pos/ref/alt with the MANE transcript.
#              Prints "PASS test_positifs".
# Usage        bash tests/test_positifs.sh

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="$HERE/../scripts"
R="$HERE/../reannot"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail() { echo "FAIL test_positifs: $*" >&2; exit 1; }
cd "$WORK"

printf 'ETFA\nRARS2\nBCS1L\nNDUFS4\n#comment\n' > genes.txt
printf 'symbol\tgtf_symbol\nOLDNAME\tNDUFS4\n' > renamed.tsv
python3 - <<'PY'
import gzip, json
h = "#NCBI_GeneID\tEnsembl_Gene\tHGNC_ID\tsymbol\tname\tRefSeq_nuc\tRefSeq_prot\tEnsembl_nuc\tEnsembl_prot\tMANE_status\tGRCh38_chr\tchr_start\tchr_end\tchr_strand\n"
rows = [("ETFA", "NM_000126.4", "ENST00000557943.6"), ("RARS2", "NM_020320.5", "ENST00000369536.10"),
        ("BCS1L", "NM_004328.5", "ENST00000359273.8"), ("NDUFS4", "NM_002495.4", "ENST00000296684.10"),
        ("TP53", "NM_000546.6", "ENST00000269305.9"),
        ("RRM2B", "NM_015713.5", "ENST00000251810.8")]
with gzip.open("mane.txt.gz", "wt") as fh:
    fh.write(h)
    for g, nm, en in rows:
        fh.write(f"1\tENSG\tHGNC\t{g}\tn\t{nm}\tNP\t{en}\tENSP\tMANE Select\tNC\t1\t2\t+\n")
cache = {
    "NM_000126.4:c.-85G>A": [["15", "76311000", "C", "T"]],
    "NM_020320.5:c.-2A>G": [["6", "88224000", "T", "C"]],
    "NM_004328.5:c.-122G>T": [],              # version inconnue -> sans version
    "NM_004328:c.-122G>T": [["2", "218658000", "G", "T"]],
    "NM_002495.4:c.*30A>G": [["5", "53560000", "A", "G"]],
    "NM_002495.4:c.-10C>T": [["5", "53460000", "C", "T"]],
    "GRCh38:8:102148058:G>T": [["8", "102148058", "G", "T"]],
}
json.dump(cache, open("cache.json", "w"))
v = ["##fileformat=VCFv4.1", "##fileDate=2026-10-01",
     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
def rec(c, p, i, r, a, sig, gene, mc, rev="criteria_provided,_single_submitter"):
    info = f"CLNSIG={sig};CLNREVSTAT={rev};GENEINFO={gene}:1;CLNDN=Leigh_syndrome"
    if mc:
        info += f";MC={mc}"
    return f"{c}\t{p}\t{i}\t{r}\t{a}\t.\t.\t{info}"
v += [rec(15, 76311000, 101, "C", "T", "Pathogenic", "ETFA", "SO:0001623|5_prime_UTR_variant"),      # = littérature
      rec(5, 53460100, 102, "G", "A", "Likely_pathogenic", "OLDNAME", "SO:0001624|3_prime_UTR_variant"),  # alias
      rec(5, 53460200, 103, "G", "A", "Pathogenic", "NDUFS4", "SO:0001583|missense_variant"),          # codant
      rec(5, 53460300, 104, "G", "A", "Pathogenic", "NDUFS4", "SO:0001623|5_prime_UTR_variant,SO:0001583|missense_variant"),
      rec(5, 53460400, 105, "G", "A", "Conflicting_classifications_of_pathogenicity", "NDUFS4", "SO:0001623|5_prime_UTR_variant"),
      rec(5, 53460500, 106, "G", "A", "Pathogenic", "TP53", "SO:0001623|5_prime_UTR_variant"),         # hors gènes
      rec(5, 53460600, 107, "G", "A", "Pathogenic", "NDUFS4", "SO:0001627|intron_variant"),            # intron
      rec(5, 53460700, 108, "G", "A", "Benign", "NDUFS4", "SO:0001623|5_prime_UTR_variant"),
      rec(5, 53460800, 109, "G", "A", "Pathogenic", "NDUFS4", ""),                                     # sans MC
      rec(5, 53460900, 110, "G", "<DEL>", "Pathogenic", "NDUFS4", "SO:0001623|5_prime_UTR_variant")]
with gzip.open("clinvar.vcf.gz", "wt") as fh:
    fh.write("\n".join(v) + "\n")
with open("supp.tsv", "w") as fh:
    fh.write("patient\tgene\tvariant\tremark\n")
    fh.write("P1\tNDUFS4\tc.*30A>G\tRNA outlier\n")
    fh.write("P2\tNDUFS4\tNM_002495.4:c.-10C>T\tMAE\n")
    fh.write("P3\tNDUFS4\tc.100A>G\tmissense\n")
    fh.write("P4\tTP53\tc.-5A>G\tnot mito\n")
    fh.write("P5\tRARS2\tc.-999A>G\tunknown to Ensembl\n")
PY

python3 "$S/positifs_mitoloop.py" --genes genes.txt --renamed renamed.tsv --mane mane.txt.gz \
    --clinvar clinvar.vcf.gz --cache cache.json --tables supp.tsv --datadir d \
    --out pos.tsv 2> log || { cat log; fail "run"; }
grep -q "variant_recoder indisponible" log || true

python3 - <<'PY' || { cat pos.tsv log; exit 1; }
import csv, sys
rows = list(csv.DictReader(open("pos.tsv"), delimiter="\t"))
k = {(r["chrom"], r["pos"]): r for r in rows}
def ok(c, m):
    if not c:
        sys.exit("FAIL test_positifs: " + m)
r = k[("chr8", "102148058")]
ok(r["gene"] == "RRM2B" and r["categorie"] == "intergenique" and r["notation_publiee"].startswith("GRCh38:"), "RRM2B")
ok(len(rows) == 7, f"7 variants attendus, {len(rows)}: {sorted(k)}")
e = k[("chr15", "76311000")]
ok("litterature" in e["source"] and "clinvar" in e["source"] and e["clinvar_id"] == "101", "fusion ETFA")
ok(e["transcrit"] == "ENST00000557943" and e["transcrit_suppose"] == "oui", "transcrit ETFA")
ok(("chr2", "218658000") in k, "BCS1L via NM sans version")
ok(("chr5", "53460100") in k and k[("chr5", "53460100")]["categorie"] == "3UTR", "alias OLDNAME")
ok(k[("chr5", "53560000")]["source"] == "supp.tsv", "table c.*30")
ok(k[("chr5", "53460000")]["transcrit_suppose"] == "non", "table NM donné")
for p in ("53460200", "53460300", "53460400", "53460500", "53460600", "53460700", "53460800"):
    ok(("chr5", p) not in k, f"{p} aurait dû être exclu")
bad = open("pos.non_convertis.tsv").read()
ok("c.-999A>G" in bad, "non converti RARS2")
PY

# avec --avec-introniques : le variant 107 entre
python3 "$S/positifs_mitoloop.py" --genes genes.txt --renamed renamed.tsv --mane mane.txt.gz \
    --clinvar clinvar.vcf.gz --cache cache.json --datadir d --out pos2.tsv \
    --avec-introniques 2> log2 || { cat log2; fail "run introns"; }
grep -q $'\t53460600\t' pos2.tsv || fail "intron avec --avec-introniques"

# lu par preparer_entree.py comme chrom/pos/ref/alt + transcrit
python3 "$R/preparer_entree.py" pos.tsv e.tsv > p.log || { cat p.log; fail "preparer"; }
grep -q "chrom=1, pos=2, ref=3, alt=4, transcrit 5" p.log || { cat p.log; fail "colonnes reconnues"; }
grep -q $'^ENST00000557943\tchr15:g.76311000C>T$' e.tsv || { cat e.tsv; fail "entrée reannot"; }
[ "$(wc -l < e.tsv)" -eq 7 ] || fail "7 lignes d'entrée"

# conversion GRCh37 -> GRCh38 et contrôle de la référence (réponses Ensembl simulées)
python3 - "$S" <<'PY' || fail "genomic"
import sys
sys.path[:0] = [sys.argv[1], sys.argv[1] + "/../reannot"]
import positifs_mitoloop as p
seq = {"G": "G"}
def fake(path):
    if path.startswith("/map/"):
        return {"mappings": [{"mapped": {"seq_region_name": "8", "start": 102148058, "end": 102148058, "strand": 1}}]}
    return {"seq": seq["G"]}
p.get_json = fake
assert p.genomic("GRCh37:8:103160286:G>T") == [["8", "102148058", "G", "T"]]
seq["G"] = "A"
try:
    p.genomic("GRCh37:8:103160286:G>T"); sys.exit("référence fausse acceptée")
except ValueError:
    pass
PY

echo "PASS test_positifs"
