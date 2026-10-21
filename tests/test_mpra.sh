#!/usr/bin/env bash
# test_mpra.sh - mpra_variants.py: sources, conversion, evaluation package,
#                evaluation
#
# Purpose      Build a small GRCh38 genome with an NRF1 site, mock files of
#              every source (Kircher 2019 saturation mutagenesis, SuRE raQTL
#              lists and full table in hg19 with a chain + and -, MPRAVarDB
#              as a BED text, Sample 2019 5'UTR library, a hand table in
#              hg19), then build, regions, mitoloop_paquet.py +
#              ajout_mitoloop.py on every variant, evaluer and entree.
#              Offline; prints "PASS test_mpra".
# Method       at the most informative degenerate position of the site (A or
#              G tolerated), C and T are significant and G is not: the loss
#              of site separates them at the same position (AUC 1) while
#              phyloP, a position annotation, cannot (0.5).
# Usage        bash tests/test_mpra.sh

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
S="$HERE/../scripts"
R="$HERE/../reannot"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail() { echo "FAIL test_mpra: $*" >&2; exit 1; }
cd "$WORK"
mkdir -p dl

python3 - <<'PY'
import gzip, random
random.seed(3)
N = 20000
seq = [random.choice("ACGT") for _ in range(N)]
s0 = 4970                                  # NRF1 site, 0-based
seq[s0:s0 + 10] = list("GCGCATGCGC")
seq = "".join(seq)
open("genome.fa", "w").write(">chr1\n" + "\n".join(seq[i:i + 60] for i in range(0, N, 60)) + "\n")
open("genome.fa.fai", "w").write(f"chr1\t{N}\t6\t60\t61\n")
comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
b38 = lambda p: seq[p - 1]                 # 1-based

# Kircher: position 5 of the site (0-based s0+4, 'A'), the matrix tolerates A and G
p = s0 + 5
rows = ["Chrom\tPos\tRef\tAlt\tBarcodes\tDNA\tRNA\tCoefficient\tpValue\tElement\tRelease"]
for alt, eff, pv in (("C", -1.2, 1e-9), ("T", -1.0, 1e-8), ("G", -0.05, 0.6)):
    rows.append(f"1\t{p}\tA\t{alt}\t50\t100\t100\t{eff}\t{pv}\tTERT-HEK\tGRCh38")
    rows.append(f"1\t{p + 7}\tA\t{alt}\t50\t100\t100\t{eff}\t{pv}\tTERT-HEK\tGRCh37")   # other build: ignored
rows.append(f"1\t{p}\tA\t-\t50\t100\t100\t-0.8\t1e-6\tTERT-HEK\tGRCh38")              # deletion
rows.append(f"1\t3000\t{b38(3000)}\t{'A' if b38(3000) != 'A' else 'C'}\t9\t9\t9\t0.01\t0.9\tTERT-HEK\tGRCh38")
rows.append(f"1\t3001\t{b38(3001)}\t{'A' if b38(3001) != 'A' else 'C'}\t9\t9\t9\t0.3\t0.001\tTERT-HEK\tGRCh38")  # intermediate
with gzip.open("dl/kircher2019_elements.tsv.gz", "wt") as fh:
    fh.write("\n".join(rows) + "\n")

# chain hg19 -> hg38: hg19 chr1 [0,10000) -> chr1 [100,10100) (+);
# hg19 chr2 [0,1000) -> chr1 [15000,16000) (-)
with gzip.open("dl/hg19ToHg38.over.chain.gz", "wt") as fh:
    fh.write("chain 1000 chr1 249250621 + 0 10000 chr1 20000 + 100 10100 1\n10000\n\n")
    fh.write("chain 900 chr2 243199373 + 0 1000 chr1 20000 - 4000 5000 2\n1000\n\n")
alt_of = lambda b: "A" if b != "A" else "C"
def plus(x):                                # hg19 chr1 1-based -> ref, alt (hg19 = hg38 strand)
    b = b38(x + 100)
    return b, alt_of(b)
def minus(y):                               # hg19 chr2 1-based -> hg38 1-based, ref hg19
    q = 4000 + y - 1
    f = 20000 - q
    return f, comp[b38(f)]
H = ["chr", "SNP_ID", "SNPabspos", "ref.element.count", "alt.element.count",
     "k562.ref.mean", "k562.alt.mean", "k562.wilcox.p.value", "k562.wilcox.p.value.random",
     "hepg2.ref.mean", "hepg2.alt.mean", "hepg2.wilcox.p.value", "hepg2.wilcox.p.value.random",
     "ref", "alt"]
def sure(chrom, sid, x, ref, alt, k=(5.0, 20.0, 1e-5), h=(5.0, 5.0, 0.9), n=(50, 50)):
    return [chrom, sid, str(x), str(n[0]), str(n[1]), *map(str, k[:3]), "0.5", *map(str, h[:3]), "0.5", ref, alt]
k1 = sure("chr1", "rs1", 2000, *plus(2000))
k2 = sure("chr1", "rs2", 2100, *plus(2100))
f3, r3 = minus(500)
h3 = sure("chr2", "rs3", 500, r3, alt_of(r3), k=(5, 5, 0.9), h=(2.0, 30.0, 1e-6))
def sign(cell, rs):
    keep = [0, 1, 2, 3, 4] + ([5, 6, 7, 8] if cell == "k562" else [9, 10, 11, 12]) + [13, 14]
    with gzip.open(f"dl/sure_{cell}.sign.txt.gz", "wt") as fh:
        fh.write("\t".join(H[i] for i in keep) + "\n")
        for r in rs:
            fh.write("\t".join(r[i] for i in keep) + "\n")
sign("k562", [k1, k2])
sign("hepg2", [h3])
negs = [sure("chr1", f"rsn{i}", 3000 + 10 * i, *plus(3000 + 10 * i), k=(18.0, 22.0, 0.8), h=(25.0, 31.0, 0.9))
        for i in range(12)]
negs.append(sure("chr1", "rslow", 4000, *plus(4000), k=(6.0, 6.5, 0.8), n=(5, 50)))     # too few fragments
negs.append(sure("chr1", "rsoff", 4100, *plus(4100), k=(0.1, 0.1, 0.8), h=(0.1, 0.1, 0.9)))  # inactive element
negs.append(sure("chr1", "rssig", 4200, *plus(4200), k=(6.0, 6.5, 0.2), h=(0.1, 0.1, 0.3)))  # p <= 0.5
with gzip.open("dl/SuRE_SNP_table_LP190708.txt.gz", "wt") as fh:
    fh.write("\t".join(H) + "\n")
    for r in [k1, k2, h3] + negs:
        fh.write("\t".join(r) + "\n")

# MPRAVarDB as BED text (autoSql names); 0-based start
with open("dl/mpravardb.bb", "w") as fh:
    fh.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\t"
             "ref\talt\trsid\tdisease\tcellline\tdesc\tlog2fc\tpvalue\tfdr\tstudy\n")
    for i, (pos, fdr) in enumerate(((6000, 0.01), (6100, 0.8), (6200, 0.2))):
        b = b38(pos)
        fh.write(f"chr1\t{pos - 1}\t{pos}\trs9{i}\t0\t.\t{pos - 1}\t{pos}\t0\t{b}\t{alt_of(b)}\t"
                 f"rs9{i}\tx\tGM12878\t.\t0.5\t0.001\t{fdr}\tTewhey2016\n")
    fh.write(f"chr1\t6299\t6300\trs99\t0\t.\t6299\t6300\t0\tN\tA\trs99\tx\tK562\t.\t0.1\t0.5\t0.9\tX\n")

# Sample 2019: a human 5'UTR and one SNV
b = b38(7000)
with gzip.open("dl/GSM3130443_designed_library.csv.gz", "wt") as fh:
    fh.write(",utr,rl,id,info1,library,mother\n")
    fh.write("0,ACGTACGT,5.0,NM_1,.,human_utrs,\n")
    fh.write(f"1,ACGTTCGT,4.2,NC_000001.11:g.7000{b}>{alt_of(b)},rs77,snv,ACGTACGT\n")
    fh.write("2,ACGTTCGA,4.0,NC_000001.11:g.7001C>T,rs78,snv,UNKNOWN\n")

# hand table in hg19 (Griesemer-like)
with open("griesemer.tsv", "w") as fh:
    r, a = plus(8000)
    fh.write(f"chrom\tpos\tref\talt\teffet\tp\tsignif\tcellule\nchr1\t8000\t{r}\t{a}\t0.7\t0.001\t1\tHEK293FT\n")
    fh.write("chr1\t9999999\tA\tC\t0.1\t0.9\t0\tHEK293FT\n")                   # outside the chain

# MitoLoop package inputs: NRF1 matrix with a degenerate position 5 (A or G)
with open("jaspar.txt", "w") as fh:
    m = {"A": [0,0,0,0,10,0,0,0,0,0], "C": [0,20,0,20,0,0,0,20,0,20],
         "G": [20,0,20,0,10,0,20,0,20,0], "T": [0,0,0,0,0,20,0,0,0,0]}
    fh.write(">MA0506.3\tNRF1\n" + "".join(f"{b}  [ {' '.join(map(str, m[b]))} ]\n" for b in "ACGT"))
with gzip.open("remap_NRF1.bed.gz", "wt") as fh:
    fh.write(f"chr1\t{s0 - 100}\t{s0 + 100}\tNRF1:K-562,HepG2,GM12878,HEK293,MCF-7,A549,IMR-90\t0\t.\n")
with open("phylop.bedgraph", "w") as fh:
    fh.write(f"chr1\t0\t{s0}\t-1\nchr1\t{s0}\t{s0 + 10}\t3.5\nchr1\t{s0 + 10}\t{N}\t-1\n")
with open("regions_mito.bed", "w") as fh:
    fh.write("#chrom\tstart\tend\tgene\tclass\tsource\ttss\tstrand\tdetail\n")
    fh.write(f"chr1\t{s0 - 200}\t{s0 + 200}\tNDUFA1\tpromoter\tgencode\t{s0 + 30}\t+\t\n")
    fh.write("chr1\t5800\t6150\tNDUFA2\tgene_body\tgencode\t6000\t+\t\n")         # gene body: ignored
open("site.txt", "w").write(str(s0))
PY

python3 "$S/mpra_variants.py" build --datadir dl --fasta genome.fa --out mpra.tsv \
    --table-hg19 griesemer.tsv 2> build.log || { cat build.log; fail "build"; }
python3 - <<'PY' || { cat build.log; fail "build content"; }
import csv
s0 = int(open("site.txt").read())
rows = list(csv.DictReader(open("mpra.tsv"), delimiter="\t"))
src = {}
for r in rows:
    src.setdefault(r["source"], []).append(r)
k = src["kircher2019"]
assert len(k) == 6, k                                 # 3 SNV + deletion + 2 others, GRCh37 ignored
d = [r for r in k if len(r["ref"]) == 2][0]
seq = "".join(l.strip() for l in open("genome.fa") if not l.startswith(">"))
assert d["pos"] == str(s0 + 4) and d["ref"] == seq[s0 + 3] + "A" and d["alt"] == seq[s0 + 3], d
sig = {r["alt"]: r["signif"] for r in k if r["pos"] == str(s0 + 5)}
assert sig == {"C": "1", "T": "1", "G": "0"}, sig
assert {r["signif"] for r in k if r["pos"] == "3001"} == {"."}
assert all(r["cellule"] == "HEK293T" for r in k)
s = src["sure2019"]
pos = {r["element"]: r for r in s}
assert pos["rs1"]["pos"] == "2100" and pos["rs1"]["signif"] == "1" and pos["rs1"]["cellule"] == "K562", pos["rs1"]
assert abs(float(pos["rs1"]["effet"]) - 1.807) < 0.01, pos["rs1"]
r3 = pos["rs3"]                                       # chain -: position and complemented alleles
assert r3["pos"] == str(20000 - (4000 + 499)) and r3["ref"] == seq[int(r3["pos"]) - 1], r3
assert r3["cellule"] == "HepG2" and r3["signif"] == "1"
neg = [r for r in s if r["signif"] == "0"]
cells = sorted(r["cellule"] for r in neg)               # 3 per raQTL: 6 K562, 3 HepG2
assert cells == ["HepG2"] * 3 + ["K562"] * 6, [(r["element"], r["cellule"]) for r in neg]
assert not {"rslow", "rsoff", "rssig"} & {r["element"] for r in neg}
m = src["mpravardb:Tewhey2016"]
assert [r["signif"] for r in sorted(m, key=lambda r: int(r["pos"]))] == ["1", "0", "."], m
assert "mpravardb:X" not in src
sa = src["sample2019"]
assert len(sa) == 1 and sa[0]["pos"] == "7000" and abs(float(sa[0]["effet"]) + 0.8) < 1e-6, sa
g = src["griesemer"]
assert len(g) == 1 and g[0]["pos"] == "8100" and g[0]["signif"] == "1", g
PY
grep -q "sure : SNPabspos + 0" build.log || { cat build.log; fail "SuRE coordinate base"; }

python3 "$S/mpra_variants.py" regions --variants mpra.tsv --out eval 2> regions.log \
    || { cat regions.log; fail "regions"; }
python3 "$R/mitoloop_paquet.py" --genes eval/genes_eval.txt --l1 eval/l1_vide.tsv \
    --regions eval/regions_eval.bed --remap NRF1=remap_NRF1.bed.gz --jaspar jaspar.txt \
    --tfs NRF1 --phylop-bedgraph phylop.bedgraph --out pkg > paquet.log 2>&1 \
    || { cat paquet.log; fail "evaluation package"; }
python3 "$R/ajout_mitoloop.py" --data pkg --fasta genome.fa --samtools no-samtools-here \
    --min-ic 0.5 mpra.tsv 2> ajout.log || { cat ajout.log; fail "ajout"; }
grep -q "hors_regions_mito" ajout.log && { cat ajout.log; fail "every variant should be annotated"; }

python3 "$S/mpra_variants.py" evaluer --variants mpra.tsv --annote mpra.mitoloop.tsv \
    --out evaluation.tsv 2> eval.log || { cat eval.log; fail "evaluer"; }
python3 - <<'PY' || { cat eval.log; fail "evaluation content"; }
import csv
rows = list(csv.DictReader(open("evaluation.tsv"), delimiter="\t"))
e = {(r["source"], r["cellule"], r["annotation"]): r for r in rows}
r = e[("kircher2019", "toutes", "perte_site")]
assert r["classe"] == "C" and r["n_paires_meme_position"] == "2" and r["auc_meme_position"] == "1.000", r
r = e[("kircher2019", "toutes", "phyloP")]
assert r["classe"] == "B" and r["auc_meme_position"] == "0.500", r
assert e[("kircher2019", "toutes", "site_TF")]["auc_meme_position"] == "1.000"
assert e[("kircher2019", "toutes", "delta_site_signe")]["classe"] == "C"
assert ("sure2019", "K562", "phyloP") in e and ("sure2019", "toutes", "phyloP") in e
assert e[("sure2019", "toutes", "phyloP")]["n_sans_effet"] == "9"
PY

python3 "$S/mpra_variants.py" entree --variants mpra.tsv --regions regions_mito.bed \
    --out mpra_mito.tsv 2> entree.log || { cat entree.log; fail "entree"; }
python3 - <<'PY' || fail "entree content"
import csv
rows = list(csv.DictReader(open("mpra_mito.tsv"), delimiter="\t"))
assert rows and all(r["gene"] == "NDUFA1" and r["classe"] == "promoter" for r in rows), rows
assert {r["source"] for r in rows} == {"kircher2019"}, rows
PY
echo "PASS test_mpra"
