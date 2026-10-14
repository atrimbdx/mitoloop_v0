#!/usr/bin/env bash
# test_wgs.sh - end-to-end test of vcf_annotate.py
#
# Purpose      Build two trios and a singleton with planted regulatory,
#              coding and structural variants (tests/make_synthetic_vcf.py),
#              run `regions` then `annotate`, and check every planted case
#              against its expected tier and category. Offline; prints
#              "PASS test_wgs" on success and exits non-zero on the first
#              failed assertion.
# Inputs       none
# Outputs      none kept; everything is written to a temporary directory
# Method       1. regions: MANE transcript, promoter, UTRs and padded coding
#                 exons on both strands, elements above the prevalence
#                 threshold, loop anchors called in muscle, capture regions
#                 with and without a gene column, alias resolution (TAZ).
#              2. annotate: trans and cis second hits from parental origin,
#                 trans from a PS block, homozygous, hemizygous, de novo in
#                 a dominant gene, X-linked female carrier, multi-allelic
#                 record, frequencies from INFO, CSQ and a gnomAD VCF, CADD
#                 (plain), ReMM (tabix when pysam is installed) and an INFO
#                 score, FILTER and GQ filters, inversion separating an
#                 element from its TSS, deletion of an element.
#              3. --require-score keeps only the variants reaching a
#                 threshold.
# Usage        bash tests/test_wgs.sh
# References   see the header of scripts/vcf_annotate.py

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/../scripts"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

fail() { echo "FAIL test_wgs: $*" >&2; exit 1; }

python3 "$HERE/make_synthetic_vcf.py" --out "$WORK/s"
cd "$WORK/s"

# ------------------------------------------------------------- 1. regions
python3 "$SCRIPTS/vcf_annotate.py" regions --gtf gencode.gtf.gz \
    --genes genes.txt --elements re2g_elements.tsv --min-prevalence 3 \
    --loops loops.tsv --called-in "psoas" --bed chic=chic.bed \
    --out regions.bed > regions.log
python3 - <<'PY' || fail "regions"
import pandas as pd
r = pd.read_csv("regions.bed", sep="\t")
r.columns = [c.lstrip("#") for c in r.columns]
def has(gene, cls, s, e):
    return ((r.gene == gene) & (r["class"] == cls) & (r.start == s)
            & (r.end == e)).any()
assert has("GENEA", "promoter", 999_000, 1_000_500)
assert has("GENEB", "promoter", 2_999_500, 3_001_000)          # - strand
assert has("GENEA", "utr5", 1_000_000, 1_000_200)
assert has("GENEB", "utr5", 2_999_800, 3_000_000)
assert has("GENEB", "utr3", 2_989_000, 2_989_700)
assert has("GENEA", "cds", 1_004_990, 1_005_310)                # +/- 10 bp
assert (r[(r.gene == "GENEA")]["class"] == "cds").sum() == 3    # MANE only
assert has("GENEA", "element", 1_049_800, 1_050_200)
assert not ((r.start == 1_079_800)).any()                       # prevalence
assert has("GENEA", "anchor", 795_000, 805_000)
assert not (r.start == 695_000).any()                           # liver only
assert has("GENEB", "chic", 3_119_500, 3_120_500)
assert has("GENEC", "chic", 5_200_000, 5_201_000)               # by distance
assert has("TAFAZZIN", "element", 50_029_800, 50_030_200)       # TAZ alias
assert (r.tss[r.gene == "GENEB"] == 3_000_000).all()
PY

# ------------------------------------------------------------ 2. annotate
python3 "$SCRIPTS/vcf_annotate.py" annotate --vcf snv.vcf.gz --vcf sv.vcf \
    --ped cohort.ped --regions regions.bed --inheritance panelapp.tsv \
    --gnomad gnomad.vcf --score CADD=cadd.tsv.gz:PHRED \
    --score ReMM=remm.tsv.gz --info-score AlphaGenome=AG \
    --min-score CADD=15 --min-score ReMM=0.9 --pass-only \
    --out annot > annotate.log
for f in candidates.tsv genes.tsv cohort_genes.tsv log.txt; do
    [ -s "annot/$f" ] || fail "missing annot/$f"
done
grep -q "FILTER 1, GQ 1" annot/log.txt || fail "FILTER and GQ filters"
python3 - <<'PY' || fail "candidates"
import pandas as pd
c = pd.read_csv("annot/candidates.tsv", sep="\t", dtype=str,
                keep_default_na=False)
truth = pd.read_csv("truth.tsv", sep="\t", dtype=str)
reg = c[c.role == "regulatory"]
for _, t in truth.iterrows():
    r = reg[(reg.proband == t.proband) & (reg.gene == t.gene)
            & (reg.pos == t.pos)]
    assert len(r) == 1, (t.to_dict(), len(r))
    r = r.iloc[0]
    assert (r.tier, r.category) == (t.tier, t.category), (t.to_dict(),
                                                          r.tier, r.category)
assert len(reg) == len(truth), sorted(zip(reg.proband, reg.pos))

def row(proband, pos):
    return c[(c.proband == proband) & (c.pos == str(pos))].iloc[0]
# frequencies: Number=A INFO, CSQ, gnomAD VCF; removed when common
assert row("P3", 999801).af == "0.0001"
assert row("P1", 1005101).af_source == "gnomADg_AF"
assert row("P3", 3120001).af_source == "gnomad"
for p, pos in (("P1", 999701), ("P1", 4999901), ("P3", 4999951),
               ("P2", 700001), ("P1", 1080001), ("P1", 5040101)):
    assert not ((c.proband == p) & (c.pos == str(pos))).any(), (p, pos)
# scores and flags
r = row("P1", 1050001)
assert (r.CADD, r.ReMM, r.score_flags) == ("18", "0.95", "CADD,ReMM")
assert r.cohort_carriers == "2" and r.origin == "maternal"
assert row("P2", 5040001).AlphaGenome == "0.8"
assert row("P3", 3120001).CADD == "12"
# second hits and phase
assert row("P1", 1005101).category == "coding_partner"
assert "(coding,trans)" in row("P1", 1050001).partners
assert row("P3", 3120001).phase_block == "chr1:100#1"
assert row("P2", 3000301).origin == "paternal"
# structural variants
r = row("P2", 1020001)
assert r.svtype == "INV" and "separated" in r.effect and r.origin == "de_novo"
r = row("P3", 5039001)
assert r.svtype == "DEL" and r.effect.startswith("deleted")
assert "chr1:1020001-1070000:INV" in row("P2", 800001).partners
PY
python3 - <<'PY' || fail "gene tables"
import pandas as pd
g = pd.read_csv("annot/cohort_genes.tsv", sep="\t").set_index("gene")
assert g.loc["GENEA", "probands_tier1"] == 2
assert g.loc["GENEA", "probands_tier1_2"] == 3
assert g.loc["GENEC", "inheritance"] == "AD"
PY

# ------------------------------------------------------- 3. --require-score
python3 "$SCRIPTS/vcf_annotate.py" annotate --vcf snv.vcf.gz \
    --ped cohort.ped --regions regions.bed --inheritance panelapp.tsv \
    --score CADD=cadd.tsv.gz --score ReMM=remm.tsv.gz \
    --min-score CADD=15 --min-score ReMM=0.9 --require-score \
    --samples P1,P2,P3 --out annot_req > req.log
python3 - <<'PY' || fail "--require-score"
import pandas as pd
c = pd.read_csv("annot_req/candidates.tsv", sep="\t", dtype=str,
                keep_default_na=False)
reg = c[c.role == "regulatory"]
assert (reg.score_flags != "").all()
assert set(reg.pos) == {"800001", "1050001", "3120001", "5040001"}, set(reg.pos)
PY

echo "PASS test_wgs"
