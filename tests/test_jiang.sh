#!/usr/bin/env bash
# test_jiang.sh - end-to-end test of jiang2020_tables.py
#
# Purpose      Write the synthetic atlas of make_synthetic_l1.py as
#              supplementary spreadsheets laid out in two plausible ways,
#              convert them with jiang2020_tables.py and check that
#              l1_score.py gives on the converted tables the scores it gives
#              on the atlas itself. Offline; prints "PASS test_jiang".
# Inputs       none
# Outputs      none kept; everything is written to a temporary directory
# Method       layout 1: sheet names 'A. ...' to 'D. ...', a title row above
#              the header, protein keyed by an accession with a gene-name
#              column, samples named <tissue abbreviation>.<donor> (tissue
#              from Table S2 sheet A), RNA as log2(TPM + 1) with versioned
#              Ensembl ids. Layout 2: a readme sheet first, 'Table S2D'
#              names, GTEx-like sample ids that differ between RNA and
#              protein (paired by donor and tissue site, tissue from a GTEx
#              attributes file), RNA as ln(TPM + 1); the per-tissue protein
#              table under the tissue names of a GTEx-like TPM table.
#              Layout 3, as published (Tables S2 and S3 of the article):
#              readme and empty sheets first, title rows, a tissue row above
#              the GTEx sample ids, 'NA' strings, RNA and protein aliquots
#              with different SM- suffixes, RNA as log2(TPM) with negative
#              values.
# Usage        bash tests/test_jiang.sh

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/../scripts"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

fail() { echo "FAIL test_jiang: $*" >&2; exit 1; }

python3 "$HERE/make_synthetic_l1.py" --out "$WORK/s" > /dev/null
cd "$WORK/s"

python3 - <<'PY'
import math
import openpyxl
import pandas as pd

rna = pd.read_csv("atlas_rna.gct", sep="\t", skiprows=2)
prot = pd.read_csv("atlas_protein.tsv", sep="\t")
smp = pd.read_csv("atlas_samples.tsv", sep="\t")
tissues = sorted(smp.tissue.unique())
abbr = {t: f"T{chr(65 + i // 26)}{chr(65 + i % 26)}" for i, t in enumerate(tissues)}
counts = rna.iloc[:, 2:]
tpm = counts / counts.sum() * 1e6


def book(path, sheets):
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    for name, rows in sheets:
        ws = wb.create_sheet(name)
        for r in rows:
            ws.append(r)
    wb.save(path)


for lay in (1, 2, 3):
    d = f"lay{lay}"
    import os
    os.makedirs(d, exist_ok=True)
    if lay == 1:
        pname = {r.protein_sample: f"{abbr[r.tissue]}.{r.donor}" for r in smp.itertuples()}
        rname = {r.rna_sample: f"{abbr[r.tissue]}.{r.donor}" for r in smp.itertuples()}
        tab = [["Abbreviation", "Tissue", "GTEx tissue"]] + \
            [[abbr[t], t.replace("_", " ") + " tissue", t] for t in tissues]
        tab = [[a, t] for a, t in ((abbr[t], t) for t in tissues)]
        tab = [["Abbreviation", "Tissue"]] + tab
        names = ["A. Tissue names", "B. Raw abundances", "C. Normalized", "D. Cleaned log2"]
        ph = ["Protein", "Gene.name"] + [pname[s] for s in prot.columns[1:]]
        prows = [["Table S2D. cleaned relative protein abundances"], ph] + \
            [[f"P{i:05d}", g] + list(v) for i, (g, v) in
             enumerate(zip(prot.gene, prot.iloc[:, 1:].values.tolist()))]
        f = lambda x: math.log2(x + 1)
        rh = ["gene.id", "gene.name"] + [rname[s] for s in tpm.columns]
        rrows = [rh] + [[e, g] + [f(x) for x in v] for e, g, v in
                        zip(rna.Name, rna.Description, tpm.values.tolist())]
        book(f"{d}/1-s2.0-X-mmc3.xlsx", [(names[0], tab), (names[1], [["x"]]),
                                         (names[2], [["x"]]), (names[3], prows)])
        book(f"{d}/1-s2.0-X-mmc4.xlsx", [("A. RNA TPM", rrows), ("B. medians", [["x"]])])
    elif lay == 3:
        # as published: readme and empty sheets first, title rows, a tissue
        # row above the GTEx ids, 'NA' strings, RNA as log2(TPM) with a floor
        site = {t: f"{1000 + i:04d}" for i, t in enumerate(tissues)}
        tname = {t: t.replace("_", " - ").title() for t in tissues}
        pid = {r.protein_sample: f"GTEX-{r.donor}-{site[r.tissue]}-SM-9K{r.Index:03d}"
               for r in smp.itertuples()}
        rid = {r.rna_sample: f"GTEX-{r.donor}-{site[r.tissue]}-SM-5H{r.Index:03d}"
               for r in smp.itertuples()}
        tis_p = {r.protein_sample: tname[r.tissue] for r in smp.itertuples()}
        tis_r = {r.rna_sample: tname[r.tissue] for r in smp.itertuples()}
        tab = [["Table S2. A. Tissue name abbreviations"], [], [],
               ["1. Tissue name and abbreviation"], [],
               ["Tissue full name", "Tissue name abbreviations"]] + \
            [[tname[t], f"X{i}"] for i, t in enumerate(tissues)]
        P = prot.iloc[:, 1:].round(2).astype(object)
        P = P.mask(P.index.to_series().apply(lambda i: i % 7 == 0).values[:, None]
                   & (pd.DataFrame(range(P.shape[1])).T.values % 3 == 0), "NA")
        prows = [[None, None, "Table S2. D. Cleaned relative protein abundances"], [],
                 [None, None] + [tis_p[s] for s in prot.columns[1:]],
                 ["gene.id.full", "gene.id"] + [pid[s] for s in prot.columns[1:]]] + \
            [[f"ENSG{i:011d}.1", f"ENSG{i:011d}"] + list(v) for i, v in
             zip(range(len(prot)), P.values.tolist())]
        # the synthetic GTF ids are ENSG00000000000-...: map them by name
        gid = dict(zip(rna.Description, rna.Name))
        prows[4:] = [[gid.get(g, r[0]), gid.get(g, r[0]).split(".")[0]] + r[2:]
                     for g, r in zip(prot.gene, prows[4:])]
        rrows = [["Table S3. A. Cleaned and normalized RNA TPM in log scale"], [],
                 [None, None] + [tis_r[s] for s in tpm.columns],
                 ["gene.id.full", "gene.id"] + [rid[s] for s in tpm.columns]] + \
            [[e, e.split(".")[0]] + [round(math.log2(max(x, 1e-3)), 3) for x in v]
             for e, v in zip(rna.Name, tpm.values.tolist())]
        rrows[-1][2:] = [-9.97] * (len(rrows[-1]) - 2)     # a gene not detected
        book(f"{d}/1-s2.0-X-mmc3.xlsx", [("readme", [["Table S2. A."], ["Table S2. D."]]),
                                         ("Sheet1", []), ("A tisse names", tab),
                                         ("B protein raw abundance", [["x"]]),
                                         ("C protein normalized abundance", [["x"]]),
                                         ("D protein relative abundance", prows)])
        book(f"{d}/1-s2.0-X-mmc4.xlsx", [("ReadMe", [["Table S3. A."]]),
                                         ("A RNA cleaned log TPM", rrows)])
        with open(f"{d}/median_tpm.gct", "w") as fh:
            fh.write("#1.2\n1\t16\nName\tDescription\t" + "\t".join(tname[t] for t in tissues) + "\n")
            fh.write("ENSG00000000000.2\tMITOG1" + "\t1" * 16 + "\n")
    else:
        site = {t: f"{1000 + i:04d}" for i, t in enumerate(tissues)}
        pid = lambda r: f"GTEX-{r.donor}-{site[r.tissue]}-SM-P{r.Index:03d}"
        rid = lambda r: f"GTEX-{r.donor}-{site[r.tissue]}-SM-R{r.Index:03d}"
        pname = {r.protein_sample: pid(r) for r in smp.itertuples()}
        rname = {r.rna_sample: rid(r) for r in smp.itertuples()}
        with open(f"{d}/median_tpm.gct", "w") as fh:     # GTEx-like names
            fh.write("#1.2\n1\t17\nName\tDescription\t" + "\t".join(
                t.replace("_", " - ").title() for t in tissues) + "\tOther\n")
            fh.write("ENSG00000000000.2\tMITOG1" + "\t1" * 17 + "\n")
        with open(f"{d}/attributes.txt", "w") as fh:
            fh.write("SAMPID\tSMTS\tSMTSD\tSMTSISCH\n")
            for r in smp.itertuples():
                fh.write(f"{rid(r)}\tX\t{r.tissue}\t{r.SMTSISCH}\n")
        ph = ["Gene"] + [pname[s] for s in prot.columns[1:]]
        prows = [ph] + [[g] + list(v) for g, v in
                        zip(prot.gene, prot.iloc[:, 1:].values.tolist())]
        rh = ["Gene"] + [rname[s] for s in tpm.columns]
        rrows = [rh] + [[e.split(".")[0]] + [math.log(x + 1) for x in v]
                        for e, v in zip(rna.Name, tpm.values.tolist())]
        book(f"{d}/1-s2.0-X-mmc3.xlsx", [("README", [["tables"]]), ("Table S2A", [["x"]]),
                                         ("Table S2B", [["x"]]), ("Table S2C", [["x"]]),
                                         ("Table S2D", prows)])
        book(f"{d}/1-s2.0-X-mmc4.xlsx", [("README", [["tables"]]), ("Table S3A", rrows)])
PY

python3 "$SCRIPTS/l1_score.py" score --rna atlas_rna.gct --protein atlas_protein.tsv \
    --samples atlas_samples.tsv --gtf gencode.gtf.gz --genes genes.txt \
    --boot 50 --out ref > /dev/null
for lay in 1 2 3; do
    att=()
    [ "$lay" = 2 ] && att=(--attributes lay$lay/attributes.txt --tpm lay$lay/median_tpm.gct)
    [ "$lay" = 3 ] && att=(--tpm lay$lay/median_tpm.gct)
    python3 "$SCRIPTS/jiang2020_tables.py" --dir lay$lay --gtf gencode.gtf.gz \
        --out lay$lay/out --report lay$lay/rep "${att[@]}" > lay$lay/conv.log \
        || { cat lay$lay/conv.log; fail "conversion of layout $lay"; }
    grep -q "64 pairs in 16 tissues" lay$lay/conv.log || { cat lay$lay/conv.log; fail "pairs, layout $lay"; }
    [ -s lay$lay/rep/layout.txt ] || fail "layout.txt, layout $lay"
    python3 "$SCRIPTS/l1_score.py" score --rna lay$lay/out/jiang2020_rna_samples.tsv \
        --protein lay$lay/out/jiang2020_protein_samples.tsv --samples lay$lay/out/jiang2020_samples.tsv \
        --gtf gencode.gtf.gz --genes genes.txt --rna-units tpm --boot 50 \
        --out lay$lay/l1 > lay$lay/score.log
    grep -q "64 matched samples" lay$lay/score.log || fail "matched samples, layout $lay"
    python3 - "$lay" <<'PY' || fail "scores, layout $1"
import sys
import pandas as pd
lay = sys.argv[1]
ref = pd.read_csv("ref/l1.tsv", sep="\t").set_index("gene")
new = pd.read_csv(f"lay{lay}/l1/l1.tsv", sep="\t").set_index("gene")
t = pd.read_csv("truth.tsv", sep="\t").set_index("gene")
m = ref[["score"]].join(new[["score"]], rsuffix="_new").join(t).dropna()
assert len(m) > 50, len(m)
assert m.score.corr(m.score_new) > 0.95, m.score.corr(m.score_new)
med = m.groupby("class").score_new.median()
assert med["transcriptional"] > 0.9 and abs(med["buffered"]) < 0.3, med
PY
done
grep -q "log base 2 taken" lay1/rep/conversion.txt || fail "log base 2"
grep -q "log base 2.72 taken" lay2/rep/conversion.txt || fail "log base e"
grep -q "pairs (same GTEx tissue sample)" lay2/rep/conversion.txt || fail "pairing by GTEx tissue sample"
grep -q "tissue row above the header" lay3/rep/conversion.txt || fail "tissue row, layout 3"
grep -q "log base 2 taken, values log(TPM + 0)" lay3/rep/conversion.txt || fail "log2 TPM without offset"
[ "$(head -1 lay3/out/jiang2020_protein.tsv | tr '\t' '\n' | wc -l)" = 17 ] || fail "protein per tissue, layout 3"
python3 - <<'PY' || fail "protein per tissue"
import pandas as pd
t = pd.read_csv("lay2/out/jiang2020_protein.tsv", sep="\t", index_col=0)
s = pd.read_csv("atlas_protein.tsv", sep="\t", index_col=0)
assert list(t.columns[:2]) == ["Tissue - 00", "Tissue - 01"], list(t.columns)
assert t.shape[1] == 16 and abs(t.loc["MITOG1", "Tissue - 00"]
                                - s.loc["MITOG1", [f"PA00_{i}" for i in range(4)]].median()) < 1e-3
PY
[ ! -e lay1/out/jiang2020_protein.tsv ] || fail "per-tissue table without --tpm"
echo "PASS test_jiang"
