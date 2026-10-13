#!/usr/bin/env bash
# test_v4c.sh - end-to-end test of v4c_pipeline.py and v4c_loops.R
#
# Purpose
#   Build a synthetic set of Hi-C libraries (make_synthetic_v4c.py: 8
#   ENCODE-like terms in 5 GTEx-like tissues, 2 libraries each, 20
#   mitochondrial and 8 control viewpoints on two chromosomes, one planted
#   loop per mitochondrial viewpoint, loop strength rising with mtDNA copy
#   number for half of them) and run the loop pipeline on it: extract from
#   .cool, .mcool and .hic, select, export, scan, call, partner extraction and
#   export, reciprocal, assoc and tissue.
#
# Inputs
#   None. The .hic part needs Juicer Tools (jar path in the environment
#   variable JUICER_TOOLS), java and the hicstraw Python module; without them
#   it is skipped with a message. R needs data.table, limma, peakC, isotone.
#
# Outputs
#   Prints PASS v4c on success; exits non-zero on the
#   first failed assertion. Works in a temporary directory removed on exit.
#
# Method
#   Assertions: npz keys and viewpoint sets; identical counts from a library
#   read as .cool, .mcool (through a custom location column) and .hic;
#   libraries.tsv and export columns; every planted loop present in loops.tsv
#   within 20 kb of its planted partner, with a non-empty `called` column and
#   no loop elsewhere; in assoc, loops of the trend genes gaining more contact
#   per doubling of mtCN than those of the flat genes; in tissue, loop
#   strength increasing with mtCN across tissues. The deseq subcommand needs
#   DESeq2, which may be absent: the test then checks that it stops with a
#   message naming the package. Option checks: an unknown option and a
#   missing --wsize are rejected with an explanation.
#
# Usage
#   JUICER_TOOLS=tools/juicer_tools.jar bash tests/test_v4c.sh
#
# References
#   Durand et al. 2016. Juicer provides a one-click system for analyzing
#       loop-resolution Hi-C experiments. Cell Systems.
#   Abdennur and Mirny 2020. Cooler: scalable storage for Hi-C data and other
#       genomically labeled arrays. Bioinformatics.
#   Geeven et al. 2018. peakC: a flexible, non-parametric peak calling
#       package for 4C and Capture-C data. Nucleic Acids Research.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="python3 $ROOT/scripts/v4c_pipeline.py"
RS="Rscript $ROOT/scripts/v4c_loops.R"
CORES=2

W="$(mktemp -d)"
trap 'rm -rf "$W"' EXIT

fail() { echo "FAIL: $*" >&2; exit 1; }

echo "== synthetic data"
python3 "$HERE/make_synthetic_v4c.py" --outdir "$W" --hic-pairs SYN03 \
    --mcool SYN05 >/dev/null
cd "$W"

GENES="--gtf gencode.gtf.gz --genes mito.txt --controls controls.txt"

echo "== extract (.cool, location column 'path')"
$PY extract --tissues tissues.tsv $GENES --outdir raw >extract.log
[ "$(ls raw/*.npz | wc -l)" -eq 16 ] || fail "extract did not write 16 npz"
python3 - <<'EOF' || fail "npz content"
import numpy as np
d = np.load("raw/SYN01.npz")
need = {"counts", "raw_total", "names", "chrom", "pos", "set", "offsets",
        "accession", "term", "gtex", "mtcn", "resolution", "anchor_half", "via"}
assert need <= set(d.files), need - set(d.files)
assert set(d["set"]) == {"mito", "control"}, set(d["set"])
assert str(d["via"]) == "path"
assert d["counts"].shape == (28, 201)
assert np.isfinite(d["counts"]).any(1).all()
EOF

echo "== extract (.mcool, --location-columns local_file)"
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","term_name","gtex_tissue","mtcn","local_file"}
            $1=="SYN05"{print $1,$2,$3,$4,"libs/SYN05.mcool"}' tissues.tsv > t_mcool.tsv
$PY extract --tissues t_mcool.tsv $GENES --outdir raw_mcool \
    --location-columns local_file >/dev/null
python3 - <<'EOF' || fail ".mcool and .cool differ"
import numpy as np
a, b = np.load("raw/SYN05.npz"), np.load("raw_mcool/SYN05.npz")
assert np.array_equal(a["counts"], b["counts"], equal_nan=True)
assert str(b["via"]) == "local_file"
EOF

echo "== extract (.hic)"
if [ -n "${JUICER_TOOLS:-}" ] && [ -f "${JUICER_TOOLS}" ] && command -v java >/dev/null \
        && python3 -c "import hicstraw" 2>/dev/null; then
    java -Xmx2g -jar "$JUICER_TOOLS" pre -r 10000 -n libs/SYN03.pairs \
        libs/SYN03.hic chrom.sizes >pre.log 2>&1 || fail "juicer_tools pre"
    awk -F'\t' 'BEGIN{OFS="\t"} NR==1 || $1=="SYN03"' tissues.tsv \
        | sed 's#libs/SYN03.cool#libs/SYN03.hic#' > t_hic.tsv
    $PY extract --tissues t_hic.tsv $GENES --outdir raw_hic >/dev/null
    python3 - <<'EOF' || fail ".hic and .cool differ"
import numpy as np
a, b = np.load("raw/SYN03.npz"), np.load("raw_hic/SYN03.npz")
assert np.array_equal(a["counts"], b["counts"], equal_nan=True)
assert np.nansum(b["counts"]) > 0
EOF
    echo "   .hic and .cool profiles identical"
else
    echo "   skipped: needs JUICER_TOOLS, java and hicstraw"
fi

echo "== select, export"
$PY select --tissues tissues.tsv --indir raw --out libraries.tsv \
    --k-sd 3 --min-keep 4 >/dev/null
python3 - <<'EOF' || fail "libraries.tsv"
import pandas as pd
d = pd.read_csv("libraries.tsv", sep="\t")
assert {"accession", "term", "gtex", "mtcn", "depth", "kept"} <= set(d.columns)
assert len(d) == 16 and d["kept"].sum() == 16
assert d["term"].nunique() == 8 and d["gtex"].nunique() == 5
EOF
$PY export --indir raw --libraries libraries.tsv --set mito \
    --outdir export/mito >/dev/null
[ -f export/mito/decay.tsv ] || fail "decay.tsv missing"
python3 - <<'EOF' || fail "export table"
import pandas as pd
d = pd.read_csv("export/mito/SYN01.tsv.gz", sep="\t")
assert list(d.columns[:3]) == ["name", "chrom", "pos"]
assert d.columns[3] == "-1000000" and d.columns[-1] == "1000000"
assert len(d) == 20 and d["name"].str.startswith("MITO").all()
EOF

echo "== v4c_loops.R option handling"
if $RS call --libraries libraries.tsv --export export/mito --wsize 3 \
        --alpha_fdrr 0.1 >opt.log 2>&1; then fail "unknown option accepted"; fi
grep -q "unknown option --alpha-fdrr" opt.log || fail "unknown option message"
if $RS call --libraries libraries.tsv --export export/mito >opt.log 2>&1; then
    fail "missing --wsize accepted"; fi
grep -q "wsize is required" opt.log && grep -q "v4c_loops.R scan" opt.log \
    || fail "missing --wsize message"
$RS --help >opt.log 2>&1 || fail "--help"
grep -q -- "--min-ref-loops" opt.log || fail "--help does not list options"

echo "== scan"
$RS scan --libraries libraries.tsv --export export/mito --wsizes 1,3,5 \
    --cores $CORES --out loops >scan.log
grep -q "proposed window" scan.log || fail "scan proposal"
[ "$(tail -n +2 loops/scan.tsv | wc -l)" -eq 3 ] || fail "scan.tsv rows"

echo "== call"
$RS call --libraries libraries.tsv --export export/mito --wsize 3 \
    --cores $CORES --out loops >call.log
[ -s loops/peaks_mito.tsv ] && [ -s loops/partners.tsv ] || fail "call outputs"
head -1 loops/partners.tsv | grep -q $'^name\tchrom\tpos\tset$' \
    || fail "partners.tsv columns"

echo "== extract and export partners"
$PY extract --tissues tissues.tsv --anchors loops/partners.tsv \
    --outdir raw_partners >/dev/null
$PY export --indir raw_partners --libraries libraries.tsv --set partner \
    --outdir export/partners >/dev/null
[ "$(ls export/partners/*.tsv.gz | wc -l)" -eq 16 ] || fail "partner export"

echo "== reciprocal"
$RS reciprocal --libraries libraries.tsv --export export/partners --wsize 3 \
    --cores $CORES --out loops >recip.log
python3 - <<'EOF' || fail "loops.tsv against planted loops"
import pandas as pd
loops = pd.read_csv("loops/loops.tsv", sep="\t")
cols = ["gene", "chrom", "vp_pos", "summit", "start", "end", "n_terms",
        "called", "loop_id"]
assert list(loops.columns) == cols, list(loops.columns)
planted = pd.read_csv("planted.tsv", sep="\t")
for r in planted.itertuples():
    hit = loops[(loops.gene == r.gene) & ((loops.summit - r.partner).abs() <= 20000)]
    assert len(hit) == 1, f"planted loop of {r.gene} not recovered"
    called = hit.called.iloc[0].split("|")
    assert len(called) == hit.n_terms.iloc[0] >= 2, r.gene
assert len(loops) == len(planted), f"{len(loops) - len(planted)} extra loops"
print(f"   {len(planted)} planted loops recovered, median "
      f"{loops.n_terms.median():g} terms each")
EOF

echo "== assoc"
$RS assoc --libraries libraries.tsv --export export/mito --wsize 3 \
    --cores $CORES --out loops --perms 200 --sets sets.tsv --min-set 5 \
    --min-loops 10 --fdr 0.1 >assoc.log
for f in assoc_loops.tsv assoc_global.tsv assoc_null.tsv assoc_sets.tsv \
         assoc_volcano.svg; do
    [ -s "loops/$f" ] || fail "assoc: $f missing"
done
python3 - <<'EOF' || fail "assoc on the planted trend"
import pandas as pd
s = pd.read_csv("loops/assoc_sets.tsv", sep="\t").set_index("set")
tr = s.loc["synthetic::trend", "mean_log2fc_per_doubling"]
fl = s.loc["synthetic::flat", "mean_log2fc_per_doubling"]
assert tr > 0.1 and abs(fl) < 0.1, (tr, fl)
g = pd.read_csv("loops/assoc_global.tsv", sep="\t").set_index("statistic")
assert g.loc["mean_t", "observed"] > 0
print(f"   log2 change per doubling of mtCN: trend {tr:+.3f}, flat {fl:+.3f}")
EOF

echo "== tissue"
$RS tissue --libraries libraries.tsv --export export/mito --wsize 3 \
    --cores $CORES --out loops --min-loops 10 >tissue.log
for f in tissue_strength.tsv tissue_strength_libraries.tsv tissue_strength.svg; do
    [ -s "loops/$f" ] || fail "tissue: $f missing"
done
python3 - <<'EOF' || fail "tissue strength against mtCN"
import pandas as pd
from scipy.stats import spearmanr
t = pd.read_csv("loops/tissue_strength.tsv", sep="\t")
assert len(t) == 5
rho = spearmanr(t["mtcn"], t["strength"])[0]
assert rho > 0.8, rho
print(f"   Spearman rho over 5 tissues {rho:+.2f}")
EOF

echo "== deseq"
if Rscript -e 'quit(status = !requireNamespace("DESeq2", quietly = TRUE))' \
        >/dev/null 2>&1; then
    echo "   DESeq2 installed; deseq is not exercised by this test"
else
    if $RS deseq --libraries libraries.tsv --export export/mito --wsize 3 \
            --out loops >deseq.log 2>&1; then fail "deseq ran without DESeq2"; fi
    grep -q "needs the Bioconductor package(s) DESeq2" deseq.log \
        || fail "deseq message: $(cat deseq.log)"
    echo "   DESeq2 absent: deseq stops with a clear message"
fi

echo "PASS v4c"
