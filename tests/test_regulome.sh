#!/usr/bin/env bash
# test_regulome.sh - end-to-end test of regulome.py
#
# Purpose      Build MitoCarta-like sheets, a small Gene Ontology, GTEx-like
#              expression, protein, eQTL and Vg tables and per-gene features
#              with tests/make_synthetic_regulome.py, run the three
#              subcommands and check them against the planted truth.
#              Offline; prints "PASS test_regulome" on success and exits
#              non-zero on the first failed assertion.
# Inputs       none
# Outputs      none kept; everything is written to a temporary directory
# Method       1. pathways: MitoPathways read as full hierarchies and as leaf
#                 names placed through sheet C, a gene with two sections, a
#                 former symbol found through the MitoCarta synonyms, the
#                 central dogma section reported as the mtDNA layer, and the
#                 section of the 11 genes without MitoPathway inferred from
#                 GO (NOT and ND annotations ignored); leave-one-out
#                 accuracy reported.
#              2. modes: the transcription score separates the genes planted
#                 with strong and with no transcriptional control.
#              3. compare: rE2G elements and loops follow the score, promoter
#                 and methylation features do not.
# Usage        bash tests/test_regulome.sh
# References   see the header of scripts/regulome.py

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/../scripts"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

fail() { echo "FAIL test_regulome: $*" >&2; exit 1; }

python3 "$HERE/make_synthetic_regulome.py" --out "$WORK/s"
cd "$WORK/s"

# ------------------------------------------------------------- 1. pathways
python3 "$SCRIPTS/regulome.py" pathways --mitocarta mitocarta_A.tsv \
    --mitopathways mitocarta_C.tsv --genes genes.txt --gaf goa.gaf.gz \
    --obo go.obo --out r > pathways.log
grep -q "section right for 100.0% of 69" pathways.log || fail "leave-one-out"
python3 - <<'PY' || fail "pathways"
import pandas as pd
p = pd.read_csv("r/pathways.tsv", sep="\t", dtype=str,
                keep_default_na=False).set_index("gene")
t = pd.read_csv("truth.tsv", sep="\t", dtype=str).set_index("gene")
assert len(p) == 60, len(p)
assert "OLDMITOG3" not in p.index and p.loc["MITOG3", "in_mitocarta"] == "1"
for g in p.index:
    secs = p.loc[g, "sections"].split("|")
    assert t.loc[g, "section"] in secs, (g, secs)
    want = "go_inferred" if int(g[5:]) >= 50 else "mitocarta"
    assert p.loc[g, "source"] == want, (g, p.loc[g, "source"])
assert p.loc["MITOG7", "layers"] == "Mitochondrial dynamics and surveillance|mtDNA"
assert p.loc["MITOG4", "layers"] == "mtDNA"
assert p.loc["MITOG4", "pathways"] == \
    "Mitochondrial central dogma > mtDNA maintenance > mtDNA replication"
assert p.loc["MITOG6", "pathways"].endswith("> TIM23 presequence pathway")
assert p.loc["MITOG3", "pathways"] == \
    "Mitochondrial dynamics and surveillance > Fusion"
for g in ("MITOG52", "MITOG53"):       # NOT and ND rows ignored
    assert p.loc[g, "sections"].count("|") == 0, p.loc[g, "sections"]
assert float(p.loc["MITOG55", "expected_precision"]) > 0.9
m = pd.read_csv("r/layer_matrix.tsv", sep="\t").set_index("gene")
assert m.loc["MITOG52", "mtDNA"] == 1 and m["mtDNA"].sum() == 16
PY

# ---------------------------------------------------------------- 2. modes
python3 "$SCRIPTS/regulome.py" modes --gtf gencode.gtf.gz \
    --genes r/pathways.tsv --tpm tpm.gct --protein protein.tsv \
    --tissue-map tissue_map.tsv --eqtl eqtl.txt.gz --vg vg.tsv \
    --out r > modes.log
grep -q "19 tissues" modes.log || fail "cultured cells not excluded"
python3 - <<'PY' || fail "modes"
import pandas as pd
m = pd.read_csv("r/modes.tsv", sep="\t").set_index("gene")
assert (m.set == "mito").sum() == 60
assert "LNC1" not in m.index
assert (m.loc[[f"MITOG{i}" for i in range(1, 61)], "n_transcription"] == 5).all()
hi = m.loc[[f"MITOG{i}" for i in range(1, 21)], "score_transcription"]
lo = m.loc[[f"MITOG{i}" for i in range(41, 61)], "score_transcription"]
assert hi.min() > lo.max(), (hi.min(), lo.max())
assert m.loc["MITOG1", "mrna_protein_rho"] > 0.8
assert m.loc["MITOG41", "eqtl_signals"] == 0
PY

# -------------------------------------------------------------- 3. compare
python3 "$SCRIPTS/regulome.py" compare --modes r/modes.tsv \
    --pathways r/pathways.tsv --promoters promoter_features.tsv \
    --methylation methylation_by_group.tsv --re2g re2g_per_gene.tsv \
    --loops loops.tsv --assoc assoc_loops.tsv --perms 300 --seed 2 \
    --out r > compare.log
python3 - <<'PY' || fail "compare"
import pandas as pd
c = pd.read_csv("r/compare.tsv", sep="\t")
c["name"] = c.source + "." + c.feature
c = c.set_index("name")
for f in ("re2g.n_union", "loops.n_loops"):
    assert c.loc[f, "partial_rho"] > 0.6 and c.loc[f, "q"] < 1e-6, f
    assert c.loc[f, "p_within_sections"] < 0.01, f
    assert c.loc[f, "high"] > c.loc[f, "low"], f
for f in ("promoter.gc", "promoter.cpg_oe", "methylation.mean_beta",
          "re2g.n_ccre_window"):
    assert c.loc[f, "q"] > 0.1, f
assert c.loc["promoter.motif_NRF1", "kind"] == "binary"
g = pd.read_csv("r/groups.tsv", sep="\t").set_index("gene")
assert (g.group == "high").sum() == 20 and (g.group == "low").sum() == 20
assert g.loc["MITOG52", "section"] == "Mitochondrial central dogma"
PY

echo "PASS test_regulome"
