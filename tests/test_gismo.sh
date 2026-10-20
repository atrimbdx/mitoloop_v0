#!/usr/bin/env bash
# test_gismo.sh - test of gismo_scores.py on a small archive laid out like the
#                 Zenodo record of GISMO (Liao et al. 2026); offline.
#                 Prints "PASS test_gismo".
# Checks       statistic of GenerateScore/01_GISMOscores.R on a hand case,
#              deciles (1 = fewest ortholog losses), one genome per species
#              (largest N50), GISMO-mis on the MANE Select transcript, and the
#              reading of the table by mitoloop_paquet.py (deciles kept).
# Usage        bash tests/test_gismo.sh

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail() { echo "FAIL test_gismo: $*" >&2; exit 1; }
cd "$WORK"

python3 - <<'PY'
import gzip, zipfile
P = "cliao5-GISMO-b109d22/"
sp = [("Sp_a", "aaa__A__HLa1", 1, 10, 100), ("Sp_b", "bbb__B__HLb1", 1, 20, 100),
      ("Sp_c", "ccc__C__HLc1", 2, 40, 100), ("Sp_d", "ddd__D__HLd1", 2, 80, 100)]
dist = "species\tcommon_name\tassembly_name\tcontig_N50_bp\tlin.num\tdiv.time\n"
for s, col, lin, t, n50 in sp:
    dist += f"{s}\tx\t{col.split('__')[-1]}\t{n50}\t{lin}\t{t}\n"
dist += "Sp_d\tx\tHLd0\t5\t2\t80\n"                     # second, poorer genome of Sp_d
# 20 genes : G00 has no loss, G19 is lost everywhere
genes = [f"ENSG{n:011d}" for n in range(20)]
rows = []
for i, g in enumerate(genes):
    lost = {0: [], 19: [0, 1, 2, 3]}.get(i, [j for j in range(4) if (i + j) % 5 == 0])
    rows.append(g + "\t" + "\t".join("one2zero" if j in lost else "one2one" for j in range(4)))
mat = "\t".join(c for _s, c, *_ in sp) + "\n" + "\n".join(rows) + "\n"
mis = ("genes\tmean.mis\tmean.syn\tmis.cal\tmean.comb\tdec\tgenename\n"
       "ENST00000000001\t1\t10\t0.1\t0.1\t1\tGENE0\n"
       "ENST00000000002\t5\t10\t0.5\t0.5\t7\tGENE0\n"     # MANE : kept
       "ENST00000000003\t2\t10\t0.2\t0.2\t3\tGENE1\n")
with zipfile.ZipFile("GISMO-GISMO.zip", "w") as z:
    z.writestr(P + "GenerateScore/inputs/gen.dist.withassembly.txt", dist)
    z.writestr(P + "GenerateScore/inputs/unmerged-species_combined_matrix_2023-07-26.tsv", mat)
    z.writestr(P + "partitioned-heritability/inputs/gismo-mis_consensus_genenames_01-08-2023.txt", mis)
with gzip.open("mini.gtf.gz", "wt") as fh:
    for i, g in enumerate(genes):
        fh.write(f'chr1\tX\tgene\t1\t2\t.\t+\t.\tgene_id "{g}.1"; gene_name "GENE{i}";\n')
    fh.write('chr1\tX\ttranscript\t1\t2\t.\t+\t.\tgene_id "ENSG00000000000.1"; '
             'transcript_id "ENST00000000002.1"; gene_name "GENE0"; tag "MANE_Select";\n')
PY

python3 "$HERE/../scripts/gismo_scores.py" --zip GISMO-GISMO.zip --gtf mini.gtf.gz \
    --out gismo_genes.tsv --iterations 200 > log || { cat log; fail "run"; }
grep -q "GISMO : 20 gènes, 4 espèces, 2 lignées" log || { cat log; fail "species / lineages"; }

python3 - "$HERE/../reannot" <<'PY' || exit 1
import sys
import pandas as pd
t = pd.read_csv("gismo_genes.tsv", sep="\t")
def ok(c, m):
    if not c:
        sys.exit("FAIL test_gismo: " + m)
g = t.set_index("gene")
# observed statistic: mean over species of (div.time / max) x loss
w = {0: 10 / 80, 1: 20 / 80, 2: 40 / 80, 3: 80 / 80}
lost5 = [j for j in range(4) if (5 + j) % 5 == 0]
ok(abs(g.at["GENE5", "gismo_obs"] - sum(w[j] for j in lost5) / 4) < 1e-6, "statistique observée")
ok(g.at["GENE0", "gismo"] == 0 and g.at["GENE0", "gismo_decile"] == 1, "gène sans perte = décile 1")
ok(g.at["GENE19", "gismo_decile"] == 10, "gène perdu partout = décile 10")
ok(sorted(t["gismo_decile"].dropna().astype(int).value_counts()) == [2] * 10, "déciles de 2 gènes")
ok(g.at["GENE0", "gismo_mis"] == 0.5 and g.at["GENE0", "gismo_mis_decile"] == 7, "GISMO-mis sur le MANE")
ok(g.at["GENE1", "gismo_mis_decile"] == 3, "GISMO-mis, transcrit unique")
sys.path.insert(0, sys.argv[1])
import mitoloop_paquet as mp
from pathlib import Path
r = mp.gismo_table(Path("gismo_genes.tsv"))
ok(r["ENSG00000000000"]["gismo_mis_decile"] == "7" and r["GENE0"]["gismo_decile"] == "1",
   "lecture par mitoloop_paquet (déciles gardés)")
ok(r["GENE19"]["gismo_decile"] == "10", "lecture par symbole")
PY

echo "PASS test_gismo"
