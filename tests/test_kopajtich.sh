#!/usr/bin/env bash
# test_kopajtich.sh - end-to-end test of kopajtich_tables.py
#
# Purpose      Write the synthetic patient cohort of make_synthetic_l1.py in a
#              layout like that of Kopajtich et al. (raw_data.zip), convert
#              it with kopajtich_tables.py and check that `l1_score.py
#              validate` finds on the converted tables the patient events it
#              finds on the cohort itself. Offline; prints "PASS
#              test_kopajtich".
# Inputs       none
# Outputs      none kept; everything is written to a temporary directory
# Method       counts keyed by versioned Ensembl ids with a symbol column;
#              proteomics keyed by protein groups ('MITOG1;MITOG1B') with
#              'NA' strings and an extra text column; an annotation table
#              with proteome id, RNA id (other names than the matrices'
#              columns' case and separators), TMT batch and instrument.
# Usage        bash tests/test_kopajtich.sh

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/../scripts"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

fail() { echo "FAIL test_kopajtich: $*" >&2; exit 1; }

python3 "$HERE/make_synthetic_l1.py" --out "$WORK/s" > /dev/null
cd "$WORK/s"
mkdir -p raw

python3 - <<'PY'
import pandas as pd
rna = pd.read_csv("patients_rna.tsv", sep="\t")
prot = pd.read_csv("patients_protein.tsv", sep="\t")
smp = pd.read_csv("patients_samples.tsv", sep="\t")
gtf = pd.read_csv("atlas_rna.gct", sep="\t", skiprows=2)[["Name", "Description"]]
ver = dict(zip(gtf.Name.str.split(".").str[0], gtf.Name))
sym = dict(zip(gtf.Name.str.split(".").str[0], gtf.Description))
rid = {s: f"MUC{1000 + i}" for i, s in enumerate(smp.rna_sample)}
pid = {s: f"P{i:03d}_TMT" for i, s in enumerate(smp.protein_sample)}
c = rna.rename(columns=rid)
c.insert(1, "gene_name", [sym.get(g, "") for g in rna.gene_id])
c["gene_id"] = [ver.get(g, g) for g in rna.gene_id]
c.rename(columns={"gene_id": "geneID"}).to_csv("raw/raw_counts.tsv", sep="\t", index=False)
p = prot.rename(columns=pid).astype(object)
p.iloc[::5, 3] = "NA"
p.insert(0, "Protein.IDs", [f"Q{i:05d}" for i in range(len(p))])
p["gene"] = [f"{g};{g}B" for g in prot.gene]
p.rename(columns={"gene": "geneID"}).to_csv("raw/proteomics_not_normalized.tsv", sep="\t", index=False)
pd.DataFrame({"SAMPLE_ID": [pid[s].lower() for s in smp.protein_sample],
              "RNA_ID": [rid[s].replace("MUC", "muc-") for s in smp.rna_sample],
              "BATCH_RUN": smp.batch, "INSTRUMENT": "QE"}).to_csv(
    "raw/proteomics_annotation.tsv", sep="\t", index=False)
PY

python3 "$SCRIPTS/l1_score.py" score --rna atlas_rna.gct --protein atlas_protein.tsv \
    --samples atlas_samples.tsv --gtf gencode.gtf.gz --genes genes.txt \
    --boot 20 --out l1 > /dev/null
python3 "$SCRIPTS/kopajtich_tables.py" --dir raw --gtf gencode.gtf.gz \
    --out conv --report rep > conv.log || { cat conv.log; fail "conversion"; }
grep -q "batch column 'BATCH_RUN'" conv.log || { cat conv.log; fail "batch column"; }
grep -q "60 individuals with RNA and protein (3 batches)" conv.log || { cat conv.log; fail "pairs"; }
[ -s rep/layout.txt ] || fail "layout.txt"
for d in direct conv; do
    if [ $d = direct ]; then
        p=(--patients-rna patients_rna.tsv --patients-protein patients_protein.tsv
           --patients-samples patients_samples.tsv)
    else
        p=(--patients-rna conv/patients_rna.tsv --patients-protein conv/patients_protein.tsv
           --patients-samples conv/patients_samples.tsv)
    fi
    python3 "$SCRIPTS/l1_score.py" validate --l1 l1/l1.tsv "${p[@]}" \
        --gtf gencode.gtf.gz --out val_$d > val_$d.log
    grep -q "protein centred within 3 batches" val_$d.log || fail "batches, $d"
done
python3 - <<'PY' || fail "events"
import pandas as pd
a = pd.read_csv("val_direct/patient_events.tsv", sep="\t")
b = pd.read_csv("val_conv/patient_events.tsv", sep="\t")
planted = {f"MITOG{i}" for i in range(1, 31)}
assert planted <= set(b.gene), planted - set(b.gene)
hi = b[b.gene.isin([f"MITOG{i}" for i in range(1, 16)])]
lo = b[b.gene.isin([f"MITOG{i}" for i in range(16, 31)]) & (b.rna_z < -6)]
assert hi.propagated.mean() > 0.8 and lo.propagated.mean() < 0.2
assert abs(len(a) - len(b)) <= 0.1 * len(a), (len(a), len(b))
PY
echo "PASS test_kopajtich"
