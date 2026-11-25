#!/usr/bin/env bash
# fetch_hic_encode.py : filtre des fichiers .hic ENCODE (biosample, exclusion, assay, type de
# sortie) et table des librairies, sur une réponse de recherche enregistrée (sans réseau)
set -euo pipefail
cd "$(dirname "$0")/.."
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
fail() { echo "FAIL test_fetch_hic_encode: $1"; exit 1; }
python3 - "$T/r.json" <<'PY'
import json, sys
f = lambda acc, term, assay, ot, fmt="hic", st="released": {"accession": acc, "file_format": fmt, "status": st,
    "assembly": "GRCh38", "output_type": ot, "href": f"/files/{acc}/@@download/{acc}.hic", "file_size": 10,
    "cloud_metadata": {"url": f"https://s3/{acc}.hic"}, "assay_title": assay, "dataset": "/experiments/E1/",
    "biosample_ontology": {"term_name": term, "classification": "cell line"}}
g = [f("A1", "GM12878", "intact Hi-C", "mapq"), f("A2", "GM12878", "intact Hi-C", "all"),
     f("A3", "K562", "intact Hi-C", "mapq"), f("A4", "IMR-90", "in situ Hi-C", "mapq"),
     f("A5", "GM12878", "Micro-C", "mapq"), f("A6", "GM12878", "intact Hi-C", "mapq", fmt="bed"),
     f("A7", "GM12878", "intact Hi-C", "mapq", st="archived")]
json.dump({"@graph": g}, open(sys.argv[1], "w"))
PY
python3 scripts/fetch_hic_encode.py --json "$T/r.json" --biosamples "GM12878|K562" --exclure K562 --liste > "$T/l.out" 2> "$T/l.log" || { cat "$T/l.log"; fail liste; }
grep -qP "^GM12878\tcell line\tintact Hi-C\tall\t1$" "$T/l.out" && grep -qP "^GM12878\tcell line\tintact Hi-C\tmapq\t1$" "$T/l.out" \
    && ! grep -q K562 "$T/l.out" && grep -q "2 fichiers .hic ; 1 biosamples" "$T/l.log" || { cat "$T/l.out" "$T/l.log"; fail "liste content"; }
python3 scripts/fetch_hic_encode.py --json "$T/r.json" --biosamples "GM12878" --output-type mapq --gtex-tissu "Whole Blood" --mtcn 2 --out "$T/t.tsv" 2> "$T/t.log" || { cat "$T/t.log"; fail table; }
[ "$(tail -n +2 "$T/t.tsv" | cut -f1,2,8,9,10,11,12 | tr '\t' ' ')" = "A1 GM12878 https://s3/A1.hic https://www.encodeproject.org/files/A1/@@download/A1.hic Whole Blood 2 1" ] \
    || { cat "$T/t.tsv"; fail "table content"; }
echo "PASS test_fetch_hic_encode"
