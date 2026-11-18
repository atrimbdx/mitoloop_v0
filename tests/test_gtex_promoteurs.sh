#!/usr/bin/env bash
# gtex_promoteurs.py sans réseau : GTF synthétique et réponses GTEx en cache
set -euo pipefail
cd "$(dirname "$0")/.."
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
tx() { printf 'chr1\tHAVANA\ttranscript\t%s\t%s\t.\t%s\t.\tgene_id "%s.1"; transcript_id "%s.1"; gene_name "%s"; transcript_name "%s";%s\n' "$@"; }
{
  # G1 brin - : TSS MANE 10000 ; TSS en aval à 9000 (+1000), même groupe de 1 kb, majoritaire
  tx 5000 10000 - ENSG1 ENST11 G1 G1-201 ' tag "MANE_Select";'
  tx 5000 9000  - ENSG1 ENST12 G1 G1-202 ''
  # G2 brin + : TSS MANE 5000 ; promoteur alternatif à 20000 (+15000), majoritaire
  tx 5000 30000  + ENSG2 ENST21 G2 G2-201 ' tag "MANE_Select";'
  tx 20000 30000 + ENSG2 ENST22 G2 G2-202 ''
  # G3 brin + : un seul promoteur
  tx 40000 50000 + ENSG3 ENST31 G3 G3-201 ' tag "MANE_Select";'
} > "$T/g.gtf"
mkdir -p "$T/cache"
j() { python3 -c 'import json,sys; a=sys.argv[1:]; print(json.dumps({"data":[{"transcriptId":a[i]+".1","tissueSiteDetailId":"Muscle_Skeletal","median":float(a[i+1])} for i in range(0,len(a),2)]}))' "$@"; }
j ENST11 2 ENST12 18 > "$T/cache/G1.json"
j ENST21 3 ENST22 27 > "$T/cache/G2.json"
j ENST31 40 > "$T/cache/G3.json"
printf 'G1 G2 G3\n' > "$T/genes.txt"
python3 scripts/gtex_promoteurs.py --genes "$T/genes.txt" --gtf "$T/g.gtf" --tissus muscle=Muscle_Skeletal \
  --cache "$T/cache" --out "$T/o.tsv" 2> "$T/log"
cat "$T/o.tsv"
ok() { echo "OK  $1"; }; ko() { echo "ECHEC $1"; exit 1; }
l() { awk -F'\t' -v g="$1" '$1==g' "$T/o.tsv"; }
[ "$(head -1 "$T/o.tsv" | cut -f11,12)" = "$(printf 'alerte\tpart_couverte')" ] && ok "colonnes" || ko "colonnes"
# G1 : le groupe MANE porte tout, mais le transcrit dominant démarre à +1000, hors de -2500/+100
[ "$(l G1 | cut -f5,7,12)" = "$(printf '1.000\t+1000\t0.100')" ] && ok "G1 part_mane 1, TSS +1000, couverte 0.1" || ko "G1 $(l G1)"
l G1 | cut -f11 | grep -qx 'hors_deletion' && ok "G1 hors_deletion seul" || ko "G1 alerte"
[ "$(l G2 | cut -f7,12)" = "$(printf '+15000\t0.100')" ] && ok "G2 TSS +15000" || ko "G2 $(l G2)"
l G2 | cut -f11 | grep -qx 'promoteur_dominant_autre,hors_deletion' && ok "G2 deux alertes" || ko "G2 alerte"
[ "$(l G3 | cut -f11,12)" = "$(printf '.\t1.000')" ] && ok "G3 sans alerte" || ko "G3 $(l G3)"
