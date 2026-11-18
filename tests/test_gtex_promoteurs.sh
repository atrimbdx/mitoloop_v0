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
# GENCODE v26 of GTEx (--gtf-gtex), tags, CAGE peaks (FANTOM5), dominant TSS out of the deletion
ex() { printf 'chr1\tHAVANA\texon\t%s\t%s\t.\t%s\t.\tgene_id "%s.1"; transcript_id "%s.1"; gene_name "%s"; exon_number 1;%s\n' "$@"; }
{
  # G4 (VDAC1-like), + : MANE TSS 60000, exon 1 of 300 bp ; G5 the same at 80000
  tx 60000 70000 + ENSG4 ENST41 G4 G4-201 ' tag "MANE_Select";'
  ex 60000 60299 + ENSG4 ENST41 G4 ' tag "MANE_Select";'
  tx 80000 90000 + ENSG5 ENST51 G5 G5-201 ' tag "MANE_Select";'
  ex 80000 80299 + ENSG5 ENST51 G5 ' tag "MANE_Select";'
} > "$T/p.gtf"
{
  # v26: G4-202 starts at +233 (exon 1 of 167 bp), untagged ; G5-202 at +233 tagged mRNA_start_NF
  tx 60000 70000 + ENSG4 ENST41 G4 G4-201 ' tag "basic";'
  tx 60233 70000 + ENSG4 ENST42 G4 G4-202 ' tag "basic"; transcript_support_level "1";'
  ex 60233 60399 + ENSG4 ENST42 G4 ''
  tx 80000 90000 + ENSG5 ENST51 G5 G5-201 ' tag "basic";'
  tx 80233 90000 + ENSG5 ENST52 G5 G5-202 ' tag "mRNA_start_NF";'
  ex 80233 80399 + ENSG5 ENST52 G5 ''
} > "$T/v26.gtf"
printf 'chr1\t60220\t60240\tp1@G4\t9\t+\nchr1\t80220\t80240\tp1@G5\t9\t+\nchr1\t59990\t60010\tp2@G4\t9\t+\n' > "$T/cage.bed"
python3 -c 'import json; json.dump({"data":[{"transcriptId":"ENST41.1","tissueSiteDetailId":"Muscle_Skeletal","median":1.0},{"transcriptId":"ENST42.1","tissueSiteDetailId":"Muscle_Skeletal","median":30.0}],"tissus":["Muscle_Skeletal"]}, open("'"$T"'/cache/G4.json","w"))'
python3 -c 'import json; json.dump({"data":[{"transcriptId":"ENST51.1","tissueSiteDetailId":"Muscle_Skeletal","median":1.0},{"transcriptId":"ENST52.1","tissueSiteDetailId":"Muscle_Skeletal","median":30.0}],"tissus":["Muscle_Skeletal"]}, open("'"$T"'/cache/G5.json","w"))'
printf 'G4 G5\n' > "$T/genes2.txt"
python3 scripts/gtex_promoteurs.py --genes "$T/genes2.txt" --gtf "$T/p.gtf" --gtf-gtex "$T/v26.gtf" --cage "$T/cage.bed" \
  --tissus muscle=Muscle_Skeletal --cache "$T/cache" --out "$T/o2.tsv" --tss-dominants "$T/tdom.tsv" 2> "$T/log2" || { cat "$T/log2"; ko "v26"; }
cat "$T/o2.tsv" "$T/tdom.tsv"
c() { awk -F'\t' -v g="$1" -v k="$2" 'NR==1{for(i=1;i<=NF;i++) h[$i]=i} $1==g{print $h[k]}' "$T/o2.tsv"; }
[ "$(c G4 statut_tss) $(c G4 part_couverte_calib) $(c G4 promoteur_dominant) $(c G4 cage_dominant) $(c G4 meme_pic_cage)" = "dominant_confirme 0.032 +233 oui non" ] \
  && ok "G4 TSS dominant +233 confirmé par la CAGE" || ko "G4 $(grep G4 "$T/o2.tsv")"
[ "$(c G5 statut_tss) $(c G5 part_mrna_start_nf)" = "dominant_non_confirme 0.968" ] && c G5 alerte | grep -q dominant_tronque_5p \
  && ok "G5 transcrit tronqué en 5' : non confirmé" || ko "G5 $(grep G5 "$T/o2.tsv")"
[ "$(tail -n +2 "$T/tdom.tsv")" = "$(printf 'G4\tmuscle\tchr1\t60232\t+\t60399\t233')" ] \
  && ok "tss_dominants : G4 seul" || ko "tss_dominants $(cat "$T/tdom.tsv")"
