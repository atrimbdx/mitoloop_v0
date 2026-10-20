#!/usr/bin/env bash
# reannot_alphagenome.sh
# (version MitoLoop : étape 5, AlphaGenome lu par son index sur les seules
#  positions des variants ; le reste est identique à la version d'origine)
# Conversion hors ligne de variants (g. / c. / n.) en CHROM/POS/REF/ALT (GRCh38),
# puis annotation : VEP (conséquence, HGVSc), SpliceAI (scores précalculés, plugin VEP),
# MORFEE (uORF / dORF, SNV des UTR, via ANNOVAR) et AlphaGenome (fichier tabix).
#
# Entrée : TSV sans en-tête
#   col1 = transcrit Ensembl (ENST00000379370, version facultative)
#   col2 = chr2:g.10543284T>A | c.541-6T>G | n.190+24517G>A
#
# Usage : ./reannot_alphagenome.sh variants.tsv [prefixe_sortie]
set -euo pipefail

IN=${1:?Usage: $0 variants.tsv [prefixe_sortie]}
OUT=${2:-reannot}
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

# ------------------------- Paramètres -------------------------
# Arborescence identique sur Curta et sur la VM (voir 00_config.sh)
if [ -z "${REANNOT_BASE:-}" ]; then
  if [ -d "/scratch/$USER" ]; then REANNOT_BASE=/scratch/$USER/reannot; else REANNOT_BASE=$HOME/reannot; fi
fi
BASE=$REANNOT_BASE
export PATH="$BASE/envs/vep/bin:$PATH"
export BCFTOOLS_PLUGINS="$BASE/envs/vep/libexec/bcftools"
THREADS=${THREADS:-4}

VEP=${VEP:-$BASE/envs/vep/bin/vep}
VEP_CACHE=${VEP_CACHE:-$BASE/data/vep_cache}
VEP_PLUGINS=${VEP_PLUGINS:-$BASE/data/vep_plugins}
FASTA=${FASTA:-$BASE/data/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz}
VEP_REL=$(cat "$BASE/VEP_RELEASE")
GTF=${GTF:-$BASE/data/ref/Homo_sapiens.GRCh38.${VEP_REL}.gtf.gz}   # même release que le cache VEP

SPLICEAI_SNV=${SPLICEAI_SNV:-$BASE/data/spliceai/spliceai_scores.masked.snv.hg38.vcf.gz}
SPLICEAI_INDEL=${SPLICEAI_INDEL:-$BASE/data/spliceai/spliceai_scores.masked.indel.hg38.vcf.gz}

ANNOVAR_DIR=${ANNOVAR_DIR:-$BASE/data/annovar}
ANNOVAR_DB=${ANNOVAR_DB:-$ANNOVAR_DIR/humandb}
ANNOVAR_PROTOCOL=${ANNOVAR_PROTOCOL:-refGene,ensGene}
ANNOVAR_OPERATION=${ANNOVAR_OPERATION:-g,g}
RSCRIPT=${RSCRIPT:-$BASE/envs/morfee/bin/Rscript}
MORFEE_RDATA=${MORFEE_RDATA:-$BASE/data/morfee/MORFEE_DATA_GRCh38.RData}
MORFEE_SCRIPT=${MORFEE_SCRIPT:-$SCRIPT_DIR/run_morfee.R}
MORFEE_FIELDS=${MORFEE_FIELDS:-MORFEE_uTIS,MORFEE_uSTOP,MORFEE_uSTOP_creation,MORFEE_dTIS,MORFEE_dSTOP,MORFEE_dSTOP_creation}

# AlphaGenome : à ajuster à la structure réelle du fichier
AG_FILE=${AG_FILE:-$BASE/data/alphagenome/alphagenome.tsv.gz}   # bgzip + index tabix
AG_FORMAT=${AG_FORMAT:-tsv}                                     # "tsv" ou "vcf"
AG_COLS=${AG_COLS:-CHROM,POS,REF,ALT,AG_SCORE}                  # si tsv : colonnes pour bcftools annotate -c
AG_TAGS=${AG_TAGS:-AG_SCORE}                                    # scores AlphaGenome à reporter
AG_CHR_PREFIX=${AG_CHR_PREFIX-chr}   # "chr" si AlphaGenome utilise chr1..chrM ; "" s'il utilise 1..MT
# ---------------------------------------------------------------

IFS=',' read -ra TAGS <<< "$AG_TAGS"
IFS=',' read -ra MF <<< "$MORFEE_FIELDS"
NVAR=$(( ${#MF[@]} + ${#TAGS[@]} ))

# 1. Conversion des notations g. / c. / n. en VCF
#    (structures de transcrits issues du GTF Ensembl, bases issues du FASTA)
python3 "$SCRIPT_DIR/hgvs_to_vcf.py" "$IN" "$GTF" "$FASTA" "$(command -v samtools)" \
        "$OUT.conv.vcf" "$OUT.map.tsv" "$OUT.conv_warnings.txt"

# 2. VEP hors ligne + plugin SpliceAI : conséquences, HGVSc, scores SpliceAI
$VEP --offline --cache --dir_cache "$VEP_CACHE" --assembly GRCh38 \
     --fasta "$FASTA" --format vcf -i "$OUT.conv.vcf" \
     --vcf -o "$OUT.vep.vcf" --force_overwrite --no_stats \
     --symbol --hgvs --mane --canonical --fork "$THREADS" \
     --dir_plugins "$VEP_PLUGINS" \
     --plugin "SpliceAI,snv=$SPLICEAI_SNV,indel=$SPLICEAI_INDEL" \
     --warning_file "$OUT.vep_warnings.txt"

# 3. Contigs, tri, normalisation, passage en nomenclature chr1..chrM
{ for c in $(seq 1 22) X Y; do printf '%s\tchr%s\n' "$c" "$c"; done
  printf 'MT\tchrM\n'; } > "$OUT.num2chr.txt"
awk -v OFS='\t' '{print $2, $1}' "$OUT.num2chr.txt" > "$OUT.chr2num.txt"

bcftools reheader -f "$FASTA.fai" "$OUT.vep.vcf" \
  | bcftools sort -Ou - \
  | bcftools norm -f "$FASTA" -c w -m -any -Ou - \
  | bcftools annotate --rename-chrs "$OUT.num2chr.txt" -Oz -o "$OUT.norm.vcf.gz" -
bcftools index -t -f "$OUT.norm.vcf.gz"

# 4. MORFEE sur les SNV (ANNOVAR puis R), résultat en TSV indexé
bcftools view -v snps "$OUT.norm.vcf.gz" -Ou | bcftools annotate -x INFO -Ov -o "$OUT.snv.vcf"
if [ "$(bcftools view -H "$OUT.snv.vcf" | head -n 1 | wc -l)" -gt 0 ]; then
  perl "$ANNOVAR_DIR/table_annovar.pl" "$OUT.snv.vcf" "$ANNOVAR_DB" -buildver hg38 \
       -out "$OUT.annovar" -remove -protocol "$ANNOVAR_PROTOCOL" \
       -operation "$ANNOVAR_OPERATION" -nastring . -vcfinput
  $RSCRIPT "$MORFEE_SCRIPT" "$OUT.annovar.hg38_multianno.vcf" "$MORFEE_RDATA" \
           "$OUT.morfee.raw.tsv" "$MORFEE_FIELDS"
else
  : > "$OUT.morfee.raw.tsv"
fi
sort -k1,1 -k2,2n "$OUT.morfee.raw.tsv" | bgzip -c > "$OUT.morfee.tsv.gz"
tabix -f -s1 -b2 -e2 "$OUT.morfee.tsv.gz"

for t in "${MF[@]}"; do
  printf '##INFO=<ID=%s,Number=.,Type=String,Description="MORFEE %s">\n' "$t" "$t"
done > "$OUT.morfee.hdr"
bcftools annotate -a "$OUT.morfee.tsv.gz" -c "CHROM,POS,REF,ALT,$MORFEE_FIELDS" \
  -h "$OUT.morfee.hdr" -Oz -o "$OUT.morfee.vcf.gz" "$OUT.norm.vcf.gz"
bcftools index -t -f "$OUT.morfee.vcf.gz"

# 5. AlphaGenome (correspondance exacte CHROM/POS/REF/ALT)
if [ "$AG_FORMAT" = "tsv" ]; then
  for t in "${TAGS[@]}"; do
    printf '##INFO=<ID=%s,Number=.,Type=String,Description="AlphaGenome %s">\n' "$t" "$t"
  done > "$OUT.ag.hdr"
  AGARGS=(-c "$AG_COLS" -h "$OUT.ag.hdr")
else
  cols=$(printf 'INFO/%s,' "${TAGS[@]}"); AGARGS=(-c "${cols%,}")
fi
#    Le fichier AlphaGenome couvre tout le génome (dizaines de Go) : on en extrait
#    d'abord, par son index, les seules lignes aux positions des variants, puis on
#    annote avec cet extrait. Sans cela, bcftools annotate parcourt tout le fichier.
#    Correspondance et colonnes inchangées (CHROM/POS/REF/ALT exacts, AG_COLS).
bcftools query -f '%CHROM\t%POS\n' "$OUT.morfee.vcf.gz" \
  | { if [ -n "$AG_CHR_PREFIX" ]; then cat; else
        awk -F'\t' -v OFS='\t' 'NR == FNR { m[$1] = $2; next } { print ($1 in m ? m[$1] : $1), $2 }' \
          "$OUT.chr2num.txt" -; fi; } \
  | awk '!vu[$0]++' > "$OUT.ag_pos.txt"
if [ "$AG_FORMAT" = "tsv" ]; then
  AG_SUB="$OUT.ag_sub.tsv.gz"
  if [ -s "$OUT.ag_pos.txt" ]; then
    tabix -h -R "$OUT.ag_pos.txt" "$AG_FILE" | bgzip -c > "$AG_SUB"
  else
    tabix -H "$AG_FILE" | bgzip -c > "$AG_SUB"
  fi
  tabix -f -s1 -b2 -e2 "$AG_SUB"
else
  AG_SUB="$OUT.ag_sub.vcf.gz"
  if [ -s "$OUT.ag_pos.txt" ]; then
    bcftools view -R "$OUT.ag_pos.txt" -Oz -o "$AG_SUB" "$AG_FILE"
  else
    bcftools view -h -Oz -o "$AG_SUB" "$AG_FILE"
  fi
  bcftools index -t -f "$AG_SUB"
fi
echo "AlphaGenome : $(zcat "$AG_SUB" | grep -vc '^#' || true) lignes extraites pour $(wc -l < "$OUT.ag_pos.txt") positions" >&2
if [ -n "$AG_CHR_PREFIX" ]; then
  bcftools annotate -a "$AG_SUB" "${AGARGS[@]}" -Oz -o "$OUT.annot.vcf.gz" "$OUT.morfee.vcf.gz"
else
  bcftools annotate --rename-chrs "$OUT.chr2num.txt" -Ou "$OUT.morfee.vcf.gz" \
    | bcftools annotate -a "$AG_SUB" "${AGARGS[@]}" -Ou - \
    | bcftools annotate --rename-chrs "$OUT.num2chr.txt" -Oz -o "$OUT.annot.vcf.gz" -
fi
bcftools index -t -f "$OUT.annot.vcf.gz"

# 6. Une ligne par transcrit VEP, puis sélection du transcrit donné en entrée
SAI="%SpliceAI_pred_DS_AG\t%SpliceAI_pred_DS_AL\t%SpliceAI_pred_DS_DG\t%SpliceAI_pred_DS_DL"
SAI="$SAI\t%SpliceAI_pred_DP_AG\t%SpliceAI_pred_DP_AL\t%SpliceAI_pred_DP_DG\t%SpliceAI_pred_DP_DL"
QV=$(printf '\t%%INFO/%s' "${MF[@]}" "${TAGS[@]}")
bcftools +split-vep "$OUT.annot.vcf.gz" -d \
  -f "%ID\t%CHROM\t%POS\t%REF\t%ALT\t%Feature\t%SYMBOL\t%Consequence\t%HGVSc\t${SAI}${QV}\n" \
  > "$OUT.split.tsv"

{
  printf 'transcrit\tvariant_entree\tchrom\tpos\tref\talt\tgene\tconsequence\tHGVSc_VEP'
  printf '\tSpliceAI_DS_AG\tSpliceAI_DS_AL\tSpliceAI_DS_DG\tSpliceAI_DS_DL\tSpliceAI_DS_max'
  printf '\tSpliceAI_DP_AG\tSpliceAI_DP_AL\tSpliceAI_DP_DG\tSpliceAI_DP_DL'
  for t in "${MF[@]}" "${TAGS[@]}"; do printf '\t%s' "$t"; done
  printf '\tstatut\tcontrole_hgvs\n'

  awk -F'\t' -v OFS='\t' -v nvar="$NVAR" '
  function ctrl(entree, hgvsc) {
    if (entree !~ /^[cn]\./) return "-"
    if (hgvsc == "." || hgvsc == "") return "absent"
    sub(/^[^:]*:/, "", hgvsc)
    gsub(/ /, "", entree)
    return (hgvsc == entree) ? "identique" : "different (" hgvsc ")"
  }
  FILENAME == ARGV[1] {
    id = $1; f = $6; sub(/\.[0-9]+$/, "", f)
    if (!(id in loc)) {
      loc[id] = $2 OFS $3 OFS $4 OFS $5
      a = ""; for (i = 18; i < 18 + nvar; i++) a = a OFS $i
      vl[id] = a
    }
    m = "."
    for (i = 10; i <= 13; i++) if ($i != "." && $i != "" && (m == "." || $i + 0 > m + 0)) m = $i
    row[id, f] = $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12 OFS $13 OFS m OFS $14 OFS $15 OFS $16 OFS $17
    hg[id, f] = $9
    next
  }
  {
    h = $1; tx = $2; sub(/\.[0-9]+$/, "", tx)
    if (!(h in loc)) {
      na = "NA"; for (i = 2; i <= 16 + nvar; i++) na = na OFS "NA"
      print $2, $3, na, "non_converti", "NA"
    } else if ((h, tx) in row) {
      print $2, $3, loc[h], row[h, tx] vl[h], "ok", ctrl($3, hg[h, tx])
    } else {
      na = "NA"; for (i = 2; i <= 12; i++) na = na OFS "NA"
      print $2, $3, loc[h], na vl[h], "transcrit_absent_du_cache", "NA"
    }
  }' "$OUT.split.tsv" "$OUT.map.tsv"
} > "$OUT.final.tsv"

echo "Résultat : $OUT.final.tsv"
echo "Non convertis : $(awk -F'\t' '$(NF-1)=="non_converti"' "$OUT.final.tsv" | wc -l) (détail : $OUT.conv_warnings.txt)"
echo "Contrôle HGVS non identique : $(awk -F'\t' 'NR>1 && $NF != "identique" && $NF != "-"' "$OUT.final.tsv" | wc -l)"
