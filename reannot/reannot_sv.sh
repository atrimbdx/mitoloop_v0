#!/usr/bin/env bash
# reannot_sv.sh — VM. Les variants structuraux en une commande, depuis un tableur
# Excel, un TSV ou un CSV (type, chromosome, début, fin ; GRCh38, base 1) ou un VCF de SV :
#   1. preparer_entree.py --sv : VCF de SV (<DEL>, <DUP>, <INV>...)
#   2. AnnotSV, s'il est installé : annotation et classe ACMG/ClinGen de chaque SV
#      (ANNOTSV = dossier d'installation d'AnnotSV, défaut $REANNOT_BASE/AnnotSV)
#   3. ajout_mitoloop_sv.py : gènes mitochondriaux et régions touchés, effet,
#      fréquence gnomAD-SV, second allèle (--snv), priorité
#
# Usage : bash reannot_sv.sh sv.xlsx [prefixe_sortie] [--feuille NOM] [--snv run1.final.spliceai2.mitoloop.tsv]
# Résultat : PREFIXE.annotsv.mitoloop.tsv (avec AnnotSV) ou PREFIXE.sv.mitoloop.tsv
set -euo pipefail
IN=${1:?Usage: $0 sv.xlsx|.tsv|.csv|.vcf [prefixe_sortie] [--feuille NOM] [--snv TABLE]}
OUT=${2:-$(basename "${IN%.*}")}
shift $(( $# >= 2 ? 2 : 1 ))
FEUILLE=(); SNV=()
while [ $# -gt 0 ]; do
  case "$1" in
    --feuille) FEUILLE=(--feuille "$2"); shift 2 ;;
    --snv) SNV=(--snv "$2"); shift 2 ;;
    *) echo "option inconnue : $1"; exit 2 ;;
  esac
done
BASE=${REANNOT_BASE:-$HOME/reannot}
export REANNOT_BASE=$BASE
S="$BASE/scripts"
export PATH="$BASE/envs/vep/bin:$PATH"          # bcftools, bedtools éventuel
log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

# 1. VCF de SV
if [[ "$IN" =~ \.vcf(\.gz)?$ ]]; then
  V="$IN"
else
  log "1. Préparation : $IN"
  V="$OUT.sv.vcf"
  python3 "$S/preparer_entree.py" "$IN" "$V" --sv ${FEUILLE[@]+"${FEUILLE[@]}"}
fi

# 2. AnnotSV
ANNOTSV=${ANNOTSV:-$BASE/AnnotSV}
T="$V"
if [ -x "$ANNOTSV/bin/AnnotSV" ]; then
  log "2. AnnotSV ($ANNOTSV)"
  export ANNOTSV
  "$ANNOTSV/bin/AnnotSV" -SVinputFile "$V" -outputFile "$OUT.annotsv.tsv" \
      -genomeBuild GRCh38 >"$OUT.annotsv.log" 2>&1 \
    || { log "AnnotSV en échec (voir $OUT.annotsv.log) : MitoLoop sur le VCF seul"; }
  [ -s "$OUT.annotsv.tsv" ] && T="$OUT.annotsv.tsv"
else
  log "2. AnnotSV non installé ($ANNOTSV/bin/AnnotSV) : étape sautée"
fi

# 3. MitoLoop
log "3. MitoLoop (SV)"
python3 "$S/ajout_mitoloop_sv.py" "$T" ${SNV[@]+"${SNV[@]}"}
log "Résultat : ${T%.*}.mitoloop.tsv"
