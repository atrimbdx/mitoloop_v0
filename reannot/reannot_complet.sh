#!/usr/bin/env bash
# reannot_complet.sh — VM. Le pipeline reannot en une commande, depuis un tableur
# Excel, un TSV, un CSV ou un VCF :
#   1. preparer_entree.py : entrée de reannot (transcrit, variant g./c./n.) ;
#      sans transcrit, le MANE Select voisin (GTF Ensembl du pipeline)
#   2. reannot_alphagenome.sh : VEP, SpliceAI, MORFEE, AlphaGenome -> PREFIXE.final.tsv
#   3. ajout_spliceai2.sh, quand SpliceAI2 est installé -> PREFIXE.final.spliceai2.tsv
#   4. ajout_mitoloop.py, quand MitoLoop est installé -> ...mitoloop.tsv
#
# Usage : bash reannot_complet.sh variants.xlsx [prefixe_sortie] [--feuille NOM]
#   (un TSV déjà au format de reannot, sans en-tête, est accepté tel quel)
# Variables : REANNOT_BASE (défaut $HOME/reannot), THREADS, et celles de
# reannot_alphagenome.sh (AG_FILE, ...).
set -euo pipefail
IN=${1:?Usage: $0 variants.xlsx|.tsv|.csv|.vcf.gz [prefixe_sortie] [--feuille NOM]}
OUT=${2:-$(basename "${IN%.*}")}
shift $(( $# >= 2 ? 2 : 1 ))
BASE=${REANNOT_BASE:-$HOME/reannot}
export REANNOT_BASE=$BASE
S="$BASE/scripts"
export PATH="$BASE/envs/vep/bin:$PATH"
VEP_REL=$(cat "$BASE/VEP_RELEASE")
GTF=${GTF:-$BASE/data/ref/Homo_sapiens.GRCh38.${VEP_REL}.gtf.gz}
log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

# 1. Entrée
log "1. Préparation de l'entrée : $IN"
python3 "$S/preparer_entree.py" "$IN" "$OUT.entree.tsv" --gtf "$GTF" "$@"

# 2. Pipeline reannot
log "2. VEP, SpliceAI, MORFEE, AlphaGenome"
bash "$S/reannot_alphagenome.sh" "$OUT.entree.tsv" "$OUT"
T="$OUT.final.tsv"

# 3. SpliceAI2
if [ -f "$BASE/data/spliceai2/spliceai2.tsv.gz" ] && [ -f "$S/ajout_spliceai2.sh" ]; then
  log "3. SpliceAI2"
  bash "$S/ajout_spliceai2.sh" "$T"
  T="${T%.tsv}.spliceai2.tsv"
else
  log "3. SpliceAI2 non installé : étape sautée"
fi

# 4. MitoLoop
if [ -f "$BASE/data/mitoloop/genes.tsv" ] && [ -f "$S/ajout_mitoloop.py" ]; then
  log "4. MitoLoop"
  python3 "$S/ajout_mitoloop.py" "$T"
  T="${T%.tsv}.mitoloop.tsv"
else
  log "4. MitoLoop non installé : étape sautée"
fi
log "Résultat : $T"
