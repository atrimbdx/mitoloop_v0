#!/bin/bash
#SBATCH --job-name=simpol_traj
#SBATCH --time=12:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_logs/traj_%A_%a.out
# Une tache du job array = une voie, un jeu de parametres, un replicat.
# Aucune communication entre taches.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh
mkdir -p "$RESULTATS" slurm_logs
export OPENMM_CPU_THREADS=1

# PLAN permet de relancer un sous-ensemble : reprendre.py ecrit la liste
# des taches manquantes, et l'array lit ce fichier-la.
PLAN=${PLAN:-plan.tsv}
L=$((SLURM_ARRAY_TASK_ID + 1))        # la ligne 1 est l'en-tete
read COND VOIE JEU REP < <(sed -n "${L}p" "$PLAN")
[ -n "${COND:-}" ] || { echo "ligne $L absente de $PLAN"; exit 1; }
echo "  tâche $SLURM_ARRAY_TASK_ID de $PLAN : $COND voie $VOIE jeu $JEU réplicat $REP"

if [ "$COND" = "MEL" ]; then MAT="$MAT_MEL"; CTCF="$CTCF_MEL";
else MAT="$MAT_FIBRO"; CTCF="$CTCF_FIBRO"; fi
OPT=""
[ -s "$CTCF" ] && OPT="--ctcf $CTCF"

$PY une_trajectoire.py --locus "$LOCUS" --condition "$COND" --voie "$VOIE" \
    --jeu "$JEU" --replicat "$REP" --matrice "$MAT" $OPT \
    --blocs "$BLOCS" --pas-par-bloc "$PAS_PAR_BLOC" --sortie "$RESULTATS"
