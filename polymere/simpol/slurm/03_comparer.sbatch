#!/bin/bash
#SBATCH --job-name=simpol_comp
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_logs/comp_%j.out
# Agrege les replicats, confronte chaque voie a la mesure et les deux
# voies entre elles. A lancer quand le job array est terminé.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh
mkdir -p slurm_logs

for C in MEL FIBRO; do
  if [ "$C" = "MEL" ]; then MAT="$MAT_MEL"; else MAT="$MAT_FIBRO"; fi
  $PY comparer.py --locus "$LOCUS" --condition "$C" --matrice "$MAT" \
      --resultats "$RESULTATS" --sortie comparaison
done
$PY bilan.py --locus "$LOCUS" --sortie bilan
