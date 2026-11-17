#!/bin/bash
# trajectoires.sbatch - simpol path A trajectories of one MitoLoop locus, one
# task of the job array per line of the plan (condition, parameter set,
# replicate). The CTCF bed of each condition is <LOCUSDIR>/ctcf_<condition>.bed
# (preparer_locus.py); trajectories go to <LOCUSDIR>/resultats.
#
# Usage (from mitoloop_v0, after preparer_locus.py preparer):
#   N=$(($(wc -l < polymere/abat/plan_ajustement.tsv) - 1))
#   sbatch --array=1-$N --export=ALL,LOCUSDIR=polymere/abat,LOCUS=abat,PLAN=plan_ajustement.tsv \
#       polymere/trajectoires.sbatch
# Variables: LOCUSDIR, LOCUS, PLAN (in LOCUSDIR), PY (python of the OpenMM
#   environment, default /scratch/$USER/envs/polymere/bin/python),
#   BLOCS (2000), PAS_PAR_BLOC (10000)
#SBATCH --job-name=polymere
#SBATCH --time=12:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/polymere_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
: "${LOCUSDIR:?LOCUSDIR}" "${LOCUS:?LOCUS}" "${PLAN:?PLAN}"
PY=${PY:-/scratch/$USER/envs/polymere/bin/python}
export SIMPOL_LOCI="$LOCUSDIR/loci.json" OPENMM_CPU_THREADS=1
L=$((SLURM_ARRAY_TASK_ID + 1))                 # line 1 is the header
read -r COND VOIE JEU REP < <(sed -n "${L}p" "$LOCUSDIR/$PLAN")
[ -n "${COND:-}" ] || { echo "line $L missing from $LOCUSDIR/$PLAN"; exit 1; }
CTCF="$LOCUSDIR/ctcf_${COND}.bed"
[ -s "$CTCF" ] || { echo "$CTCF missing or empty"; exit 1; }
ETIQ=$(printf "%s_%s_%s_jeu%03d_rep%02d" "$LOCUS" "$COND" "$VOIE" "$JEU" "$REP")
if [ -s "$LOCUSDIR/resultats/$ETIQ.npy" ]; then echo "$ETIQ already done"; exit 0; fi
echo "== $(date -Is) $ETIQ on $(hostname)"
"$PY" polymere/simpol/une_trajectoire.py --locus "$LOCUS" --condition "$COND" --voie "$VOIE" \
    --jeu "$JEU" --replicat "$REP" --ctcf "$CTCF" --matrice "$LOCUSDIR/matrice.tsv" \
    --blocs "${BLOCS:-2000}" --pas-par-bloc "${PAS_PAR_BLOC:-10000}" --sortie "$LOCUSDIR/resultats"
