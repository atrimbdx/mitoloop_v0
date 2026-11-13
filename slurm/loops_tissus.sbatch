#!/bin/bash
#SBATCH --job-name=loops_tissus
#SBATCH --array=0-3
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
# Partition, account and QOS directives are site specific: add them above.
#
# Loops of each tissue on its own Hi-C libraries pooled as replicates (peakC
# combined.analysis, single.analysis for a tissue with one library): one array
# task per tissue (0 muscle, 1 cardiomyocyte, 2 neurone, 3 rein), each in its
# own folder, so that loop coordinates never merge across tissues.
#
# Environment variables
#   CONDA_PROFILE   path to <conda>/etc/profile.d/conda.sh (required)
#   MITOLOOP_DIR    working directory holding scripts/ (required)
#   CONDA_ENV_R     R environment (default mitoloop_v0-r)
#   TISSUS_DIR      folder written by loops_par_tissu.py librairies
#                   (default results/loops_tissus)
#   EXPORT          export folder: profiles of the mitochondrial viewpoints for
#                   call, of the partner anchors for reciprocal (required)
#   WSIZE           peakC window in bins (default 3, as for the ANR loops)
#
# Usage
#   EXPORT=/scratch/.../MITODARK/export/both sbatch slurm/loops_tissus.sbatch call
#   EXPORT=<merged partner export> sbatch slurm/loops_tissus.sbatch reciprocal

set -euo pipefail
STEP=${1:?call ou reciprocal}
TISSUS=(muscle cardiomyocyte neurone rein)
T=${TISSUS[${SLURM_ARRAY_TASK_ID:-0}]}
source "${CONDA_PROFILE:?set CONDA_PROFILE to <conda>/etc/profile.d/conda.sh}"
conda activate "${CONDA_ENV_R:-mitoloop_v0-r}"
cd "${MITOLOOP_DIR:?set MITOLOOP_DIR to the working directory}"
D="${TISSUS_DIR:-results/loops_tissus}/$T"
[ -s "$D/libraries.tsv" ] || { echo "$T : $D/libraries.tsv absent (aucune librairie Hi-C ?)"; exit 0; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "== $(date -Is) $STEP $T ($(($(wc -l < "$D/libraries.tsv") - 1)) librairies)"
Rscript scripts/v4c_loops.R "$STEP" --libraries "$D/libraries.tsv" \
    --export "${EXPORT:?set EXPORT}" --wsize "${WSIZE:-3}" --out "$D" \
    --cores "${SLURM_CPUS_PER_TASK:-1}"
echo "== $(date -Is) fin $STEP $T"
