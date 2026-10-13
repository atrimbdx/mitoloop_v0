#!/bin/bash
#SBATCH --job-name=v4c_extract
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%A_%a.out
# Partition, account and QOS directives are site specific: add them above.
#
# One array task per library, that is per data row of the library table.
#
# Environment variables
#   CONDA_PROFILE   path to <conda>/etc/profile.d/conda.sh (required)
#   MITOLOOP_DIR    working directory holding scripts/ and the inputs (required)
#   CONDA_ENV       Python environment (default mitoloop_v0)
#   LIBRARY_TABLE   library table (default tissues.tsv)
#
# Usage
#   mkdir -p logs
#   N=$(( $(wc -l < tissues.tsv) - 1 ))
#   sbatch --array=1-${N}%10 slurm/extract.sbatch --gtf data/gencode.gtf.gz \
#       --genes genes.txt --controls controls.txt --outdir results/raw
#   sbatch --array=1-${N}%10 slurm/extract.sbatch \
#       --anchors results/loops/partners.tsv --outdir results/raw_partners

set -euo pipefail
source "${CONDA_PROFILE:?set CONDA_PROFILE to <conda>/etc/profile.d/conda.sh}"
conda activate "${CONDA_ENV:-mitoloop_v0}"
cd "${MITOLOOP_DIR:?set MITOLOOP_DIR to the working directory}"

python3 scripts/v4c_pipeline.py extract \
    --tissues "${LIBRARY_TABLE:-tissues.tsv}" \
    --index "$SLURM_ARRAY_TASK_ID" "$@"
