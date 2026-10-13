#!/bin/bash
#SBATCH --job-name=v4c_r
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
# Partition, account and QOS directives are site specific: add them above.
#
# Runs one subcommand of v4c_loops.R with as many processes as allocated CPUs.
#
# Environment variables
#   CONDA_PROFILE   path to <conda>/etc/profile.d/conda.sh (required)
#   MITOLOOP_DIR    working directory holding scripts/ and the inputs (required)
#   CONDA_ENV_R     R environment (default mitoloop_v0-r)
#
# Usage
#   mkdir -p logs
#   sbatch slurm/run_r.sbatch call --libraries results/libraries.tsv \
#       --export results/export/both --wsize 3 --out results/loops

set -euo pipefail
source "${CONDA_PROFILE:?set CONDA_PROFILE to <conda>/etc/profile.d/conda.sh}"
conda activate "${CONDA_ENV_R:-mitoloop_v0-r}"
cd "${MITOLOOP_DIR:?set MITOLOOP_DIR to the working directory}"

# one BLAS thread per R process, so that parallel processes do not
# oversubscribe the allocated cores
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

Rscript scripts/v4c_loops.R "$@" --cores "${SLURM_CPUS_PER_TASK:-1}"
