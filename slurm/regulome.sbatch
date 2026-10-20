#!/bin/bash
#SBATCH --job-name=regulome
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
# Partition, account and QOS directives are site specific: add them above.
#
# Runs the computing steps of run_mcia.sh (run, reannot, then bundle) on a compute
# node, after `bash run_mcia.sh setup` and `bash run_mcia.sh fetch` on the
# login node (compute nodes have no internet access).
#
# Usage
#   mkdir -p logs
#   sbatch slurm/regulome.sbatch

cd "${SLURM_SUBMIT_DIR:-.}"
bash run_mcia.sh run
rc=$?
bash run_mcia.sh reannot || rc=1
bash run_mcia.sh bundle
exit $rc
