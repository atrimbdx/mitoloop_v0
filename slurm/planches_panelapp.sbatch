#!/bin/bash
#SBATCH --job-name=planches
#SBATCH --partition=ubx2
#SBATCH --array=0-19%10
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#
# Plates of the PanelApp mitochondrial-disorder genes (scripts/planches_panelapp.py),
# one array task per (tissue, share of the genes): with TISSUS="muscle cardiomyocyte"
# and TRANCHES=10, tasks 0-9 draw the muscle, 10-19 the heart (--array=0-19).
# %10 keeps at most 10 tasks at a time (each reads its .hic from ENCODE).
# The compute nodes must reach the internet (the AG2 campaign does from them).
#
# Submit from ~/mitoloop_v0 with the mitodark environment active: sbatch passes
# the environment on, so python3 is the one of mitodark (hicstraw, Pillow).
#
# Environment variables
#   TISSUS      tissues, space-separated (default "muscle cardiomyocyte")
#   TRANCHES    shares per tissue (default 10); --array must be 0-(tissues x shares - 1)
#   REANNOT     reannot clone (default ~/reannot_git)
#   NIVEAUX     PanelApp levels (default 3)
#
# Usage
#   sbatch slurm/planches_panelapp.sbatch
#   then, once all tasks are done, for each tissue:
#   python3 scripts/planches_panelapp.py --tissu muscle --out results/planches_panelapp/muscle --fusionner

set -euo pipefail
read -r -a TS <<< "${TISSUS:-muscle cardiomyocyte}"
N=${TRANCHES:-10}
I=${SLURM_ARRAY_TASK_ID:-0}
T=${TS[$((I / N))]:?tâche $I au-delà de ${#TS[@]} tissus x $N tranches}
K=$((I % N))
R=${REANNOT:-$HOME/reannot_git}
python3 -c "import hicstraw, PIL" || { echo "hicstraw ou Pillow absent : activer mitodark avant sbatch"; exit 1; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "== $(date -Is) $T tranche $K/$N sur $(hostname)"
python3 scripts/planches_panelapp.py --panelapp data/panelapp_genes_all.tsv --renommes data/genes_renamed.tsv \
    --niveaux "${NIVEAUX:-3}" --tissu "$T" --dir results/loops_tissus --gtf data/gencode.gtf.gz \
    --elements "$R/results/paquet/re2g_tous/re2g_elements.tsv" --ag2-dir "$R/results/ag2_v2" --reannot "$R" \
    --tranche "$K/$N" --reprendre --out "results/planches_panelapp/$T"
echo "== $(date -Is) fin"
