#!/bin/bash
#SBATCH --job-name=simpol_verif
#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_logs/verif_%j.out
# Les trois controles que la comparaison par profils ne fait pas, puis le
# rendu 3D. A lancer apres 03_comparer.sbatch, dont le bilan.tsv sert a
# choisir le jeu retenu pour chaque voie.
#
# Le nul par permutation du differentiel est le poste couteux : chaque
# tirage recalcule deux matrices de contact sur l'ensemble des
# conformations. 24 tirages sur 790 billes tiennent dans le temps demande.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh
mkdir -p slurm_logs figures

$PY verifier.py --locus "$LOCUS" \
    --mel "$MAT_MEL" --fibro "$MAT_FIBRO" \
    --resultats "$RESULTATS" --bilan bilan.tsv \
    --tirages-nul "${TIRAGES_NUL:-24}" \
    --sortie figures/verification

$PY rendu3d.py --locus "$LOCUS" --resultats "$RESULTATS" \
    --bilan bilan.tsv --sortie figures/conformation

echo
echo "  figures/verification.png          les trois contrôles"
echo "  figures/verification_*.tsv        les tableaux"
echo "  figures/conformation_*.png        conformations représentatives"
echo "  figures/conformation_*.pdb        à ouvrir dans PyMOL ou ChimeraX"
