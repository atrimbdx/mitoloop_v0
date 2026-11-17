#!/bin/bash
#SBATCH --job-name=simpol_essai
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_logs/essai_%j.out
# Une trajectoire courte, sur un noeud de calcul, avant d'en lancer 1080.
# Elle verifie ce que le noeud de login ne peut pas verifier : qu'OpenMM
# demarre la ou le calcul aura lieu, que les chemins de config.sh sont
# lisibles depuis un noeud, et que les deux voies construisent bien leurs
# entrees. Le job array en depend : s'il echoue, rien n'est soumis.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh
mkdir -p essai slurm_logs
export OPENMM_CPU_THREADS=1

echo "=== noeud : $(hostname)"
$PY -c "import openmm, numpy, scipy, pandas, matplotlib
print('  OpenMM', openmm.version.version)
print('  plateformes', [openmm.Platform.getPlatform(i).getName()
                        for i in range(openmm.Platform.getNumPlatforms())])"

for V in A B; do
  if [ "$V" = "A" ]; then OPT="--ctcf $CTCF_MEL"; else OPT=""; fi
  echo "=== essai voie $V"
  $PY une_trajectoire.py --locus "$LOCUS" --condition MEL --voie "$V" \
      --jeu 0 --replicat 0 --matrice "$MAT_MEL" $OPT \
      --blocs 40 --pas-par-bloc 1000 --sortie essai
done

echo "=== essai concluant : les deux voies démarrent sur un nœud de calcul"
