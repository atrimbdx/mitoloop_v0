#!/bin/bash
#SBATCH --job-name=simpol_ctcf
#SBATCH --time=00:40:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_logs/ctcf_%j.out
# Voie A, entree : sites CTCF orientes depuis ENCODE.
# Une seule tache, et une fois pour toutes.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh
mkdir -p ctcf slurm_logs

python3 - <<'PY' > /tmp/bornes_$$.txt
import config as cfg, os
c = cfg.LOCI[os.environ.get("LOCUS", "tyrp1")]
n, d0, f0 = cfg.geometrie(os.environ.get("LOCUS", "tyrp1"))
print(c["chrom"], d0, f0)
PY
read CHROM DEB FIN < /tmp/bornes_$$.txt
rm -f /tmp/bornes_$$.txt
echo "  fenêtre simulée : $CHROM:$DEB-$FIN"

for C in MEL FIBRO; do
  BIO=$($PY -c "import config,os;print(config.LOCI[os.environ['LOCUS']]['encode']['$C'])")
  echo "=== $C : ENCODE « $BIO »"
  $PY ctcf_encode.py --biosample "$BIO" --genome "$GENOME" \
      --chrom "$CHROM" --debut "$DEB" --fin "$FIN" \
      --sortie "ctcf/${C}_ctcf_oriente.bed"
done
