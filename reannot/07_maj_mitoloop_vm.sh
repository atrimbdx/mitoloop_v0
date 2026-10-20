#!/usr/bin/env bash
# 07_maj_mitoloop_vm.sh — VM. À lancer depuis le dossier transféré (transfert_mitoloop) :
#   bash 07_maj_mitoloop_vm.sh
# Installe les données MitoLoop et le script ajout_mitoloop.py.
# Le pipeline existant (VEP, SpliceAI, SpliceAI2, MORFEE, AlphaGenome) n'est pas modifié ;
# une version précédente des données MitoLoop est remplacée.
set -euo pipefail
SRC=$(cd "$(dirname "$0")" && pwd)
BASE=${REANNOT_BASE:-$HOME/reannot}
export PATH="$BASE/envs/vep/bin:$PATH"
cd "$SRC"

echo "Vérification des sommes de contrôle..."
sha256sum -c SHA256SUMS

M="$BASE/data/mitoloop"
rm -rf "$M.nouveau"; mkdir -p "$M.nouveau/annotations"
cp genes.tsv regions.bed.gz MANIFEST.tsv VERSION "$M.nouveau/"
cp annotations/* "$M.nouveau/annotations/"
touch "$M.nouveau"/annotations/*.tbi 2>/dev/null || true   # index plus récent que les données
[ -f motifs.jaspar ] && cp motifs.jaspar "$M.nouveau/"
[ -f README_MitoLoop.md ] && cp README_MitoLoop.md "$M.nouveau/"
rm -rf "$M.ancien"; [ -d "$M" ] && mv "$M" "$M.ancien"
mv "$M.nouveau" "$M"
cp ajout_mitoloop.py "$BASE/scripts/"

echo "Contrôle de lecture :"
python3 - "$M" <<'PY'
import gzip, sys
from pathlib import Path
m = Path(sys.argv[1])
n = sum(1 for _ in open(m / "genes.tsv")) - 1
r = sum(1 for _ in gzip.open(m / "regions.bed.gz", "rt"))
a = sorted(p.name[:-7] for p in (m / "annotations").glob("*.bed.gz"))
print(f"  {n} gènes, {r} régions, annotations : {', '.join(a)}")
PY
echo "Version : $(cat "$M/VERSION")"
echo "Mise à jour MitoLoop installée dans $M"
echo "Usage : python3 $BASE/scripts/ajout_mitoloop.py run.final.tsv  (ou run.final.spliceai2.tsv)"
