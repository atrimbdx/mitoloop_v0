#!/usr/bin/env bash
# 07_maj_mitoloop_vm.sh — VM. À lancer depuis le dossier transféré (transfert_mitoloop) :
#   export REANNOT_BASE=/home/jovyan/DURAAK/reannot
#   bash 07_maj_mitoloop_vm.sh
#
# Installe ou met à jour MitoLoop :
#   - données dans  $REANNOT_BASE/data/mitoloop
#   - scripts       $REANNOT_BASE/scripts/{ajout_mitoloop.py, preparer_entree.py, reannot_complet.sh}
# Remplace reannot_alphagenome.sh par sa version accélérée (étape AlphaGenome)
# seulement s'il est exactement la version d'origine, gardée dans archives/.
# Ne touche à rien d'autre : les données déjà installées (alphagenome, annovar,
# morfee, ref, spliceai, spliceai2, vep_cache, vep_plugins), les environnements
# et les autres scripts du pipeline restent tels quels.
#
# La version précédente de MitoLoop (données et ces trois scripts) est archivée
# dans $REANNOT_BASE/archives/mitoloop_<version>.tar.gz ; seules les
# GARDER (défaut 2) archives les plus récentes sont conservées.
set -euo pipefail
SRC=$(cd "$(dirname "$0")" && pwd)
BASE=${REANNOT_BASE:-$HOME/reannot}
GARDER=${GARDER:-2}
MINE="ajout_mitoloop.py preparer_entree.py reannot_complet.sh"
export PATH="$BASE/envs/vep/bin:$PATH"
cd "$SRC"

[ -d "$BASE/data" ] && [ -d "$BASE/scripts" ] || {
  echo "Erreur : $BASE ne ressemble pas au pipeline reannot (data/ et scripts/ attendus)."
  echo "Indiquer le bon dossier : export REANNOT_BASE=/home/jovyan/DURAAK/reannot"; exit 1; }

echo "Vérification des sommes de contrôle..."
sha256sum -c --quiet SHA256SUMS && echo "  OK"

# 1. Archive de la version en place
M="$BASE/data/mitoloop"
if [ -d "$M" ]; then
  mkdir -p "$BASE/archives"
  V=$(head -n 1 "$M/VERSION" 2>/dev/null | tr -c 'A-Za-z0-9._\n-' '_' | tr -d '\n' | head -c 60 || true)
  A="$BASE/archives/mitoloop_${V:-$(date +%Y%m%d_%H%M)}.tar.gz"
  S=(); for f in $MINE; do [ -f "$BASE/scripts/$f" ] && S+=("scripts/$f"); done
  tar -czf "$A" -C "$BASE" data/mitoloop ${S[@]+"${S[@]}"}
  echo "Version précédente archivée : $A"
  # anciennes archives au-delà de GARDER
  ls -1t "$BASE"/archives/mitoloop_*.tar.gz 2>/dev/null | tail -n +$((GARDER + 1)) | while read -r old; do
    rm -f "$old"; echo "  archive supprimée : $(basename "$old")"; done
fi
rm -rf "$BASE/data/mitoloop.ancien"          # ancien format de sauvegarde, remplacé par archives/

# 2. Nouvelle version (copiée à côté, puis mise en place d'un coup)
rm -rf "$M.nouveau"; mkdir -p "$M.nouveau/annotations"
cp genes.tsv regions.bed.gz MANIFEST.tsv VERSION "$M.nouveau/"
cp annotations/* "$M.nouveau/annotations/"
[ -f motifs.jaspar ] && cp motifs.jaspar "$M.nouveau/"
[ -f README_MitoLoop.md ] && cp README_MitoLoop.md "$M.nouveau/"
touch "$M.nouveau"/annotations/*.tbi 2>/dev/null || true   # index plus récent que les données
rm -rf "$M"; mv "$M.nouveau" "$M"
for f in $MINE; do cp "$f" "$BASE/scripts/"; done
chmod +x "$BASE/scripts/reannot_complet.sh"

# 3. reannot_alphagenome.sh : étape 5 (AlphaGenome) lue par l'index sur les seules
#    positions des variants, résultat identique. Appliqué seulement si le script en
#    place est exactement la version d'origine connue ; l'original est gardé dans
#    archives/reannot_alphagenome.sh.origine (jamais supprimé).
AG_ORIG=9b28a24a508e096da42cd5e09a98b5d083d89fb2f03366fad1fd6d90c80ae4f4
AG_NEW=dd3eba4f20ee8990b4de98df5a5052477826e13ec5525076b769e71d9897cc33
if [ -f reannot_alphagenome.sh ] && [ -f "$BASE/scripts/reannot_alphagenome.sh" ]; then
  cur=$(sha256sum "$BASE/scripts/reannot_alphagenome.sh" | cut -d' ' -f1)
  if [ "$cur" = "$AG_ORIG" ]; then
    mkdir -p "$BASE/archives"
    cp -p "$BASE/scripts/reannot_alphagenome.sh" "$BASE/archives/reannot_alphagenome.sh.origine"
    cp reannot_alphagenome.sh "$BASE/scripts/reannot_alphagenome.sh"
    chmod +x "$BASE/scripts/reannot_alphagenome.sh"
    echo "reannot_alphagenome.sh : étape AlphaGenome accélérée (original : archives/reannot_alphagenome.sh.origine)"
  elif [ "$cur" = "$AG_NEW" ]; then
    echo "reannot_alphagenome.sh : déjà la version accélérée"
  else
    echo "reannot_alphagenome.sh : version locale modifiée, non remplacée."
    echo "  Pour accélérer AlphaGenome, voir l'étape 5 de $SRC/reannot_alphagenome.sh"
  fi
fi

# 4. Contrôle
echo "Contrôle de lecture :"
python3 - "$M" <<'PY'
import gzip, sys
from pathlib import Path
m = Path(sys.argv[1])
n = sum(1 for _ in open(m / "genes.tsv")) - 1
r = sum(1 for _ in gzip.open(m / "regions.bed.gz", "rt"))
a = sorted(p.name for p in (m / "annotations").iterdir() if not p.name.endswith(".tbi"))
print(f"  {n} gènes, {r} régions, annotations : {', '.join(a)}")
PY
echo "Version installée : $(cat "$M/VERSION")"
echo "Données du pipeline non modifiées : $(ls "$BASE/data" | grep -v '^mitoloop' | tr '\n' ' ')"
echo
echo "Utilisation :"
echo "  bash $BASE/scripts/reannot_complet.sh variants.xlsx run1      (tout le pipeline)"
echo "  python3 $BASE/scripts/ajout_mitoloop.py run1.final.tsv          (MitoLoop seul)"
echo "Le dossier transféré ($SRC) peut être supprimé."
