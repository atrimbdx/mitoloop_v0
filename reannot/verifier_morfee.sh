#!/usr/bin/env bash
# verifier_morfee.sh — VM. Compare l'ancien et le nouveau run_morfee.R sur la
# sortie ANNOVAR d'un run déjà fait : même résultat attendu, et les deux durées.
#
# Usage : bash verifier_morfee.sh PREFIXE      (PREFIXE.annovar.hg38_multianno.vcf
#                                               doit exister, laissé par un run)
# THREADS : nombre de processus de la nouvelle version (défaut 4).
set -euo pipefail
P=${1:?Usage: $0 PREFIXE_DU_RUN}
BASE=${REANNOT_BASE:-$HOME/reannot}
RS=$BASE/envs/morfee/bin/Rscript
RD=$BASE/data/morfee/MORFEE_DATA_GRCh38.RData
F=MORFEE_uTIS,MORFEE_uSTOP,MORFEE_uSTOP_creation,MORFEE_dTIS,MORFEE_dSTOP,MORFEE_dSTOP_creation
IN="$P.annovar.hg38_multianno.vcf"
ORIG="$BASE/archives/run_morfee.R.origine"
NEW="$BASE/scripts/run_morfee.R"
[ -f "$IN" ] || { echo "Erreur : $IN introuvable (lancer d'abord un run avec ce préfixe)"; exit 1; }
[ -f "$ORIG" ] || { echo "Erreur : $ORIG introuvable (installé par 07_maj_mitoloop_vm.sh)"; exit 1; }
for v in ancien nouveau; do
  s=$([ $v = ancien ] && echo "$ORIG" || echo "$NEW")
  t0=$(date +%s)
  "$RS" "$s" "$IN" "$RD" "$P.verif_$v.tsv" "$F" > "$P.verif_$v.log" 2>&1
  echo "$v : $(( $(date +%s) - t0 )) s, $(wc -l < "$P.verif_$v.tsv") variants annotés"
done
if diff <(sort "$P.verif_ancien.tsv") <(sort "$P.verif_nouveau.tsv") > "$P.verif_diff.txt"; then
  echo "Résultats identiques."
else
  echo "DIFFÉRENCES : voir $P.verif_diff.txt (et $P.verif_*.log)"; exit 1
fi
