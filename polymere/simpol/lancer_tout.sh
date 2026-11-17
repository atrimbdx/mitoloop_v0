#!/bin/bash
# =====================================================================
# Lance la campagne entiere, du depart a la figure finale.
#
# Les etapes qui ont besoin du reseau sont faites ici, sur le noeud de
# login : creation de l'environnement, telechargement du chromosome et
# des pics ENCODE. Beaucoup de clusters coupent la sortie internet sur
# les noeuds de calcul, et une campagne qui echoue a l'etape 1 apres
# plusieurs heures de file est une campagne perdue.
#
# Les etapes de calcul partent dans Slurm, chainees par dependance : si
# l'une echoue, les suivantes ne demarrent pas et restent visibles dans
# squeue avec la raison DependencyNeverSatisfied.
#
#   00_essai        une trajectoire courte, sur un noeud       ~5 min
#   02_trajectoires le job array                               ~3,3 h par tache
#   03_comparer     agregation et metriques                    ~1 h
#   04_verifier     controles et rendu 3D                      ~1 h
#
# Le script est relancable : chaque etape deja faite est sautee.
#
# USAGE
#   bash lancer_tout.sh              tout faire
#   bash lancer_tout.sh --preparer   s'arreter avant la soumission
#   bash lancer_tout.sh --replicats 5 --throttle 100
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"

REPLICATS=20
THROTTLE=200
PREPARER=0
while [ $# -gt 0 ]; do
  case "$1" in
    --replicats) REPLICATS="$2"; shift 2 ;;
    --throttle)  THROTTLE="$2";  shift 2 ;;
    --preparer)  PREPARER=1;     shift ;;
    *) echo "option inconnue : $1"; exit 2 ;;
  esac
done

source config.sh
mkdir -p slurm_logs ctcf figures "$RESULTATS"
ENV_DIR="$(dirname "$(dirname "$PY")")"

note() { printf '\n\033[1m== %s\033[0m\n' "$*"; }
echec() { printf '\n  ÉCHEC : %s\n' "$*" >&2; exit 1; }

# Filet : sous set -e, une commande qui echoue sans traitement explicite
# arrete le script sans rien afficher. Ce piege dit au moins ou.
trap 'c=$?; printf "\n  ARRÊT ligne %s, code %s. Rien de plus n'"'"'a été lancé.\n" \
      "$LINENO" "$c" >&2' ERR

# ---------------------------------------------------------------------
note "1/6  environnement Python"
if [ -x "$PY" ] && "$PY" -c "import openmm, scipy, pandas, matplotlib" 2>/dev/null; then
  echo "  déjà présent : $PY"
else
  CREER=$(command -v mamba || command -v conda || true)
  [ -n "$CREER" ] || echec "ni mamba ni conda dans le PATH"
  echo "  création avec $CREER (plusieurs minutes)"
  "$CREER" create -y -p "$ENV_DIR" -c conda-forge --override-channels \
      python=3.11 openmm numpy scipy pandas matplotlib
fi
"$PY" -c "import openmm; print('  OpenMM', openmm.version.version)"

# ---------------------------------------------------------------------
note "2/6  matrices mesurées"
for M in "$MAT_MEL" "$MAT_FIBRO"; do
  [ -s "$M" ] || echec "matrice absente : $M  (produire avec extraire_matrice_locus.sh)"
  echo "  $(wc -l < "$M") lignes  $M"
done

# ---------------------------------------------------------------------
note "3/6  séquence du chromosome"
read -r CHROM DEB FIN <<< "$("$PY" - <<PY
import config as cfg
c = cfg.LOCI["$LOCUS"]; n, d0, f0 = cfg.geometrie("$LOCUS")
print(c["chrom"], d0, f0)
PY
)"
echo "  locus $LOCUS → $CHROM:$DEB-$FIN"

GENOME_DIR="$(dirname "$GENOME")"
mkdir -p "$GENOME_DIR"
if [ -s "$GENOME" ]; then
  echo "  déjà présent : $GENOME"
else
  # Un seul chromosome suffit au balayage du motif : 40 Mo au lieu des
  # 3 Go du génome entier. Le nom de séquence d'UCSC (>chr9) correspond
  # à celui des matrices.
  NU="${CHROM#chr}"
  UCSC="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/${CHROM}.fa.gz"
  ENS="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${NU}.fa.gz"
  TMP="$GENOME_DIR/.${CHROM}.fa.gz.partiel"
  aspirer() {
    if command -v curl >/dev/null; then curl -fsSL --retry 3 -o "$TMP" "$1"
    elif command -v wget >/dev/null; then wget -q -O "$TMP" "$1"
    else echec "ni curl ni wget"; fi
  }
  rm -f "$TMP"
  if aspirer "$UCSC"; then echo "  UCSC : $UCSC"
  elif aspirer "$ENS"; then echo "  Ensembl : $ENS"
  else
    rm -f "$TMP"
    echec "aucune des deux sources n'est joignable depuis ce nœud.
  Récupérer le fichier depuis une machine qui a accès, puis le déposer :
    curl -L -o ${CHROM}.fa.gz '$UCSC'
    gunzip -c ${CHROM}.fa.gz > $GENOME"
  fi
  # Ensembl nomme ses sequences « 15 », UCSC « chr15 ». Les matrices et
  # config.py utilisent la forme UCSC : on normalise la premiere ligne.
  gunzip -c "$TMP" | sed "1s/^>.*/>${CHROM}/" > "$GENOME"
  rm -f "$TMP"
  echo "  $(du -h "$GENOME" | cut -f1)  $GENOME"
fi
"$PY" -c "
from ctcf_encode import lire_fenetre
s = lire_fenetre('$GENOME', '$CHROM', $DEB, $FIN)
assert len(s) > 0, 'fenêtre vide'
assert set(s) - {'N'}, 'la fenêtre ne contient que des N : assemblage incohérent'
print('  séquence lisible : %s:%d-%d, %d pb, %.1f%% de N'
      % ('$CHROM', $DEB, $FIN, len(s), 100 * s.count('N') / len(s)))"

# ---------------------------------------------------------------------
note "4/6  sites CTCF orientés (ENCODE)"
for C in MEL FIBRO; do
  BED="ctcf/${C}_ctcf_oriente.bed"
  if [ -s "$BED" ]; then echo "  déjà présent : $BED ($(wc -l < "$BED") sites)"; continue; fi

  # Un fichier de pics designe dans config.sh court-circuite l'API.
  eval "URL=\${PICS_$C:-}"
  SRC=()
  if [ -n "$URL" ]; then
    PICS="ctcf/${C}_pics_encode.bed.gz"
    if [ ! -s "$PICS" ]; then
      echo "  $C ← $URL"
      if command -v curl >/dev/null; then curl -fsSL --retry 3 -o "$PICS" "$URL"
      else wget -q -O "$PICS" "$URL"; fi || echec "téléchargement impossible : $URL"
    fi
    # Le balayage se cale sur le sommet du pic, colonne 10 du narrowPeak.
    # Un fichier a 9 colonnes retombe sur le milieu du pic, ce qui elargit
    # la fenetre balayee : autant le savoir.
    #
    # Un seul parcours, et surtout aucune sortie anticipee dans le tube :
    # « zcat | head -1 » faisait recevoir un SIGPIPE a zcat des que le
    # fichier depassait quelques milliers de lignes, et pipefail propageait
    # ce 141 jusqu'a set -e, qui tuait le script sans rien afficher.
    read -r NPICS NCOL < <(zcat -f "$PICS" | awk 'NR==1{c=NF} END{print NR, c+0}')
    [ "${NPICS:-0}" -gt 0 ] || echec "fichier de pics vide ou illisible : $PICS"
    echo "  $NPICS pics, $NCOL colonnes"
    [ "$NCOL" -ge 10 ] || echo "  ATTENTION moins de 10 colonnes : pas un narrowPeak,
  le sommet du pic est inconnu et le milieu sera utilisé à la place."
    SRC=(--pics "$PICS")
  else
    BIO=$("$PY" -c "import config;print(config.LOCI['$LOCUS']['encode']['$C'])")
    echo "  $C ← recherche ENCODE « $BIO »"
    # tableau, et non chaine : le libelle contient une espace
    SRC=(--biosample "$BIO")
  fi

  "$PY" ctcf_encode.py "${SRC[@]}" --genome "$GENOME" \
      --chrom "$CHROM" --debut "$DEB" --fin "$FIN" --sortie "$BED" \
    || echec "ctcf_encode.py a échoué pour $C.
  Si ENCODE n'a pas ce type cellulaire, voir la liste de ce qui existe :
    $PY ctcf_encode.py --lister
  puis corriger le champ 'encode' du locus dans config.py.
  Si c'est le nœud de login qui n'a pas d'accès sortant :
    sbatch slurm/01_ctcf.sbatch"
done
# Les deux conditions doivent recevoir des sites differents, sans quoi la
# voie A n'a aucune entree qui les distingue et son differentiel MEL/FIBRO
# est un nul par construction. verifier.py le signalera aussi.
if [ "$(md5sum < ctcf/MEL_ctcf_oriente.bed)" = "$(md5sum < ctcf/FIBRO_ctcf_oriente.bed)" ]; then
  echo "  ATTENTION les deux BED CTCF sont identiques : le différentiel de"
  echo "  la voie A sera un nul par construction. La comparaison des voies"
  echo "  et les contrôles par condition restent valables."
fi

# ---------------------------------------------------------------------
note "5/6  plan de campagne"
# la ligne « soumettre » de plan.py est filtree : ici c'est ce script qui
# soumet, et la recopier ferait partir le job array une seconde fois
# sortie capturee puis filtree : « | grep -v » renvoie 1 quand il ne
# laisse passer aucune ligne, ce qui suffirait a arreter le script
PLAN=$("$PY" plan.py --locus "$LOCUS" --replicats "$REPLICATS" --sortie plan.tsv)
printf '%s\n' "$PLAN" | grep -v "soumettre" || true
N=$(($(wc -l < plan.tsv) - 1))
echo "  $N tâches, $REPLICATS réplicats, au plus $THROTTLE simultanées"

if [ "$PREPARER" = 1 ]; then
  note "préparation terminée, rien n'a été soumis"
  echo "  pour soumettre : bash lancer_tout.sh"
  exit 0
fi

# ---------------------------------------------------------------------
note "6/6  soumission"
J0=$(sbatch --parsable slurm/00_essai.sbatch)
echo "  essai                 $J0"
J1=$(sbatch --parsable --dependency=afterok:"$J0" \
     --array=1-"$N"%"$THROTTLE" slurm/02_trajectoires.sbatch)
echo "  trajectoires          $J1   (array 1-$N%$THROTTLE)"
J2=$(sbatch --parsable --dependency=afterok:"$J1" slurm/03_comparer.sbatch)
echo "  comparaison           $J2"
J3=$(sbatch --parsable --dependency=afterok:"$J2" slurm/04_verifier.sbatch)
echo "  contrôles et rendu    $J3"

cat <<FIN

  Suivi
    squeue -u \$USER
    tail -f slurm_logs/essai_$J0.out

  Si l'essai échoue, les trois autres restent en attente avec la raison
  DependencyNeverSatisfied. Corriger, puis les annuler et relancer :
    scancel $J1 $J2 $J3 && bash lancer_tout.sh

  À l'arrivée
    figures/verification.png       les trois contrôles
    figures/conformation_*.png     conformations représentatives
    figures/conformation_*.pdb     à ouvrir dans PyMOL ou ChimeraX
    bilan.png, bilan.tsv           voie A contre voie B
FIN
