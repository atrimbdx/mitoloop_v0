#!/usr/bin/env bash
# Chemin 3 (micropeptides) : téléchargements sur Curta, machine connectée.
# Aucune donnée patient n'est utilisée ni envoyée. Usage : bash telecharger_curta.sh [dossier]
# Les URL marquées « à vérifier » n'ont pas pu être testées depuis le conteneur de
# développement (ftp.ebi.ac.uk, broadinstitute.org et openprot.org y étaient bloqués).
set -uo pipefail

DEST="${1:-$PWD/donnees_chemin3}"
GENCODE_REL="${GENCODE_REL:-46}"   # aligner sur la version GENCODE/VEP de la VM
mkdir -p "$DEST/riboseq" "$DEST/openprot" "$DEST/gencode" "$DEST/mitocarta" "$DEST/outils"
cd "$DEST" || exit 1
date -u +"%Y-%m-%dT%H:%MZ" > DATE_TELECHARGEMENT.txt

echo "== 1. Fichiers d'entrée de la phase I des Ribo-seq ORFs (GitHub, licence MIT ; testé)"
if [ ! -d riboseq/gencode-riboseqORFs ]; then
  git clone --depth 1 https://github.com/jorruior/gencode-riboseqORFs.git riboseq/gencode-riboseqORFs
fi
git -C riboseq/gencode-riboseqORFs rev-parse HEAD >> DATE_TELECHARGEMENT.txt

echo "== 2. Catalogue GENCODE final des Ribo-seq ORFs (Mudge et al. 2022) : répertoire EBI (à vérifier)"
wget -q -r -np -nd -e robots=off -A "*.bed,*.bed.gz,*.gtf,*.gtf.gz,*.fa,*.fa.gz,*.fasta,*.txt,*.tsv,*.tsv.gz" -P riboseq/gencode_final "https://ftp.ebi.ac.uk/pub/databases/gencode/riboseq_orfs/" || echo "ÉCHEC : chercher le lien sur https://www.gencodegenes.org (rubrique Ribo-seq ORFs)"
ls -la riboseq/gencode_final 2>/dev/null

echo "== 3. GENCODE : annotation GTF et génome hg38"
wget -q -c -P gencode "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_REL}/gencode.v${GENCODE_REL}.annotation.gtf.gz"
wget -q -c -P gencode "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_REL}/GRCh38.primary_assembly.genome.fa.gz"
[ -f gencode/GRCh38.primary_assembly.genome.fa ] || gunzip -k gencode/GRCh38.primary_assembly.genome.fa.gz

echo "== 4. MitoCarta3.0 humaine (URL à vérifier)"
wget -q -c -P mitocarta "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls" || echo "ÉCHEC : télécharger depuis la page MitoCarta3.0 du Broad Institute"

echo "== 5. OpenProt : téléchargement manuel"
echo "Sur https://www.openprot.org (rubrique Downloads), prendre pour Homo sapiens, GRCh38 :"
echo "  les AltProts en BED (ou GTF) et leurs séquences protéiques FASTA ; les déposer dans $DEST/openprot/."
echo "  Le nom exact des fichiers dépend de la version d'OpenProt : noter la version dans DATE_TELECHARGEMENT.txt."

echo "== 6. Outils Python"
python3 -m pip install --user pyfaidx xlrd pandas
echo "TMbed (licence Apache-2.0, modèle ProtT5 téléchargé au premier lancement, plusieurs Go) :"
echo "  python3 -m pip install --user git+https://github.com/BernhoferM/TMbed.git"
echo "TargetP 2.0 et DeepTMHMM : licence académique DTU, téléchargement après inscription (services.healthtech.dtu.dk)."
echo "MitoFates : archive Perl à télécharger depuis le site des auteurs (voir Fukasawa et al. 2015)."
echo "Fin. Contenu : $DEST"
