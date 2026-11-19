#!/usr/bin/env bash
# Chemin 2 - preparation sur Curta (machine connectee). Aucune donnee patient.
# Produit le paquet a transferer sur la VM : donnees/arn_coordonnees_hg38.tsv, donnees/ref_regions_hg38.fa
# (+ fichiers de verification Rfam et ClinVar a relire a la main).
# Usage : bash curta_preparer_chemin2.sh DOSSIER_TRAVAIL
set -euo pipefail
W="${1:?dossier de travail}"
ICI="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$W/dl" "$W/paquet"
cd "$W/dl"

# 1. Annotation GENCODE a jour (la version du depot a ete construite avec GENCODE v27, faute d'acces a l'EBI).
wget -nc https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz

# 2. Genome hg38 (meme FASTA que celui utilise pour construire le paquet du depot).
wget -nc https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
wget -nc https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai

# 3. Coordonnees et sequences de reference.
python3 "$ICI/preparer_references.py" coordonnees --gtf gencode.v49.primary_assembly.annotation.gtf.gz --genes "$ICI/../donnees/genes_chemin2.txt" --out "$W/paquet/arn_coordonnees_hg38.tsv" --source "GENCODE v49 primary_assembly"
python3 "$ICI/preparer_references.py" sequences --coord "$W/paquet/arn_coordonnees_hg38.tsv" --fasta Homo_sapiens_assembly38.fasta --out "$W/paquet/ref_regions_hg38.fa"

# 4. Rfam + Infernal : positions structurales (helice I et boucle E de l'ARNr 5S, RF00001 ; RMRP ; RPPH1)
#    pour remplacer les coordonnees derivees par ViennaRNA dans elements_mito.tsv (relecture manuelle).
wget -nc https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip -kf Rfam.cm.gz
cmpress -F Rfam.cm
python3 -c "import sys; [print(l.strip()) for l in open(sys.argv[1]) if not l.startswith(';')]" "$W/paquet/ref_regions_hg38.fa" > regions_pour_cmscan.fa
cmscan --cut_ga --rfam --nohmmonly --tblout "$W/paquet/rfam_regions_chemin2.tblout" Rfam.cm regions_pour_cmscan.fa > /dev/null

# 5. ClinVar : coordonnee hg38 exacte du fondateur RMRP n.71A>G / n.72A>G (temoin T03) et autres variants RMRP publics.
wget -nc https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
zcat clinvar.vcf.gz | awk -F'\t' '$1=="9" && $2>=35657700 && $2<=35658100' > "$W/paquet/clinvar_RMRP_hg38.vcf"

# 6. ViennaRNA pour la VM (roue a transferer si la VM n'a pas pip en ligne).
python3 -m pip download ViennaRNA==2.8.0 -d "$W/paquet/roues"

# 7. Donnees APEX matrice (Fazal 2019, Kaewsapsak 2017) : relever les numeros GEO dans la section
#    "Data and code availability" des articles (non recopies ici car non verifies), puis telecharger les
#    tables de comptes et verifier l'enrichissement de RMRP, RPPH1, RNA5S*, TERC, MALAT1 dans la matrice.
#    Rappel : APEX-RIP utilise les lncRNA nucleaires comme faux positifs (circularite) et sous-echantillonne les petits ARN.

echo "Paquet pret dans $W/paquet : copier arn_coordonnees_hg38.tsv et ref_regions_hg38.fa dans chemin2/donnees/ apres relecture."
