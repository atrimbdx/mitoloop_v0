#!/usr/bin/env bash
# Chemin 1 : téléchargements et préparation du paquet sur Curta (machine connectée). Aucune donnée patient.
# Usage : bash telecharger_curta.sh <dossier_donnees> <dossier_paquet>
# Réutilise les fichiers déjà présents sur Curta (variables ci-dessous) ; ne télécharge que ce qui manque.
#   FASTA          génome hg38, noms chr1..chrM, avec .fai (défaut : copie Broad téléchargée)
#   MITOCARTA_XLS  Human.MitoCarta3.0.xls (défaut : téléchargé)
#   PANELAPP       liste PanelApp « Mitochondrial disorders » (défaut : /scratch/$USER/MITODARK/data/panelapp_genes.txt)
#   GTEX_ATTR      GTEx_v8 SampleAttributesDS (défaut : téléchargé)
#   GTEX_GCT       TPM GTEx v8 par échantillon (défaut : téléchargé, environ 3 Go ; SAUTER_GTEX=1 pour sauter l'étape 8)
set -euo pipefail
D="${1:?dossier_donnees}"
P="${2:?dossier_paquet}"
ICI="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$D/hgnc" "$D/ref" "$D/mirbase" "$D/mitocarta" "$D/intronic" "$D/gtrnadb" "$D/gtex" "$P"
FASTA="${FASTA:-$D/ref/Homo_sapiens_assembly38.fasta}"
MITOCARTA_XLS="${MITOCARTA_XLS:-$D/mitocarta/Human.MitoCarta3.0.xls}"
PANELAPP="${PANELAPP:-/scratch/$USER/MITODARK/data/panelapp_genes.txt}"
GTEX_ATTR="${GTEX_ATTR:-$D/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt}"
GTEX_GCT="${GTEX_GCT:-$D/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz}"
# télécharge seulement si le fichier est absent ou vide
tele() { if [ -s "$2" ]; then echo "déjà présent : $2"; else echo "téléchargement : $1"; curl -fsSL -o "$2" "$1"; fi; }
[ -s "$PANELAPP" ] || { echo "PanelApp absent : $PANELAPP (variable PANELAPP)"; exit 1; }

# 1. HGNC (symboles à jour ; testé le 19/11/2026 depuis la session)
tele https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt "$D/hgnc/hgnc_complete_set.txt"

# 2. Génome hg38 et annotation GENCODE v27 au format refFlat (seau public Broad ; testé)
tele https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta "$FASTA"
tele https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai "$FASTA.fai"
grep -q -P "^chr1\t" "$FASTA.fai" || { echo "FASTA sans noms chr1..chrM : $FASTA (le refFlat GENCODE utilise chr1)"; exit 1; }
tele https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/GRCh38_gencode.v27.refFlat.txt "$D/ref/GRCh38_gencode.v27.refFlat.txt"

# 3. miRBase (officiel, v22.1). A VERIFIER : mirbase.org refusé depuis la session ; la session a utilisé la copie
#    nf-core/test-datasets (branche smrnaseq, reference/{mature,hairpin}.fa, humain seulement, 1881 épingles).
tele https://www.mirbase.org/download/mature.fa "$D/mirbase/mature.fa"
tele https://www.mirbase.org/download/hairpin.fa "$D/mirbase/hairpin.fa"

# 4. MitoCarta3.0 (officiel) ; PanelApp « Mitochondrial disorders » : reprendre la liste déjà présente sur Curta
#    (/scratch/atrimouille/MITODARK/data/panelapp_genes.txt et sa provenance JSON). A VERIFIER : URL Broad et numéro de la feuille « A Human MitoCarta3.0 » (sheet_name=1 supposé).
tele https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls "$MITOCARTA_XLS"
python3 -c "import pandas as pd, sys; pd.read_excel(sys.argv[1], sheet_name="A Human MitoCarta3.0").to_csv(sys.argv[2], sep='\t', index=False)" "$MITOCARTA_XLS" "$D/mitocarta/Human.MitoCarta3.0.tsv"
python3 "$ICI/prep_chemin1.py" mito --hgnc "$D/hgnc/hgnc_complete_set.txt" --mitocarta "$D/mitocarta/Human.MitoCarta3.0.tsv" --panelapp "$PANELAPP" --out "$P/genes_mito.tsv"

# 5. Tables du paquet
python3 "$ICI/prep_chemin1.py" utr3 --refflat "$D/ref/GRCh38_gencode.v27.refFlat.txt" --fasta "$FASTA" --hgnc "$D/hgnc/hgnc_complete_set.txt" --out "$P/utr3.tsv.gz"
python3 "$ICI/prep_chemin1.py" mirna --refflat "$D/ref/GRCh38_gencode.v27.refFlat.txt" --fasta "$FASTA" --hgnc "$D/hgnc/hgnc_complete_set.txt" --hairpin "$D/mirbase/hairpin.fa" --mature "$D/mirbase/mature.fa" --out "$P/mirna.tsv.gz"
python3 "$ICI/prep_chemin1.py" epissage --refflat "$D/ref/GRCh38_gencode.v27.refFlat.txt" --fasta "$FASTA" --hgnc "$D/hgnc/hgnc_complete_set.txt" --out "$P/sites_epissage.tsv.gz"
python3 "$ICI/prep_chemin1.py" elements --refflat "$D/ref/GRCh38_gencode.v27.refFlat.txt" --fasta "$FASTA" --hgnc "$D/hgnc/hgnc_complete_set.txt" --out "$P/elements_ncrna.tsv.gz"
python3 "$ICI/prep_chemin1.py" partenaires --hgnc "$D/hgnc/hgnc_complete_set.txt" --mito "$P/genes_mito.tsv" --out "$P/partenaires_pol3_arnt_arnr.tsv"

# 6. Introns U12 génome entier (intronIC 3.0.0, PyPI ; scikit-learn 1.7.2 pour charger le modèle pré-entraîné sans avertissement)
python3 -m pip install --user intronIC==3.0.0 scikit-learn==1.7.2
python3 "$ICI/prep_chemin1.py" gtf --refflat "$D/ref/GRCh38_gencode.v27.refFlat.txt" --hgnc "$D/hgnc/hgnc_complete_set.txt" --out "$D/ref/gencode27_codant.gtf"
intronIC classify -g "$FASTA" -a "$D/ref/gencode27_codant.gtf" -n homo_sapiens -p 4 -o "$D/intronic"
python3 "$ICI/prep_chemin1.py" u12 --hgnc "$D/hgnc/hgnc_complete_set.txt" --table "$ICI/../../../introns_U12_genes_mito_hg38.tsv" --mito "$P/genes_mito.tsv" --intronic "$D/intronic/homo_sapiens.meta.iic" --gtf "$D/ref/gencode27_codant.gtf" --outdir "$P"

# 7. ARNt nucléaires (GtRNAdb hg38). A VERIFIER : URL ; convertir le BED en lignes d'elements_ncrna (classe ARNt) avant usage.
tele https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz "$D/gtrnadb/hg38-tRNAs.tar.gz"

# 8. GTEx v8 (co-expression lncRNA / OXPHOS) : seau public adult-gtex, lu en HTTPS (testé le 20/11/2026 ; gsutil inutile)
if [ "${SAUTER_GTEX:-0}" = 1 ]; then echo "étape 8 (GTEx) sautée"; else
tele https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz "$GTEX_GCT"
tele https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt "$GTEX_ATTR"
python3 "$ICI/prep_chemin1.py" coexpression --hgnc "$D/hgnc/hgnc_complete_set.txt" --gct "$GTEX_GCT" --attributs "$GTEX_ATTR" --mitocarta "$D/mitocarta/Human.MitoCarta3.0.tsv" --out "$P/coexpression_lncrna_oxphos.tsv"
fi

# 9. snoDB et boucles mitoloop : formats à inspecter avant conversion (voir README, section « squelettes »).
#    snoDB : télécharger l'export des interactions depuis le site de snoDB (URL exacte à relever), puis :
#    python3 "$ICI/prep_chemin1.py" snodb --hgnc "$D/hgnc/hgnc_complete_set.txt" --table <export_snodb.tsv> --col-snorna <colonne> --col-cible <colonne> --out "$P/snodb_cibles.tsv"
#    Boucles : /scratch/atrimouille/MITODARK/loops/loops.tsv (ou assoc_loops.tsv), puis :
#    python3 "$ICI/prep_chemin1.py" boucles --hgnc "$D/hgnc/hgnc_complete_set.txt" --table /scratch/atrimouille/MITODARK/loops/loops.tsv --col-gene <colonne> --col-chrom <colonne> --col-debut <colonne> --col-fin <colonne> --out "$P/boucles_mito.tsv"
echo "Paquet prêt dans $P"
