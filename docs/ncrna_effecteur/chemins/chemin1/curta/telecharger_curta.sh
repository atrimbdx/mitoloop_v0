#!/usr/bin/env bash
# Chemin 1 : téléchargements et préparation du paquet sur Curta (machine connectée). Aucune donnée patient.
# Usage : bash telecharger_curta.sh <dossier_donnees> <dossier_paquet>
# Chaque commande tient sur une ligne. Les étapes marquées « A VERIFIER » ont une URL non testée depuis la session
# de développement (hôte refusé par le proxy) : vérifier l'URL avant de lancer.
set -euo pipefail
D="${1:?dossier_donnees}"
P="${2:?dossier_paquet}"
ICI="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$D/hgnc" "$D/ref" "$D/mirbase" "$D/mitocarta" "$D/intronic" "$D/gtrnadb" "$D/gtex" "$P"

# 1. HGNC (symboles à jour ; testé le 19/11/2026 depuis la session)
curl -fsSL -o "$D/hgnc/hgnc_complete_set.txt" https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt

# 2. Génome hg38 et annotation GENCODE v27 au format refFlat (seau public Broad ; testé)
curl -fsSL -o "$D/ref/Homo_sapiens_assembly38.fasta" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
curl -fsSL -o "$D/ref/Homo_sapiens_assembly38.fasta.fai" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
curl -fsSL -o "$D/ref/GRCh38_gencode.v27.refFlat.txt" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/GRCh38_gencode.v27.refFlat.txt

# 3. miRBase (officiel, v22.1). A VERIFIER : mirbase.org refusé depuis la session ; la session a utilisé la copie
#    nf-core/test-datasets (branche smrnaseq, reference/{mature,hairpin}.fa, humain seulement, 1881 épingles).
curl -fsSL -o "$D/mirbase/mature.fa" https://www.mirbase.org/download/mature.fa
curl -fsSL -o "$D/mirbase/hairpin.fa" https://www.mirbase.org/download/hairpin.fa

# 4. MitoCarta3.0 (officiel) ; PanelApp « Mitochondrial disorders » : reprendre la liste déjà présente sur Curta
#    (/scratch/atrimouille/MITODARK/data/panelapp_genes.txt et sa provenance JSON). A VERIFIER : URL Broad et numéro de la feuille « A Human MitoCarta3.0 » (sheet_name=1 supposé).
curl -fsSL -o "$D/mitocarta/Human.MitoCarta3.0.xls" https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls
python3 -c "import pandas as pd, sys; pd.read_excel(sys.argv[1], sheet_name=1).to_csv(sys.argv[2], sep='\t', index=False)" "$D/mitocarta/Human.MitoCarta3.0.xls" "$D/mitocarta/Human.MitoCarta3.0.tsv"
python3 "$ICI/prep_chemin1.py" mito --hgnc "$D/hgnc/hgnc_complete_set.txt" --mitocarta "$D/mitocarta/Human.MitoCarta3.0.tsv" --panelapp /scratch/atrimouille/MITODARK/data/panelapp_genes.txt --out "$P/genes_mito.tsv"

# 5. Tables du paquet
python3 "$ICI/prep_chemin1.py" utr3 --refflat "$D/ref/GRCh38_gencode.v27.refFlat.txt" --fasta "$D/ref/Homo_sapiens_assembly38.fasta" --hgnc "$D/hgnc/hgnc_complete_set.txt" --out "$P/utr3.tsv.gz"
python3 "$ICI/prep_chemin1.py" mirna --refflat "$D/ref/GRCh38_gencode.v27.refFlat.txt" --fasta "$D/ref/Homo_sapiens_assembly38.fasta" --hgnc "$D/hgnc/hgnc_complete_set.txt" --hairpin "$D/mirbase/hairpin.fa" --mature "$D/mirbase/mature.fa" --out "$P/mirna.tsv.gz"
python3 "$ICI/prep_chemin1.py" epissage --refflat "$D/ref/GRCh38_gencode.v27.refFlat.txt" --fasta "$D/ref/Homo_sapiens_assembly38.fasta" --hgnc "$D/hgnc/hgnc_complete_set.txt" --out "$P/sites_epissage.tsv.gz"
python3 "$ICI/prep_chemin1.py" elements --refflat "$D/ref/GRCh38_gencode.v27.refFlat.txt" --fasta "$D/ref/Homo_sapiens_assembly38.fasta" --hgnc "$D/hgnc/hgnc_complete_set.txt" --out "$P/elements_ncrna.tsv.gz"
python3 "$ICI/prep_chemin1.py" partenaires --hgnc "$D/hgnc/hgnc_complete_set.txt" --mito "$P/genes_mito.tsv" --out "$P/partenaires_pol3_arnt_arnr.tsv"

# 6. Introns U12 génome entier (intronIC 3.0.0, PyPI ; scikit-learn 1.7.2 pour charger le modèle pré-entraîné sans avertissement)
python3 -m pip install --user intronIC==3.0.0 scikit-learn==1.7.2
python3 "$ICI/prep_chemin1.py" gtf --refflat "$D/ref/GRCh38_gencode.v27.refFlat.txt" --hgnc "$D/hgnc/hgnc_complete_set.txt" --out "$D/ref/gencode27_codant.gtf"
intronIC classify -g "$D/ref/Homo_sapiens_assembly38.fasta" -a "$D/ref/gencode27_codant.gtf" -n homo_sapiens -p 4 -o "$D/intronic"
python3 "$ICI/prep_chemin1.py" u12 --hgnc "$D/hgnc/hgnc_complete_set.txt" --table "$ICI/../../../introns_U12_genes_mito_hg38.tsv" --mito "$P/genes_mito.tsv" --intronic "$D/intronic/homo_sapiens.meta.iic" --gtf "$D/ref/gencode27_codant.gtf" --outdir "$P"

# 7. ARNt nucléaires (GtRNAdb hg38). A VERIFIER : URL ; convertir le BED en lignes d'elements_ncrna (classe ARNt) avant usage.
curl -fsSL -o "$D/gtrnadb/hg38-tRNAs.tar.gz" https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz

# 8. GTEx v8 (co-expression lncRNA / OXPHOS) : seau public adult-gtex. A VERIFIER : noms exacts des fichiers (gsutil ls gs://adult-gtex/bulk-gex/v8/rna-seq/)
gsutil cp gs://adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz "$D/gtex/"
gsutil cp gs://adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt "$D/gtex/"
python3 "$ICI/prep_chemin1.py" coexpression --hgnc "$D/hgnc/hgnc_complete_set.txt" --gct "$D/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz" --attributs "$D/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt" --mitocarta "$D/mitocarta/Human.MitoCarta3.0.tsv" --out "$P/coexpression_lncrna_oxphos.tsv"

# 9. snoDB et boucles mitoloop : formats à inspecter avant conversion (voir README, section « squelettes »).
#    snoDB : télécharger l'export des interactions depuis le site de snoDB (URL exacte à relever), puis :
#    python3 "$ICI/prep_chemin1.py" snodb --hgnc "$D/hgnc/hgnc_complete_set.txt" --table <export_snodb.tsv> --col-snorna <colonne> --col-cible <colonne> --out "$P/snodb_cibles.tsv"
#    Boucles : /scratch/atrimouille/MITODARK/loops/loops.tsv (ou assoc_loops.tsv), puis :
#    python3 "$ICI/prep_chemin1.py" boucles --hgnc "$D/hgnc/hgnc_complete_set.txt" --table /scratch/atrimouille/MITODARK/loops/loops.tsv --col-gene <colonne> --col-chrom <colonne> --col-debut <colonne> --col-fin <colonne> --out "$P/boucles_mito.tsv"
echo "Paquet prêt dans $P"
