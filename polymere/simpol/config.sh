# Chemins et reglages de la campagne. A sourcer avant les sbatch.
# Tout ce qui est propre a votre installation est ici, et nulle part ailleurs.

PROJ=${PROJ:-/scratch/atrimouille/chic_albinisme}
SIMPOL=${SIMPOL:-$PROJ/simpol}
PY=${PY:-/scratch/atrimouille/envs/polymere/bin/python}

LOCUS=${LOCUS:-tyrp1}
GENOME=${GENOME:-$PROJ/reference/genome.fa}

# matrices de contacts, produites par extraire_matrice_locus.sh
MAT_MEL=${MAT_MEL:-$PROJ/TYRP1_matrice_MEL_UNSTIM.tsv}
MAT_FIBRO=${MAT_FIBRO:-$PROJ/TYRP1_matrice_FIBRO.tsv}

# Pics CTCF ENCODE, choisis a la main. Renseignes, ils court-circuitent la
# recherche par l'API : c'est le mode sur lequel compter, le libelle d'un
# type cellulaire dans ENCODE n'etant pas stable d'un jeu a l'autre.
#   MEL   : ENCFF004AXJ, epithelium pigmentaire retinien
#   FIBRO : ENCFF204WXC, fibroblaste
# Laisser vide pour revenir a la recherche automatique.
PICS_MEL=${PICS_MEL:-https://www.encodeproject.org/files/ENCFF004AXJ/@@download/ENCFF004AXJ.bed.gz}
PICS_FIBRO=${PICS_FIBRO:-https://www.encodeproject.org/files/ENCFF204WXC/@@download/ENCFF204WXC.bed.gz}

# sites CTCF orientes, produits par 01_ctcf.sbatch ou par lancer_tout.sh
CTCF_MEL=${CTCF_MEL:-$SIMPOL/ctcf/MEL_ctcf_oriente.bed}
CTCF_FIBRO=${CTCF_FIBRO:-$SIMPOL/ctcf/FIBRO_ctcf_oriente.bed}

RESULTATS=${RESULTATS:-$SIMPOL/resultats}
BLOCS=${BLOCS:-2000}
PAS_PAR_BLOC=${PAS_PAR_BLOC:-10000}
REPLICATS=${REPLICATS:-20}
