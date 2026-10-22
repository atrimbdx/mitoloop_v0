# Module MitoLoop du pipeline reannot

Ajoute aux tables du pipeline (`*.final.tsv`, ou `*.final.spliceai2.tsv`) des
arguments sur l'effet possible d'un variant non codant sur l'expression des
gènes mitochondriaux. Ce n'est pas un score : chaque argument vient d'une
donnée publiée, et la liste dit pourquoi un variant remonte.

## Préparation (Curta, dans mitoloop_v0)

    bash run_mcia.sh all        # régulome, L1, annotations, puis le paquet
    # ou, après un run complet : bash run_mcia.sh reannot

Le paquet est écrit dans `mitoloop_v0/transfert_mitoloop/` : à copier sur la VM.

## Installation et mise à jour (VM)

    export REANNOT_BASE=/home/jovyan/DURAAK/reannot     # dossier du pipeline sur la VM
    cd transfert_mitoloop && bash 07_maj_mitoloop_vm.sh

Les données vont dans `$REANNOT_BASE/data/mitoloop`, les scripts
(`ajout_mitoloop.py`, `preparer_entree.py`, `reannot_complet.sh`) dans
`$REANNOT_BASE/scripts`. Une mise à jour remplace la version précédente
(gardée dans `data/mitoloop.ancien`).

## Utilisation (VM)

Tout le pipeline en une commande, depuis un Excel (.xlsx), un TSV, un CSV ou
un VCF :

    export REANNOT_BASE=/home/jovyan/DURAAK/reannot
    bash $REANNOT_BASE/scripts/reannot_complet.sh variants.xlsx run1 [--feuille NOM]

Étapes : préparation de l'entrée (colonnes transcrit et variant c./n./g.
reconnues par leur contenu, ou chrom/pos/ref/alt ; sans transcrit, le MANE
Select voisin), VEP + SpliceAI + MORFEE + AlphaGenome
(`reannot_alphagenome.sh`), SpliceAI2 s'il est installé, MitoLoop. Résultat :
`run1.final.spliceai2.mitoloop.tsv` (ou `run1.final.mitoloop.tsv`).

Notation c./n. sans transcrit, avec une colonne chromosome (1..22, X, Y, MT) et
une colonne position (GRCh38) : le transcrit pris est le MANE Select à cette
position. Comme la notation a pu être écrite sur un autre transcrit, la
position obtenue après conversion est comparée à la position saisie : colonne
`controle_position` (concordante, ou discordante avec la position saisie).
Une discordance veut dire que la notation porte sur un autre transcrit (le
donner dans une colonne ENST) ou que les positions ne sont pas en GRCh38.

MitoLoop seul, sur une table déjà produite ou sur un VCF :

    python3 $REANNOT_BASE/scripts/ajout_mitoloop.py run1.final.tsv
    python3 $REANNOT_BASE/scripts/ajout_mitoloop.py variants.vcf.gz

Seuls les variants des régions régulatrices des gènes mitochondriaux
(promoteur, UTR, éléments ENCODE-rE2G, ancres de boucles) sont annotés ; les
autres ont le statut `hors_regions_mito`. Trier sur `MITO_n_arguments` à
l'intérieur des priorités habituelles.

## Colonnes

| Colonne | Contenu |
|---|---|
| MITO_gene | gène(s) mitochondrial(aux) et classe de région |
| MITO_bidirectionnel | partenaire du promoteur bidirectionnel (autre gène, TSS sur le brin opposé à moins de 1 kb : distance, type) |
| MITO_contact_promoteur | ancre de boucle du gène qui porte le promoteur d'un autre gène (promoteur à fonction d'enhancer possible) |
| MITO_type_promoteur | type du TSS le plus proche : menage, tissu, developpement, indetermine, avec ses critères (classe CpG, TATA, Inr, tau GTEx, largeur du pic CAGE, Polycomb) |
| MITO_L1_classe | transmitted / uncertain / not_transmitted / not_estimable |
| MITO_L1, MITO_L1_IC | L1 et son intervalle à 95 % |
| MITO_L1_note | sous-unité tenue par son complexe ; corrélation ARN-protéine dans les fibroblastes |
| MITO_dist_TSS | distance au TSS EPDnew le plus proche (négatif = en amont) |
| MITO_CAGE | sur un pic CAGE FANTOM5 |
| MITO_fixation | facteurs fixés (ReMap 2022) et nombre de types cellulaires |
| MITO_motifs | sites JASPAR touchés : position, information (bits), variation du score |
| MITO_CpG | CpG détruit ou créé ; o/e CpG, GC et classe (HCP/ICP/LCP) de la fenêtre de 500 pb, référence > variant ; sites touchés de facteurs sensibles à la méthylation (NRF1, ETS, CREB, E2F...) |
| MITO_chromatine | état ChromHMM (muscle E108, cerveau E073) |
| MITO_phyloP | conservation de la base sur 241 mammifères (Zoonomia) ; argument au-dessus de 2,27 |
| MITO_gnomAD | gène : LoF observés/attendus, o/e LoF, LOEUF, pLI, o/e et z faux-sens (gnomAD v4.1, MANE Select) |
| MITO_GISMO | gène : GISMO et GISMO-mis (Liao et al. 2026) et leur décile ; déciles 1-4 de GISMO pour les gènes récessifs |
| MITO_arguments, MITO_n_arguments | arguments réunis (les contraintes du gène, gnomAD et GISMO, n'en font pas partie) |
| MITO_evenements | événements portés par l'allèle, tous outils réunis, avec leur force : SpliceAI/SpliceAI2, MORFEE, AlphaGenome, site_TF, CpG_site_TF, bascule_ilot_CpG, TATA_detruite, Inr_detruit, G4_detruit |
| MITO_priorite | 1_fort, 2_modere, 3_contexte, 4_aucun |

## Priorité

Les annotations du gène et de la position sont les mêmes pour un variant
pathogène et un variant bénin au même endroit ; seules celles de l'allèle les
séparent. `MITO_evenements` réunit donc les événements de l'allèle, de tous les
outils de la table, et `MITO_priorite` les classe :

| Priorité | Condition |
|---|---|
| 1_fort | au moins un événement fort, gène dont la classe L1 n'est pas not_transmitted |
| 2_modere | un événement modéré, ou fort dans un gène not_transmitted |
| 3_contexte | pas d'événement de l'allèle, au moins deux arguments de position (promoteur minimal, zone de pause, CAGE, fixation, chromatine active, conservation) |
| 4_aucun | le reste |

Fort : boîte TATA détruite ; SpliceAI ≥ 0,2 ; SpliceAI2 ≥ 0,1 ; |AlphaGenome| ≥ 20 ; un champ
MORFEE du 5'UTR renseigné (uTIS, uSTOP, uSTOP_creation : création ou perte
d'une uORF, qui freine la traduction) ; site_TF d'au moins 2 (log2) sur une
position d'au moins 1,5 bits ; un événement modéré du promoteur minimal, ou
d'un élément distal (rE2G, ancre) d'un gène dont le promoteur est de type tissu
ou développement (les promoteurs de ménage dépendent peu des éléments distaux).
Modéré : un champ MORFEE du 3'UTR renseigné (dTIS, dSTOP, dSTOP_creation) ;
site_TF (ESRRA/ESRRG notés « ERR:effet_tissulaire ») ; CpG_site_TF ;
bascule_ilot_CpG ; initiateur détruit ; G-quadruplex cassé. `zone_pause` (+20
à +60 du TSS) compte parmi les arguments de position. Les seuils sont des options de
`ajout_mitoloop.py` (`--spliceai-seuil`, `--spliceai2-seuil`,
`--alphagenome-seuil`...). Ceux de SpliceAI, SpliceAI2 et AlphaGenome sont
fixés à la main : leurs sorties ne servent jamais à entraîner ni à ajuster un
paramètre, comme l'exigent leurs conditions d'utilisation ; la priorisation
qui les utilise se fait sur la VM. Les variants hors des régions MitoLoop
reçoivent aussi les événements d'épissage, MORFEE et AlphaGenome. Trier sur `MITO_priorite`, puis sur `MITO_n_arguments`.

## Sens et limites de L1

L1 dit si une baisse de l'ARNm d'un gène se retrouve sur sa protéine
(protéome GTEx, Jiang et al. 2020), calculé ici sans l'ajustement du contenu
mitochondrial : un événement chez un patient est une baisse absolue de la
protéine. Classe `transmitted` : l'intervalle à 95 % est au-dessus de 0,5 ;
`not_transmitted` : au-dessous de 0,3. La validation sur les fibroblastes de
patients (Kopajtich et al.) est dans `results/l1_raw/validation/validation.tsv`
(ligne `outlier_propagation_by_l1_class`). `not_estimable` ne veut pas dire
non transmis (POLG, TFAM, OPA1, SURF1 sont dans ce cas) : leurs variants
restent annotés. Une annotation oriente, elle ne permet pas d'exclure.
