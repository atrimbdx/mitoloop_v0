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

    cd transfert_mitoloop && bash 07_maj_mitoloop_vm.sh

Les données vont dans `$REANNOT_BASE/data/mitoloop`, le script dans
`$REANNOT_BASE/scripts`. Une mise à jour remplace la version précédente
(gardée dans `data/mitoloop.ancien`).

## Utilisation (VM)

    python3 ~/reannot/scripts/ajout_mitoloop.py run1.final.tsv [run2.final.tsv ...]

Résultat à côté de chaque fichier : `run1.final.mitoloop.tsv`. Seuls les
variants des régions régulatrices des gènes mitochondriaux (promoteur, UTR,
éléments ENCODE-rE2G, ancres de boucles) sont annotés ; les autres ont le
statut `hors_regions_mito`. Trier sur `MITO_n_arguments` à l'intérieur des
priorités habituelles.

## Colonnes

| Colonne | Contenu |
|---|---|
| MITO_gene | gène(s) mitochondrial(aux) et classe de région |
| MITO_L1_classe | transmitted / uncertain / not_transmitted / not_estimable |
| MITO_L1, MITO_L1_IC | L1 et son intervalle à 95 % |
| MITO_L1_note | sous-unité tenue par son complexe ; corrélation ARN-protéine dans les fibroblastes |
| MITO_dist_TSS | distance au TSS EPDnew le plus proche (négatif = en amont) |
| MITO_CAGE | sur un pic CAGE FANTOM5 |
| MITO_fixation | facteurs fixés (ReMap 2022) et nombre de types cellulaires |
| MITO_motifs | sites JASPAR touchés : position, information (bits), variation du score |
| MITO_chromatine | état ChromHMM (muscle E108, cerveau E073) |
| MITO_arguments, MITO_n_arguments | arguments réunis |

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
