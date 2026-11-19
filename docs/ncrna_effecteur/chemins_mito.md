# Trois chemins pour relier un variant non codant « effecteur » à la mitochondrie

*Version 0, 19 novembre 2026. Méthode validée oralement par l'utilisateur ; première implémentation en cours.*

## Objectif

Interpréter un variant rare de patient (WGS, VM hors ligne) situé dans un élément non codant effecteur ou un
micropeptide, pour **découvrir** un mécanisme mitochondrial nouveau. Le catalogue `catalogue_ncrna_mito.tsv`
sert de témoins positifs (surtout des phénocopies). Hors périmètre : ADNmt (déjà exploré chez les patients),
éléments régulateurs (promoteurs, enhancers, UTR, introns profonds, expansions : pris en charge par mitoloop),
ARN-seq et tests OXPHOS (non traités ici).

## Les 12 types d'éléments

1. snRNA du spliceosome majeur (U1, U2, U4, U5, U6)
2. snRNA du spliceosome mineur (U11, U12, U4atac, U6atac)
3. snRNA U7
4. snoRNA C/D, H/ACA, scaRNA
5. miARN
6. lncRNA (cas général)
7. lncRNA soumis à empreinte
8. lncRNA porteurs d'expansions
9. ARN de RNP transcrits par la Pol III (RMRP, RPPH1, 7SK, TERC)
10. ARNt nucléaires
11. ARNr 5S et 45S, ARN de voûte, ARN Y, 7SL
12. Micropeptides (sORF traduits dans des « lncRNA », UTR, uORF ; nucléaires uniquement)

## Principes communs

- Annotations préparées sur une machine connectée (Curta) à partir de données publiques, transférées en paquet ;
  calcul pour le variant sur la VM, hors ligne.
- Comparaison allèle de référence contre allèle alternatif chaque fois que c'est possible.
- Ensemble « mito » de référence : MitoCarta3.0 ∪ PanelApp « Mitochondrial disorders » (≈ 1 220 gènes, ≈ 6 %
  des gènes) ; tout score est comparé à ce bruit de fond.
- Seuils fixés a priori et écrits dans le code. Aucun entraînement ni ajustement de seuil sur des sorties
  d'AlphaGenome, SpliceAI, SpliceAI2 ou AlphaFold 3, ni sur des variants de patients. AlphaFold 3 n'est pas utilisé.
- Sortie : une colonne par chemin (`NC_chemin1`, `NC_chemin2`, `NC_chemin3`) avec score, niveau
  (`fort` / `modere` / `faible` / `hypothese` / `aucun`) et justification en clair.
- Contrôles : témoins phénocopies du catalogue (spécificité) et variants rares tirés de gnomAD (bruit de fond,
  mesure seulement).

## Chemin 1 — l'élément agit sur des gènes mitochondriaux (cibles ou épissage)

- miARN : sites cibles REF contre ALT dans les UTR 3' (graine 7-8mer, type TargetScan), gains et pertes dans
  les gènes mito, test hypergéométrique ; hors graine, ΔΔG ViennaRNA et position par rapport aux sites
  Drosha et Dicer (perte de maturation = perte de toutes les cibles, part de cibles mito).
- Spliceosome majeur : si le variant touche une région d'appariement (U1 5' / site 5', boîte ACAGA de U6,
  U2 / point de branchement), recalcul thermodynamique de l'appariement avec chaque site des gènes mito.
- Spliceosome mineur : gènes mito à intron U12 (`introns_U12_genes_mito_hg38.tsv`, symboles HGNC à jour).
- snoRNA : guide antisens touché ? cibles connues (snoDB) ou prédites (type PLEXY / RNAsnoop) sur les ARNm mito.
- lncRNA : cis (antisens chevauchant, même TAD ou boucle via mitoloop) ; trans (co-expression GTEx avec les
  modules OXPHOS, seuil fixé a priori).
- Pol III, ARNt, ARNr : table manuelle des partenaires protéiques à fonction mitochondriale.

## Chemin 2 — l'ARN est présent dans la mitochondrie

- Table de localisation mitochondriale des ARN non codants nucléaires, avec source et niveau de preuve
  (APEX-seq matrice > mitochondries purifiées + RNase > imagerie > fractionnement seul).
- Variant dans un ARN de cette table : touche-t-il l'élément d'import (tige-boucle PNPT1 de RMRP / RPPH1), un
  domaine de liaison, ou le repliement (ΔΔG ViennaRNA, position Rfam) ?
- Recherche exploratoire du motif d'import dans les autres ARN non codants : hypothèse uniquement, jamais de
  points ACMG.

## Chemin 3 — l'élément porte un micropeptide mitochondrial

- sORF traduits : catalogue GENCODE des Ribo-seq ORFs, OpenProt ; conservation (PhyloCSF).
- Lien mito du peptide : MitoCarta3.0, protéomes mitochondriaux publiés, séquence d'adressage (TargetP,
  MitoFates), hélice transmembranaire (DeepTMHMM).
- Variant traité comme codant via un GTF personnalisé dans VEP : perte d'initiation, stop, décalage, faux-sens
  conservé.

## Implémentation

Dossier `chemins/` : un sous-dossier par chemin (`chemin1/`, `chemin2/`, `chemin3/`), chacun avec un README,
les scripts Curta (préparation) et VM (annotation), les tables produites quand les données étaient
accessibles, et un test sur les témoins.
