# Chemin 1 : l'ARN non codant agit sur des gènes mitochondriaux (cibles ou épissage)

*Version 0, 19 novembre 2026. Première implémentation du chemin 1 de `../../chemins_mito.md`.*

Question posée pour chaque variant : le variant touche-t-il un ARN non codant nucléaire (miARN, snRNA du
spliceosome majeur ou mineur, snoRNA, lncRNA, ARN Pol III, ARNt, ARNr) dont les **cibles** ou les **sites
d'épissage** sont dans des gènes mitochondriaux, plus souvent que ne le voudrait le hasard ?

Sortie : une colonne `NC_chemin1` = `score|niveau|justification`, niveau dans
`fort` / `modere` / `faible` / `hypothese` / `aucun`. Si le variant touche plusieurs éléments, le meilleur
niveau est retenu et les justifications sont concaténées (séparateur ` ; `).

## État d'avancement

| Volet | État | Testé sur |
|---|---|---|
| Spliceosome mineur : symboles HGNC à jour, table gène mito → introns U12, liste U12 génome entier (intronIC) | **fait** | données réelles (GENCODE v27, hg38) |
| miARN : graine REF/ALT, sites 8mer / 7mer-m8 / 7mer-A1 dans les UTR 3', gains / pertes, hypergéométrique, ajustement longueur (info) | **fait** | réel + jeu construit |
| miARN hors graine : ΔΔG ViennaRNA du précurseur, distance aux sites Drosha / Dicer | **fait** | réel (MIR96) |
| Spliceosome majeur : ΔΔG d'appariement U1 (nt 1-11) / site 5', boîte ACAGAGA de U6 / intron +1..+12, U2 (GUAGUA) / 3' de l'intron | **fait (v0)** | réel (RNU1-1, RNU2-2) |
| Pol III / ARNt / ARNr : table des partenaires protéiques et de leur appartenance à l'ensemble mito | **fait** | réel + jeu construit |
| lncRNA : antisens chevauchant un gène mito (cis) | **fait** | jeu construit |
| lncRNA : boucles mitoloop, co-expression GTEx ; snoRNA : cibles snoDB | **squelette** (code écrit, données sur Curta) | jeu synthétique (co-expression) |
| ARNt nucléaires (GtRNAdb), guide antisens des snoRNA (type PLEXY / RNAsnoop) | **non fait** | — |

## Fichiers

```
chemin1/
  README.md
  commun.py                    lecture VCF/TSV, FASTA indexé, refFlat, HGNC, hypergéométrique, format de sortie
  curta/prep_chemin1.py        préparation des tables (machine connectée) : sous-commandes ci-dessous
  curta/telecharger_curta.sh   téléchargements + préparation du paquet, une commande par ligne
  vm/annoter_chemin1.py        annotation hors ligne (VM) : VCF ou TSV chr/pos/ref/alt -> colonne NC_chemin1
  tests/mini_paquet.py         JEU DE TEST CONSTRUIT À LA MAIN (gènes et séquences fictifs)
  tests/construire_temoins.py  témoins du catalogue convertis en chr/pos/ref/alt hg38 (contrôle de la base REF)
  tests/test_chemin1.sh        tests (jeu construit toujours ; témoins si PAQUET=<dossier>)
  tables/                      tables produites dans cette session (petites)
```

Tables produites (`tables/`) :

| Fichier | Contenu |
|---|---|
| `introns_U12_genes_mito_hg38_hgnc.tsv` | la table d'origine avec `gene_hgnc` et `statut_hgnc` (GARS → GARS1, C21orf33 → GATD3 ; les 12 autres inchangés) |
| `gene_mito_introns_U12.tsv` | table d'origine regroupée par gène : 14 gènes mito, 1 intron U12 chacun |
| `genes_U12_genome_entier.tsv` | intronIC sur tous les transcrits codants GENCODE v27 : 647 gènes à intron U12 |
| `gene_mito_introns_U12_intronIC.tsv` | les 19 gènes de l'ensemble mito parmi eux, avec coordonnées, dinucléotides et score |
| `resume_U12_genome_entier.tsv` | effectifs et test hypergéométrique |
| `partenaires_pol3_arnt_arnr.tsv` | partenaires protéiques (groupes HGNC + PNPT1), appartenance à l'ensemble mito, action sur l'ARN nucléaire |
| `temoins_annotes.tsv` | sortie de `annoter_chemin1.py` sur les 32 témoins |

Le paquet complet (UTR 3', sites d'épissage, éléments, miARN : 18 Mo) n'est pas dans le dépôt ; il se
régénère avec `curta/telecharger_curta.sh`.

## Méthode

### Ensemble mito et fond

Ensemble mito = MitoCarta3.0 ∪ PanelApp « Mitochondrial disorders », symboles mis à jour par HGNC
(symbole approuvé, sinon ancien symbole non ambigu, sinon alias non ambigu). **Dans cette session, seule
MitoCarta3.0 (1 136 gènes) a pu être chargée** ; PanelApp est inaccessible d'ici. Sur Curta, la liste
PanelApp existe déjà (`/scratch/atrimouille/MITODARK/data/panelapp_genes.txt`, avec sa provenance JSON) ;
le paramètre `--panelapp` de `prep_chemin1.py mito` l'ajoute, et `--mito` de `annoter_chemin1.py` permet
de passer n'importe quelle autre liste (un symbole par ligne, ou TSV avec une colonne `gene` / `Symbol`).

Le fond de chaque test est l'ensemble des gènes codants qui portent l'objet testé : gènes ayant une UTR 3'
(18 635) pour les miARN, gènes ayant au moins un intron codant (≈ 19 000) pour l'épissage. La part mito
du fond est donc ≈ 6 %.

### Spliceosome mineur (U11, U12, U4atac, U6atac)

1. Symboles HGNC de `introns_U12_genes_mito_hg38.tsv` mis à jour (sous-commande `u12`).
2. Introns U12 du génome entier : intronIC 3.0.0 (Moyer et al. 2020, méthode de l'IAOD) sur un GTF tiré du
   refFlat GENCODE v27 (transcrits codants), seuil par défaut de l'outil (probabilité U12 ≥ 90 %). Les
   introns identiques partagés entre transcrits sont fusionnés par intronIC ; on les rattache à tous les
   gènes de leurs transcrits (fichier `dupe_map`), faute de quoi MTG1 était perdu (son intron AT-AC est
   représenté par un transcrit de lecture continue).
3. Annotation : tout variant dans RNU11, RNU12, RNU4ATAC, RNU6ATAC renvoie la liste des gènes mito à intron
   U12 et le test d'enrichissement génome entier. Niveau `faible` seulement si cet enrichissement est
   significatif (P < 0,05), sinon `hypothese` : l'effet attendu est global, pas spécifique des gènes mito.

### miARN

- Précurseurs : transcrits GENCODE de type « RNA, micro » (HGNC), séquence lue dans hg38 avec 15 nt de flanc ;
  bras matures miRBase localisés par recherche exacte de leur séquence dans le précurseur.
- Sites (Bartel 2009 ; Agarwal et al. 2015, TargetScan 7) : 8mer, 7mer-m8 et 7mer-A1 dans l'UTR 3'. Un gène
  est cible s'il porte au moins un de ces sites. Les 6mer et les sites non canoniques sont ignorés.
- UTR 3' : pour chaque gène codant, transcrit codant à l'UTR 3' la plus longue (règle a priori, simple).
- Variant dans la graine (nt 2-8 d'un bras mature) : cibles REF et ALT, gènes perdus et gagnés, nombre de
  gènes mito parmi eux, P hypergéométrique (queue supérieure) contre le fond. On donne aussi, à titre
  d'information seulement, un nombre attendu de gènes mito stratifié par déciles de longueur d'UTR 3' et
  son z-score : les UTR des gènes mito sont plus courtes, ce qui rend le test simple conservateur.
- Variant hors graine (reste du mature, ou précurseur hors bras) : ΔΔG = MFE(ALT) − MFE(REF) du précurseur
  (ViennaRNA 2.8, RNAfold, paramètres Turner 2004 par défaut ; Lorenz et al. 2011) et distance aux
  extrémités des bras matures (sites Drosha et Dicer). Si l'un des deux seuils est franchi : perte de
  maturation possible, donc perte de toutes les cibles des deux bras, avec leur part mito et le test.

Niveaux (miARN, graine) :

| Niveau | Règle (fixée avant tout test) |
|---|---|
| fort | P < 0,001 et ≥ 5 gènes mito perdus (ou gagnés) |
| modere | P < 0,05 et ≥ 3 gènes mito |
| hypothese | au moins un gène mito perdu ou gagné, sans enrichissement |
| aucun | aucun gène mito touché, ou REF discordante |

Hors graine : `faible` si la perte de maturation toucherait des cibles enrichies en gènes mito (P < 0,05,
≥ 3 gènes), sinon `hypothese` ; `aucun` sous les deux seuils.

### Spliceosome majeur (U1, U2, U6 et paralogues)

Un variant n'est évalué que s'il tombe dans une région d'appariement aux sites, repérée par son motif dans
la séquence du paralogue (pas de numérotation codée en dur) :

| snRNA | Région | Comparée à |
|---|---|---|
| U1 | nt 1-11 (`AUACUUACCUG`) | site 5' de chaque intron codant, −3..+6 |
| U6 | boîte `ACAGAGA` ± 2 nt | intron +1..+12 |
| U2 | `GUAGUA` ± 3 nt | 42 premiers des 50 derniers nt de l'intron (zone du point de branchement) |

Énergie d'hybridation : `RNA.duplexfold` (ViennaRNA), ΔΔG = ALT − REF pour chaque site ; un gène est touché
si l'un de ses sites change de |ΔΔG| ≥ 1,5 kcal/mol. Puis hypergéométrique gènes mito contre fond. Niveau
`modere` si P < 0,001 et ≥ 5 gènes mito, `hypothese` si des gènes mito sont touchés sans enrichissement,
`aucun` sinon. U4, U5, U7 et les variants hors région d'appariement : `hypothese` (effet global possible,
par exemple ReNU pour RNU4-2, sans spécificité mito attendue).

### snoRNA (squelette)

Si `--snodb` est fourni (TSV `snorna`, `cible_gene`, produit par `prep_chemin1.py snodb`) : cibles connues
dans des gènes mito → `faible`. Sans table : `aucun`, « non évalué ». Restent à faire : position du variant
dans la boîte antisens (guide) et cibles prédites sur les ARNm mito (type PLEXY / RNAsnoop).

### lncRNA

- cis, fait : lncRNA antisens chevauchant un gène codant mito (calculé à la préparation) → `faible`.
- cis, squelette : ancres de boucles mitoloop reliées à un promoteur mito (`--boucles`, TSV `gene_mito`,
  `chrom`, `debut0`, `fin`, produit par `prep_chemin1.py boucles` à partir de
  `/scratch/atrimouille/MITODARK/loops/loops.tsv` ; les noms de colonnes de ce fichier sont à passer en
  option, ils n'ont pas pu être lus d'ici) → `hypothese`. Les TAD ne sont pas encore utilisés.
- trans, squelette : co-expression GTEx v8 (`--coexpression`, produit par `prep_chemin1.py coexpression`) :
  par tissu (muscle squelettique, ventricule gauche, cortex cérébral, cortex rénal), score du module
  OXPHOS = médiane des z-scores log2(TPM+1) des gènes « OXPHOS subunits » de MitoCarta3.0, puis r de
  Spearman avec chaque lncRNA (TPM médian ≥ 1, ≥ 30 échantillons). |r| ≥ 0,5 → `hypothese`.

### Pol III, ARNt, ARNr

Table `partenaires_pol3_arnt_arnr.tsv`, construite uniquement à partir de faits vérifiables :
appartenance à un groupe HGNC (RNase MRP complex subunits pour RMRP et RPPH1 ; Telomerase holoenzyme
subunits et H/ACA ribonucleoprotein complex pour TERC ; 7SK snRNP complex ; Signal recognition particle
pour 7SL ; tous les groupes HGNC dont le nom contient « tRNA » pour les ARNt), plus PNPT1, dont dépend
l'import mitochondrial de l'ARN de la RNase P, de RMRP et de l'ARNr 5S (Wang et al. 2010).
La colonne `agit_sur_arn_nucleaire` écarte les paralogues mitochondriaux (nom HGNC contenant
« mitochondrial » : ils agissent sur les ARNt mitochondriaux, pas sur les ARNt nucléaires) et marque
`a_verifier` les gènes de l'ensemble mito sans cette mention (GARS1, KARS1, YARS2, GATB, GATC, QRSL1) :
leur action sur les ARNt nucléaires n'est pas établie par cette table.

Résultat : aucune sous-unité de la RNase MRP / P, de la télomérase, de 7SK ou de la SRP n'est dans
MitoCarta3.0 ; le seul partenaire mito retenu pour RMRP, RPPH1 et l'ARNr 5S est PNPT1. Niveau
`hypothese` au plus (une interaction ne dit pas que le variant la modifie).

## Seuils a priori et justification

Aucun seuil n'a été ajusté sur des variants de patients ni sur des sorties d'AlphaGenome, SpliceAI,
SpliceAI2 ou AlphaFold 3. Les témoins n'ont servi qu'à mesurer, après coup ; aucun seuil n'a été changé
après leur passage.

| Seuil | Valeur | Justification |
|---|---|---|
| types de sites miARN | 8mer, 7mer-m8, 7mer-A1 | sites canoniques efficaces (Bartel 2009 ; Agarwal et al. 2015) ; 6mer et non canoniques exclus |
| P « modéré » | 0,05 | convention, un test par direction (pertes, gains) |
| P « fort » | 0,001 | plusieurs tests par variant (pertes, gains, bras, régions) et des centaines de variants de petits ARN par génome : seuil plus strict choisi a priori, non calibré |
| effectif minimal | 3 (modéré), 5 (fort) gènes mito | éviter qu'un ou deux gènes donnent seuls un P faible |
| ΔΔG précurseur | 2,0 kcal/mol | choix a priori, nettement au-dessus de l'arrondi des énergies de plus proches voisins ; **non calibré** |
| distance au clivage | 2 nt | couvrir un décalage d'1 à 2 nt des extrémités ; **hypothèse de travail, non sourcée** |
| ΔΔG appariement snRNA | 1,5 kcal/mol | choix a priori ; **non calibré** (voir points douteux) |
| intron U12 | probabilité ≥ 90 % | seuil par défaut d'intronIC (Moyer et al. 2020) |
| co-expression | \|r\| ≥ 0,5 | convention « corrélation forte », a priori |

Ces seuils sont dans `SEUILS` en tête de `vm/annoter_chemin1.py`. Les recalibrer exige des données
indépendantes (par exemple des mesures expérimentales publiées de maturation des miARN ou d'appariement),
jamais des sorties des modèles cités ni des variants de patients.

## Résultats sur les témoins (cette session)

Données : GENCODE v27 refFlat et hg38 (seau public Broad), HGNC du 19/11/2026, copie GitHub de
MitoCarta3.0, copie nf-core de miRBase (voir Sources). **MitoCarta seule, sans PanelApp.**

Spliceosome mineur :
- Table d'origine : 14 gènes, tous dans MitoCarta3.0 ; deux symboles mis à jour (GARS → GARS1,
  C21orf33 → GATD3).
- intronIC : 697 introns U12 (dont 169 AT-AC) dans 647 gènes. 19 sont des gènes mito : les 14 de la table
  d'origine plus ARL2, COA3, FHIT, HADH et HDHD5. Part mito : 19/647 = 2,9 % contre 5,8 % dans le fond
  (1 110 / 19 093 gènes à introns) ; attendu 37,6. Pas d'enrichissement (P = 1) ; appauvrissement,
  P(X ≤ 19) = 3,9·10⁻⁴. Un variant du spliceosome mineur n'a donc pas de raison de cibler préférentiellement
  les gènes mito ; il touche 19 gènes mito, dont HARS2, GARS1, PHB2, DNA2.

miARN (cibles REF, UTR 3' la plus longue, sans filtre de conservation) :
- miR-96-5p : 4 081 gènes cibles, dont 204 mito (5,0 % ; attendu ajusté sur la longueur 216, z = −1,0).
- miR-204-5p : 5 133 gènes, dont 259 mito (5,0 % ; attendu ajusté 271,7, z = −1,0).
- miR-184 : 994 gènes, dont 61 mito (6,1 % ; z = +1,0).
- MIR204 n.37C>T (graine, position 5 du mature, REF contrôlée) : 2 153 cibles perdues dont 116 mito
  (P = 0,91), 3 753 gagnées dont 211 mito (P = 0,88) → `hypothese`.
- MIR184 c.57C>T : perdues 587 dont 40 mito (P = 0,22) → `hypothese`.
- MIR96 : les 21 substitutions de la graine de miR-96-5p (jeu construit, les positions publiées ne sont
  pas dans le résumé de Mencía et al. 2009) donnent toutes `hypothese` ; le plus petit P vaut 0,081
  (gains, pos 2 A>C).

Spliceosome majeur et autres :
- RNU1-1 r.3A>G (somatique) : 11 769 gènes avec |ΔΔG| ≥ 1,5 au site 5', dont 730 mito (P = 0,33).
- RNU2-2 n.35A>G : 16 499 gènes touchés, dont 1 010 mito (P = 0,66) ; n.4G>A hors région → `hypothese`.
- RNU4-2 n.64_65insT, RNU12 n.84C>T, RNU4ATAC n.51G>A, RNU6ATAC n.43G>A, RMRP n.71A>G / n.72A>G :
  `hypothese`.

Bilan : les 32 témoins (tous des phénocopies ou des éléments sans lien mito établi) sortent au niveau
`hypothese` ; aucun n'atteint `faible` ou plus. C'est le résultat attendu pour la spécificité.

## Sources et provenance

| Donnée | Source utilisée ici | Empreinte md5 | À utiliser sur Curta |
|---|---|---|---|
| HGNC complete set | `storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt`, 19/11/2026 | 2a30aa650f51d6fed261f79057463d5a | idem |
| hg38 + .fai | `gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta` | — | idem |
| GENCODE v27 (refFlat) | `gcp-public-data--broad-references/hg38/v0/GRCh38_gencode.v27.refFlat.txt` | — | idem, ou GENCODE récent (GTF avec biotypes) |
| MitoCarta3.0 | **copie non officielle** : `raw.githubusercontent.com/hnolCol/instantclue/master/src/main/python/annotations/intern/Human.MitoCarta3.0.txt` (1 136 gènes, comme la publication) | 587bccd450d1ec6e736edefa44edd2e9 | fichier officiel du Broad |
| miRBase | **copie** : `raw.githubusercontent.com/nf-core/test-datasets/smrnaseq/reference/{mature,hairpin}.fa` (humain seulement, 2 588 matures, 1 881 épingles ; version non indiquée dans le fichier) | mature 46db10d71667b8cd40334fa0a969bbbd ; hairpin 466e1b446da10657c919db1bc2df6320 | miRBase v22.1 officiel |
| intronIC | PyPI `intronIC` 3.0.0, modèle pré-entraîné fourni, scikit-learn 1.7.2 | — | idem |
| ViennaRNA | PyPI `ViennaRNA` 2.8.0 | — | idem (à installer sur la VM) |

Inaccessibles depuis la session (refusés par le proxy) : ftp.ebi.ac.uk, mirbase.org, Broad (MitoCarta),
PanelApp, snoDB, introns.info (WtMTA), UCSC, Ensembl, NCBI.

## Ce qui reste à télécharger ou à produire sur Curta

1. MitoCarta3.0 officiel et la liste PanelApp déjà présente sur Curta → `genes_mito.tsv` (relancer tout ce qui
   en dépend : U12, partenaires, et rien d'autre, l'annotation lisant la liste à chaque fois).
2. miRBase v22.1 officiel (`mature.fa`, `hairpin.fa`) → `mirna.tsv.gz`.
3. ARNt nucléaires : GtRNAdb hg38 (GENCODE n'en contient pas) ; un convertisseur BED → `elements_ncrna` reste à écrire.
4. snoDB (export des interactions) → `snodb_cibles.tsv` ; formats de colonnes à relever.
5. Boucles mitoloop (`/scratch/atrimouille/MITODARK/loops/loops.tsv` ou `assoc_loops.tsv`) → `boucles_mito.tsv` ;
   TAD : à ajouter.
6. GTEx v8 TPM + attributs d'échantillons → `coexpression_lncrna_oxphos.tsv`.
7. Optionnel : GENCODE récent en GTF, pour classer les lncRNA par biotype (les noms de clones de GENCODE v27,
   absents de HGNC, sont aujourd'hui ignorés : 3 883 lncRNA seulement).

Commandes (sur Curta, depuis le dépôt) :

```
cd docs/ncrna_effecteur/chemins/chemin1
bash curta/telecharger_curta.sh /scratch/$USER/chemin1_donnees /scratch/$USER/chemin1_paquet
PAQUET=/scratch/$USER/chemin1_paquet bash tests/test_chemin1.sh
```

Sur la VM (paquet transféré, ViennaRNA installé) :

```
python3 vm/annoter_chemin1.py --variants variants.vcf.gz --paquet chemin1_paquet --out variants_chemin1.tsv
python3 vm/annoter_chemin1.py --variants variants.tsv --paquet chemin1_paquet --mito genes_mito.txt --snodb chemin1_paquet/snodb_cibles.tsv --boucles chemin1_paquet/boucles_mito.tsv --coexpression chemin1_paquet/coexpression_lncrna_oxphos.tsv --out variants_chemin1.tsv
bash tests/test_chemin1.sh
```

L'annotation charge les tables à la demande ; un variant dans U1, U2 ou U6 coûte environ 40 s (≈ 260 000
sites recalculés), un variant de graine de miARN quelques secondes (lecture des UTR 3' au premier appel).

## Limites

- Variants pris en compte : SNV et petits indels à l'intérieur de l'élément GENCODE (± 15 nt pour les
  miARN). Les grandes délétions (MIR17HG), les épimutations et les variants du pri-miARN hors flancs ne sont
  pas traités.
- Les miARN : pas de filtre de conservation, pas de contexte de site (TargetScan context++), UTR 3' unique
  par gène (la plus longue) ; les nombres de cibles sont donc plus élevés que ceux de TargetScan. Isoformes
  de miARN (isomiR) et changement de bras dominant ignorés.
- Spliceosome majeur : seuls les sites annotés sont recalculés, pas les sites cryptiques (alors que l'effet
  connu de U1 r.3A>G passe par des sites 5' cryptiques) ; le point de branchement n'est pas localisé (meilleur
  duplex sur 42 nt).
- Nombreuses copies : GENCODE v27 donne le nom RNU1-1 à 13 loci ; tous sont évalués comme U1.
- HGNC sert à classer les éléments (« RNA, micro », « RNA, small nuclear »…) : les gènes de GENCODE v27 sans
  symbole HGNC sont ignorés.
- Les niveaux traduisent une plausibilité mécanistique, pas une pathogénicité ; ils ne donnent aucun point
  ACMG.

## Points douteux

- MitoCarta et miRBase viennent de copies tierces (provenance et md5 ci-dessus) ; à remplacer par les fichiers
  officiels avant tout usage.
- Le seuil de 1,5 kcal/mol pour U2 est peu discriminant : avec la recherche du meilleur duplex sur 42 nt,
  RNU2-2 n.35A>G touche 16 499 gènes sur ≈ 19 000. Le constat vient d'un témoin ; ne pas recalibrer sur ces
  témoins, mais décider a priori d'une autre mesure (point de branchement localisé, ou seuil sur le duplex
  attendu U2/BPS seulement).
- Les numérotations n. du catalogue supposent que le transcrit GENCODE commence comme le NR_ ou l'épingle
  miRBase ; la base REF a été vérifiée pour les 10 témoins convertis (toutes concordantes), mais RMRP n.71 et
  n.72 portent tous deux un A : la conversion ne permet pas de trancher entre les deux numérotations.
- `GATB`, `GATC`, `QRSL1`, `YARS2`, `GARS1`, `KARS1` : action sur les ARNt nucléaires marquée `a_verifier`.
- La version de miRBase de la copie nf-core n'est pas indiquée (1 881 épingles humaines).

## Références (vérifiées dans PubMed)

- Agarwal V, Bell GW, Nam JW, Bartel DP. Predicting effective microRNA target sites in mammalian mRNAs.
  eLife 2015;4:e05005. doi:10.7554/eLife.05005
- Bartel DP. MicroRNAs: target recognition and regulatory functions. Cell 2009;136:215-233.
  doi:10.1016/j.cell.2009.01.002
- Moyer DC, Larue GE, Hershberger CE, Roy SW, Padgett RA. Comprehensive database and evolutionary dynamics of
  U12-type introns. Nucleic Acids Res 2020;48:7066-7078. doi:10.1093/nar/gkaa464
- Rath S et al. MitoCarta3.0: an updated mitochondrial proteome now with sub-organelle localization and
  pathway annotations. Nucleic Acids Res 2021;49:D1541-D1547. doi:10.1093/nar/gkaa1011
- Lorenz R et al. ViennaRNA Package 2.0. Algorithms Mol Biol 2011;6:26. doi:10.1186/1748-7188-6-26
- Wang G et al. PNPASE regulates RNA import into mitochondria. Cell 2010;142:456-467.
  doi:10.1016/j.cell.2010.06.035
- Wu J et al. Cryo-EM structure of the human ribonuclease P holoenzyme. Cell 2018;175:1393-1404.
  doi:10.1016/j.cell.2018.10.003 (composition protéique de la RNase P humaine)
- Mencía A et al. Mutations in the seed region of human miR-96 are responsible for nonsyndromic progressive
  hearing loss. Nat Genet 2009;41:609-613. doi:10.1038/ng.355
- Conte I et al. MiR-204 is responsible for inherited retinal dystrophy associated with ocular coloboma.
  Proc Natl Acad Sci USA 2015;112:E3236-E3245. doi:10.1073/pnas.1401464112
