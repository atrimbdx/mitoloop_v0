# Cibles des éléments non codants effecteurs : sources, conversion, résumé

*Version 0, 20 novembre 2026. Interface machine : `cibles_ncrna.tsv` (10 colonnes fixées avec reannot + colonne facultative finale `region_guide`).*

`region_guide` : région antisens du snoRNA/scaRNA au format `guide:A-B` (positions 1-based dans l'ARN mature du snoRNA), lue par reannot. Elle vaut `.` partout : aucune source vérifiée accessible (snoDB ne la fournit pas) ; à remplir depuis snoRNABase (appariements guide-cible) sur Curta. `position_cible` reste la position dans la cible.

## Règle d'usage par reannot

Un variant touche un élément (`type_id` + `element` ou `element_id`, ou `element = *` pour toute la classe).
On lit toutes ses cibles ; si une `cible` appartient à l'ensemble mito (MitoCarta3.0 ∪ PanelApp), on annote.
Cibles en cascade : la table donne la cible **directe** (ex. `RNU12`, `RNU4ATAC`) ; reannot joint ensuite
`chemins/chemin1/tables/gene_mito_introns_U12.tsv` (19 gènes mito à intron U12) pour U11/U12/U4atac/U6atac.
La table n'est **pas** filtrée sur la mitochondrie.

Conventions :

- ARNr : `18S` → `RNA18SN1`, `28S` → `RNA28SN1`, `5.8S` → `RNA5-8SN1`, `5S` → `RNA5S1` (copie 1 choisie
  par convention ; les 5 unités 45S sont identiques). snRNA : `U1` → `RNU1-1`, `U2` → `RNU2-1`, `U4` → `RNU4-1`,
  `U5` → `RNU5A-1`, `U6` → `RNU6-1`, `U7` → `RNU7-1`, `U11` → `RNU11`, `U12` → `RNU12`, `U4atac` → `RNU4ATAC`,
  `U6atac` → `RNU6ATAC`, `U3` → `SNORD3A`. Le jeton brut de snoDB est gardé dans `source`.
- `type_cible` hors liste n'existe pas : les cibles lncRNA et pseudogènes des interactions snoDB sont mises en
  `arnm` (transcrit Pol II), les petits ARN non codants (7SK, 7SL, ARN Y, RPPH1, scaRNA) en `snrna`, `MT-T*`
  en `arnt`. Les interactions snoRNA-snoRNA et snoRNA-miARN de snoDB (602 jetons) sont **ignorées**.
- `preuve = base_de_donnees` pour tout ce qui vient de snoDB, HGNC ou GENCODE ; `prediction` pour les voisins
  cis calculés ; `experimentale_humaine` / `experimentale_animale` pour les articles lus (PubMed).
- Les ENSG manquants sont complétés depuis HGNC à la fusion ; les anciens symboles de cible sont renommés
  (ex. ATP5S → DMAC2L), avec mention dans `source`.

## Lignes par classe (cibles_ncrna.tsv, 22 994 lignes)

| type_id | lignes | source |
|---|---|---|
| snorna_cd | 671 | snoDB v2.0 + Elliott 2019 / 2025 |
| snorna_haca | 492 | snoDB v2.0 (dont AluACA) |
| scarna | 98 | snoDB v2.0 |
| snorna_orphelin_15q11_14q32 | 201 | snoDB v2.0 (interactions haut débit seulement, pas de site de modification) |
| snorna_u3_u8_u13 | 104 | snoDB v2.0 |
| sdrna | 1 | Ender 2008 (ACA45 = SCARNA15 → CDK19) |
| mirna_canonique | 46 | catalogue (articles vérifiés) ; miRTarBase 11.0 ajouté sur Curta par reannot `nc_sources` (NC_MIRTARBASE ; « Functional MTI » seulement par défaut, paires dédoublonnées, preuve experimentale_humaine si au moins un Functional MTI non Weak) |
| svrna | 1 | Persson 2009 (svRNAb → CYP3A4) |
| vault_rna | 5 | Horos 2019, Kong 2000, Kickhoefer 2002 |
| vault_nc886 | 1 | Lee 2011 (PKR) |
| y_rna | 1 | HGNC groupe Ro60 |
| rn7sk | 5 | HGNC + Nguyen 2001, Yik 2003 |
| rmrp_rpph1, terc, srp_7sl, trna_nucleaire, rrna_5s | 24, 6, 24, 61, 1 | table chemin1 `partenaires_pol3_arnt_arnr.tsv` (HGNC + Wang 2010) |
| rdna_pRNA_papas | 2 | Feng 2023 (RIEP → C1QBP, CHCHD2) |
| lncrna_mito | 7 | catalogue (GAS5, SAMMSON, TUG1) |
| lncrna_general | 20 685 | GENCODE v27 (cis calculé) + 3 lignes d'articles |
| pseudogene_eponge | 544 | HGNC (nom « X pseudogene N ») × GENCODE v27 (biotypes transcrits) |
| micropeptide_lncrna_sorf | 2 | Ge 2021 (LINC00467/ASAP) ; le reste est dans chemin3 |
| micropeptide_numt_humanin | 13 | HGNC (MTRNR2L1-13 → MT-RNR2) |

Classes **sans** ligne (source à télécharger sur Curta) : mirtron, agotron, pirna_cluster, tdr_trf_tirna,
trna_lookalike_mt, ysrna, circrna, sno_lncrna_spa. Classes sans cible au sens de cette table : snRNA (les cibles
sont les introns, gérés par chemin1), ARNr 45S, ribozymes, TERRA, satellites, eRNA/PROMPT.

## Sources par classe

### snoRNA C/D, H/ACA, scaRNA, orphelins, U3/U8/U13 — snoDB v2.0 (utilisée)

- Référence : Bergeron 2023 NAR, doi:10.1093/nar/gkac835. Site : https://bioinfo-scottgroup.med.usherbrooke.ca/snoDB/
  (inaccessible depuis la session). Licence : non indiquée sur la copie ; usage académique, citer l'article.
- **Copie utilisée** (provenance) : dépôt GitHub `Thomson-RNA-Lab/Integrated-rRNA-snoRNA-seq` (licence MIT du
  dépôt ; Chan, Barbaric, Thomson 2026 bioRxiv doi:10.64898/2026.08.10.743918), commit
  5660d1faa5924afa80c778265f33dd988712ce8f, fichier `fig/snoRNA/snoDB_All_V2.0.tsv`, 2 123 entrées,
  md5 f8dc6b0648296a396bfca52f6f4cf677. Colonnes identiques à l'export « All » de snoDB 2.0.
- Contenu : colonnes `rrna_targets`, `snrna_targets` (sites guidés, ex. `U12-19`), et cibles d'interactions
  ARN-ARN haut débit sans position (`protein_coding_targets`, `lncrna_targets`…).
- Sur Curta, revalider avec l'export officiel puis relancer la conversion :

```
wget -O snoDB_All_V2.0.tsv https://raw.githubusercontent.com/Thomson-RNA-Lab/Integrated-rRNA-snoRNA-seq/5660d1faa5924afa80c778265f33dd988712ce8f/fig/snoRNA/snoDB_All_V2.0.tsv
md5sum snoDB_All_V2.0.tsv
python3 convertir_cibles.py snodb snoDB_All_V2.0.tsv > t_snodb.tsv
```

  (export officiel : page « Download » de snoDB, fichier « All snoRNAs » en TSV ; l'URL exacte n'a pas pu être
  vérifiée.)
- Compléments à télécharger sur Curta :
  - snoRNA-LBME-db / snoRNABase (Lestrade & Weber 2006 NAR, doi:10.1093/nar/gkj002), http://www-snorna.biotoul.fr/ :
    type de modification (Nm ou Ψ) de chaque site snRNA, pour remplacer l'heuristique scaRNA (voir points douteux).
    Pas de licence explicite (base académique). Export manuel des pages « Find guide RNA » puis
    `python3 convertir_cibles.py generique export.tsv --type-id scarna --col-element guide --col-cible cible --col-position position --type-cible snrna --nature pseudouridylation --preuve base_de_donnees --source "snoRNABase" --reference "snoRNABase (date d'export)"`.
  - RMBase v3.0 (Xuan 2024 NAR, doi:10.1093/nar/gkad1070), http://bioinformaticsscience.cn/rmbase/ : module
    « Mechanisms » (23 382 sites guidés par snoRNA ou enzyme), dont sites Nm/Ψ dans des ARNm. Licence : usage
    académique. Télécharger le tableau « snoRNA-guided » humain hg38 depuis la page Download, puis
    `python3 convertir_cibles.py generique rmbase_sno_hg38.tsv --type-id snorna_cd --col-element snoRNA --col-cible geneName --col-position modStart --type-cible arnm --nature modification_2OMe --preuve base_de_donnees --source "RMBase v3.0 Mechanisms" --reference "10.1093/nar/gkad1070"`
    (noms de colonnes à relever sur le fichier réel ; un appel par type de modification).
  - Prédiction sur les ARNm et pré-ARNm mito : snoGloBe (Deschamps-Francoeur 2022 NAR, doi:10.1093/nar/gkac475,
    GPL-3, https://github.com/scottgroup/snoGloBe ; C/D), PLEXY (Kehr 2011 Bioinformatics, doi:10.1093/bioinformatics/btq642,
    http://www.bioinf.uni-leipzig.de/Software/PLEXY ; C/D) et RNAsnoop (Tafer 2010 Bioinformatics,
    doi:10.1093/bioinformatics/btp680 ; H/ACA, distribué avec ViennaRNA). Sortie → `generique` avec
    `--preuve prediction`.

```
git clone https://github.com/scottgroup/snoGloBe.git
conda install -c bioconda viennarna
```

### miARN canoniques, mirtrons, agotrons

- Utilisé : 46 lignes du catalogue `catalogue_ncrna_mito.tsv` (articles déjà vérifiés). Le calcul par graine
  (TargetScan-like, REF/ALT) reste celui de chemin1.
- À télécharger : miRTarBase 2025 (Cui 2025 NAR, doi:10.1093/nar/gkae1072), https://mirtarbase.cuhk.edu.cn/
  (inaccessible ici). Licence : usage académique libre, citer. Fichier humain des MTI (nom exact à relever sur
  la page Download ; colonnes attendues : `miRTarBase ID, miRNA, Species (miRNA), Target Gene, …, Support Type,
  References (PMID)` ; convertir en CSV si fourni en xlsx).

```
python3 convertir_cibles.py mirtarbase hsa_MTI.csv --hairpin hairpin.fa --hgnc hgnc_complete_set.txt --version "miRTarBase 2025" > t_mtb.tsv
```

  `Functional MTI` → `experimentale_humaine` ; `Functional MTI (Weak)` et `Non-Functional` → `base_de_donnees`.
  Le miARN mature est rattaché à ses précurseurs miRBase puis au symbole HGNC (`hsa-mir-181a-1` → `MIR181A1`).
- Mirtrons : mêmes cibles que les miARN (graine) ; agotrons : liste de Hansen 2016 (doi:10.1038/ncomms11538,
  tableaux supplémentaires) puis graine, comme chemin1.

### sdRNA, svRNA, YsRNA, tRF/tiRNA, piARN

- sdRNA : 1 ligne (SCARNA15/ACA45 → CDK19, Ender 2008). Pas de base dédiée : liste à compléter depuis la
  littérature.
- svRNA : 1 ligne (`*` → CYP3A4, Persson 2009) : l'élément exact (vtRNA d'origine de svRNAb) n'a pas pu être
  confirmé dans le résumé, d'où `*`.
- YsRNA : aucune cible publiée retrouvée (inventaire : « aucun lien mito »).
- tRF/tiRNA : tRFtarget (cibles prédites et validées) et MINTbase v2 (https://cm.jefferson.edu/MINTbase/,
  licence académique) sur Curta, puis `generique --type-id tdr_trf_tirna --type-cible arnm --nature repression_graine`.
  L'extrait synthétique `tests/generique_synthetique.tsv` montre le format.
- piARN : piRBase v3.0 (Wang 2022 NAR, doi:10.1093/nar/gkab1012), http://bigdata.ibp.ac.cn/piRBase (« piRNA
  targets »). Licence : accès libre, citer. Conversion `generique --type-id pirna_cluster --nature repression_graine`.

### lncRNA

- Cis : calcul GENCODE v27 (`lncrna_cis`) : antisens chevauchant, même brin chevauchant, bidirectionnel
  (TSS opposés à ≤ 1 kb), et gène codant le plus proche de chaque côté dans une fenêtre de 50 kb ;
  `preuve = prediction`. Les boucles mitoloop (TAD, contacts) restent à ajouter par reannot.
- Trans : partenaires d'articles (GAS5, SAMMSON, TUG1 ; `lncrna_mito`). La co-expression GTEx (chemin1)
  n'est pas une cible au sens de cette table.

### Pseudogènes éponges

- Gène parent selon la nomenclature HGNC (nom approuvé « <parent> pseudogene N »), restreint aux biotypes
  GENCODE v27 transcrits de l'inventaire. Les pseudogènes nommés autrement (ex. « … like ») ne sont pas pris.
  Complément sur Curta : psiCube / pseudogene.org (parent par alignement), même conversion `generique
  --type-cible gene_parent --nature origine`.

### Sosies nucléaires d'ARNt mitochondriaux (NUMT)

- Telonis 2014 Front Genet (doi:10.3389/fgene.2014.00344) : liste des tRNA-lookalikes avec l'ARNt d'origine
  (tableaux supplémentaires, coordonnées hg19 à convertir en hg38). Conversion : `generique --type-id
  trna_lookalike_mt --type-cible arnt_mito_origine --nature origine --preuve base_de_donnees`.

### Micropeptides

- Déjà traités par chemin3 (`micropeptides_mito_etablis.tsv`). Ici : LINC00467/ASAP (partenaires ATP5F1A,
  ATP5F1C) et la famille humanine MTRNR2L1-13 (origine MT-RNR2, nom HGNC).

### ARN de voûte, nc886, RMRP, RPPH1, TERC, 7SK, 7SL, ARN Y

- Table chemin1 (groupes HGNC + PNPT1, Wang 2010 Cell doi:10.1016/j.cell.2010.06.035) convertie par `pol3`.
- Ajouts d'articles vérifiés : voûte (MVP, TEP1, PARP4 : Kong 2000 ; SSB : Kickhoefer 2002 ; SQSTM1 : Horos 2019),
  nc886 → EIF2AK2 (Lee 2011), 7SK → CDK9, CCNT1 (Nguyen 2001), HEXIM1 (Yik 2003), ARN Y → RO60 (groupe HGNC).

## Reproduire la table

```
python3 convertir_cibles.py snodb snoDB_All_V2.0.tsv > t_snodb.tsv
python3 convertir_cibles.py curees cibles_curees.tsv > t_curees.tsv
python3 convertir_cibles.py pol3 ../chemins/chemin1/tables/partenaires_pol3_arnt_arnr.tsv > t_pol3.tsv
python3 convertir_cibles.py humanin --hgnc hgnc_complete_set.txt > t_humanin.tsv
python3 convertir_cibles.py pseudogenes --hgnc hgnc_complete_set.txt --gtf gencode.v27.primary_assembly.annotation.gtf > t_pseudo.tsv
python3 convertir_cibles.py lncrna_cis --gtf gencode.v27.primary_assembly.annotation.gtf > t_lnc.tsv
python3 convertir_cibles.py fusion --hgnc hgnc_complete_set.txt --out cibles_ncrna.tsv t_snodb.tsv t_curees.tsv t_pol3.tsv t_humanin.tsv t_pseudo.tsv t_lnc.tsv
sh tests/test_convertir.sh hairpin.fa hgnc_complete_set.txt
```

`cibles_curees.tsv` contient les 74 lignes tirées d'articles (format de sortie sans `element_id`).
`tests/` contient des extraits **synthétiques étiquetés** (miRTarBase, format générique) : ils ne doivent
jamais être fusionnés dans la table réelle.

## Résumé : cibles mitochondriales (jointure MitoCarta3.0 seule ; PanelApp non disponible ici)

### Directes

- **SNORD32A → COX4I2** : confirmé, mais seulement dans un préprint (Elliott 2025 bioRxiv,
  doi:10.1101/2025.07.23.666475) : la délétion de U32A dans des HEK293T réduit la 2'-O-méthylation de l'ARNm de
  COX4I2 (RTL-P) et augmente la protéine sans changer l'ARNm. Position non publiée dans le résumé.
- snoDB, interactions haut débit sans position (`liaison`, à prendre comme faibles) : SNORD13 → PRODH2,
  SNORD118 → ABCB7, GARS1, MT-CO2, MT-CYB, MT-RNR1 ; SNORA26 → SLC25A3 ; SNORA80E → ATP5F1B, MT-RNR1 ;
  SNORA59B → MT-RNR1 ; ENSG00000252050 (H/ACA) → DMAC2L ; SNORD35A → FKBP10 ; SNORD10 → MT-CO1 ;
  SNORD54 → MT-ND5 ; SNORD3D → MT-CO3, MT-ND2 ; SNORD12B → MT-TF. Les cibles codées par l'ADNmt sont
  probablement des chimères artéfactuelles (ARN mitochondriaux très abondants).
- miARN (articles) : MIR181C → MT-CO1 ; MIR21 → MT-CYB ; MIR378A → MT-ATP6, CRAT ; MIR762 → MT-ND2 ;
  MIR4485 → MT-RNR2 ; MIR210 → ISCU, COX10 ; MIR338 → COX4I1, ATP5MC1 ; MIR33A/B → HADHB, CPT1A ;
  MIR484 → FIS1 ; MIR494 → TFAM ; MIR761 → MFF ; MIR9-1/2/3 → TRMU, GTPBP3, MTO1 ; MIR181A1/A2/B1/B2 →
  COX11, COQ10B, PRDX3.
- Partenaires : RMRP, RPPH1, ARNr 5S → PNPT1 ; GAS5 → FH, MDH2, CS ; SAMMSON → C1QBP ; PAPAS/RIEP → C1QBP,
  CHCHD2 ; LINC00467 → ATP5F1A, ATP5F1C ; ARNt → aaRS mitochondriales (table chemin1).
- Origine : MTRNR2L1-13 → MT-RNR2 ; 43 pseudogènes transcrits à parent mito (dont SDHA, CYCS, VDAC1,
  HSPD1, HSPE1, COX5B, NDUFB4, NDUFAF1, SLC25A3, GLUD1, GOT2, MT-ND1, MT-ND6).
- lncRNA cis : 982 lncRNA ont au moins un voisin cis mito (629 gènes mito distincts) : prédiction pure.

### En cascade par le spliceosome mineur (cible directe = snRNA)

| snRNA | Guides à site de modification (snoDB) | Interactions sans position |
|---|---|---|
| RNU12 (U12) | SCARNA17 (21, 22 ; C/D), SCARNA20 (27/28 ; H/ACA), SCARNA21 (17-19 ; H/ACA), SCARNA21B (18), SNORD140 (83), SNORD162 (81), 6 SNORD-like non nommés | SNORD13, SNORD20, SNORD22, SNORD25, SNORD27, SNORD30, SNORD41, SNORD48, SNORD49A, SNORD51, SNORD53, SNORD74L3, SNORD104, SNORA75 |
| RNU4ATAC (U4atac) | SNORA7A/7B (14), SNORA47 (5), SNORA58/58B (9), SNORA62 (28), SNORA70F (19), SNORA71E (19), SNORA77 (41), SNORA3B (43), SNORA79B (90), SNORA92 (90), SNORA113, SNORA115, SNORA116, SNORD109A/B, SNORD124 (7), SNORD28B (102), 15 SNORD-like | SNORD10, SNORD30, SNORD44L1 |
| RNU6ATAC (U6atac) | SCARNA21 (83), SNORA119 (101), 3 SNORD-like | SNORA51L7, SNORD63, SNORD83B |
| RNU11 (U11) | SNORA69, SNORA120 (18), SNORA60 (97), SNORA74C-1/-2 (100), SNORD175 (45), 5 snoRNA non nommés | SNORD13, SNORD22, SNORD27, SNORD30, SNORD31B, SNORD48, SNORD49A, SNORD75L2, SNORA75 |

Spliceosome majeur : RNU2-1 (U2) est guidé par SCARNA1, 2, 4, 8, 9, 9L, 13, 14, 15, 17, 28, SNORA9, 11C, 21B,
63E, 94, 101A/B, 119, SNORD132, 137 ; RNU6-1 (U6) par SCARNA3, SCARNA23, SNORA3C, 70E, 79, 79B, SNORD7, 8, 9,
10, 67, 94 (et une dizaine d'interactions sans position).

Ces snoRNA ne visent les 19 gènes mito à intron U12 (DAP3, DNA2, GARS1, HARS2, PHB2, TXNRD2…) qu'**en
cascade** et sans aucune preuve directe : la perte d'un site de modification d'un snRNA mineur n'a jamais été
reliée à un défaut d'épissage d'un gène mito.

## Points douteux

1. snoDB mélange sites de modification (souvent issus de snoRNABase ou de prédictions de l'atlas des snoRNA)
   et interactions haut débit (PARIS, LIGR, SPLASH…). Les deux sont en `base_de_donnees` ; la `source` dit
   lequel. Les cibles codées par l'ADNmt sont suspectes (chimères).
2. Numérotation des snRNA : snoDB donne parfois deux positions pour un même site (`U12-19` et `U12.1-18`,
   `U5-45` et `U5.1-43`), décalées d'un ou deux nucléotides selon la séquence de référence. Les deux sont
   gardées (positions différentes) ; vérifier avec snoRNABase avant tout calcul positionnel.
3. Type de guide des scaRNA déduit des boîtes (ACA terminal → Ψ ; sinon boîtes C/D → Nm). Les scaRNA
   composites (SCARNA10/U85, éventuellement SCARNA12, SCARNA23) peuvent porter les deux activités :
   ex. SCARNA10 → U5 45/46 est noté 2OMe pour les deux positions. À corriger avec snoRNABase.
4. SNORD51 (homologue redondant de SNORD32A sur PXDN et peut-être COX4I2) n'est pas ajouté : les noms de
   snoRNA sont effacés du texte intégral PMC ; à vérifier dans le PDF d'Elliott 2019 et 2025.
5. miR-1 → transcrits de l'ADNmt (Zhang 2014 Cell, doi:10.1016/j.cell.2014.05.047) non ajouté : les ARNm cibles
   ne sont pas nommés dans le résumé.
6. Preuve des cibles de miR-181a/b mise en `experimentale_animale` par prudence (espèce des tests luciférase
   non vérifiée).
7. Copie snoDB tierce (dépôt GitHub) : à revalider contre l'export officiel sur Curta (md5 ci-dessus).
8. PanelApp « Mitochondrial disorders » non joint ici : la liste des cibles mito ci-dessus est un minimum.
