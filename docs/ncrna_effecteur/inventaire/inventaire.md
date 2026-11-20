# Inventaire des classes d'éléments non codants « effecteurs » chez l'homme

*Version 0, 20 novembre 2026. Interface machine : `inventaire_ncrna_effecteurs.tsv` (47 classes ou sous-classes, 16 colonnes, schéma fixé). Ce document en est la lecture humaine.*

## Méthode et conventions

- **Périmètre.** Éléments dont le produit fonctionnel est un ARN non traduit, ou un micropeptide porté par un transcrit dit non codant. Hors périmètre : ADNmt, promoteurs, enhancers, UTR (sauf uORF codantes), introns profonds et expansions. Les eRNA et PROMPT figurent dans la table pour être classés, mais ils sont renvoyés à l'outil des éléments régulateurs.
- **Réutilisation.** Les 329 entrées de `catalogue_ncrna_mito.tsv` et les tables des chemins 2 et 3 n'ont pas été refaites. Leurs DOI déjà vérifiés sont repris tels quels.
- **Nouvelles vérifications.** Les nouveaux DOI ont été vérifiés dans PubMed pendant cette session (titre, revue, résumé lus).
- **Comptes de loci.** Ils viennent du GTF GENCODE v27 local (copie du chemin 3) quand c'est indiqué « v27 ». Ailleurs, la mention est « à vérifier ».
- **Niveau du lien mitochondrial.**
  - `direct` : l'ARN ou le peptide agit dans la mitochondrie ou sur elle, preuve expérimentale publiée (souvent en lignées ou chez l'animal).
  - `indirect` : voie publiée qui passe par un intermédiaire.
  - `hypothese` : rien n'est démontré.
  - `aucun_connu` : aucune publication trouvée.

  Le niveau qualifie le meilleur exemple de la classe, pas chaque locus.
- **Priorité utilisateur.**
  - `+++` : classe novatrice ou peu connue ayant au moins une piste mitochondriale, ou classe déjà à fort enjeu (snRNA mineurs, ARNt, miARN, micropeptides, RMRP/RPPH1). Les ARN de voûte sont en `+++` à la demande de l'utilisateur.
  - `++` : classe nouvelle sans lien mitochondrial, ou lien faible.
  - `+` : classe à faible rendement attendu.

## Tableau par classe

Abréviations : « v27 » = GENCODE v27 ; c1/c2/c3 = chemins 1, 2 et 3.

### Spliceosome et snRNA (types 1 à 3)

| Classe (type_id) | Ce que c'est | Combien | Lien mito | Ce qu'un variant pourrait faire | Comment l'annoter |
|---|---|---|---|---|---|
| snrna_u1 | U1 : site 5' d'épissage et protection contre la polyadénylation prématurée (telescripting) | 4 RNU1 + 12 RNVU1 (v27) | indirect (c1) | Modifier l'appariement aux sites 5' des gènes mito | biotype snRNA, nom RNU1-1..4 et RNVU1 |
| snrna_majeur_u2_u4_u5_u6 | Paralogues fonctionnels pathogènes (ReNU, RNU2-2, RNU5B-1, RNU6) | Quelques fonctionnels parmi des centaines de pseudogènes | indirect (phénocopies) | Usage anormal des sites 5' | liste des paralogues du catalogue |
| snrna_mineur | U11, U12, U4atac, U6atac | 4 | indirect (c1, gènes mito à intron U12) | Rétention d'introns U12 | nom exact |
| snrna_u7 | Maturation 3' des ARNm d'histones | 1 | aucun connu | Interféronopathie | RNU7-1 |

### snoRNA et ARN apparentés (type 4)

| Classe (type_id) | Ce que c'est | Combien | Lien mito | Ce qu'un variant pourrait faire | Comment l'annoter |
|---|---|---|---|---|---|
| snorna_cd, snorna_haca | Guides de 2'-O-méthylation et de pseudouridylation | 375 SNORD, 421 SNORA (v27, copies incluses) | indirect : SNORD32A–35A (ROS ; ARNm de COX4I2), SNORA73 (respiration) | Perte du snoRNP (boîtes) ou d'une cible (guide) | biotype snoRNA + snoDB |
| scarna | Guides des snRNA dans les corps de Cajal | 49 (v27) | hypothèse | Effet indirect sur l'épissage | biotype scaRNA |
| snorna_orphelin_15q11_14q32 | Amas SNORD115/116 (paternels) et SNORD113/114 (maternels) | ~30 à 47 copies par amas | hypothèse (Prader-Willi) | Surtout des SV ou CNV | nom ^SNORD11[3-6] |
| snorna_u3_u8_u13 | Maturation des ARNr ; U8 = SNORD118 (LCC) | Quelques gènes | aucun connu | Ribosomopathie | nom exact |
| sno_lncrna_spa | lncRNA bornés par des snoRNA ; deux SPA en 15q11 qui séquestrent TDP43, RBFOX2, hnRNP M | à vérifier | hypothèse | Modification de l'épissage alternatif | BED dédié |

### miARN et petits ARN guides (type 5 et nouveautés)

| Classe (type_id) | Ce que c'est | Combien | Lien mito | Ce qu'un variant pourrait faire | Comment l'annoter |
|---|---|---|---|---|---|
| mirna_canonique | Graine de 7-8 nt | 1 881 (v27) ; MirGeneDB plus strict | **direct** : miR-1, miR-2392 (c1+c2) | Gain ou perte de cibles mito ; perte de maturation | biotype miRNA + MirGeneDB |
| mirtron | Pré-miARN = intron court, sans Drosha | à vérifier | hypothèse | Perte du miARN (sites d'épissage) | miARN chevauchant un intron (calcul) |
| agotron | Intron court chargé sur AGO, sans Dicer | à vérifier | hypothèse | Graine | BED dédié |
| sdrna | Petits ARN dérivés de snoRNA à fonction de miARN (ACA45) | à vérifier | hypothèse | Graine | snoDB + liste publiée |
| pirna_cluster | piARN germinaux en clusters ; les « piARN » somatiques des bases sont surtout des fragments d'autres ARN | Centaines de loci précurseurs (à vérifier) | indirect : biogenèse à la surface mitochondriale (MITOPLD, souris) | Peu prédictible | BED des clusters testiculaires |

### ARNt et dérivés (type 10 et nouveautés)

| Classe (type_id) | Ce que c'est | Combien | Lien mito | Ce qu'un variant pourrait faire | Comment l'annoter |
|---|---|---|---|---|---|
| trna_nucleaire | ARNt Pol III | ~420 de haute confiance (GtRNAdb, à vérifier) | indirect (import de l'ARNt-Gln, ISR) (c2) | Changement d'identité (anticodon), perte d'expression | GtRNAdb score ≥ 50 |
| tdr_trf_tirna | tRF et moitiés d'ARNt produites par l'angiogénine | Dérivés des ARNt | **direct** : tiRNA-5 His/Glu enrichis dans la mitochondrie (rat, Musashi-1, Mpc1) ; tDR induits par l'angiogénine sous inhibiteurs de la chaîne respiratoire (c1+c2) | Modification du site de coupure ou de la graine | région du fragment dans l'ARNt |
| trna_lookalike_mt | « Sosies » nucléaires d'ARNt mitochondriaux, parfois transcrits | Centaines (à vérifier) | hypothèse | Inconnu ; risque d'artefact NUMT | BED Telonis 2014 |

### ARN ribosomiques et ADNr (type 11)

| Classe (type_id) | Ce que c'est | Combien | Lien mito | Ce qu'un variant pourrait faire | Comment l'annoter |
|---|---|---|---|---|---|
| rrna_5s | ARNr 5S | 17 + 496 pseudogènes (v27) | **controversé** (import rapporté, absent du mitoribosome) | Effet nul attendu (redondance) | nom ^RNA5S[0-9]+$ |
| rrna_45s | Unités d'ADNr 45S | Centaines de copies | indirect (association nombre de copies d'ADNr / ADNmt) | CNV seulement | non interprétable en SNV |
| rdna_pRNA_papas | ARN antisens du promoteur de l'ADNr (PAPAS) codant la protéine RIEP | Répété | **direct** (c3) : RIEP mitochondriale, partenaires C1QBP et CHCHD2 | Inconnu (multicopie) | BED dédié |

### ARN de voûte (type 11)

| Classe (type_id) | Ce que c'est | Combien | Lien mito | Ce qu'un variant pourrait faire | Comment l'annoter |
|---|---|---|---|---|---|
| vault_rna | vtRNA1-1/1-2/1-3 ; vtRNA1-1 freine l'autophagie via p62 | 3 | hypothèse (mitophagie) | ΔΔG, perte d'expression (boîtes A/B) | nom ^VTRNA1-[123]$ ; biotype vault_RNA/vaultRNA/misc_RNA |
| vault_nc886 | VTRNA2-1, inhibiteur de PKR, empreinte métastable | 1 | hypothèse (PKR et ARNdb mitochondrial) | Perte d'inhibition de PKR ; la méthylation est invisible en WGS | VTRNA2-1 / MIR886 |
| svrna | Petits ARN dérivés des vtRNA, chargés sur AGO (svRNAb → CYP3A4) | Dérivés de 4 vtRNA | hypothèse | Graine | BED de la région svRNA |

### Autres ARN Pol III et ARN de RNP (types 9 et 11)

| Classe (type_id) | Ce que c'est | Combien | Lien mito | Ce qu'un variant pourrait faire | Comment l'annoter |
|---|---|---|---|---|---|
| y_rna / ysrna | RNY1/3/4/5 (Ro60, réplication) et leurs fragments | 4 + 758 rétrocopies (v27 ; RNY5 absent du v27) | aucun connu | ΔΔG | noms exacts |
| srp_7sl | RN7SL1/2/3 | 3 + ~850 copies | aucun connu | ΔΔG | noms exacts |
| alu_rna | Transcrits Alu | > 10⁶ copies | indirect : l'ARN Alu provoque la sortie de l'ADNmt vers cGAS (DMLA) | Non interprétable en SNV | signalement de SV seulement |
| bc200 | BCYRN1, ARN dendritique | 1 | aucun connu | ΔΔG | BCYRN1 |
| snar | snaR, liés à NF90, propres aux grands singes | Dizaines (à vérifier) | aucun connu | Inconnu | BED |
| rn7sk | 7SK / P-TEFb | 1 | aucun connu | ΔΔG, sites HEXIM/LARP7 | RN7SK |
| rmrp_rpph1 | RNase MRP, RNase P | 2 | hypothèse (import PNPase, rôle contesté) (c2) | Élément d'import, cœur catalytique | noms exacts ; biotype ribozyme en v27 |
| terc | ARN de la télomérase | 1 | indirect (télomères → p53 → PGC-1) + import (TERC-53) (c1+c2) | Positions pathogènes connues | TERC |
| mascrna | ARN de type ARNt coupé en 3' de MALAT1 et NEAT1 | 2 | hypothèse | Défaut de maturation de MALAT1 | biotype sRNA, Rfam RF01684 |

### Répétitions, ribozymes, lncRNA, circRNA

| Classe (type_id) | Ce que c'est | Combien | Lien mito | Ce qu'un variant pourrait faire | Comment l'annoter |
|---|---|---|---|---|---|
| terra | ARN des télomères | Subtélomères | indirect : TERRA-ZBP1 sur la membrane externe mitochondriale → MAVS | Non interprétable | signalement |
| satellite_rna | HSATII, alpha-satellites | Mb | aucun connu | — | — |
| ribozyme_humain | CPEB3, CoTC, hammerhead | 4-5 | aucun connu | Perte d'autoclivage | biotype ribozyme |
| lncrna_architectural | NEAT1, MALAT1 | 2 | indirect : NEAT1 retient des ARNm mito dans les paraspeckles (c1+c2) | Inconnu sauf domaines | noms exacts |
| lncrna_mito | SAMMSON, GAS5, CEROX1, TUG1, LINC00467, PINK1-AS, SRA | ~10 | **direct** (c1+c2) | Sites d'éponge, domaines | liste du catalogue |
| lncrna_general | Cas général et lncRNA soumis à empreinte | ~13 000 (v27) | hypothèse | Contrainte seulement | biotypes lncRNA |
| pseudogene_eponge | Pseudogènes transcrits jouant l'éponge (PTENP1) | ~1 400 transcrits (v27) | hypothèse | Sites miARN partagés avec un parent mito | biotypes transcribed_* |
| circrna | circRNA ; plusieurs localisés dans la mitochondrie (SCAR, circSamd4, circBRAKE) | Dizaines de milliers (circAtlas) | **direct** (c2 ; lignées et souris) | Perte de la jonction de rétro-épissage | circAtlas, ±10 nt des jonctions |
| ciRNA_sisRNA | Lariats introniques stables | à vérifier | hypothèse | Perte des éléments GU-riche ou C-riche | BED |
| erna_prompt | Transcrits d'enhancers et de promoteurs | Dizaines de milliers | — | Hors périmètre | autre outil |

### Micropeptides (type 12)

| Classe (type_id) | Ce que c'est | Combien | Lien mito | Ce qu'un variant pourrait faire | Comment l'annoter |
|---|---|---|---|---|---|
| micropeptide_lncrna_sorf | sORF de « lncRNA » : MTLN, SMIM26, ASAP, PIGBOS | ≥ 7 000 ORF non canoniques traduites ; une dizaine mito | **direct** (c3) | Comme une CDS | Ribo-seq ORFs GENCODE, OpenProt |
| micropeptide_uorf_altorf | uORF et altORF : SLC35A4-MP, altMiD51, MP31 | Milliers | **direct** (c3) | Effet sur le peptide et sur la CDS | idem |
| micropeptide_numt_humanin | MTRNR2L1-13 (peptides de type humanine) | 13 | indirect (c3) | Comme une CDS ; attention aux NUMT | nom ^MTRNR2L |

## Classes nouvelles par rapport aux 12 types

- **ARN de voûte.** Sous-classes séparées : vtRNA1, nc886 et **svRNA**.
- **Fragments et dérivés.**
  - **tRF/tiRNA** ;
  - **sosies d'ARNt mitochondriaux** ;
  - **YsRNA** ;
  - **sdRNA** ;
  - **mirtrons** et **agotrons** ;
  - **mascRNA/menRNA**.
- **piARN** et clusters.
- **sno-lncRNA/SPA**, et les snoRNA orphelins comme sous-classe à part.
- **Répétitions et ARN apparentés.**
  - **ARN Alu** ;
  - **BC200** ;
  - **snaR** ;
  - **TERRA** ;
  - **ARN satellites**.
- **Ribozymes humains** (CPEB3, CoTC, hammerhead).
- **ARN de l'ADNr** : **pRNA/PAPAS**, porteur de la protéine RIEP.
- **circRNA**, **ciRNA/sisRNA** et **pseudogènes éponges**.
- **Micropeptides de type humanine (NUMT)**.
- Classement explicite des **eRNA/PROMPT**, hors périmètre.

## Classes les plus novatrices à surveiller

1. **Fragments d'ARNt (tiRNA-5, i-tRF).** C'est la seule classe de petits ARN nucléaires avec un enrichissement mitochondrial et un effet sur la respiration publiés :
   - Bayazit 2022, îlots de rat ;
   - Al-Mesitef 2026, HEK293T sous inhibiteurs des complexes I à V.

   Un variant d'ARNt touchant la région 5' d'un isodécodeur His ou Glu mérite d'être signalé, avec le niveau hypothèse.
2. **PAPAS/RIEP.** Une protéine de 25 kDa propre aux primates, codée par un ARN de l'ADNr dit non codant. Elle est localisée dans la mitochondrie et se lie à C1QBP et CHCHD2, deux gènes de maladie mitochondriale (Feng 2023). La détection WGS est presque impossible (multicopie) ; c'est une piste de recherche.
3. **ARN de voûte et svRNA.** Les pistes p62/mitophagie (vtRNA1-1) et PKR/ARNdb mitochondrial (nc886) restent des hypothèses. Les svRNA (svRNAb) sont de vrais petits ARN guides : leur graine peut être testée sur les cibles mito comme pour un miARN.
4. **circRNA localisés dans la mitochondrie.** Trois exemples :
   - SCAR lie ATP5B et ferme le mPTP ;
   - circSamd4 limite les ROS (souris) ;
   - circBRAKE, dérivé de Birc6, désassemble l'ATP synthase (souris, 2026).

   Un variant de jonction de rétro-épissage est prédictible ; l'ARN-seq est nécessaire pour confirmer.
5. **TERRA → ZBP1 → MAVS.** Signalisation télomère-mitochondrie (Nassour 2023), voie indirecte. Elle n'est pas exploitable en SNV, mais elle est à garder pour l'interprétation de phénotypes inflammatoires.
6. **ARN Alu → ADNmt cytosolique → cGAS** (Kerur 2018). Même statut que TERRA.
7. **mascRNA** (3' de MALAT1/NEAT1) et **sosies d'ARNt mitochondriaux.** Ce sont des ARN de type ARNt mal annotés, sources possibles de confusion avec des fragments mitochondriaux.

## Limites et points douteux

- **Accessions Rfam.** Seules RF00622 (ribozyme CPEB3) et RF01684 (mascRNA-menRNA) ont été vérifiées en ligne. Les autres (RF00001-RF00030, RF00066, RF00096, RF00100, RF00548, RF00618, RF00619) sont des accessions usuelles reprises de mémoire : **à vérifier** sur `family.txt` de Rfam sur Curta avant usage. Les accessions de snaR, BC200 et CoTC n'ont pas été trouvées. Ces classes s'appuient sur `gtf_nom`.
- **Comptes de loci.** Ils viennent de GENCODE v27 et incluent les pseudogènes. Les orthographes de biotypes évoluent : vaultRNA en v27, vault_RNA ensuite ; lincRNA et antisense regroupés sous lncRNA dans les versions récentes. RNY5 et les gènes SNAR sont absents du GTF v27. Les ARNt sont absents du GTF principal (fichier tRNA séparé ou GtRNAdb).
- **Liens mitochondriaux.** Ils sont majoritairement expérimentaux, en lignées ou chez l'animal :
  - tiRNA chez le rat ;
  - circSamd4 et circBRAKE chez la souris ;
  - MITOPLD chez la souris ;
  - SCAR dans des fibroblastes et chez la souris.

  L'origine génomique de SCAR (nucléaire ou ADNmt) n'a pas été vérifiée. Si SCAR provient de l'ADNmt, il sort du périmètre.
- **Classes contestées.**
  - L'import du 5S, de RMRP et de RPPH1 reste discuté.
  - L'APEX-seq de la matrice ne trouve pas d'ARN nucléaire fortement enrichi (voir la table du chemin 2).
  - Les « piARN somatiques » des bases sont surtout des fragments d'autres ARN : ne jamais annoter un variant comme « piARN » hors cluster germinal.
- **BED dédiés à construire.** Mirtrons, agotrons, sdRNA, svRNA, sno-lncRNA/SPA, ciRNA, snaR, sosies d'ARNt et clusters de piARN : aucune liste n'a été téléchargée (réseau limité). Les règles de sélection sont écrites, les fichiers sont à produire sur Curta à partir des publications citées.
- **Répétitions** (Alu, satellites, TERRA, ADNr, 5S, amas SNORD116) : les SNV ne sont pas interprétables en lectures courtes. Seuls les CNV et les SV le sont, et seulement en partie.
- **Méthylation.** Le niveau de nc886 dépend surtout de la méthylation, invisible en WGS.
- **Niveau de preuve.** Aucune de ces classes n'apporte de points ACMG pour un lien mitochondrial. Les niveaux `hypothese` servent à signaler, pas à classer.
