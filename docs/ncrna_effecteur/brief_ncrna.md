# Brief commun : le non codant « effecteur » et les phénotypes mitochondriaux

## La question de l'utilisateur
Jusqu'ici, le projet traite le non codant **régulateur** : promoteurs, enhancers (éléments rE2G), ancres de
boucles des gènes nucléaires mitochondriaux (campagne AlphaGenome 2). Le nouveau chantier est le non codant
**effecteur** : des gènes dont le produit fonctionnel est un ARN non traduit (snRNA RNU, snoRNA/scaRNA,
miRNA, lncRNA, RMRP, RPPH1, 7SK, 7SL, ARN Y, ARN de voûte, ARNr 5S/45S, ARNt nucléaires, ARN mitochondriaux
non codants…), dont une variation peut donner un **phénotype mitochondrial** chez un patient : soit par un
mécanisme mitochondrial (fonction dans la mitochondrie, contrôle de l'OXPHOS, de l'ADNmt, de la traduction
mitochondriale, import d'ARN…), soit par une phénocopie clinique d'une maladie mitochondriale (leucoencéphalopathie,
calcifications, surdité, rétinopathie, myopathie, cardiomyopathie, acidose lactique, tableau de Leigh,
épilepsie et retard de développement…). Il faut proposer à l'utilisateur une **stratégie** (demain matin).

## Contexte clinique et infrastructure
- Équipe de diagnostic et de recherche des maladies mitochondriales (CHU / INSERM). Patients atteints de maladie
  mitochondriale suspectée, non résolus après génome (WGS courtes lectures). Le laboratoire dispose des tests
  fonctionnels mitochondriaux : western blot, BN-PAGE, activités des complexes de la chaîne respiratoire,
  respiration (oxygraphie / Seahorse), sur biopsie musculaire et fibroblastes.
- Pipeline « reannot » (/home/user/reannot) : réannotation de WGS sur une **VM hors ligne** (aucune donnée
  patient ne sort ; les annotations sont préparées sur une machine connectée, le cluster « Curta », puis
  transférées en « paquet »). Ensemble de gènes : MitoCarta3.0 + panel PanelApp « Mitochondrial disorders »
  (≈ 1 220 gènes, tous protéiques).
- **Ce qui existe déjà pour les petits ARN** (à lire : /home/user/reannot/mcia/petits_arn.py,
  /home/user/reannot/vm/ajout_petits_arn.py, /home/user/reannot/vm/trier_variants.py, SOURCES.md section 4 bis,
  vm/README_MitoLoop.md, tests/test_petits_arn.sh) : petits ARN de GENCODE (snRNA, snoRNA, scaRNA, miRNA,
  misc_RNA, ribozyme, sRNA, scRNA, vault_RNA, Y_RNA ; pseudogènes et chrM exclus), contrainte nucléotide par
  nucléotide dans gnomAD v4.1 (taux Roulette), position dans la structure Rfam (Infernal), repliement ViennaRNA,
  graines de miRNA (miRBase) et sites cibles TargetScan ; colonnes ARN_* sur la VM, critère ARN_signal dans le
  tri. **Ce module ne relie pas encore ces ARN aux phénotypes mitochondriaux**, ne couvre pas les lncRNA, ni
  l'ADNmt, ni les micropeptides, ni la confirmation fonctionnelle.
- Projet « mitoloop_v0 » (/home/user/mitoloop_v0) : architecture régulatrice (Hi-C, rE2G, GTEx, polymère).
- Tissus d'intérêt de la campagne : muscle squelettique, cœur, cortex cérébral, rein ; fibroblastes et sang
  pour la confirmation chez le patient.

## Contraintes impératives
- Interdiction d'entraîner ou d'ajuster quoi que ce soit (modèle, seuil) sur des sorties AlphaGenome,
  SpliceAI, SpliceAI2 ou AlphaFold 3. Seuils fixés a priori ou tirés de données indépendantes.
- Aucune donnée patient vers un service externe ; requêtes bibliographiques génériques uniquement (aucun
  nom de personne, d'institution ou de patient du demandeur dans une requête).
- Le fichier OMIM est sous licence : on peut s'en servir localement, jamais le redistribuer ni le publier.
- Ne modifier aucun fichier des dépôts, ne rien commiter, ne lancer aucun appel AlphaGenome.
- Citer uniquement des références vérifiées (outils PubMed via ToolSearch :
  "select:mcp__PubMed__search_articles,mcp__PubMed__get_article_metadata,mcp__PubMed__find_related_articles",
  bioRxiv : "select:mcp__bioRxiv__search_preprints,mcp__bioRxiv__get_preprint", et WebSearch), avec leur DOI.
  N'inventer aucune référence, aucun gène, aucune association. Distinguer clairement établi / probable /
  hypothèse.
- Réponses et rapports en français.
