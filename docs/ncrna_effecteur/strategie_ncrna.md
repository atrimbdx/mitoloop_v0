# Non codant « effecteur » et phénotypes mitochondriaux : catalogue et stratégie

*Rapport du 19 novembre 2026, version 2 après contre-revue (références toutes vérifiées existantes ; barème des cas types, ADNmt et RNU2-2 corrigés ; 12 gènes restés sans verdict). Rien n'a été modifié dans les dépôts `reannot` et `mitoloop_v0`, rien n'a été commité, aucun appel AlphaGenome n'a été lancé. Aucune donnée patient n'a servi aux recherches bibliographiques. Les références de la section (h) ont été vérifiées dans PubMed (ou sur medRxiv et bioRxiv pour les préprints, signalés comme tels). Les points que personne n'a encore relus dans le texte intégral sont marqués « à relire ».*

Vocabulaire de preuve utilisé partout :
- **établi** : relation gène-maladie mendélienne reconnue (plusieurs familles indépendantes, preuve fonctionnelle) ;
- **probable** : quelques familles, ou une seule avec un appui fonctionnel fort ;
- **cas isolés** : une famille, ou des associations ;
- **expérimental** : données cellulaires ou animales, sans patient ;
- **hypothèse** : raisonnement plausible, sans donnée directe ;
- **aucun** : rien de publié.

---

## (a) Synthèse en 15 lignes

1. Le catalogue vérifié compte environ 150 entrées confirmées. Le seul mécanisme mitochondrial **direct et établi chez l'homme** porté par un ARN non codant est celui des ARNt et ARNr de l'**ADNmt**.
2. Aucun ARN nucléaire n'a de mécanisme mitochondrial primaire établi chez l'homme. Les mitomiR, snoRNA de RPL13A, lncRNA (CEROX1, TUG1, SAMMSON…) restent **expérimentaux**.
3. Les ARN nucléaires utiles au diagnostic sont des **phénocopies** : RNU4-2, RNU2-2, RNU5B-1, RNU4ATAC, RNU6ATAC, RNU12, RNU7-1, SNORD118, RNU6-x, MIR96, MIR204, plus RMRP et TERC.
4. reannot annote déjà les petits ARN (contrainte, structure, repliement, miRNA), mais ne les relie à aucun phénotype. Il ignore l'ADNmt, les lncRNA effecteurs et la preuve fonctionnelle.
5. Priorité 1 : une **route chrM** (réappel depuis le CRAM, hétéroplasmie, bases dédiées). Elle dépend d'une décision d'infrastructure préalable (**D0** : CRAM et conteneurs sur la VM).
6. Priorité 2 : un **catalogue versionné** branché sur `petits_arn`, avec un critère de tri réservé aux gènes établis à phénocopie (N2) et une **voie « TND »** explicite pour RNU4-2, RNU2-2 et RNU5B-1.
7. Les conclusions **D1 et D2** suivent désormais la classification **ACMG/AMP en points**. La chaîne de preuve M0 à M4 n'est plus qu'une liste de vérifications.
8. Pour un ARN nucléaire, le bilan OXPHOS **qualifie le phénotype**, sans prouver ni innocenter le variant. La preuve vient de l'ARN (M1) et de sa signature en aval (M2), qui comptent pour **un seul** code PS3, supporting tant qu'il n'existe pas 11 variants témoins validés.
9. Pour l'ADNmt, les activités, le BN-PAGE et le Seahorse **décrivent** le phénotype. Cybrides et PCR en fibre unique valent PS3_supporting (McCormick 2020). La preuve forte vient de la ségrégation et du caractère de novo.
10. Une phénocopie confirmée se rend comme un diagnostic **« non mitochondrial »**, avec orientation vers l'équipe concernée.
11. Aucun seuil n'est ajusté sur AlphaGenome, SpliceAI, SpliceAI2 ou AlphaFold 3, ni sur des variants de patients. Les positifs publiés mesurent le rappel, des génomes déjà résolus mesurent la spécificité.
12. Rendement attendu modeste, et c'est une hypothèse : environ 0,6 % dans la seule série mitochondriale publiée (2 RNU4-2 sur 357), probablement davantage avec RNU2-2 récessif et l'ADNmt.
13. Effort : environ 4 à 5 semaines de bio-informatique avec 1 personne (route chrM en plus, 2 à 4 semaines), ou 4 semaines avec 2 personnes en parallèle.
14. Les diagnostics confirmés par l'ARN sont attendus à **6 à 9 mois**. Ils dépendent d'une référence locale d'ARN-seq d'ARN total et de témoins positifs obtenus par collaboration.
15. Treize décisions vous attendent (§ g). Les plus urgentes sont D0 (infrastructure), D1 (qui porte chrM), D2 (génotypes dans l'entrée) et D12 (découvertes secondaires, dont MT-RNR1).

---

## (b) La question et le périmètre

**La question.** Jusqu'ici, le projet traite le non codant **régulateur** : promoteurs, éléments rE2G, ancres de boucles des gènes nucléaires mitochondriaux. Le nouveau chantier est le non codant **effecteur** : des gènes dont le produit fonctionnel est un ARN non traduit, et dont une variation peut donner un phénotype mitochondrial chez un patient :
- soit par un **mécanisme mitochondrial** (fonction dans la mitochondrie, contrôle de l'OXPHOS, de l'ADNmt, de la traduction mitochondriale) ;
- soit par une **phénocopie clinique** (leucoencéphalopathie, calcifications, surdité, rétinopathie, myopathie, cardiomyopathie, acidose lactique, tableau de Leigh, épilepsie et retard).

**Les patients.** Maladie mitochondriale suspectée, non résolus après WGS en lectures courtes. Le laboratoire dispose du western blot, du BN-PAGE, des activités des complexes, de l'oxygraphie et du Seahorse, sur muscle et fibroblastes.

**L'infrastructure.** reannot réannote les WGS sur une **VM hors ligne**. Les annotations sont préparées sur le cluster connecté Curta puis transférées en « paquet ». L'ensemble de gènes (MitoCarta3.0 et PanelApp « Mitochondrial disorders », environ 1 220 à 1 258 symboles) est entièrement protéique.

**Classes d'ARN couvertes.** snRNA (spliceosomes majeur et mineur, U7), snoRNA et scaRNA, miRNA, lncRNA, ARN Pol III divers (RMRP, RPPH1, 7SK, 7SL, ARN Y, ARN de voûte, ARNt-Sec, TERC), ARNr 5S et 45S, ARNt nucléaires, ARN de l'ADNmt (ARNt, ARNr, région de contrôle, lncRNA antisens).

**Hors périmètre, mais signalés.** Les gènes protéiques de la machinerie des ARN (POLR3A/B, EXOSC*, TSEN*, LARP7…) et les micropeptides (MTLN, UQCC5/6…) relèvent de l'ensemble de gènes protéiques (décision D5). Les loci à expansion (FXN, ATXN8OS, NUTM2B-AS1…) relèvent d'une route STR.

**Contraintes.**
- Aucun seuil ni modèle entraîné ou ajusté sur des sorties AlphaGenome, SpliceAI, SpliceAI2 ou AlphaFold 3.
- Aucune donnée patient vers un service externe.
- OMIM utilisé localement seulement, jamais redistribué.

---

## (c) Ce qui existe déjà dans reannot

Chemins relatifs à `/home/user/reannot` (R) et `/home/user/mitoloop_v0` (M). Lecture seule.

### c.1 Module petits ARN

**Préparation sur Curta** (`R/mcia/petits_arn.py`, étape `do_petits_arn` de `R/mcia/run_paquet.sh`) :
- **gènes** : types GENCODE snRNA, snoRNA, scaRNA, miRNA, misc_RNA, ribozyme, sRNA, scRNA, vault_RNA, Y_RNA. Exclus : chrM (`:126`), pseudogènes, lncRNA, rRNA, gènes de plus de 500 nt, copies `_PAR_Y`. Toutes les copies sont gardées, sans aucune donnée d'expression ;
- **contrainte** : gnomAD v4.1 génomes, SNV PASS, attendu Roulette (sinon trinucléotidique), calibration sur des flancs de ±500 pb, fenêtre de ±10 pb, O/E et borne haute de l'IC à 95 % ;
- **structure** : Infernal et Rfam (base appariée, partenaire) ;
- **repliement** : ViennaRNA, ΔMFE de tous les SNV, gènes de 400 nt au plus ;
- **miRNA** : graines miRBase et sites TargetScan 8 ;
- l'étape **n'est pas comprise dans `all`** ; Rfam et miRBase sont téléchargés sans numéro de version.

**Annotation sur la VM** (`R/vm/ajout_petits_arn.py`) : colonnes `ARN_gene`, `ARN_contrainte`, `ARN_structure`, `ARN_repliement`, `ARN_miRNA` et drapeau `ARN_signal`. Ce drapeau est `oui` si la position est contrainte (borne haute < 1, **en dur**), si un appariement est cassé ou si une graine est touchée.

**Tri** (`R/vm/trier_variants.py:88-90`) : un variant passe si `ARN_signal = oui`, quel que soit le gène, et il est classé après les priorités MitoLoop.

**Test** : `R/tests/test_petits_arn.sh`, synthétique. Il ne couvre ni le mode Roulette, ni les indels, ni l'étape 4b.

### c.2 lncRNA

Ils ne sont vus que comme **régulateurs en cis** d'un gène mitochondrial (`R/mcia/mitoloop_paquet.py:223-316`, colonne `MITO_lncRNA`) : lncRNA divergent ou antisens, promoteur et jonctions d'épissage. Une délétion d'un lncRNA compte comme événement SV modéré. Les exons d'un lncRNA effecteur ne sont jamais annotés.

### c.3 ADNmt

**Rien n'est fait** :
- reannot part de listes de variants déjà appelées ;
- `R/vm/preparer_entree.py:200-212` ne garde que chrom, pos, ref et alt. FORMAT, AF et l'hétéroplasmie sont perdus ;
- VEP tourne sans `--af` et sans base mitochondriale ;
- chrM est exclu des petits ARN (`petits_arn.py:126`) et des régions MitoLoop (`R/mcia/vcf_annotate.py:335-336`).

### c.4 Défauts qui gênent le nouveau chantier

- **OMIM du voisin.** Le gène attribué à un variant est le transcrit MANE codant le plus proche dans un rayon de 1 Mb (`preparer_entree.py:215-243`). Un variant de RNU4-2 reçoit donc l'OMIM d'un gène codant voisin.
- **Clé par nom.** La table de contrainte par gène est indexée par nom (`ajout_petits_arn.py:95-104`, `:174`), ce qui mélange les copies portant le même nom.
- **Pas de fréquence** de population dans la table finale.
- **Pas de lien** vers l'ARN-seq, le petit ARN-seq ou les tests fonctionnels.
- **Aucun positif** ARN non codant ni ADNmt pour l'évaluation (`R/donnees/positifs/positifs_mitoloop.tsv` : 30 variants, tous codants).
- **Traçabilité.** `MANIFEST.tsv` est écrit avant la copie de `petits_arn/`, donc sans les versions de Rfam, miRBase et TargetScan.
- **Documentation.** `ajout_petits_arn.py` manque dans `README_MitoLoop.md:20-24` et dans l'en-tête de `07_maj_mitoloop_vm.sh`.
- **Infrastructure.** Aucun script de `vm/` ne lit de BAM ni de CRAM, et le dépôt ne mentionne aucun moteur de conteneurs (docker, singularity, apptainer).

### c.5 Conventions à respecter

- Sous-commandes en français ; docstring « — Curta. » ou « — VM. » ; colonnes préfixées (`ARN_*`) ; valeurs `oui`, `non`, `.` ; arguments plutôt que score.
- Seuils exposés en options et fixés a priori.
- Sur la VM, bibliothèque standard plus pysam ou tabix.
- Tests `tests/test_<module>.sh` hors ligne, avec `PASS test_x`.
- OMIM jamais dans le dépôt ni dans le paquet. reannot n'a pas de fichier LICENSE.

---

## (d) Le catalogue

**Catalogue complet** (toutes les entrées, y compris les statuts `contexte_non_verifie` et `non_verifie_priorite_basse`, avec toutes les références et corrections) : `/tmp/claude-0/-home-user/04d773a8-9dce-5a49-9259-4d11be1339d9/scratchpad/ncrna/catalogue_ncrna_mito.tsv`.

**Colonnes des tableaux** :
- **Lien mito** : direct, indirect (mécanisme mitochondrial indirect), phénocopie, hypothèse, aucun. C'est le lien vérifié, après corrections.
- **Niveau** : niveau de preuve gène-maladie (vocabulaire défini en tête).
- **Test attendu** : ce qui devrait être anormal chez un porteur. « OXPHOS normal » signifie activités, BN-PAGE, respiration et ADNmt attendus normaux.

### d.1 ARN de l'ADNmt (route chrM)

**ARNt mitochondriaux.** Pour tous : muscle (COX-/SDH+, PCR en fibre unique), activités des complexes et BN-PAGE (déficit combiné CI, CIII et CIV en général), northern de l'ARNt, hétéroplasmie dans plusieurs tissus. Seules les particularités sont indiquées.

| Gène | Maladie | Lien mito | Niveau | Particularité du test |
|---|---|---|---|---|
| MT-TA | Myopathie proximale ; PEO | direct | établi | Charge plus élevée dans les fibres COX- |
| MT-TR | Encéphalomyopathie avec déficit en COX (cas unique) | direct | probable | m.10437G>A de novo, absent du sang et des fibroblastes : **muscle ou urine obligatoires** |
| MT-TN | Encéphalomyopathie, PEO, épilepsie | direct | établi | CI et CIV diminués ; WB MT-CO1 et NDUFB8 |
| MT-TD | Myopathie (rapports isolés) | direct | probable | m.7526A>G probablement de novo |
| MT-TC | Surdité ; encéphalopathie (ClinGen Limited) | direct | probable | ARNt-Cys et ARNt-Tyr diminués ; notation brin L à harmoniser |
| MT-TQ | Cardiomyopathie hypertrophique (possible) | direct | cas isolés | m.4336T>C est un variant d'haplogroupe, à ne pas traiter comme pathogène |
| MT-TE | Myopathie infantile réversible (m.14674T>C homoplasmique) ; myopathie et diabète | direct | établi | Biopsie précoce ; la fibre unique n'apporte rien si homoplasmique ; **vérifier que l'annotateur inclut la borne 14674** |
| MT-TG | Cardiomyopathie hypertrophique | direct | probable | Traduction diminuée surtout pour les protéines riches en glycine ; « diabète » non vérifié |
| MT-TH | Encéphalomyopathie (MELAS ou MERRF), rétinopathie | direct | établi | Les délétions qui emportent MT-TH sont un objet distinct |
| MT-TI | Cardiomyopathie hypertrophique maternelle (m.4300A>G) ; myopathie, ataxie, surdité ; tubulopathie | direct | établi | Muscle squelettique possiblement normal : **tissu cardiaque** ou cybrides. m.4267A>G est sporadique |
| MT-TL2 | PEO, myopathie, diabète, surdité | direct | établi | m.12315G>A restreint au muscle, invisible sur un WGS sanguin |
| MT-TL1 | MELAS, MIDD, PEO, myopathie, cardiomyopathie, atteinte rénale | direct | établi | Sang avec correction d'âge, urine, muscle ; déficit CI prédominant |
| MT-TK | MERRF ; Leigh, cardiomyopathie, surdité | direct | établi | Protéines riches en lysine plus touchées |
| MT-TM | Retard, intolérance à l'effort, myopathie ; Leigh rapporté | direct | probable | Cas Born 2015 détecté seulement dans le muscle et les fibroblastes cutanés, absent des autres tissus et chez la mère (de novo probable) : invisible sur un WGS sanguin ; muscle ou fibroblastes |
| MT-TF | Néphropathie tubulo-interstitielle avec hyperuricémie ; épilepsie ; MELAS | direct | établi | Muscle souvent normal ; m.591C>T pathogène (7 familles), m.643A>G VUS ; couverture faible près de la boucle D |
| MT-TP | Myopathie, encéphalomyopathie | direct | établi | m.16023G>A « Reported » seulement dans MITOMAP |
| MT-TS2 | Encéphalomyopathie avec déficit CI | direct | probable | m.12207G>A (Wong 2006), m.12258C>A (Whittaker 2007) |
| MT-TS1 | Surdité ± kératodermie ; encéphalomyopathie | direct | établi | m.7445A>G est **hors des coordonnées de l'ARNt** (site de maturation 3') ; m.7471dupC = m.7472insC |
| MT-TT | Associations contestées (LHON, coronaropathie) | direct | probable | m.15923A>G LP (ClinGen, non revérifié directement) ; haplogroupe obligatoire |
| MT-TW | Encéphalomyopathie, Leigh, myopathie, cardiomyopathie | direct | établi | m.5545C>T pathogène à **faible hétéroplasmie** : appel à faible VAF |
| MT-TY | Intolérance à l'effort, déficit CIII | direct | établi | Variant musculaire invisible dans le sang : ADNmt du muscle obligatoire |
| MT-TV | Leigh (y compris à l'âge adulte), encéphalomyopathie | direct | établi | Hétéroplasmie sanguine possible à 2 % : appel à faible VAF et contrôle urinaire ; phénocopie rénale m.1642G>A |

**ARNr et région de contrôle**

| Gène ou élément | Maladie | Lien mito | Niveau | Test attendu |
|---|---|---|---|---|
| MT-RNR1 | m.1555A>G : surdité non syndromique ou induite par les aminosides ; peptide MOTS-c (souris) | direct | établi | Variant homoplasmique, identique dans tous les tissus ; OXPHOS normal. **Pharmacogénétique** : m.1555A>G, m.1494C>T, m.1095T>C (CPIC) |
| MT-RNR2 | Une famille avec cardiomyopathie hypertrophique (m.2336T>C) ; humanine et SHLP (expérimental) | direct | cas isolés | ARNr 16S diminué (rapport MT-RNR2/MT-RNR1), MRPL en WB, traduction mitochondriale |
| ARN 7S | Régulateur négatif de la transcription mitochondriale ; aucune maladie | direct | expérimental | RNA-seq brin-spécifique ; région HVS-II et CSB très polymorphe |
| MT-OHR | Variant de longueur à chrM:302 (dans CSB2) qui contrôle l'abondance de l'ADNmt | direct | cas isolés (association) | Nombre de copies d'ADNmt, avec m.302 en covariable |
| MT-OLR | Aucune maladie propre | direct | expérimental | Nombre de copies, délétions |
| LSP2 | Néphropathie tubulo-interstitielle maternelle (**préprint medRxiv 2026**) | direct | probable, à confirmer | ARNt du brin L et MT-ND6 en RNA-seq orientée |
| MT-LSP | Aucune maladie mendélienne ; mutations somatiques dans le cerveau Alzheimer (seul m.414 est dans le LSP ; études négatives) | direct | cas isolés (association) | MT-ND6, ARN 7S, nombre de copies |
| MT-HSP1 | Aucune | direct | expérimental | ARNt du brin H, traduction |
| CSB3 | Aucune ; rôle modulateur in vitro | direct | hypothèse | Seulement si variant rare |
| lncND5, lncND6, lncCytb | Aucune ; biogenèse directe, fonction régulatrice non démontrée ; LIPCAR est un biomarqueur distinct | direct (biogenèse) | expérimental | RNA-seq orientée ; explorer le gène chevauchant (CI ou CIII) ; ARN double brin et ISG chez les patients PNPT1 |

### d.2 snRNA

| Gène | Maladie | Lien mito | Niveau | Test attendu |
|---|---|---|---|---|
| RNU4-2 | Syndrome ReNU (dominant de novo) ; forme récessive avec atteinte de la substance blanche (Rius 2026) ; rétinite pigmentaire | phénocopie | établi | **WGS** (l'exome manque n.64_65insT) ; ARN-seq d'ARN total : sites 5' alternatifs (A5SS) ; épisignature ; en récessif, **niveau de RNU4-2 diminué**. 2 cas sur 357 avec OCR diminuée dans une série mitochondriale |
| RNU2-2 | Encéphalopathie développementale et épileptique (dominante) ; TND récessif | phénocopie | établi | Récessif : rapport U2-2/U2-1 abaissé dans le sang (ARN total déplété en ARNr) ; épissage visible dans les **fibroblastes**, pas dans le sang |
| RNU5B-1 | TND dominant de novo | phénocopie | établi | Lymphocytes cultivés : anomalies d'épissage discrètes, propres au variant |
| RNU5A-1 | Candidat TND (2 à 3 probands, une étude) | aucun | cas isolés | Exploratoire seulement |
| RNU4ATAC | MOPD1/Taybi-Linder, Roifman, Lowry-Wood ; formes Joubert-like, diabète auto-immun | phénocopie | établi | Sang : excès d'outliers de rétention d'introns dans les gènes à intron mineur (MIG) ; OXPHOS non spécifique |
| RNU6ATAC | Diabète auto-immun très précoce avec hypogammaglobulinémie ; spliceopathie mineure (environ 10 individus, 7 familles) | phénocopie | établi | Même signature MIG dans le sang ; immunophénotype B |
| RNU12 | SCAR33 (ataxie précoce) ; syndrome CDAGS | phénocopie | probable | Rétention U12 dans le sang ; pour l'ataxie seule, une famille |
| RNU7-1 | Syndrome d'Aicardi-Goutières (récessif) | phénocopie | établi | Fibroblastes en prolifération : ARNm d'histones lus au-delà de la tige-boucle et polyadénylés ; score ISG |
| RNU6-1, RNU6-2, RNU6-8, RNU6-9 | Rétinite pigmentaire autosomique dominante | phénocopie | établi | OXPHOS normal ; ERG et OCT ; contrôle du paralogue (n.55_56insG) |
| RNU4-1 | Aucune propre ; dénominateur d'expression ; forme AR proposée dans un préprint (2 patients) | aucun | cas isolés | Augmenté dans la forme récessive de RNU4-2 |
| RNU2-1 | Aucune ; réseau multicopie ; dénominateur du rapport U2-2/U2-1 | aucun | aucun | Agréger les copies |
| RNU6-7 | Aucune ; porte n.55_56insG chez des témoins | aucun | aucun | — |
| RNU5D-1 | Aucune | aucun | aucun | — |
| RNU5E-1, RNU5F-1 | Aucune ; 2 et 3 variants de novo de signification inconnue (cohorte PFMG2025) | aucun | hypothèse | Feuille « recherche » |
| RNU1-1 à RNU1-4 | Aucune maladie germinale ; mutation somatique r.3A>G dans le médulloblastome SHH | aucun | aucun | **Zone aveugle** en lectures courtes |
| RNVU1-8 | Aucune | aucun | expérimental | — |
| RNU11 | Aucune maladie humaine ; souris : microcéphalie | aucun | expérimental | — |

### d.3 snoRNA

| Gène | Maladie | Lien mito | Niveau | Test attendu |
|---|---|---|---|---|
| SNORD118 | Leucoencéphalopathie avec calcifications et kystes (LCC, Labrune) | phénocopie | établi | Panel mitochondrial **normal**, qui n'exclut pas la LCC ; IRM ; U8 et précurseur prolongé en 3' |
| SNORD116-1 (amas) | Prader-Willi (microdélétions de l'amas ; la plus courte publiée fait 71 kb) | phénocopie | établi | Analyse de CNV : la méthylation SNRPN est normale dans ces microdélétions |
| SNORD114-1 | Syndrome de Temple (porté par le locus, DLK1) | aucun (corrigé) | expérimental | Méthylation 14q32 |
| SNORA31 | Encéphalite herpétique (HSV-1) | aucun | probable | OXPHOS normal |
| SNORA73A | Sa perte protège du stress lipotoxique (avec SNORA73B) | indirect | expérimental | Respiration attendue **augmentée** |
| SNORD32A, SNORD33, SNORD34, SNORD35A | Stress lipotoxique et oxydant (knock-out du locus RPL13A) | indirect au niveau du locus ; hypothèse pour chaque snoRNA seul | expérimental | Aucun déficit OXPHOS attendu |
| SNORD115-1 | Aucune maladie démontrée | aucun | expérimental | — |
| SNORD3A | Aucune | aucun | expérimental | — |

### d.4 miRNA

| Gène | Maladie | Lien mito | Niveau | Test attendu |
|---|---|---|---|---|
| MIR96 | Surdité DFNA50 (ClinGen modéré, 3 familles) | phénocopie | établi | Expression limitée à l'oreille interne : audiogramme, OXPHOS normal |
| MIR204 | Dystrophie rétinienne avec colobome (3 familles) | phénocopie | établi | Examen ophtalmologique ; OXPHOS normal |
| MIR184 | Kératocône sévère avec cataracte ; EDICT | aucun | établi | Ophtalmologie |
| MIR140 | Dysplasie spondylo-épiphysaire type Nishimura | aucun | établi | Spécifique du chondrocyte |
| MIR17HG | Syndrome de Feingold type 2 (délétion) | aucun | établi | CNV |
| MIR145 | Dysfonction multisystémique du muscle lisse (un cas) | phénocopie | cas isolés | Phénotypage digestif et vasculaire |
| MIR208A | Variants rares du pré-miRNA dans la cardiomyopathie dilatée (défaut de maturation in vitro) | phénocopie | cas isolés | Tissu cardiaque seulement |
| MIR182, MIR183 | Aucun variant humain ; surdité chez la souris | hypothèse | expérimental | OXPHOS normal |
| MIR181A1, MIR181B1, MIR181A2, MIR181B2 | Renouvellement mitochondrial dans la rétine (souris) | indirect | expérimental | Seul un **gain** d'expression aurait un effet |
| MIR181C | Réprime mt-CO1 (rat) | direct | expérimental | CIV, MT-CO1 |
| MIR1-1, MIR1-2 | Stimulent la traduction mitochondriale dans le muscle (souris) | direct | expérimental | Muscle seulement ; fibroblastes normaux |
| MIR378A | Réprime ATP6 ; oxydation des acides gras | direct | expérimental | La perte augmente la capacité oxydative |
| MIR21, MIR762, MIR2392 | mt-CYB (rat), ND2 (ischémie), transcription de l'ADNmt (tumeur) | direct | expérimental | — |
| MIR4485 | Séquence d'origine mitochondriale (ARNr 16S), incluse dans MTRNR2L8 | direct | expérimental | **Alerte NUMT** : un « variant » peut venir de l'ADNmt |
| MIR210, MIR338, MIR33A, MIR33B, MIR499A, MIR484, MIR494, MIR761, MIR1229 | Hypoxie (ISCU, COX10), OXPHOS axonale, β-oxydation, dynamique mitochondriale | indirect | expérimental | Une perte de fonction n'explique pas un déficit OXPHOS |
| MIR9-1, MIR9-2, MIR9-3 | Médiateur de la physiopathologie du MELAS (cybrides) | indirect | expérimental | Modulateur secondaire seulement |

### d.5 lncRNA (y compris empreinte et expansions)

| Gène | Maladie | Lien mito | Niveau | Test attendu |
|---|---|---|---|---|
| ATXN8OS | Ataxie spinocérébelleuse SCA8 (expansion CTA/CTG) | phénocopie | établi | OXPHOS normal ; ExpansionHunter puis confirmation |
| NUTM2B-AS1 | Myopathies oculopharyngodistale et oculopharyngée avec leucoencéphalopathie (expansion CGG) | phénocopie | établi | Biopsie : vacuoles bordées, inclusions p62+ ; FISH des foyers d'ARN CGG dans les fibroblastes |
| GNAS-AS1 | Pseudohypoparathyroïdie 1B dominante (calcifications des noyaux gris possibles) | phénocopie | établi | PTH, calcémie, méthylation du locus GNAS |
| H19 | Beckwith-Wiedemann, Silver-Russell (ICR1) | aucun | établi (locus) | MS-MLPA |
| KCNQ1OT1 | Beckwith-Wiedemann (ICR2) | aucun | établi (locus) | Méthylation |
| MEG3 | Temple, Kagami-Ogata (domaine 14q32) | aucun | établi (locus) | MS-MLPA 14q32, UPD, CNV |
| HYMAI | Diabète néonatal transitoire type 1 (locus PLAGL1/HYMAI) | aucun | établi (locus) | MS-MLPA, UPD6 |
| LINC01081 | Dysplasie alvéolo-capillaire (amplificateur de FOXF1) | aucun | cas isolés | — |
| FAST-1 (antisens de FXN) | Augmenté dans l'ataxie de Friedreich (2 lignées) | indirect | expérimental | **D'abord** génotyper l'expansion GAA de FXN |
| CEROX1 (HGNC:53928) | Augmente l'activité du CI (cellules) | indirect | expérimental | Seulement en cas de perte de fonction |
| TUG1 | PGC-1α (souris, podocytes) | indirect | expérimental | Hypothétique |
| LINC00467 (micropeptide ASAP) | ATP synthase (cellules tumorales) | direct | expérimental | BN-PAGE CV, hypothétique |
| GAS5 | Flux du cycle de Krebs (tumeur, une étude) | direct, fragile | expérimental | OXPHOS normal |
| SAMMSON | Traduction mitochondriale dans le mélanome ; non exprimé ailleurs | direct (tumoral) | expérimental | Sans objet en germinal |
| SLC16A1-AS1 | MCT1 (cellules tumorales) ; le crible Cas13 initial est **rétracté** | indirect | expérimental | Transport du lactate, hypothétique |

### d.6 Autres ARN Pol III, ARNr et télomérase

| Gène | Maladie | Lien mito | Niveau | Test attendu |
|---|---|---|---|---|
| RMRP | Hypoplasie cartilage-cheveux, dysplasie anauxétique, déficit immunitaire isolé | hypothèse (corrigé) | établi | ARN RMRP effondré (ARN total), maturation du pré-ARNr, cycline B2 |
| TERC | Télomeropathies : dyskératose congénitale, aplasie, fibrose pulmonaire | indirect | établi | Flow-FISH des télomères < 1er percentile |
| TRU-TCA1-1 (ARNt-Sec) | Déficit en sélénoprotéines « de stress » ; 2 familles homozygotes, la seconde peu symptomatique | phénocopie | probable | Sélénium, GPX3, SELENOP, bilan thyroïdien |
| RPPH1 | Aucune ; la RNase P mitochondriale est protéique | hypothèse | expérimental | Aucun test mitochondrial anormal attendu |
| RN7SK | Aucune ; partenaire LARP7 | hypothèse | aucun | — |
| RN7SL1 | Aucune ; partenaire SRP54 | indirect, fragile | expérimental | — |
| RN7SL2 | Aucune | hypothèse | aucun | — |
| VTRNA1-1 | Autophagie (p62) | hypothèse | expérimental | — |
| VTRNA2-1 (nc886) | Épiallèle métastable ; inhibe PKR | hypothèse | cas isolés (association) | — |
| VTRNA1-2, VTRNA1-3 | Aucune | aucun | expérimental, aucun | — |
| RNY1, RNY3, RNY4, RNY5 | Aucune | aucun | aucun | — |
| RNA45SN1 à RNA45SN5 | Aucune maladie par variant de l'ADNr | hypothèse | cas isolés | Zone aveugle ; nombre de copies en covariable |

### d.7 Gènes protéiques et micropeptides présents dans le catalogue confirmé

Ils sont retirés du catalogue d'ARN et proposés pour la liste « machinerie ARN et micropeptides » (décision D5).

| Gène | Maladie | Lien mito | Niveau | Remarque |
|---|---|---|---|---|
| LARP7 | Syndrome d'Alazami | indirect (souris) | établi | Tests mitochondriaux humains attendus normaux |
| MEPCE | Retard et épilepsie (un cas de novo) | phénocopie | cas isolés | 7SK diminué dans les fibroblastes |
| DIS3L2 | Syndrome de Perlman | aucun | établi | Témoin positif « machinerie » |
| INTS13 | Ciliopathie de type oro-facio-digital | aucun | probable | Sous-unité de l'Integrator |
| MTLN (ex-LINC00116) | Micropeptide mitochondrial (respiration et lipides) | direct | expérimental | Localisation controversée |

### d.8 Entrées contestées

| Gène | Raison |
|---|---|
| MT-HSP2 | Second promoteur du brin lourd, situé dans MT-TF près de MT-RNR1. Tout variant s'interprète d'abord comme variant de l'ARNt-Phe ou de l'ARNr 12S |
| HELLPAR | Syndrome HELLP d'origine placentaire, sans réplication. Le drapeau anamnestique pertinent est HADHA/HADHB (LCHAD), pas HELLPAR |
| HNF1A-AS1 | Aucun patient. La fonction est portée par le promoteur (HASTER) et passe par HNF1A (MODY3) |
| PINK1-AS | Une étude ; effet sur la fonction mitochondriale jamais testé |

### d.9 Entrées réfutées

- **TRR-TCT3-1** (ARNt nucléaire Arg-UCU), proposé comme orthologue humain de la famille murine n-Tr20 : réfuté par la vérification. L'exemple murin n-Tr20 correspond à tRNA-Arg-TCT-4-1, classé en priorité basse. Le motif détaillé est dans le TSV.

### d.10 Entrées non vérifiées ou sans résultat

- **Contexte non vérifié** (statut `contexte_non_verifie`) : environ 140 gènes protéiques ou éléments hors classe ARN. On y trouve la machinerie Pol III et TSEN, l'exosome, l'Integrator, le complexe SMN, les loci à expansion (FXN, RFC1, FMR1, DMPK, NOTCH2NLC, C9orf72…), les gènes de maintenance de l'ADNmt, les micropeptides, les ribonucléases et les gènes porteurs d'un intron U12. Ils servent de contexte, pas de catalogue.
- **Priorité basse non vérifiée** : RNA5S1, la famille des ARNt nucléaires, MTRNR2L8, MTRNR2L1, CRNDE, LIPCAR, XIST, FENDRR, MEG8, LINC00299, HOTAIRM1, LINC-PINT, RMST, MIR143, MIR133A1/A2, MIR206, MIR208B, MIR451A, MIR1974, MIR16-1, MIR218-2, MIR137, MIR127, MIR2861, SNORA55.
- **Soumis à vérification, sans résultat retourné** : CHASERR, SNHG14, CARMN, Maenli, MDRL, MALAT1, MIR23A, MIR29B1, MIR195, MIR15B, MIR141, MIR30A. Leur statut est dans le TSV.
  - Pour **CHASERR**, la référence primaire (Ganesh 2024, NEJM) est vérifiée. Il est gardé comme « piste sans variant » (délétion), pas comme gène confirmé.

---

## (e) La stratégie finale

### e.0 Ce qui a changé après la critique

**Problèmes majeurs, tous corrigés**
1. **Force de la preuve pour l'ADNmt.** Le « PS3 » attribué aux tests OXPHOS est supprimé. Cybrides et PCR en fibre unique valent PS3_supporting, dans les conditions de McCormick 2020 : déficit biochimique documenté dans la lignée du patient (moins de 20 % des témoins ou plus de 2 ET sous leur moyenne) et charge mutante d'au moins 60 % dans les cybrides. Les activités, le BN-PAGE et le Seahorse décrivent le phénotype (§ e.4).
2. **B2 contredisait N2.** J'ai retenu l'option (b) : une **voie « TND » fixée a priori** (§ e.3). RNU2-2 et RNU5B-1 restent au tri par cette voie, avec un motif visible dans le catalogue.
3. **Infrastructure de la route chrM.** Ajout de la décision bloquante **D0**. E1 devient conditionnelle, avec une solution de repli : un VCF chrM de la plateforme, lu par `--garder-format`.
4. **Règle de rendu.** D1 et D2 sont définis par la **classification ACMG/AMP en points** (Tavtigian 2020). M1 et M2 comptent pour un seul code PS3. Une table de cas types est ajoutée (§ e.5).

**Problèmes moyens, tous traités**
5. **Motif OUTRIDER (Jackson 2026).** Corrigé. OUTRIDER et le rapport U2-2/U2-1 sont utilisés ensemble. La distribution du rapport est vérifiée sur des témoins locaux. La sensibilité est signalée comme inconnue.
6. **Seuil de Brnich.** « 11 variants témoins (P/LP et B/LB) », distinct de la référence de cohorte (individus).
7. **Filtre récessif.** Il porte sur la **combinaison génotypique** observée chez les témoins (méthode de Jackson 2026). *Non retenu* : la borne chiffrée de 0,5 % proposée par le critique pour l'allèle hypomorphe. Elle n'est pas vérifiée ; la borne sera fixée en E2, à partir des fréquences publiées des allèles hypomorphes connus, et documentée.
8. **TRU-TCA1-1 passe en N2c.** *Motif différent de celui du critique.* La seconde famille existe (Geslot 2021, homozygote non apparenté), mais elle est peu symptomatique et les auteurs concluent que le variant n'est pas toujours causal. Le niveau gène-maladie reste « probable » dans le catalogue, mais le gène entre au tri comme candidat.
9. **Spécificité.** Ajout de témoins négatifs (génomes déjà résolus) en E5 et d'une validation de la route chrM sur des porteurs connus en E1.
10. **Calendrier et effectifs.** Effectifs indiqués ; diagnostics confirmés par l'ARN reportés à 6 à 9 mois et conditionnés à D4 et D6.
11. **Découvertes secondaires et consentement.** Ajout de la décision D12.

**Problèmes mineurs, tous corrigés**
- XIST retiré des gènes d'empreinte. Ajout d'une piste « biais d'inactivation de l'X ».
- Seuil −0,2 et règle « 50 % à 20× » attribués à Jackson 2025 (fenêtres de 10 pb). Ce point n'est plus « à relire ».
- Dawes 2026 reformulé (signature de la boucle T, sang et neurones iPSC).
- Seuils d'hétéroplasmie reformulés selon le mtCN ; la borne « 75 à 100 », non sourcée, est retirée.
- BS3_supporting limité aux classes qui ont un test M1 validé dans le tissu prélevé.
- Régions récessives de RNU2-2 alignées sur Leitão 2026 (contre-revue).
- Ordres de grandeur d'hypermutabilité séparés.

*Retenu en partie* : le critique propose de calculer le drapeau `hypermutable` « sur les données locales ». Je le calcule sur **gnomAD v4.1** (données indépendantes), et non sur la cohorte de patients, pour respecter la règle « aucun seuil tiré des patients ».

**Autres corrections issues de la vérification du catalogue**
- RNU5A-1 sort du tri (lien mito « aucun », cas isolés) et va dans la feuille « recherche ».
- MIR208A entre en N2c (cardiomyopathie, cas isolés).
- RMRP garde le niveau N2, mais son lien mitochondrial est « hypothèse ». Il n'entre au tri que sur une cytopénie ou un déficit immunitaire.
- SLC16A1-AS1 et SAMMSON passent en N4 (annotation seulement).
- FXN (expansion GAA) rejoint les pistes STR, à cause de FAST-1.

### e.1 Les quatre principes

1. **ADNmt en priorité 1.** C'est la seule classe d'ARN à mécanisme mitochondrial direct établi, et probablement la plus rentable. Le lot chrM démarre dès que D0 est tranchée.
2. **Pour un ARN nucléaire, le bilan OXPHOS ne prouve ni n'innocente le variant.** Il qualifie le phénotype : déficit primaire, déficit secondaire ou phénocopie. La preuve du variant vient de l'ARN lui-même (M1) et de sa signature en aval (M2). **Pour l'ADNmt**, la preuve vient surtout de la ségrégation hétéroplasmie-phénotype, du caractère de novo, de la PCR en fibre unique et des cybrides (Wong 2020 ; Yarham 2011 ; McCormick 2020). Le bilan OXPHOS seul reste descriptif.
3. **Une phénocopie confirmée se rend comme un diagnostic « non mitochondrial »** (ReNU, RNU2-2, RNU7-1/Aicardi-Goutières, SNORD118/LCC…).
4. **Aucun seuil n'est ajusté** sur AlphaGenome, SpliceAI, SpliceAI2 ou AlphaFold 3, ni sur des variants de patients. Les positifs publiés mesurent le rappel ; des génomes déjà résolus mesurent la spécificité.

### e.2 Les quatre conclusions possibles

| Conclusion | Définition | Rendu |
|---|---|---|
| **D1** : diagnostic, mécanisme mitochondrial direct | Variant de l'ADNmt (ARNt, ARNr, région de contrôle) classé **P ou LP** selon les spécifications ClinGen pour l'ADNmt (McCormick 2020) | Filière mitochondriale |
| **D2** : diagnostic, phénocopie | Variant d'un ARN nucléaire du catalogue (N2 ou N2c) classé **P ou LP** selon l'ACMG/AMP en points (§ e.5) | **« Maladie non mitochondriale »** ; bilan OXPHOS décrit ; orientation vers l'équipe concernée |
| **R** : recherche | Variant plausible ou récurrent, gène non établi, ou absence de marqueur validé | Feuille « recherche » séparée, sous protocole (D12) ; pas de rendu |
| **VUS** | Classe VUS | Réévaluation programmée |

La chaîne M0 à M4 (§ e.4) est une **liste de vérifications** qui alimente les codes. Elle n'est pas une règle de rendu.

### e.3 Niveaux du catalogue et règles d'entrée au tri

**Note A, lien gène-maladie**
- A1 : GenCC Definitive ou Strong, ou au moins 2 familles avec preuve fonctionnelle.
- A2 : GenCC Moderate, ou 2 familles.
- A3 : cas isolé.
- A4 : aucun.

**Note B, lien mitochondrial**
- B1 : mécanisme direct chez l'homme.
- B2 : phénocopie. Il faut au moins **1 terme HPO spécifique** publié pour le gène (acidose lactique, tableau de Leigh, leucoencéphalopathie ou hypomyélinisation, calcifications des noyaux gris, surdité neurosensorielle, rétinopathie pigmentaire, atrophie optique, cardiomyopathie, myopathie, ophtalmoplégie, ataxie, diabète précoce, cytopénie), **ou** des patients publiés d'abord comme suspects mitochondriaux.
- **B2-TND (nouvelle voie, fixée a priori)** : gène A1 dont le phénotype principal publié est un trouble du neurodéveloppement (TND) et qui figure parmi les causes non codantes les plus fréquentes des TND non résolus. Il n'entre au tri que si le patient a lui-même un TND (retard global, déficience intellectuelle ou épilepsie).
  - Gènes concernés aujourd'hui : RNU4-2 (qui satisfait aussi B2 par Nakamura 2025 et par la leucoencéphalopathie récessive), RNU2-2 et RNU5B-1.
  - La liste est fermée et révisée seulement à la curation trimestrielle.
- B3 : mécanisme indirect expérimental.
- B4 : aucun.

| Niveau | Règle | Usage | Gènes |
|---|---|---|---|
| **N1** | B1 | Route chrM ; classement par variant | 22 MT-T, MT-RNR1, MT-RNR2 (cas isolés au niveau du gène), LSP2 (préprint), éléments de la région de contrôle en annotation |
| **N2** | A1 ou A2, avec B2 ou B2-TND | **Critère de tri** | RNU4-2 (dominant et récessif), RNU2-2 (dominant et récessif ; B2-TND), RNU5B-1 (B2-TND), RNU4ATAC, RNU6ATAC, RNU12, RNU7-1, SNORD118, RMRP (cytopénie ou déficit immunitaire), TERC (cytopénie), RNU6-1/-2/-8/-9 et MIR96/MIR204 (`test=clinique`, au tri seulement si le HPO spécifique concorde) |
| **N2c** | A3 (ou A2 fragile), avec B2 | Tri, marqué « candidat » | TRU-TCA1-1, MIR145, MIR208A |
| **N3** | A1 ou A2, avec B4 | Annoté ; au tri seulement si un HPO spécifique concorde | MIR184, MIR140, SNORA31, MIR17HG (CNV) |
| **Pistes sans variant** | Routes SV, STR, méthylation ou inactivation de l'X, absentes de reannot | Liste de vérifications, hors tri SNV | CHASERR (délétion) ; SNORD116 et 15q11 (CNV, méthylation) ; 14q32/MEG3/SNORD114 ; H19, KCNQ1OT1, GNAS-AS1, HYMAI (méthylation) ; ATXN8OS, NUTM2B-AS1, FXN (STR) ; TERC en délétion 3q26 ; **biais d'inactivation de l'X** chez les femmes hétérozygotes pour un gène mitochondrial lié à l'X (AIFM1, PDHA1…) |
| **R** | Candidats sans lien établi | Feuille « recherche » | RNU5A-1, RNU5E-1, RNU5F-1, RNU4-1 (forme AR proposée) |
| **N4** | B1 ou B3, avec A3 ou A4 | Annotation descriptive, **jamais au tri** | mitomiR, snoRNA de RPL13A, SNORA73A, CEROX1, TUG1, GAS5, LINC00467, SAMMSON, SLC16A1-AS1, FAST-1, RPPH1, RN7SK, RN7SL1/2, VTRNA*, RNY*, lncND5/6/Cytb, MIR4485 (alerte NUMT) |
| **Exclu** | Contesté, réfuté, multicopie non résoluble, gène codant | — | HELLPAR, HNF1A-AS1, PINK1-AS, MT-HSP2 (interprété comme MT-TF ou MT-RNR1), TRR-TCT3-1 ; zones aveugles RNU1-1 à -4, RNU2-1, RNA45SN*, RNA5S*, RNVU1* ; RNU4-1 et RNU2-1 restent **dénominateurs** d'expression |

Le catalogue porte une colonne `motif_tri` (`B2`, `B2-TND`, `HPO_specifique_requis`, `candidat`) pour que la règle appliquée à chaque gène soit visible.

### e.4 Chaîne de preuve et force des tests

| Maillon | Question | Tests | Code ACMG/AMP |
|---|---|---|---|
| **M0** génétique | Rare, bien cartographié, dans une région critique, transmission compatible | Trio, phase, IGV, Sanger, contrôle des paralogues | PM2_supporting ; PS2 (**moderate** dans un gène `hypermutable`) ; PM3 (barème des combinaisons en trans) ; PM1 assimilé pour une région critique publiée |
| **M1** l'ARN | ARN diminué, mal maturé, déséquilibré par rapport au paralogue | Rapport U2-2/U2-1 et outlier OUTRIDER ; niveau de RNU4-2 ; RT-qPCR ou northern (U8, RMRP, U12, U4atac) ; précurseur prolongé en 3' | **M1 et M2 = un seul code PS3** (même mécanisme ; Brnich 2019) : supporting en dessous de 11 variants témoins validés, P/LP et B/LB confondus |
| **M2** signature en aval | La fonction est atteinte dans le transcriptome | A5SS ReNU (6 sites) ; outliers de rétention U12 dans les MIG ; histones lues au-delà de la tige-boucle ; pré-ARNr et cycline B2 ; épisignature ; score ISG | (même code que M1) |
| **M3** phénotype mitochondrial | Déficit OXPHOS, et de quel type | Activités CI-IV rapportées à la citrate synthase, BN-PAGE, western blot, oxygraphie ou Seahorse (galactose), traduction mitochondriale, FGF21/GDF15, protéomique | **ARN nucléaire** : aucun code. **ADNmt** : PP4 (supporting) si le déficit est cohérent avec le variant (CI, CIII, CIV pour un ARNt), mesuré dans le muscle (fibroblastes seulement si plusieurs probands), sous 20 % de la moyenne des témoins, causes nucléaires exclues (McCormick 2020) ; une activité normale n'est jamais un critère bénin ; cybrides et fibre unique : PS3_supporting |
| **M4** causalité | L'effet vient bien de ce variant | Ségrégation (PP1) ; PCR en fibre unique et cybrides (ADNmt : PS3_supporting, conditions de McCormick 2020) ; iPSC ; sauvetage ; SGE publiée | Selon le test |

**Règle d'arrêt (BS3_supporting)**
- Elle ne s'applique qu'aux classes qui ont un **test M1 validé dans le tissu prélevé**, avec la bonne banque **et** un témoin positif dans la même série.
- Exclusions :
  - variants récessifs hypomorphes ;
  - classes sans signature dans le sang (RNU2-2 dominant, RNU5B-1, RNU6-x) ;
  - score SGE normal de RNU4-2 pour un variant récessif.
- Sans témoin positif, un résultat négatif n'a aucune valeur.

### e.5 Du maillon aux points ACMG/AMP : cas types

Barème en points (Tavtigian 2020) : supporting = 1, moderate = 2, strong = 4, very strong = 8. Classes : VUS de 0 à 5, LP de 6 à 9, P à partir de 10. **Table indicative, à valider par le binôme de curation avant usage.**

| Scénario | Codes | Points | Classe |
|---|---|---|---|
| Dominant, point chaud publié récurrent (ex. RNU4-2 n.64_65insT), de novo confirmé | Classement publié repris ; s'y ajoutent PS2_moderate (2), PM2_supporting (1), PS3_supporting (épisignature ou A5SS, 1) | — | **P ou LP** (classement existant) |
| Dominant, variant inédit dans une région critique, de novo confirmé, signature M2 positive | PS2 selon le barème SVI (supporting 1 pour un TND isolé ; moderate 2 si phénotype compatible non spécifique) + PM1_supporting 1 (Ellingford 2022) + PM2_sup 1 + PS3_sup 1 | 4 à 5 | VUS → R ; LP seulement avec une preuve en plus (second proband de novo, PS3_moderate après validation sur ≥ 11 témoins, ou PM1 moderate justifié par une spécification publiée pour ce gène) |
| Même cas, sans signature M2 | PS2 (1 ou 2) + PM1_sup 1 + PM2_sup 1 | 3 à 4 | VUS → R |
| Récessif, deux variants inédits en trans, M1 positif (rapport U2-2/U2-1 abaissé) | PM2_sup 1 + PS3_sup 1 (+ PP4 1 seulement si HPO spécifique) | 2 à 3 | VUS → R |
| Récessif, un allèle publié P en trans d'un allèle inédit dans la région critique, M1 positif | PM3 2 + PM1_supporting 1 + PM2_sup 1 (+ PS3_sup 1 seulement si l'effet est attribué à l'allèle inédit : expression allèle-spécifique, ou baisse nettement supérieure à celle des porteurs hétérozygotes du seul allèle P) | 4 à 5 | VUS ; LP seulement avec PP1 (germain atteint) ou un second proband |
| Homozygote inédit, M1 positif, sans autre atteint | PM3_supporting 1 + PM1_sup 1 + PM2_sup 1 + PS3_sup 1 | 4 | VUS (LP seulement avec ségrégation, PP1 chez d'autres atteints) |
| ARNt mitochondrial, variant inédit hétéroplasmique, absent chez la mère (NGS sur plusieurs tissus ; maternité confirmée par un génome mitochondrial identique), MitoTIP et HmtVar concordants, fibre unique positive | PS2_moderate 2 (barème SVI, phénotype « compatible non spécifique », McCormick 2020) + PM2_sup 1 + PP3 1 + PS3_sup 1 (+ PP4 1 si déficit concordant dans le muscle) | 5 (6 avec PP4) | VUS ; LP seulement avec PP4 ou un second proband de novo |
| ARNt mitochondrial, variant inédit hérité, sans fibre unique ni cybrides | PM2_sup 1 + PP3 1 (+ PP1 si ségrégation hétéroplasmie-phénotype) | 2 à 3 | VUS |

**Conséquence assumée** : beaucoup de cas récessifs à deux allèles inédits resteront en VUS ou R tant que les tests ARN ne sont pas validés localement (11 variants témoins). C'est voulu.

### e.6 Étapes ordonnées

**Effectifs** : bio-informaticien (BI), biologiste (BIO), clinicien (CLIN). Effort en jours-personne (j).

| # | Étape | Livrable | Effort | Dépend de |
|---|---|---|---|---|
| **E0** | **Constats**, sans modifier les dépôts :<br>• gene_type GENCODE v49 de TERC, RMRP, RPPH1, RN7SK et RNU2-2 ;<br>• gènes N2 présents dans `petits_arn.bed` ;<br>• contenu des VCF livrés (chrM ? GT ? AF ? filtrés ?) ;<br>• **CRAM et moteur de conteneurs sur la VM (D0)** ;<br>• bibliothèques WGS **sans PCR** ou non (STR, régions riches en GC comme RMRP) ;<br>• ARN-seq existantes (poly(A) ou ARN total) ;<br>• PAXgene et fibroblastes en banque ;<br>• patients non résolus avec HPO ;<br>• patients déjà criblés pour les 200 snRNA dans PFMG2025 (Leitão 2026) ? | Note d'une page | 2 j BI + 1 j BIO | — |
| **E1** | **Lot chrM**, conditionnel à D0 : `R/vm/reannot_mt.sh` (ou repli sur le VCF chrM de la plateforme) ; paquet d'annotation ; `--garder-format` ; seuils figés. **Validation sur les CRAM de porteurs locaux connus** (m.3243A>G, m.8344A>G, délétion unique), pour mesurer la sensibilité selon l'hétéroplasmie et le mtCN, sans toucher aux seuils | Candidats N1 et liste des « cas expliqués » | 2 à 4 semaines BI | E0, D0, D1, D2 |
| **E2** | **Catalogue v1** (niveaux, `motif_tri`, colonnes de preuve) ; liste HPO figée ; positifs publiés ; relecture des seuils dans les articles (Jackson 2026, Dawes 2026, Arriaga 2025, Mendez 2026) ; bornes des allèles hypomorphes | `catalogue_arn.tsv`, `hpo_mito.tsv`, `positifs_arn.tsv` | 4 j BI + 2 j CLIN/BIO | E0 |
| **E3** | **Curta** : liste blanche, flancs, sous-commandes `catalogue`, `frequences`, `paralogues` ; drapeaux `hypermutable` (calculé sur gnomAD) et `appelable` ; clé ENSG ; versions figées ; MANIFEST ; `all` | Paquet `petits_arn/` v2 et univers figé | 5 j BI | E2 |
| **E4** | **VM** : colonnes `ARN_*`, critère de tri, `arn_patients.py`, `recurrence_arn.py`, tests, documentation | Branche avec tests PASS (aucun commit sans votre accord) | 5 j BI | E3 |
| **E5** | **Passage sur la cohorte** après exclusion des cas expliqués (y compris par l'ADNmt) ; revue IGV ; Sanger ; trio ; staff. **Passage témoin** sur des génomes déjà résolus par un autre diagnostic, pour compter les signalements par génome | Candidats avec le test indiqué ; pistes sans variant ; feuille « recherche » ; mesure de spécificité | 3 j BI + 1 j CLIN + laboratoire | E1, E4 |
| **E6** | **Route ARN** sur la machine retenue (D3) : DROP et rMATS ; `ratios_snrna.py`, `score_signatures_arn.py`, registre des preuves ; référence locale par tissu. En attendant : ARN-seq poly(A) existantes pour les signatures M2 d'épissage (A5SS, rétention U12), mais **pas** pour M1 | Pipeline validé | 10 à 15 j BI ; coût de séquençage | D3, D4, D6 |
| **E7** | **Confirmation** des candidats d'E5 (§ e.8). PAXgene : 2 à 3 semaines ; fibroblastes : 3 à 6 semaines | Dossier par patient, conclusion D1, D2, R ou VUS | 2 à 5 j par patient | E5, E6 |
| **E8** | **Étude de calibration (recherche, axe de publication)** : bilan OXPHOS standardisé chez des patients ReNU, RNU2-2, RNU4ATAC/RNU6ATAC confirmés, contre des témoins et des maladies mitochondriales connues ; enrichissement MitoCarta parmi les événements d'épissage | Article | 3 à 6 mois ; éthique et MTA | D6 |
| **E9** | **Phase 2** : kit fédéré proposé à PFMG2025, puis à d'autres cohortes ; GeneMatcher ; `nc_contexte.py` dans mitoloop_v0 | Comptes partenaires | 2 personnes-semaines ; accords de 3 à 6 mois | D7 |

**Calendrier**

| Scénario | Semaines 1 à 5 | Semaines 6 à 9 | Mois 3 à 9 |
|---|---|---|---|
| **1 BI** | E0, E2, E3, E4 | E1 (chrM), puis E5 | E6, puis E7 dès que la référence et les témoins positifs existent |
| **2 BI** | E0, E2 à E4 en parallèle de E1 | E5 (semaines 5 et 6), E6 | E7 |

- Premiers candidats (SNV) : semaines 6 à 9.
- Premiers diagnostics confirmés par l'ARN : **6 à 9 mois**, conditionnés à D4 (référence locale) et D6 (témoins positifs, MTA de 3 à 6 mois).
- Les diagnostics de la route chrM peuvent arriver plus tôt : les tests mitochondriaux sont déjà maîtrisés.

### e.7 Modifications proposées dans les dépôts (aucune n'est faite)

**Curta** : `R/mcia/petits_arn.py`, `R/mcia/run_paquet.sh`

1. **Liste blanche** `genes --liste-blanche catalogue_arn.tsv`, indépendante du biotype et de `--max-len`. Elle couvre TERC, RN7SK, RMRP, RPPH1 et RNU2-2 selon leur `gene_type` réel, les copies annotées « pseudogène » des familles fonctionnelles (liste des 200 snRNA de Leitão 2026) et l'ARNt-Sec (GTF des ARNt de GENCODE ou GtRNAdb haute confiance).
2. **Sous-commande `catalogue`** : validation (ENSG présent dans le GTF, coordonnées, symbole HGNC, DOI non vide) ; écriture de `catalogue_arn.bed.gz`, `catalogue_genes.tsv` et `petits_arn_regions.bed`, avec des flancs fixés a priori (promoteur −300 pb, extrémité 3' +50 pb pour les gènes à promoteur propre).
3. **Sous-commande `frequences`** : AC, AF et nhomalt de gnomAD v4.1 pour chaque position des petits ARN et de leurs flancs.
4. **Sous-commande `paralogues`** : alignement de chaque ARN avec ±50 nt ; copies à au moins 95 % d'identité ; colonne `cartographie`.
5. **Drapeaux dans `contrainte`** :
   - `hypermutable` : densité de SNV gnomAD dans le corps rapportée aux flancs, seuil fixé avant tout passage sur les patients ;
   - `appelable` : `non` si moins de 50 % des individus gnomAD sont à 20× ou plus, ou si la région recouvre les régions difficiles GIAB (règle de Jackson 2025) ;
   - mesure sans modèle de Jackson 2025 (proportion de SNV normalisée sur des fenêtres de 10 pb, seuil −0,2), à côté de l'O/E Roulette ;
   - indels comptés à part.
6. **Indexer `contrainte_genes.tsv` par ENSG.**
7. **`run_paquet.sh`** :
   - figer les versions de Rfam et de miRBase ;
   - ajouter `hp.obo`, `phenotype.hpoa`, MirGeneDB 3.0, GtRNAdb, la table GenCC, les listes de signatures publiées (6 sites ReNU, introns U12 et MIG, positions discriminantes entre paralogues, histones réplication-dépendantes) ;
   - écrire le MANIFEST **après** la copie de `petits_arn/` ;
   - ajouter `petits_arn` à `all`.

**VM**

- **`preparer_entree.py --garder-format`** : transmet GT, AD et DP, pour la logique récessive et pour l'hétéroplasmie.
- **`ajout_petits_arn.py`** : nouvelles colonnes.
  - `ARN_catalogue` : `GENE:N2:AR:phenocopie:motif`.
  - `ARN_region` : `point_chaud:…`, `region_critique:…`, `promoteur`, `extension_3p`, `corps` ou `hors_region`, avec comparaison exacte REF/ALT, indels compris.
  - `ARN_gnomad`.
  - `ARN_cartographie` : `unique`, `paralogue:n` ou `multicopie`, plus `appelable` et `hypermutable`.
  - `ARN_maladie`.
  - `ARN_test` : M1, M2, tissu, banque, et « sensibilité inconnue » pour le rapport U2-2/U2-1.
  - `ARN_alerte`.
- **Arguments ajoutés** :
  - `point_chaud_connu`, `region_critique_arn`, `frequence_compatible` ;
  - `combinaison_absente_temoins` : la combinaison génotypique (homozygote, ou même paire en trans) est absente des témoins ; chaque variant rare a une fréquence < 0,1 % ; un allèle hypomorphe connu est accepté en trans sous la borne du catalogue ;
  - `deux_variants_meme_ARN`, `homozygote_arn`.
- **Rattachement du gène par les coordonnées** (`petits_arn_regions.bed`), et non par la colonne `gene` : on évite l'OMIM du voisin.
- **`trier_variants.py`** : critère `PetitARN_effecteur`, distinct d'`ARN_signal`, avec les options `--arn-catalogue` et `--arn-niveaux N2,N2c`. Ordre : point chaud, puis région critique, puis concordance HPO.
- **`arn_patients.py` (nouveau)** : lit un fichier local hors dépôt (HPO, trio, tissus en banque, cas déjà expliqués). Il produit :
  - `ARN_HPO` (`specifique`, `TND`, `appoint`, `non`), argument d'ordre et jamais d'exclusion ;
  - `rapport_arn_patients.tsv` ;
  - `pistes_sans_variant.tsv`, avec la check-list : CNV et SV, MS-MLPA, STR, inactivation de l'X.
- **`recurrence_arn.py` (nouveau)** : régions touchées chez au moins 2 patients non apparentés, après exclusion des cas expliqués. 2 porteurs dans une unité `hypermutable` = piste ; 3 porteurs, ou 2 hors zone hypermutable = signal. Sortie dans une **feuille « recherche » séparée**, jamais dans le tri.
- **Route ARN-seq (scripts légers)** : `ratios_snrna.py` (lectures uniques), `score_signatures_arn.py` (lit les sorties DROP et rMATS produites hors de la VM), `ajout_preuves_arn.py` (registre local des preuves ; il ajoute les codes ACMG et la classe proposée).
- **Enchaînement** : étape 4b de `reannot_complet.sh`, `MINE` et `cp -r` dans `07_maj_mitoloop_vm.sh`, documentation.

**Route chrM** (`R/vm/reannot_mt.sh`, conditionnelle à D0)
- Réappel depuis le CRAM (MitoHPC, ou GATK Mutect2 en mode mitochondrial avec double alignement), Haplogrep 3, eKLIPse pour les délétions.
- Annotation préparée sur Curta, licences à vérifier : gnomAD v3.1 chrM, HelixMTdb, MitImpact 4, contrainte de Lake 2024, catalogue NUMT (Wei 2022), tables MITOMAP si la licence le permet.
- MitImpact ne couvre que les variants non synonymes des gènes codants. Pour les ARNt : MitoTIP, HmtVar, PON-mt-tRNA.
- **Seuils a priori** :

| mtCN | Règle |
|---|---|
| < 50 | Appels de faible VAF non interprétables |
| 50 à 500 | Variant inédit signalé dès 5 %, avec drapeau NUMT. Seuil de Laricchia 2022, établi sur des mélanges in silico ; gnomAD n'a publié que les variants à 10 % ou plus. Entre 5 et 10 %, confirmation sur un second tissu (urine) |
| > 500 | Suspicion de lignée cellulaire |
| Tous | Site P/LP connu signalé dès 1 % (Ratnaike 2025), avec correction d'âge pour m.3243A>G (Grady 2018) ; drapeaux homopolymères (m.302-316, m.16184-16193), NUMT et contamination |

- Sortie double : candidats N1, et liste des cas expliqués.

**mitoloop_v0** : aucun changement obligatoire. Ne pas ajouter d'ARN à `genes.txt`.

**Tests synthétiques** : `test_arn_catalogue.sh`, `test_recurrence_arn.sh`, `test_ratios_snrna.sh`, `test_signatures_arn.sh`, `test_reannot_mt.sh`. Il faut aussi combler les manques de `test_petits_arn.sh` (mode Roulette, indels, étape 4b).

**Régions critiques à recalculer sur Curta avec VariantValidator**
- RNU4-2 : n.62-79 en dominant ; tige II, k-turn et site Sm en récessif ; positif de référence NR_003137.2 n.64_65insT = chr12:120291839 T>TA.
- RNU2-2 : n.3-7, n.35 et n.38 en dominant. En récessif (Leitão 2026) : PM1 pour n.1-45 et pour le site Sm n.97-107 ; PM1_supporting pour SLIIa n.47-66, SLIIb n.68-84, SLIII n.112-144 et SLIV n.147-184 ; aucun filtre de position sur le second allèle (n.100T>C et n.107_118del sont récurrents ; délétions en trans). La région n.1-67 de Jackson 2026 n'est qu'une zone d'enrichissement statistique.
- RNU5B-1 : boucle I, n.37-44.
- RNU6-x : n.55_56insG (n.56_57insG non vérifié pour RNU6-8).
- RNU4ATAC (NR_023343.1, décalage de 1 nt avec Ensembl).
- RNU6ATAC : n.4-6, n.28-43, n.64 ; n.52 fait suspecter un artefact du paralogue.
- RNU12 : extension 3'.
- SNORD118 : promoteur et extension 3'.
- RMRP : promoteur ; n.70, n.71 et n.72A>G désignent le même variant.

### e.8 Confirmation fonctionnelle par classe

**Règles transversales**
- ARN-seq d'**ARN total déplété en ARNr** (et en globine pour le sang) pour M1 ; jamais de poly(A) ni de capture.
- Double prélèvement par défaut : PAXgene et fibroblastes.
- Référence locale par tissu et par protocole. La taille minimale (« au moins 25 ») est **à relire** dans la documentation de DROP.
- Témoin positif dans chaque série.

| Classe | M1 | M2 | M3 (rôle) | Tissu | Force maximale avant validation locale |
|---|---|---|---|---|---|
| ARNt et ARNr de l'ADNmt | Hétéroplasmie multi-tissus ; northern | — | Description ; cybrides et fibre unique = PS3_supporting (McCormick 2020) | Muscle, urine, sang, fibroblastes | Selon McCormick 2020 |
| RNU4-2 dominant | Facultatif | 6 sites A5SS (Dawes 2026, **préprint** : spécifique de la boucle T, validé dans le sang et les neurones dérivés d'iPSC, non évalué dans les fibroblastes ni le muscle ; la tige III demande d'autres événements) ; épisignature ; SGE | Recherche (E8) | Sang, lymphocytes | PS3_supporting |
| RNU4-2 récessif | Niveau de RNU4-2 diminué (Rius 2026) | rMATS, FRASER2 | E8 ; IRM de la substance blanche | Sang (ARN total) | PS3_supporting ; un SGE normal n'est pas un BS3 |
| RNU2-2 récessif | Outlier OUTRIDER **et** rapport U2-2/U2-1 (Jackson 2026 : 0,352 à 0,688 chez 8 patients déjà détectés par OUTRIDER, contre 0,76 à 1,45 chez 5 448 témoins de 100kGP). Sensibilité inconnue ; distribution vérifiée sur des témoins locaux, **sans réajuster le seuil** | Épissage dans les fibroblastes (Henry 2026) | E8 | Sang et fibroblastes | PS3_supporting |
| RNU2-2 dominant, RNU5B-1 | — | Lymphocytes cultivés ; effets discrets | E8 | Lymphocytes | Supporting au mieux ; un négatif n'a aucune valeur |
| RNU4ATAC, RNU6ATAC, RNU12 | RT-qPCR de U4atac ou U12 | Excès d'outliers de rétention dans les MIG (Arriaga 2025 ; seuil à relire) | MIG mitochondriaux (liste à relire) | Sang | PS3_supporting |
| RNU7-1 | RT-qPCR de U7 | Histones lues au-delà de la tige-boucle ; score ISG | Non | Fibroblastes en prolifération, sang | PS3_supporting |
| SNORD118 | U8, précurseur prolongé en 3' | Pré-ARNr ; IRM ; exclure CTC1 | Non | Fibroblastes, sang | PS3_supporting |
| RMRP | Niveau de RMRP | Pré-ARNr 5.8S, cycline B2 | Recherche | Fibroblastes, lymphocytes | PS3_supporting |
| TERC | Niveau de hTR | Flow-FISH | Non | Sang | Laboratoire de référence |
| TRU-TCA1-1 | ARNt-Sec | Sélénium, GPX3, SELENOP, thyroïde | Non | Plasma, fibroblastes | PS3_supporting |
| RNU6-x, MIR96, MIR204, MIR145, MIR208A | — | Phénotypage ciblé (ERG et OCT, audiogramme, digestif et vasculaire, cœur) et ségrégation | Non | — | PP4 et ségrégation |
| Pistes SV, STR, méthylation, inactivation de l'X | — | Confirmation du SV, MS-MLPA, RP-PCR ou lectures longues, nanopore résolu en haplotypes (Nakamura 2026) | Non | Sang | Selon la route |

**Tests transversaux à ajouter** (D4) : protéomique non ciblée (Hock 2025), FGF21 et GDF15, épisignature EPIC v2 pour ReNU.

### e.9 Interprétation des variants : règles propres à la classe

1. **Cadre** : ACMG/AMP (Richards 2015) adapté au non codant (Ellingford 2022), en points (Tavtigian 2020). Les recommandations snRNA de D'Souza 2026 (**préprint**) servent de guide seulement.
2. **Fréquences fixées a priori**
   - Dominant de novo : absent de gnomAD.
   - Récessif : chaque variant rare à moins de 0,1 % ; **combinaison génotypique absente des témoins** ; allèle hypomorphe connu accepté en trans sous la borne du catalogue.
   - Pas de pré-filtre CADD.
3. **Hypermutabilité.** Les chiffres publiés ne mesurent pas la même chose :
   - Roulette : taux environ 10 fois plus élevé dans les gènes transcrits par Pol III (Seplyarskiy 2023) ;
   - D'Souza 2026 (préprint) : taux de novo environ 50 fois plus élevé dans les snRNA que dans l'intergénique.

   Conséquences : le drapeau est calculé sur gnomAD, sans référence à ces chiffres ; PS2 passe en moderate ; l'O/E se lit avec réserve ; deux de novo dans un même snRNA font une piste, pas une preuve.
4. **Pièges**, rappelés dans `ARN_alerte` : paralogues (U6, U5, RNU6ATAC2P), ancien nom de pseudogène (RNU2-2P), décalages de notation (RNU4ATAC, RMRP), miRNA introniques classés « intronique » par VEP (MIR204 dans TRPM3, MIR140 dans WWP2), notation brin L de l'ADNmt (MT-TC, MT-OLR), MIR4485 et NUMT, positions de bord (m.14674 de MT-TE, m.7445 de MT-TS1).
5. **SpliceAI** : colonne informative seulement, jamais utilisée pour choisir un intron U12 ni pour fixer un seuil.

---

## (f) Risques et limites

1. **Rendement faible.** Une partie a peut-être déjà été trouvée par PFMG2025. Aucun ARN nucléaire à mécanisme mitochondrial primaire n'est établi chez l'homme : ne rien promettre.
2. **Infrastructure.** Sans CRAM ni conteneurs sur la VM (D0), la route chrM et les vérifications SV et STR ne sont pas exécutables. Le repli sur un VCF chrM de la plateforme perd les faibles hétéroplasmies si elle filtre.
3. **Données d'entrée.** Sans GT, il n'y a ni récessivité complète ni hétéroplasmie. Avec des listes filtrées, la récurrence est biaisée.
4. **Spécificité inconnue** tant que le passage sur des génomes déjà résolus n'est pas fait. Classe hypermutable et riche en paralogues : principales sources de faux signaux. Réponse : drapeaux, revue IGV, Sanger.
5. **Variants invisibles dans le sang.** Plusieurs ARNt mitochondriaux portent des variants restreints au muscle ou de novo somatiques (MT-TR, MT-TL2, MT-TM, MT-TY). Un WGS sanguin négatif n'exclut rien : ADNmt du muscle ou de l'urine.
6. **Signatures partiellement validées** : boucle T seulement (préprint) ; effets discrets pour RNU2-2 dominant et RNU5B-1 ; pas d'épissage visible dans le sang pour RNU2-2 récessif ; sensibilité du rapport U2-2/U2-1 inconnue.
7. **Force des preuves.** Tant que 11 variants témoins par classe ne sont pas réunis, beaucoup de cas récessifs restent en VUS. C'est attendu, et c'est honnête.
8. **Tissu et banque.** Une ARN-seq poly(A) ne sert pas à M1. Elle reste utilisable pour les signatures d'épissage.
9. **Attribution abusive d'un mécanisme mitochondrial.** Un OXPHOS anormal chez un porteur de variant d'ARN nucléaire peut être secondaire ou lié à la culture. D'où la règle « M3 sans code » et l'étude E8.
10. **Maintenance.** La littérature bouge vite (cinq articles Nat Genet sur les snRNA en 2026). Curation trimestrielle ; catalogue versionné.
11. **Catalogue incomplet.** Douze entrées sont revenues sans résultat de vérification. Les priorités basses et le contexte ne sont pas vérifiés. Plusieurs points restent « à relire » (liste en (h)).
12. **Licences.** OMIM en local seulement. Conditions de HPO, GenCC, MITOMAP et MitImpact à vérifier. reannot n'a pas de fichier LICENSE. Prestations externes sous contrat clinique. Aucune donnée patient dans une requête.

---

## (g) Décisions à prendre

- **D0 (bloquante).** CRAM disponibles sur la VM ? Moteur de conteneurs (Apptainer, par exemple) autorisé ? Sinon, extraction de chrM et appel dans un autre environnement sécurisé, ou repli sur un VCF chrM de la plateforme avec AF.
- **D1.** Lot chrM lancé dès D0 tranchée. Qui le porte ? Avec 1 ou 2 bio-informaticiens ?
- **D2.** Ajouter `--garder-format` (GT, AD, DP). Peut-on avoir des VCF non filtrés sur la VM ?
- **D3.** Lieu d'exécution de DROP, R et rMATS. Seules des tables reviennent à la VM.
- **D4.** Standard ARN total déplété pour toute nouvelle ARN-seq ; budget de la référence locale ; protéomique, FGF21/GDF15 et EPIC en prestation ou en interne.
- **D5.** Liste séparée « machinerie ARN et micropeptides » (LARP7, MEPCE, DIS3L2, INTS13, MTLN, POLR3A/B, EXOSC*, TSEN*…) dans l'ensemble de gènes, et sous quel drapeau.
- **D6.** Collaborations pour les témoins positifs et pour E8 (PFMG2025, cohortes ReNU et RNU2-2) ; éthique et MTA.
- **D7.** Fédération en phase 2 : avis du DPO, taille minimale de cellule, partage dans GeneMatcher.
- **D8.** Fichier HPO des patients sur la VM : format, responsable, fréquence de mise à jour.
- **D9.** Appliquer les seuils publiés sans réajustement local. Si la distribution locale des témoins sort de l'intervalle publié, on ne déplace pas le seuil : on déclare le protocole non comparable et on se limite à l'outlier OUTRIDER. Cela vaut jusqu'à ce que 11 variants témoins soient réunis.
- **D10.** Rendu des phénocopies comme diagnostics « non mitochondriaux » ; circuit vers la neurologie, l'ophtalmologie, l'immunologie, l'endocrinologie.
- **D11.** Binôme de curation bio-informatique et clinique ; rythme de mise à jour ; validation de la table des cas types (§ e.5).
- **D12.** Politique de rendu des découvertes pharmacogénétiques et secondaires : au minimum MT-RNR1 m.1555A>G et m.1494C>T (recommandation CPIC, McDermott 2022), plus TERC (aplasie, fibrose) et les anomalies d'empreinte. Le consentement existant couvre-t-il la feuille « recherche » (R), les VUS et le recontact ?

---

## (h) Références (DOI vérifiés ; préprints signalés)

**Cadre d'interprétation**
- Richards S 2015 Genet Med, doi:10.1038/gim.2015.30
- Tavtigian SV 2020 Hum Mutat, doi:10.1002/humu.24088
- Brnich SE 2019 Genome Med, doi:10.1186/s13073-019-0690-2
- Ellingford JM 2022 Genome Med, doi:10.1186/s13073-022-01073-3
- McCormick EM 2020 Hum Mutat, doi:10.1002/humu.24107
- Wong LJ 2020 Genet Med, doi:10.1038/s41436-019-0746-0
- Yarham JW 2011 Hum Mutat, doi:10.1002/humu.21575
- DiStefano MT 2022 Genet Med (GenCC), doi:10.1016/j.gim.2022.04.017
- Chong JX 2024 Genet Med, doi:10.1016/j.gim.2024.101199

**snRNA et cohortes**
- Chen Y 2024 Nature, doi:10.1038/s41586-024-07773-7
- Greene D 2024 Nat Med, doi:10.1038/s41591-024-03085-5
- Nava C 2025 Nat Genet, doi:10.1038/s41588-025-02184-4
- Jackson A 2025 Nat Genet, doi:10.1038/s41588-025-02209-y
- Greene D 2025 Nat Genet, doi:10.1038/s41588-025-02159-5
- Leitão E 2026 Nat Genet, doi:10.1038/s41588-026-02547-5
- Greene D 2026 Nat Genet, doi:10.1038/s41588-026-02539-5
- Jackson A 2026 Nat Genet, doi:10.1038/s41588-026-02551-9
- Rius R 2026 Nat Genet, doi:10.1038/s41588-026-02554-6 ; correction doi:10.1038/s41588-026-02636-5
- De Jonghe J 2026 Nature, doi:10.1038/s41586-026-10334-9
- Quinodoz M 2026 Nat Genet, doi:10.1038/s41588-025-02451-4
- Henry O 2026 Epilepsia, doi:10.1002/epi.70473
- Nakano Y 2026 iScience, doi:10.1016/j.isci.2026.116814
- D'Souza 2026 medRxiv (**préprint**), doi:10.64898/2026.08.03.26359558
- Dawes R 2026 medRxiv (**préprint**), doi:10.64898/2026.09.26.26364077

**Lien mitochondrial des snRNA, inactivation de l'X**
- Nakamura K 2025 J Hum Genet, doi:10.1038/s10038-025-01356-8
- De Biase I 2026 Clin Genet, doi:10.1111/cge.70272
- Nakamura K 2026 HGG Adv, doi:10.1016/j.xhgg.2026.100632

**Spliceosome mineur**
- He H 2011 Science, doi:10.1126/science.1200587
- Edery P 2011 Science, doi:10.1126/science.1202205
- Elsaid MF 2017 Ann Neurol, doi:10.1002/ana.24826
- Cologne A 2019 RNA, doi:10.1261/rna.071423.119
- Arriaga T 2025 AJHG, doi:10.1016/j.ajhg.2025.08.018
- Johnson MB 2026 AJHG, doi:10.1016/j.ajhg.2026.02.017
- Mendez R 2026 HGG Adv, doi:10.1016/j.xhgg.2026.100588
- Matalon DR 2026 Genet Med, doi:10.1016/j.gim.2026.102633
- Cuinat S 2026 Genet Med, doi:10.1016/j.gim.2026.102632

**Autres ARN mendéliens**
- Uggenti C 2020 Nat Genet, doi:10.1038/s41588-020-00737-3
- Rice GI 2013 Lancet Neurol, doi:10.1016/S1474-4422(13)70258-8
- Jenkinson EM 2016 Nat Genet, doi:10.1038/ng.3661
- McNeill N 2017 Neurol Genet, doi:10.1212/NXG.0000000000000162
- Ridanpää M 2001 Cell, doi:10.1016/s0092-8674(01)00205-7
- Thiel CT 2005 AJHG, doi:10.1086/497708
- Vulliamy T 2001 Nature, doi:10.1038/35096585
- Alter BP 2007 Blood, doi:10.1182/blood-2007-02-075598
- Schoenmakers E 2016 J Clin Invest, doi:10.1172/JCI84747
- Geslot A 2021 Eur Thyroid J, doi:10.1159/000518275
- Mencía A 2009 Nat Genet, doi:10.1038/ng.355
- Conte I 2015 PNAS, doi:10.1073/pnas.1401464112
- Hughes AE 2011 AJHG, doi:10.1016/j.ajhg.2011.09.014
- Grigelioniene G 2019 Nat Med, doi:10.1038/s41591-019-0353-2
- de Pontual L 2011 Nat Genet, doi:10.1038/ng.915
- Lino Cardenas CL 2023 J Clin Invest, doi:10.1172/JCI166497
- Ganesh VS 2024 N Engl J Med (CHASERR), doi:10.1056/NEJMoa2400718
- Ishiura H 2019 Nat Genet (NUTM2B-AS1), doi:10.1038/s41588-019-0458-z
- Koob MD 1999 Nat Genet (SCA8), doi:10.1038/7710
- Sahoo T 2008 Nat Genet (SNORD116), doi:10.1038/ng.158
- Lafaille FG 2019 Nat Med (SNORA31), doi:10.1038/s41591-019-0672-3
- Bastepe M 2005 Nat Genet (GNAS-AS1), doi:10.1038/ng1487
- Brioude F 2018 Nat Rev Endocrinol (Beckwith-Wiedemann), doi:10.1038/nrendo.2017.166
- Kagami M 2008 Nat Genet (14q32), doi:10.1038/ng.2007.56
- Alazami AM 2012 Hum Mutat (LARP7), doi:10.1002/humu.22175
- Astuti D 2012 Nat Genet (DIS3L2), doi:10.1038/ng.1071

**ARN à lien mitochondrial expérimental (sélection)**
- Zhang X 2014 Cell (miR-1), doi:10.1016/j.cell.2014.05.047
- Das S 2012 Circ Res (miR-181c), doi:10.1161/CIRCRESAHA.112.267732
- Indrieri A 2019 EMBO Mol Med (miR-181a/b), doi:10.15252/emmm.201708734
- Meseguer S 2015 Hum Mol Genet (miR-9, MELAS), doi:10.1093/hmg/ddu427
- Sirey TM 2019 eLife (CEROX1), doi:10.7554/eLife.45051
- Michel CI 2011 Cell Metab (snoRNA de RPL13A), doi:10.1016/j.cmet.2011.04.009
- Zhu X 2022 Cell (ARN 7S), doi:10.1016/j.cell.2022.05.006
- Rackham O 2011 RNA (lncND5/6/Cytb), doi:10.1261/rna.029405.111
- Holzmann J 2008 Cell (RNase P mitochondriale protéique), doi:10.1016/j.cell.2008.09.013
- Leucci E 2016 Nature (SAMMSON), doi:10.1038/nature17161
- Liang WW 2024 Cell (**rétracté**), doi:10.1016/j.cell.2024.10.021

**ADNmt**
- Lake NJ 2024 Nature, doi:10.1038/s41586-024-08048-x
- Laricchia KM 2022 Genome Res, doi:10.1101/gr.276013.121
- Battle SL 2022 NAR Genom Bioinform, doi:10.1093/nargab/lqac034
- Gupta 2023 Nature (variant chrM:302), doi:10.1038/s41586-023-06426-5
- Schönherr S 2023 NAR, doi:10.1093/nar/gkad284
- Giovannetti 2026 NAR (MitImpact 4), doi:10.1093/nar/gkag1070
- Goudenège D 2019 Genet Med, doi:10.1038/s41436-018-0350-8
- Wei W 2022 Nature, doi:10.1038/s41586-022-05288-7
- Grady JP 2018 EMBO Mol Med, doi:10.15252/emmm.201708262
- Stenton SL 2025 HGG Adv, doi:10.1016/j.xhgg.2025.100441
- Ratnaike T 2025 AJHG, doi:10.1016/j.ajhg.2025.04.003
- Schon KR 2021 BMJ, doi:10.1136/bmj-2021-066288
- Prezant TR 1993 Nat Genet (m.1555A>G), doi:10.1038/ng0793-289
- McDermott JH 2022 Clin Pharmacol Ther (CPIC, MT-RNR1), doi:10.1002/cpt.2309
- Goto Y 1990 Nature (m.3243A>G), doi:10.1038/348651a0
- Shoffner JM 1990 Cell (m.8344A>G), doi:10.1016/0092-8674(90)90059-n
- Blakely EL 2013 Hum Mutat, doi:10.1002/humu.22358

**Preuve fonctionnelle et ARN-seq**
- Kremer LS 2017 Nat Commun, doi:10.1038/ncomms15824
- Yépez VA 2021 Nat Protoc, doi:10.1038/s41596-020-00462-5
- Yépez VA 2022 Genome Med, doi:10.1186/s13073-022-01019-9
- Brechtmann F 2018 AJHG, doi:10.1016/j.ajhg.2018.10.025
- Mertes C 2021 Nat Commun, doi:10.1038/s41467-020-20573-7
- Scheller IF 2023 AJHG, doi:10.1016/j.ajhg.2023.10.014
- Mohammadi P 2019 Science, doi:10.1126/science.aay0256
- Hock DH 2025 Genome Med, doi:10.1186/s13073-025-01467-z
- Zech M 2025 Brain, doi:10.1093/brain/awaf059
- Levy MA 2022 HGG Adv, doi:10.1016/j.xhgg.2021.100075
- Lehtonen JM 2016 Neurology, doi:10.1212/WNL.0000000000003374
- Lehtonen JM 2021 J Inherit Metab Dis, doi:10.1002/jimd.12307
- Spinazzi M 2012 Nat Protoc, doi:10.1038/nprot.2012.058
- Wittig I 2006 Nat Protoc, doi:10.1038/nprot.2006.62
- King MP & Attardi G 1989 Science, doi:10.1126/science.2814477

**Catalogues et contrainte**
- Rath S 2021 NAR (MitoCarta3.0), doi:10.1093/nar/gkaa1011
- Martin AR 2019 Nat Genet (PanelApp), doi:10.1038/s41588-019-0528-2
- Gargano MA 2024 NAR (HPO), doi:10.1093/nar/gkad1005
- Chen S 2024 Nature (Gnocchi), doi:10.1038/s41586-023-06045-0
- Seplyarskiy V 2023 Nat Genet (Roulette), doi:10.1038/s41588-023-01562-0
- Clarke AW 2025 NAR (MirGeneDB 3.0), doi:10.1093/nar/gkae1094
- Chan PP & Lowe TM 2016 NAR (GtRNAdb), doi:10.1093/nar/gkv1309
- Dolzhenko E 2020 Genome Biol (ExpansionHunter Denovo) ; Dolzhenko E 2019 Bioinformatics (ExpansionHunter), doi:10.1093/bioinformatics/btz431, doi:10.1186/s13059-020-02017-z

**Points restant « à relire » avant de figer les seuils**
- Seuil d'excès d'outliers MIG : Arriaga 2025 (texte intégral non accessible) et Mendez 2026.
- Taille minimale de cohorte pour FRASER et OUTRIDER dans DROP (le chiffre « au moins 25 » n'est pas vérifié).
- Chiffre de 1,48 % de RNUopathies (Jackson 2026, texte intégral).
- Barème exact de McCormick 2020 pour le caractère de novo, la fibre unique et la ségrégation (table des cas types, § e.5).
- Niveau GenCC de TRU-TCA1-1 et de MIR145 ; nombre de familles RNU12 (Elsaid 2017, Xing 2021).
- Biotypes GENCODE v49 de TERC, RMRP, RPPH1, RN7SK et RNU2-2.
- Liste des MIG mitochondriaux : `/tmp/claude-0/-home-user/04d773a8-9dce-5a49-9259-4d11be1339d9/scratchpad/ncrna/spliceosome_mineur/introns_U12_genes_mito_hg38.tsv` (symboles HGNC à mettre à jour, par exemple GARS devenu GARS1).

**Fichiers des dépôts concernés** (lus seulement, aucun modifié) :
- `/home/user/reannot/mcia/petits_arn.py`
- `/home/user/reannot/mcia/run_paquet.sh`
- `/home/user/reannot/vm/ajout_petits_arn.py`
- `/home/user/reannot/vm/trier_variants.py`
- `/home/user/reannot/vm/preparer_entree.py`
- `/home/user/reannot/vm/reannot_complet.sh`
- `/home/user/reannot/vm/07_maj_mitoloop_vm.sh`
- `/home/user/reannot/tests/test_petits_arn.sh`

Instantané HGNC local : `/tmp/claude-0/-home-user/04d773a8-9dce-5a49-9259-4d11be1339d9/scratchpad/ncrna/ressources_dl_hgnc/hgnc_complete_set.txt`.
