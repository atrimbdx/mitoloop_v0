# État du chantier « non codant effecteur et phénotypes mitochondriaux » — 24 novembre 2026

Document de passation, écrit pour être lu par Claude dans une autre conversation (audit croisé avec
d'autres chantiers). Il résume ce qui a été fait, où sont les fichiers, ce qui a tourné, et ce qui reste
ouvert. Les faits ci-dessous ont été vérifiés dans les dépôts et les logs Curta transmis par l'utilisateur,
sauf mention « non vérifié ».

## 1. But et périmètre (décisions de l'utilisateur)

- Un patient (WGS ; bientôt VCF complet pré-annoté, déjà filtré par fréquence gnomAD et ségrégation
  familiale) porte un variant dans un **élément non codant effecteur** (ARN non traduit ou micropeptide) :
  reannot doit le **signaler** (colonne d'alerte).
- Condition d'alerte : (1) la voie biologique de l'élément passe par la mitochondrie — règle unique : une
  cible de l'élément est un gène mito (MitoCarta3.0 ∪ PanelApp « Mitochondrial disorders »), directement ou
  en cascade — **et** (2) un effet est prédit d'après la biologie de l'ARN (région fonctionnelle touchée).
- **Pas un classifieur** : ni sensibilité ni spécificité à optimiser ; l'utilisateur approfondit lui-même.
- Hors périmètre : ADNmt (déjà exploré) ; régulateur, eRNA, PROMPT (couverts par mitoloop) ; campagne
  AlphaGenome (AG2) par délétions abandonnée pour l'instant ; ARN-seq et tests OXPHOS.
- Inventaire exhaustif des classes (47), recherche fondamentale comprise ; priorités utilisateur : ARN de
  voûte (+++), piARN (+++, biogenèse à la surface mitochondriale), sosies nucléaires d'ARNt mito (NUMT,
  hypothèse novatrice). Les ARN non codants pathologiques connus servent de témoins positifs.

### Contraintes permanentes
- Aucun entraînement ni ajustement de seuil sur des sorties d'AlphaGenome, SpliceAI, SpliceAI2 ou
  AlphaFold 3 ; seuils fixés a priori dans le code.
- Aucune donnée patient vers un service externe ; OMIM jamais redistribué.
- Commandes données à l'utilisateur : une par ligne, sans « \ » de continuation.
- Ne lancer aucun agent ni workflow sans demande explicite. Effort moyen.
- Fichiers bruts d'éditeurs (données supplémentaires d'articles) dans `reannot/donnees/` : **usage interne
  seulement** ; avant toute diffusion hors de l'équipe, vérifier la licence ou les retirer du dépôt **et de
  son historique** (consigne dans `reannot/donnees/LISEZMOI_editeurs.md`, `README.md`, `SOURCES.md`).
- Coordination avec la conversation AG2 : elle utilise `~/reannot_git` (branche `main`) sur Curta — ne
  jamais y faire checkout/pull/reset. Ce chantier utilise `~/reannot_ncrna`, branche `claude/ncrna-effecteur`.
  Pas de commit sur `main`, pas de force-push (fusion de `main` dans la branche à la place), pas de PR sans
  demande.

## 2. Dépôts, branches, derniers commits

| Dépôt | Branche | Dernier commit | Contenu du chantier |
|---|---|---|---|
| `atrimbdx/reannot` | `claude/ncrna-effecteur` (36 commits devant `origin/main`, `main` fusionné le 23/11) | `ab38f53` | intégration dans reannot (VM hors ligne) |
| `atrimbdx/mitoloop_v0` | `claude/eloquent-babbage-yz27n1` | `84b4da3` | `docs/ncrna_effecteur/` seulement |

Décision **en attente** : PR `claude/ncrna-effecteur` → `main` (recommandée, pour que la VM et la
conversation AG2 partent de la même version) ou paquet VM préparé depuis `~/reannot_ncrna`.

## 3. Architecture

### mitoloop_v0/docs/ncrna_effecteur/
- `README.md` (contraintes), `brief_ncrna.md`, `strategie_ncrna.md`, `chemins_mito.md` (méthode des 3 chemins).
- `inventaire/inventaire_ncrna_effecteurs.tsv` + `inventaire.md` : 47 classes.
- `catalogue_ncrna_mito.tsv` : témoins.
- `cibles/cibles_ncrna.tsv` (22 994 lignes, 11 colonnes dont `region_guide`), `convertir_cibles.py`
  (sous-commandes snodb, curees, pol3, pseudogenes, humanin, lncrna_cis, mirtarbase, generique, fusion),
  `cibles_sources.md`.
- `chemins/chemin1/` (cibles et épissage ; cascade snoRNA/scaRNA → snRNA mineurs → 19-21 gènes mito à intron
  U12), `chemins/chemin2/` (localisation mitochondriale des ARN), `chemins/chemin3/` (micropeptides :
  `scripts/preparer_sorf.py`, `sorf_lib.py`, `annoter_chemin3.py`, `micropeptides_mito_etablis.tsv`).

### reannot (branche claude/ncrna-effecteur)
- `mcia/elements_nc.py` : liste unifiée des éléments (GTF de VEP Ensembl 116, jointure par ENSG/HGNC, jamais
  par symbole), chemins 1-2, paquet. Colonnes VM : `NC_element`, `NC_cibles_mito`, `NC_chemin1/2/3`,
  `NC_alerte`, `NC_motif`. Copies de snRNA « fonctionnelles » d'après Leitão 2026
  (`donnees/snrna_fonctionnels_leitao2026.tsv`, 247 copies) + cCRE ENCODE ; mirtrons calculés (intron du GTF
  à ± 2 nt).
- `mcia/sources_nc.py` : convertit les fichiers d'éditeurs et calcule des régions → BED hg38 dans
  `NC_BED_DIR`. Sous-commandes : pirna (Özata 2020), numt (Telonis 2014, sosies d'ARNt mito), agotrons
  (Hansen 2016, hg19), sdrna (Ender 2008), cirna (Zhang 2013, hg19), svrna (ENA FM207500-505), guides
  (régions guides des snoRNA depuis snoDB + ARNr RefSeq), trf (fragments d'ARNt depuis GtRNAdb), ysrna
  (RNY1/3/4/5), terra (10 kb subtélomériques), rdna (2 kb en amont de RNA45SN*), micropeptides (PIGBOS,
  SLC35A4-MP, altMiD51, MP31), spa (Wu 2016 Table S1), tout. Rapport `sources_nc.rapport.tsv` (une
  sous-commande seule écrit `sources_nc.<cmd>.rapport.tsv` depuis `e330557`).
- `mcia/run_paquet.sh` : étapes `nc_sources` (téléchargements HGNC, chaîne hg19→hg38, GtRNAdb, svRNA ENA,
  ARNr RefSeq, miRBase hairpin ; puis sources_nc tout ; puis miRTarBase → fusion des cibles ; puis guides),
  `nc_effecteur` (elements → chemin1 → chemin2 → paquet), `paquet`.
- `donnees/nc_types.tsv` (47 classes, sources et règles), `nc_temoins.tsv`, fichiers d'éditeurs.
- Tests : `tests/run_all.sh` (19 tests, tous passent au 24/11), dont `test_sources_nc.sh`,
  `test_nc_effecteur.sh`, `test_mirtrons.sh` (jeux synthétiques).

## 4. Sources : décisions et état

| Source | Statut |
|---|---|
| snoDB v2.0 | utilisée (remplace snoRNABase, non téléchargeable) ; régions guides calculées : 277 trouvées / 686 paires |
| RMBase | non intégrée (redondante avec snoDB : 295/302 paires identiques) |
| MirGeneDB, miRBase GFF | utilisées (MirGeneDB téléchargé automatiquement) |
| miRTarBase 11.0 (1,6 Go, Curta) | intégrée : « Functional MTI » seulement par défaut (`NC_MIRTARBASE_WEAK=oui` pour ajouter les Weak), paires dédoublonnées, règle HGNC MIR1-1/MIR181A1 ; t_mirtarbase.tsv 2,8 Mo |
| MINTbase/MINTmap | remplacés par des régions calculées sur GtRNAdb |
| piRBase | déconseillée (FASTA sans coordonnées) ; Özata 2020 utilisé (266 précurseurs ; 1 ligne pi-TTN.4 avec séparateur fautif corrigée) |
| circAtlas 3.0 | utilisé (`donnees/circAtlas.zip`) |
| `circ_mito.bed` | abandon proposé (SCAR viendrait de l'ADNmt — non vérifié ; Samd4/BRAKE murins) |
| GENCODE Ribo-seq ORFs phase 2 « primary » | utilisé (10 127 ORF, `donnees/`) |
| OpenProt 2.2 min_2_pep | utilisé sur Curta (protéines alternatives IP_ seules ; gène hôte par chevauchement avec le GTF) |
| Wu 2016 Table S1 (SPA) | utilisée : 9 lignes (SPA1, SPA2, 7 sno-lncRNA), hg19 → hg38 |
| ADNr (pRNA, PAPAS) | aucun RNA45SN* dans le GTF Ensembl 116 : BED vide ; PAPAS seulement par nom GTF |
| ARN de voûte | littérature : liens mito seulement indirects (vtRNA1-1/p62, nc886/PKR, MVP/mitophagie) ; proposition en attente : relais SQSTM1, EIF2AK2, PRKN, MVP comme cibles |

## 5. Ce qui a tourné sur Curta (23-24/11)

Environnement `/scratch/$USER/envs/ncrna`. Variables : `NC_BED_DIR=/scratch/$USER/nc_bed`,
`NC_SNODB`, `NC_MIRTARBASE`, `NC_MIRBASE`, `NC_RFAM`, `NC_CCRE`, `SORF_TABLE`, `CHEMIN1_PAQUET`,
`CHEMIN2_DIR`, `REANNOT_BASE=/scratch/$USER/reannot`, `MITOLOOP_V0=$HOME/mitoloop_v0`.

- Run du 23/11 : `nc_eff` 21509874 en échec (ligne BED12 malformée des piARN) → corrigé (`753d60a`).
- Run du 24/11 (commits reannot `753d60a`, mitoloop `8ecef15`) :
  - `nc_src` 21527051 COMPLETED (nc_sources, nc_mirtarbase, nc_cibles_fusion, nc_guides : status 0).
    piARN 266 ; micropeptides : PIGBOS (CDS de PIGBOS1 reprise), MP31 (ORF unique), SLC35A4-MP 2 candidats
    (non écrit), altMiD51 0 candidat ; SPA : 9 convertis, snoRNA en 5' trouvé pour 8/9 (SPA2 → SNORD109A,
    sno-lncRNA → SNORD116-x), **SPA1 : aucun** ; RNY5 absent du GTF.
  - `chemin3` 21527052 COMPLETED : 571 568 sORF lus, 123 987 après dédoublonnage ; appariés : MTLN, BRAWNIN
    (UQCC6), SMIM4 (UQCC5), MOCCI, SMIM26, PIGBOS, ASAP, STMP1 ; non appariés : SLC35A4-MP, altMiD51, MP31 ;
    niveaux mito fort 5, modéré 3, hypothèse 36 766.
  - `nc_eff` 21527053 : **en cours** au moment de la rédaction (~2 h 30 attendues). Résultat à vérifier :
    `grep -E "status|Traceback|Error" /scratch/$USER/nc_eff_21527053.log`.
- Dernier `nc_effecteur` complet antérieur (20/11) : 130/132 témoins retrouvés ; 2 057 éléments à cible
  mito directe, 45 en cascade.

## 6. Corrections faites après le run du 24/11 (pas encore exécutées sur Curta)

- `e330557` : rapport de sources_nc non écrasé par `guides` ; altMiD51 cherché en uORF **ou** altORF.
- `ab38f53` : micropeptides — candidats partageant ATG et stop (même ORF, épissages différents) tous écrits
  (cas probable de SLC35A4-MP) ; coordonnées des candidats au rapport sinon. Non testé sur données réelles.
- mitoloop `84b4da3` : SLC35A4-MP apparié sur sa longueur publiée (103 aa) dans le chemin 3.

## 7. Points ouverts (à auditer)

1. Résultat de `nc_eff` 21527053 (statut, nombre de témoins retrouvés, éléments à cible mito).
2. altMiD51 introuvable (70 aa) : vérifier les ORF de MIEF1 dans la table du chemin 3 (commande donnée à
   l'utilisateur) ; position/longueur/codon d'initiation à confirmer dans Samandi 2017.
3. SPA1 sans snoRNA en 5' (attendu SNORD107, de mémoire) : vérifier le voisinage hg38 dans le GTF.
4. MP31 : longueur 31 aa déduite du nom (texte intégral non accessible) — à confirmer.
5. Table S1 de Wu 2016 : convention de coordonnées mixte (3/9 lignes 1-based inclusives), écart ≤ 1 nt, laissé tel quel.
6. Relancement complet après les corrections du §6, puis étape `paquet` et transfert VM.
7. Décision PR vers `main` (§2).
8. Proposition non demandée en attente : relais des ARN de voûte ; tige-boucle apicale de nc886.
9. Diffusion : fichiers d'éditeurs à vérifier ou retirer de l'historique avant toute diffusion externe.

## 8. Références vérifiées (PubMed) pour les micropeptides et SPA

- PIGBOS, 54 aa, membrane externe mito : Chu 2019, doi:10.1038/s41467-019-12816-z
- SLC35A4-MP, 103 aa, membrane interne : Rocha 2024, doi:10.1016/j.jmb.2024.168559 ; renommé STREMI (MICOS) : Guo 2026, doi:10.1038/s44319-026-00783-8
- altMiD51, 70 aa : Samandi 2017, doi:10.7554/eLife.27860
- MP31, uORF de PTEN : Huang 2021, doi:10.1016/j.cmet.2020.12.008
- SPA de 15q11 : Wu 2016, doi:10.1016/j.molcel.2016.10.007
