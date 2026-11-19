# Chemin 2 : l'ARN non codant nucléaire est présent dans la mitochondrie

*Version 0, 19 novembre 2026. Première implémentation de la section « Chemin 2 » de
`../../chemins_mito.md`. C'est un **chemin d'hypothèse** : le niveau rendu ne dépasse jamais
`faible` et n'apporte **aucun point ACMG**.*

## 1. Question

Un variant rare de patient est situé dans un ARN non codant nucléaire. Deux conditions doivent être
remplies pour que le chemin 2 s'applique :

1. l'ARN est rapporté dans la mitochondrie ;
2. le variant touche un élément nécessaire : signal d'import, site de liaison d'un facteur ou repliement.

L'import d'ARN nucléaires dans la mitochondrie des mammifères est **débattu**. Les témoins les plus
solides (RMRP, RPPH1, ARNr 5S) ont tous une contre-preuve publiée ou une fonction intramitochondriale
réfutée. La seule technique de proximité dans la matrice (APEX) n'a trouvé **aucun** ARN nucléaire
enrichi (voir 3.3).

## 2. Fichiers

| Fichier | Contenu |
|---|---|
| `localisation_mito_ncrna.tsv` | Table de localisation, avec une ligne par couple (ARN, publication) : 30 lignes, dont 6 contre-preuves ou résultats APEX négatifs |
| `donnees/elements_mito.tsv` | Éléments d'import et de liaison en coordonnées hg38 : publiés, dérivés ou non cartographiés |
| `donnees/arn_coordonnees_hg38.tsv` | Coordonnées GENCODE v27 des ARN de la table (gène et union des exons) |
| `donnees/ref_regions_hg38.fa` | Séquences de référence hg38 (brin +) des régions à annoter : gène ± 100 nt, ou exons ± 150 nt pour SAMMSON |
| `donnees/genes_chemin2.txt` | Liste des symboles utilisés pour construire le paquet |
| `scripts/preparer_references.py` | Préparation sur une machine connectée : extraction des coordonnées depuis un GTF et des séquences depuis un FASTA local ou une URL indexée |
| `scripts/curta_preparer_chemin2.sh` | Préparation complète sur Curta : GENCODE v49, Rfam, ClinVar, roue ViennaRNA |
| `scripts/annoter_chemin2.py` | Annotation hors ligne sur la VM, qui produit la colonne `NC_chemin2` |
| `tests/temoins_chemin2.tsv`, `tests/test_chemin2.sh`, `tests/sortie_temoins_chemin2.tsv` | Test sur 11 témoins : 10 sont **construits à la main** et étiquetés `CONSTRUIT_A_LA_MAIN`, 1 vient du catalogue |

### Provenance des données du paquet

L'EBI, Ensembl, l'UCSC et NCBI eutils sont inaccessibles depuis la session de construction. Le paquet
a donc été construit à partir du bucket public du Broad (`gs://gcp-public-data--broad-references/hg38/v0/`) :

- coordonnées : `gencode.v27.primary_assembly.annotation.gtf` ;
- séquences : `Homo_sapiens_assembly38.fasta`, lu par requêtes HTTP Range guidées par son `.fai`.

Les séquences contrôlées à l'œil commencent comme attendu : RMRP par GGUUCG…, RPPH1 par GGCGGAGGG…,
ARNr 5S par GUCUACGGCC… et TERC par GGGUUGCGG…

À refaire sur Curta avec GENCODE v49 (`curta_preparer_chemin2.sh`). Deux points d'annotation sont à
surveiller :

- GENCODE v27 contient deux gènes RMRP (lincRNA de 268 nt et ribozyme de 264 nt) et deux gènes MALAT1.
  Le script garde le plus long.
- Le RPPH1 de v27 (ribozyme) fait 333 nt, contre 340 à 341 nt pour l'ARN H1 publié (voir 4).

## 3. Méthode

### 3.1 Construction de la table (revue bibliographique)

Chaque référence a été retrouvée dans PubMed. Le DOI vient de la notice PubMed, sauf Zhao 2019, qui n'a
pas de DOI dans PubMed et est cité par son PMID. Quand le texte intégral était accessible dans PMC, il a
été lu : Wang 2010, Jeandard 2019, Fazal 2019, Kaewsapsak 2017, Mercer 2011, Noh 2016 et Silver 2024.
Pour les autres, seul le résumé a été lu. Toute méthode non vérifiée en texte intégral est signalée dans
la ligne (« non vérifié »).

Candidats examinés :

- RMRP, RPPH1, ARNr 5S, TERC ;
- mitomiR : miR-1, miR-181c, miR-2392, miARN de Barrey 2011, jeux de Bandiera 2011 et de Silver 2024 ;
- lncRNA : MALAT1, SAMMSON, jeu de Taylor 2024 ;
- ARNt nucléaires (Mercer 2011) ;
- données APEX de la matrice (Fazal 2019, Kaewsapsak 2017).

Revues de synthèse utilisées pour les contestations :

- Jeandard 2019, *Cells* 8:286, doi:10.3390/cells8030286 ;
- Kim 2017, *BMB Rep* 50:164, doi:10.5483/bmbrep.2017.50.4.013 ;
- Huang 2021, *Front Cell Dev Biol* 9:688523, doi:10.3389/fcell.2021.688523 ;
- Wang 2012, *Biochim Biophys Acta* 1819:998, doi:10.1016/j.bbagrm.2011.10.001.

Contexte génétique humain du système d'import :

- PNPT1 c.1160A>G : baisse de l'import de l'ARNr 5S et d'ARN liés à MRP dans les fibroblastes, avec
  déficit combiné de la chaîne respiratoire (Vedrenne 2012, doi:10.1016/j.ajhg.2012.09.001) ;
- PNPT1 et surdité (von Ameln 2012, doi:10.1016/j.ajhg.2012.09.002).

### 3.2 Niveaux de preuve (colonne `type_preuve` et `rang_preuve`)

| Rang | Type | Sens |
|---|---|---|
| 1 | `apex_matrice` | Marquage de proximité dans la matrice (APEX-seq, APEX-RIP) |
| 2 | `mito_purifiees_rnase` | Mitochondries ou mitoplastes purifiés et traités par une nucléase (RNase A, nucléase micrococcale) |
| 3 | `imagerie` | Hybridation in situ ou imagerie avec co-localisation mitochondriale |
| 4 | `fractionnement_seul` | Fraction mitochondriale sans contrôle nucléase vérifié |
| 5 | `import_in_vitro` | Aptitude à l'import dans des mitochondries isolées, ce qui n'est pas une localisation in vivo |
| — | `contre_preuve`, `apex_matrice_negatif` | Résultat contraire ou absence d'enrichissement |

### 3.3 Résultat de la revue (table : 12 entités nommées, 3 jeux non nominatifs)

**Meilleure preuve par ARN :**

- Aucun ARN nucléaire n'atteint le rang 1. Les deux études APEX de la matrice (HEK293T) sont négatives
  pour les ARN nucléaires. Ce résultat vaut seulement pour les ARN longs : les petits ARN sont
  sous-échantillonnés, et APEX-RIP emploie les lncRNA nucléaires comme faux positifs de son seuil.
- **Rang 2, chez l'homme (4 entités) :**
  - RMRP (Puranam 2001 ; Noh 2016 ; souris : Chang 1987) ;
  - RPPH1 (Puranam 2001 ; Wang 2010) ;
  - ARNr 5S, RNA5S1 à RNA5S17 (Magalhães 1998, chez l'homme et le rat) ;
  - TERC (Cheng 2018).
- Rang 2, rat seulement : un jeu de plus de 200 miARN dans le muscle (Silver 2024).
- **Rang 3, imagerie (5 entités) :**
  - pré-miR-302a, pré-let-7b et miR-365 (Barrey 2011, muscle humain) ;
  - MALAT1 (hépatocarcinome) ;
  - SAMMSON (mélanome).
- **Rang 4, fractionnement seul (3 entités et 2 jeux) :**
  - miR-1 (souris) ;
  - miR-181c (rat) ;
  - miR-2392 (tumeur humaine) ;
  - les 57 miARN de HeLa (Bandiera 2011) ;
  - les ARNt nucléaires (Mercer 2011).

**Contestations notées dans la table :**

- RMRP : Kiss et Filipowicz 1992 (absence dans la mitochondrie). Seule la moitié 3' serait importée,
  et le rôle d'amorçage de l'ADNmt est remis en cause (Jeandard 2019).
- RPPH1 : la RNase P mitochondriale est protéique (Holzmann 2008) ; voir aussi la lettre de Rossmanith
  et Potuschak 2001.
- ARNr 5S : absent du mitoribosome (Rorbach 2016). Sa présence n'est pas contestée, sa fonction l'est.
- TERC : la maturation dans la matrice n'est pas démontrée. Supprimer la région « tige-boucle »
  **augmente** l'import.
- miARN : 7 miARN seulement sont retrouvés dans au moins 3 études chez l'homme (Jeandard 2019).
- MALAT1 et SAMMSON : seulement dans des cellules tumorales, et la négativité APEX s'y oppose.

### 3.4 Éléments touchés (`donnees/elements_mito.tsv`)

**Publié avec coordonnées (`publie_plage`, poids 2) : RPPH1 seulement.**

- Le signal d'import dépendant de PNPASE se situe entre les nt 103 et 140 de l'ARN H1 de 340 nt
  (cartographie par délétions, Wang 2010). Une tige-boucle de 20 nt de cette région suffit à faire
  importer un ARN non importé, et elle est réutilisée en cellules humaines (Wang 2012, PNAS,
  doi:10.1073/pnas.1116792109).
- Le décalage entre la numérotation publiée (340 nt) et le transcrit v27 (333 nt) est inconnu. La plage
  a donc été élargie à 97–140 du transcrit v27, soit chr14:20343268-20343311 (brin −).
- RNAfold prédit une tige-boucle en 107–140, cohérente avec l'élément publié. La séquence exacte des
  20 nt n'a pas pu être lue : elle n'apparaît que dans une figure, et PMC ainsi que Google Patents
  étaient inaccessibles.

**Domaine publié, coordonnées dérivées (`derive_vienna`, poids 1) : ARNr 5S.** Deux éléments, d'après
Smirnov 2008 :

- partie proximale de l'hélice I, avec la paire G:U conservée ;
- région boucle E – hélice IV.

Les positions viennent de la structure d'énergie minimale (MFE) calculée par ViennaRNA 2.8.0 pour
RNA5S1 : hélice I en 1–9/110–118, boucle E – hélice IV en 73–102. Elles sont recopiées sur les 17
copies et doivent être confirmées par Rfam RF00001 sur Curta.

**Non cartographiés (`non_cartographie`) :**

- la tige-boucle de RMRP : elle existe selon Wang 2010, mais sa position n'apparaît que dans une figure ;
- la région de TERC, qui n'est pas un élément nécessaire.

**Sites des facteurs :** aucun n'est cartographié en coordonnées, qu'il s'agisse de GRSF1 et HuR sur
RMRP, de MRP-L18 et de la rhodanèse sur l'ARNr 5S, ou d'Ago2 sur les miARN. Pour la graine des miARN,
voir le chemin 1.

### 3.5 Annotation sur la VM (`scripts/annoter_chemin2.py`)

Pour chaque variant (chr, pos, ref, alt en hg38, au format TSV ou VCF) situé dans un exon d'un ARN de la
table, le script calcule trois choses :

1. **Preuve :** le meilleur rang, le meilleur rang chez l'homme, le nombre de publications humaines
   distinctes de rang ≤ 2 (comptées par DOI) et les contre-preuves.
2. **Élément :** le variant chevauche-t-il un élément de `elements_mito.tsv` ?
3. **Pliage ViennaRNA, REF contre ALT :**
   - ARN entier si sa longueur est ≤ 400 nt, sinon fenêtre de ± 75 nt sur le brin du transcrit
     (l'épissage n'est pas reconstruit) ;
   - ΔΔG = ΔG_MFE(ALT) − ΔG_MFE(REF) ;
   - distance de paires entre les deux structures MFE (non calculée pour les indels) ;
   - statut apparié ou non du nucléotide dans la structure REF ;
   - contrôle que REF concorde avec la référence.

Le résultat est écrit dans la colonne `NC_chemin2 = score|niveau|justification`. Si un variant touche
plusieurs ARN, ils sont séparés par ` // ` ; s'il a plusieurs ALT, par ` ;; `.

### 3.6 Seuils fixés a priori (écrits dans le code)

Aucun ajustement n'a été fait sur des variants de patients ni sur AlphaGenome, SpliceAI, SpliceAI2 ou
AlphaFold 3.

**Score** = points de preuve + points d'élément + points de pliage :

- points de preuve : rang 1 = 3 ; 2 = 2 ; 3 = 1 ; 4 = 0,5 ; 5 = 0,25 ;
- points d'élément : 2 pour un élément publié, 1 pour un élément dérivé ;
- pliage : +1 si |ΔΔG| ≥ 2,0 kcal/mol. Cette valeur correspond environ à un empilement de paires GC
  perdu, soit −2,4 à −3,4 kcal/mol selon les paramètres de Turner 2004 par défaut de ViennaRNA ;
- pliage : +1 si la distance de paires est ≥ max(5 ; 10 % de la longueur repliée).

**Niveaux :**

- `aucun` si le variant est hors de la table.
- `faible`, plafond du chemin, seulement si toutes ces conditions sont remplies :
  - preuve de rang ≤ 2 chez l'homme dans au moins 2 publications distinctes ;
  - variant dans un élément `publie_plage` ;
  - REF concordant ;
  - ARN non multicopie.
- `hypothese` dans tous les autres cas. Le pliage ne fait jamais passer à `faible`.

Les familles multicopies (RNA5S*) sont plafonnées à `hypothese`.

**Recherche exploratoire du motif d'import dans d'autres ARN : non implémentée.** La séquence exacte de
la tige-boucle n'a pas pu être vérifiée. Une recherche de motif reposant sur une séquence non vérifiée
serait trompeuse.

## 4. Limites et points douteux

- **La revue a pu manquer des données.** Le texte intégral était inaccessible pour Zhang 2014,
  Smirnov 2008, Smirnov 2011, Cheng 2018, Leucci 2016 et Wang 2012. La classification de leur méthode,
  notamment la présence d'un traitement nucléase, repose sur le résumé et la synthèse de Jeandard 2019.
  Elle a été classée de façon conservatrice (`fractionnement_seul` quand ce n'était pas vérifié).
- **Une seule plage d'import a des coordonnées publiées, celle de RPPH1,** et elle est approximative
  (± 6 nt) à cause de la numérotation. Le seul niveau `faible` atteignable concerne donc RPPH1, un ARN
  dont la fonction mitochondriale est réfutée (RNase P protéique). C'est cohérent avec le statut
  d'hypothèse.
- **Pliage.** La MFE d'un ARN de 250 à 350 nt est instable : la distance de paires peut être grande pour
  un seul nucléotide. Par exemple, le témoin T03, au fondateur RMRP, donne 101 paires. Un indicateur
  fondé sur l'ensemble (probabilités de paires, distance d'ensemble) serait plus robuste ; c'est une
  piste pour la v1. Une fenêtre génomique de ± 75 nt ne représente pas la structure des lncRNA épissés.
- **Coordonnées GENCODE v27.** RMRP v27 n.72 correspond au variant fondateur noté n.71A>G ou n.72A>G
  selon les nomenclatures. La position exacte du témoin T03 reste à confirmer dans ClinVar (script
  Curta, étape 5).
- **Pas de bruit de fond gnomAD** dans cette version (inaccessible ici). Il est à mesurer sur Curta, par
  la proportion de variants rares gnomAD des ARN de la table qui tombent en `faible`, sans ajuster les
  seuils.
- **Les ARN de la table ne sont pas des gènes de maladie mitochondriale.** Les phénocopies cliniques
  (RMRP, TERC) relèvent d'autres mécanismes, le plus souvent nucléaires, et doivent être interprétées
  par le chemin 1 et le catalogue, pas par ce chemin.

## 5. Données à télécharger sur Curta

Le script `scripts/curta_preparer_chemin2.sh` récupère :

1. GENCODE v49 (GTF primary assembly) ;
2. le FASTA hg38 du Broad et son `.fai` ;
3. Rfam.cm, avec Infernal `cmscan`, pour remplacer les positions dérivées de l'ARNr 5S et situer RMRP et
   RPPH1 ;
4. ClinVar GRCh38, pour la région RMRP ;
5. la roue `ViennaRNA==2.8.0`, pour la VM.

Les numéros GEO des données APEX matrice (Fazal 2019, Kaewsapsak 2017) sont à relever dans les articles.
Ils ne sont pas recopiés ici, faute de vérification. À récupérer à la main :

- figure 6 et supplément de Wang 2010 : séquence et position des tiges-boucles RP et MRP ;
- texte intégral de Smirnov 2008 : positions des mutations de l'ARNr 5S ;
- supplément de Silver 2024 : liste des miARN.

## 6. Commandes

Sur la VM (ViennaRNA installé), depuis ce dossier :

```
python3 -m pip install ViennaRNA==2.8.0
bash tests/test_chemin2.sh
python3 scripts/annoter_chemin2.py --variants variants.tsv --out variants_chemin2.tsv
```

Sur Curta :

```
bash scripts/curta_preparer_chemin2.sh /chemin/vers/travail_chemin2
```

Reconstruction du paquet à partir d'un GTF et d'un FASTA déjà présents :

```
python3 scripts/preparer_references.py coordonnees --gtf gencode.gtf.gz --genes donnees/genes_chemin2.txt --out donnees/arn_coordonnees_hg38.tsv
python3 scripts/preparer_references.py sequences --coord donnees/arn_coordonnees_hg38.tsv --fasta Homo_sapiens_assembly38.fasta --out donnees/ref_regions_hg38.fa
```
