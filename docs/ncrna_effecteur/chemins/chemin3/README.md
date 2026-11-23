# Chemin 3 : l'élément « non codant » porte un micropeptide mitochondrial

*Version 0, 19 novembre 2026. Première implémentation du chemin 3 décrit dans `../../chemins_mito.md`.*

## Objet

Un variant situé dans un lncRNA, une UTR ou un uORF peut tomber dans un **sORF traduit** (petit cadre
ouvert dont la traduction est attestée). Si le peptide produit est mitochondrial, le variant doit être lu
comme **codant** : perte d'initiation, stop gagné, décalage du cadre, stop perdu, indel en phase, faux-sens.
Le chemin 3 produit, pour chaque variant, la colonne `NC_chemin3 = score|niveau|justification`.

Périmètre : génome nucléaire seulement. Les peptides de l'ADNmt (humanine, MOTS-c, SHLP) sont exclus :
tout sORF ou variant sur chrM/MT est écarté. Les paralogues nucléaires « humanin-like » (MTRNR2L*) ne sont
pas retenus non plus, faute de preuve établie (catalogue : `non_verifie_priorite_basse`).

## Fichiers

| Fichier | Rôle | Où |
|---|---|---|
| `micropeptides_mito_etablis.tsv` | Table manuelle des micropeptides mitochondriaux nucléaires établis (11 entrées) | — |
| `scripts/telecharger_curta.sh` | Téléchargements (Ribo-seq ORFs, GENCODE, MitoCarta3.0, outils) | Curta |
| `scripts/preparer_sorf.py` | Table sORF → niveau mito, GTF et BED pour VEP, FASTA pour les prédicteurs | Curta |
| `scripts/annoter_chemin3.py` | Annotation des variants, hors ligne, bibliothèque standard seulement | VM |
| `scripts/sorf_lib.py` | Constantes a priori, lecture des sORF, règles de niveau | Curta et VM |
| `test/lancer_tests.py` | Test sur un **jeu construit à la main** (synthétique) | partout |
| `tables/phaseI_entree_sans_genome.*` | Caractéristiques des peptides des fichiers d'entrée de la phase I (voir plus bas) | — |

## Sources

- **Ribo-seq ORFs GENCODE** : Mudge et al. 2022, *Nat Biotechnol* 40:994-999,
  doi:10.1038/s41587-022-01369-0 (PMID 35831657). Le catalogue final est sur le FTP de l'EBI, **bloqué depuis le conteneur de développement**.
  Le dépôt GitHub des auteurs (`jorruior/gencode-riboseqORFs`, licence MIT, commit `b9ff0cf`) est, lui,
  accessible : son dossier `phaseI/` contient les **fichiers d'entrée** de la phase I (ORF de 7 études
  Ribo-seq, hg38, BED6 par exon en coordonnées 1-based, stop compris, et FASTA protéique). Ce n'est pas la
  liste finale unifiée : il manque l'identifiant GENCODE, le type d'ORF (uORF, lncRNA…) et le filtrage.
- **OpenProt** : Brunet et al. 2021, *Nucleic Acids Res*, doi:10.1093/nar/gkaa1036 (PMID 33179748) ;
  OpenProt 2.0 : Leblanc et al. 2024, *Nucleic Acids Res*, doi:10.1093/nar/gkad1050 (PMID 37956315).
  Site bloqué ici ; aucune copie publique trouvée sur GitHub ou bioconda. Téléchargement manuel sur Curta.
- **MitoCarta3.0** : Rath et al. 2021, *Nucleic Acids Res*, doi:10.1093/nar/gkaa1011 (PMID 33174596).
  Bloqué ici ; aucune copie locale.
- **Séquence d'adressage** : TargetP 2.0, Almagro Armenteros et al. 2019, *Life Sci Alliance*,
  doi:10.26508/lsa.201900429 (PMID 31570514) ; MitoFates, Fukasawa et al. 2015, *Mol Cell Proteomics*,
  doi:10.1074/mcp.M114.043083 (PMID 25670805). Ni l'un ni l'autre n'est sur bioconda ou PyPI (vérifié) :
  TargetP demande une licence DTU, MitoFates une archive Perl des auteurs.
- **Hélice transmembranaire** : TMbed, Bernhofer et al. 2022, *BMC Bioinformatics*,
  doi:10.1186/s12859-022-04873-x (PMID 35941534), installable par pip depuis GitHub. DeepTMHMM (licence DTU,
  ou service BioLib : n'y envoyer que les peptides de référence publics, jamais de séquence de patient).
  Repli interne : hydropathie de Kyte J et Doolittle RF 1982, *J Mol Biol*, doi:10.1016/0022-2836(82)90515-0
  (PMID 7108955).
- Contexte sur les présequences (hélices amphiphiles) : von Heijne G 1986, *EMBO J*,
  doi:10.1002/j.1460-2075.1986.tb04364.x (PMID 3015599).
- Troncature distale : arbre PVS1, Abou Tayoun et al. 2018, *Hum Mutat*, doi:10.1002/humu.23626
  (PMID 30192042).
- Conservation (prévue, **non implémentée**) : PhyloCSF, Lin et al. 2011, *Bioinformatics*,
  doi:10.1093/bioinformatics/btr209 (PMID 21685081).

Les références des micropeptides sont dans la colonne `references` de la table manuelle. Toutes les
références de ce dossier ont été vérifiées dans PubMed le 19 novembre 2026.

## Table manuelle des micropeptides mitochondriaux établis

Critère d'entrée : peptide nucléaire de moins d'environ 100 aa, issu d'un sORF (lncRNA reclassé, uORF,
ORF alternatif) ou découvert comme tel, avec **localisation mitochondriale démontrée** et **fonction
mitochondriale** testée. Les petites sous-unités OXPHOS canoniques (ATP5MK…) n'y sont pas : ce sont des
gènes MitoCarta classiques.

- `fort` : localisation démontrée et perte de fonction avec défaut mitochondrial mesuré (OXPHOS, assemblage,
  traduction mitochondriale), ou deux équipes indépendantes : MTLN, BRAWNIN (C12orf73), SMIM4, MOCCI
  (C15orf48), SMIM26, SLC35A4-MP.
- `modere` : une seule équipe, contexte tumoral, fonction non bioénergétique ou données discordantes :
  altMiD51/MIEF1-MP, PIGBOS, ASAP (LINC00467), MP31 (uORF de PTEN), STMP1.

La colonne `appariement` dit comment relier l'entrée à un sORF : `symbole:` (CDS canonique GENCODE, pour les
peptides déjà annotés codants) ou `hote:…;longueur:…` (sORF du catalogue). **S'il reste plusieurs
candidats, aucun niveau manuel n'est donné** et le cas est écrit dans le journal. Plusieurs longueurs et
symboles HGNC actuels (UQCC5, UQCC6, PIGBOS1, ancien nom de STMP1) sont marqués « à vérifier ».

## Critères a priori

Aucun seuil n'a été ajusté sur des variants de patients ni sur des sorties d'AlphaGenome, SpliceAI,
SpliceAI2 ou AlphaFold 3 (AlphaFold 3 n'est pas utilisé). Tous les seuils sont dans `scripts/sorf_lib.py`.

### Niveau mito du peptide (premier critère satisfait)

| Niveau | Règle | Poids |
|---|---|---|
| `fort` / `modere` | Entrée appariée de la table manuelle (niveau de la table) | 1,0 / 0,7 |
| `modere` | Peptide (séquence exacte ou `orf_id`) dans un protéome mitochondrial publié fourni en entrée | 0,7 |
| `faible` | Présequence prédite par un outil dédié : TargetP 2.0 classe `mTP` ou MitoFates « Possessing mitochondrial presequence » (décision par défaut de chaque outil) | 0,4 |
| `hypothese` | Hélice TM prédite (TMbed, DeepTMHMM) ; ou, **seulement si aucun outil n'a tourné**, repli : Kyte-Doolittle ≥ 1,6 sur 19 résidus, ou au moins 3 R/K et au plus 1 D/E dans les 20 premiers résidus d'un peptide d'au moins 30 aa | 0,15 |
| `aucun` | Rien | 0 |

Une hélice TM seule n'est pas spécifique de la mitochondrie (RE, membrane plasmique), d'où le niveau
`hypothese`. Si le gène hôte d'un uORF est dans MitoCarta3.0, c'est noté dans la justification sans
relever le niveau : le variant peut alors agir sur la CDS principale (voie régulatrice, mitoloop).

### Conséquence sur le sORF

| Conséquence | Poids |
|---|---|
| `start_lost` (codon d'initiation modifié, AUG ou non), `stop_gained`, `frameshift` | 1,0 |
| `stop_gained` ou `frameshift` retirant moins de 10 % du peptide (analogie PVS1) | 0,5 |
| `splice_site_orf` (±2 nt d'une jonction interne d'un sORF épissé) | 0,8 |
| `stop_lost`, `inframe_indel` | 0,6 |
| `missense` | 0,5 |
| `synonymous`, `stop_retained` | 0 |

### Score et niveau final

`score = poids_conséquence × poids_mito` ; niveau : `fort` si ≥ 0,8, `modere` si ≥ 0,45, `faible` si ≥ 0,2,
`hypothese` si > 0, sinon `aucun`. Exemples : perte d'initiation d'un peptide `fort` → 1,0 `fort` ;
faux-sens dans un peptide `fort` → 0,5 `modere` ; décalage dans un peptide `hypothese` → 0,15 `hypothese`.
Si plusieurs sORF sont touchés, le plus haut score est gardé et le nombre des autres est indiqué. Les allèles
multiples sont séparés par `;;`. Ce niveau sert au tri : il ne donne aucun point ACMG.

## Préparation sur Curta (machine connectée)

Commandes à lancer une par une, depuis le dossier `chemin3/` :

```
bash scripts/telecharger_curta.sh $HOME/donnees_chemin3
D=$HOME/donnees_chemin3
python3 -c "import pandas as pd; pd.read_excel('$D/mitocarta/Human.MitoCarta3.0.xls', sheet_name='A Human MitoCarta3.0').to_csv('$D/mitocarta/Human.MitoCarta3.0.tsv', sep='\t', index=False)"
python3 scripts/preparer_sorf.py --phaseI-bed $D/riboseq/gencode-riboseqORFs/phaseI/sORFs_genomic_hg38.prot.20210329.bed --phaseI-fa $D/riboseq/gencode-riboseqORFs/phaseI/sORFs_genomic_hg38.prot.20210329.fa --genome $D/gencode/GRCh38.primary_assembly.genome.fa --gencode-gtf $D/gencode/gencode.v46.annotation.gtf.gz --mitocarta $D/mitocarta/Human.MitoCarta3.0.tsv --fasta-peptides $D/peptides_chemin3.fa --sortie $D/chemin3_passe1
```

Quand le catalogue GENCODE final (BED12) et OpenProt (BED12 + FASTA) sont là, les ajouter avec
`--bed12 FICHIER --bed12-fa FASTA --bed12-source NOM` (une fois par source). **Le format exact de ces
fichiers n'a pas pu être vu** : vérifier qu'il s'agit bien de BED12 0-based avec le CDS en
thickStart/thickEnd, sinon adapter `lire_bed12`.

**Catalogue GENCODE phase 2 « primary » (vérifié le 23/11/2026)** : `Ribo-seq_ORFs.primary.bed` (BED12
0-based, CDS = thickStart/thickEnd, codon stop compris, hg38) et `Ribo-seq_ORFs.primary.faa` (même nom
d'ORF que la colonne 4, `c1riboseqorf2598:LINC00467:lncRNA` ; gène hôte = 2e champ), 10 127 ORF (uORF 6 437,
lncRNA 2 108, dORF 657, PT 538, uoORF 267, intORF 107, doORF 13), déposés dans `reannot/donnees/`. ASAP
(LINC00467, 94 aa) y est (c1riboseqorf2598). Commande sur Curta (chemins de la préparation du 20/11) :

```
python3 scripts/preparer_sorf.py --phaseI-bed /scratch/$USER/chemin3_donnees/riboseq/gencode-riboseqORFs/phaseI/sORFs_genomic_hg38.prot.20210329.bed --phaseI-fa /scratch/$USER/chemin3_donnees/riboseq/gencode-riboseqORFs/phaseI/sORFs_genomic_hg38.prot.20210329.fa --bed12 $HOME/reannot_ncrna/donnees/Ribo-seq_ORFs.primary.bed --bed12-fa $HOME/reannot_ncrna/donnees/Ribo-seq_ORFs.primary.faa --bed12-source gencode_riboseq_primary --genome $HOME/mitoloop_v0/data/GRCh38.fa --gencode-gtf $HOME/mitoloop_v0/data/gencode.gtf.gz --mitocarta /scratch/$USER/chemin1_donnees/mitocarta/Human.MitoCarta3.0.tsv --fasta-peptides /scratch/$USER/chemin3_donnees/peptides_chemin3.fa --sortie /scratch/$USER/chemin3_donnees/chemin3
```

**OpenProt 2.2 (« min_2_pep », au moins deux peptides uniques, GRCh38, Ensembl 106 ; format vérifié le
23/11/2026)** : déposé sur Curta (trop lourd pour GitHub) dans `/scratch/atrimouille/nc_sources/`. BED12 0-based,
noms en `chr`, CDS = thickStart/thickEnd (102 nt pour 33 aa + stop), nom = accession seule (`IP_057218` pour
une protéine alternative, `II_…` pour une nouvelle isoforme d'une protéine connue), sans gène hôte ; une même
accession peut être sur plusieurs loci (paralogues) et une même ligne répétée (une par transcrit). Traitement :
- `--bed12-garder '^IP_'` : protéines alternatives seulement (les isoformes II_ chevauchent des CDS déjà
  annotées, que VEP interprète) ; « - » pour une source sans filtre ;
- lignes répétées fusionnées par le dédoublonnage par coordonnées ;
- gène hôte attribué par chevauchement avec `--gencode-gtf` (gène contenant le 1er exon, même brin, sinon
  « antisens:GENE ») ;
- FASTA : clé = début de l'en-tête jusqu'au premier « | » ou espace.

Commande complète (Ribo-seq primary + OpenProt), après décompression :

```
mkdir -p /scratch/$USER/nc_sources/openprot
unzip -o /scratch/$USER/nc_sources/human-openprot-2_2-altprots+isoforms-min_2_pep-GRCh38.106.bed.zip -d /scratch/$USER/nc_sources/openprot
unzip -o /scratch/$USER/nc_sources/human-openprot-2_2-altprots+isoforms-min_2_pep-GRCh38.106.fasta.zip -d /scratch/$USER/nc_sources/openprot
python3 scripts/preparer_sorf.py --phaseI-bed /scratch/$USER/chemin3_donnees/riboseq/gencode-riboseqORFs/phaseI/sORFs_genomic_hg38.prot.20210329.bed --phaseI-fa /scratch/$USER/chemin3_donnees/riboseq/gencode-riboseqORFs/phaseI/sORFs_genomic_hg38.prot.20210329.fa --bed12 $HOME/reannot_ncrna/donnees/Ribo-seq_ORFs.primary.bed --bed12-fa $HOME/reannot_ncrna/donnees/Ribo-seq_ORFs.primary.faa --bed12-source gencode_riboseq_primary --bed12-garder - --bed12 /scratch/$USER/nc_sources/openprot/human-openprot-2_2-altprots+isoforms-min_2_pep-GRCh38.106.bed --bed12-fa /scratch/$USER/nc_sources/openprot/human-openprot-2_2-altprots+isoforms-min_2_pep-GRCh38.106.fasta --bed12-source openprot_min2pep --bed12-garder '^IP_' --genome $HOME/mitoloop_v0/data/GRCh38.fa --gencode-gtf $HOME/mitoloop_v0/data/gencode.gtf.gz --mitocarta /scratch/$USER/chemin1_donnees/mitocarta/Human.MitoCarta3.0.tsv --fasta-peptides /scratch/$USER/chemin3_donnees/peptides_chemin3.fa --sortie /scratch/$USER/chemin3_donnees/chemin3
```

Les noms des fichiers décompressés sont supposés identiques à ceux des archives sans `.zip` : vérifier avec
`ls /scratch/$USER/nc_sources/openprot` avant de lancer.

Prédicteurs, lancés sur `peptides_chemin3.fa` (noms de peptides = empreinte stable de la séquence) :

```
targetp -fasta $D/peptides_chemin3.fa -org non-pl -format short -prefix $D/targetp_chemin3
perl MitoFates.pl $D/peptides_chemin3.fa metazoa > $D/mitofates_chemin3.tsv
python3 -m tmbed predict -f $D/peptides_chemin3.fa -p $D/tmbed_chemin3.pred
```

(Syntaxe à confirmer avec l'aide de chaque outil installé.) Puis la passe finale :

```
python3 scripts/preparer_sorf.py --phaseI-bed $D/riboseq/gencode-riboseqORFs/phaseI/sORFs_genomic_hg38.prot.20210329.bed --phaseI-fa $D/riboseq/gencode-riboseqORFs/phaseI/sORFs_genomic_hg38.prot.20210329.fa --genome $D/gencode/GRCh38.primary_assembly.genome.fa --gencode-gtf $D/gencode/gencode.v46.annotation.gtf.gz --mitocarta $D/mitocarta/Human.MitoCarta3.0.tsv --targetp $D/targetp_chemin3_summary.targetp2 --mitofates $D/mitofates_chemin3.tsv --tm $D/tmbed_chemin3.pred --sortie $D/chemin3
sort -k1,1 -k4,4n $D/chemin3.gtf > $D/chemin3.tri.gtf
bgzip -f $D/chemin3.tri.gtf
tabix -f -p gff $D/chemin3.tri.gtf.gz
bgzip -f $D/chemin3.bed
tabix -f -p bed $D/chemin3.bed.gz
```

Lire `chemin3.log.txt` : appariements de la table manuelle, nombre de peptides discordants avec le génome,
répartition des niveaux. Les listes de peptides de protéomes mitochondriaux publiés (`--proteome-mito`,
une séquence ou un `orf_id` par ligne, ou colonne `peptide`) restent à constituer à partir des
suppléments d'articles ; aucune n'est fournie ici.

**Paquet à transférer sur la VM** : `chemin3.table_sorf.tsv.gz`, `chemin3.tri.gtf.gz(.tbi)`,
`chemin3.bed.gz(.tbi)`, `scripts/sorf_lib.py`, `scripts/annoter_chemin3.py`, le journal et
`DATE_TELECHARGEMENT.txt`.

## Annotation sur la VM (hors ligne)

```
python3 scripts/annoter_chemin3.py --variants variants.vcf --table chemin3.table_sorf.tsv.gz --sortie variants.chemin3.tsv
```

Entrée : VCF (colonnes `#CHROM POS REF ALT`) ou TSV avec `chr pos ref alt`, hg38, avec ou sans « chr ».
L'allèle de référence est contrôlé contre la séquence du sORF stockée dans la table ; en cas de désaccord, le
variant reçoit `0|aucun|allèle de référence discordant`.

Option VEP (même logique, conséquences calculées par VEP) : `--gtf chemin3.tri.gtf.gz --fasta
GRCh38.primary_assembly.genome.fa` pour les conséquences dans les sORF, et `--custom` sur `chemin3.bed.gz`
pour récupérer `orf_id|niveau_mito|gène_hôte` (la syntaxe de `--custom` dépend de la version de VEP). Le GTF ne
contient que le sORF (pas d'UTR) ; la CDS exclut le codon stop, comme dans GENCODE.

## Ce qui a été fait ici

- Test sur le **jeu construit à la main** (`test/jeu_main/`, génome `chrT1` synthétique, quatre sORF
  inventés TOYA à TOYD, variants inventés) : `python3 test/lancer_tests.py` → `OK`. Il vérifie les quatre
  niveaux mito, la lecture des sorties TargetP, MitoFates et 3 lignes (formats reconstitués d'après la
  documentation, à valider sur de vraies sorties), la traduction sur les deux brins, un sORF épissé, les
  conséquences (perte d'initiation, faux-sens, stop gagné proximal et distal, stop perdu, décalages, site
  d'épissage, synonyme), le refus d'une référence discordante, l'exclusion de chrM et la CDS du GTF. Ces
  sorties ne sont pas des résultats.
- Lecture des **fichiers d'entrée réels** de la phase I (GitHub) sans génome :
  `tables/phaseI_entree_sans_genome.table_sorf.tsv.gz` (8 762 sORF uniques par coordonnées, sur 15 141 ORF
  avec coordonnées ; 24 647 autres séquences du FASTA n'ont pas de coordonnées). Sans génome, la table n'a
  pas de séquence nucléotidique : **elle n'est pas utilisable pour annoter des variants**. Sans outil ni
  MitoCarta, seuls les replis jouent : 1 852 sORF `hypothese`, 6 910 `aucun`. Le repli N-terminal marque
  1 548 sORF (18 %) : il est peu spécifique, d'où son plafonnement à `hypothese` et son abandon dès qu'un
  outil dédié a tourné. Aucun des 11 micropeptides de la table manuelle n'est apparié dans ces fichiers
  (les CDS déjà annotées n'y figurent pas ; SLC35A4 a deux candidats, donc pas d'appariement).

## Limites et points douteux

- Catalogue GENCODE final, OpenProt, MitoCarta3.0 et le génome hg38 n'ont pas pu être téléchargés ici :
  aucun niveau `faible` ou plus n'a été calculé sur des données réelles.
- Les URL EBI (catalogue Ribo-seq et GENCODE) et Broad (MitoCarta) du script de téléchargement n'ont pas
  pu être testées. Celles d'OpenProt ne sont pas données : téléchargement manuel.
- Appariement de la table manuelle : par symbole ou par gène hôte et longueur. Les uORF (MIEF1, SLC35A4,
  PTEN) demandent les coordonnées exactes, à reprendre dans les articles ou dans le catalogue final.
- Un variant d'uORF peut aussi modifier la traduction de la CDS principale (régulation, hors chemin 3) ;
  pour MP31, PTEN est lui-même un gène de maladie.
- Les ORF alternatifs chevauchant une CDS canonique dans un autre cadre donnent deux lectures (VEP standard
  et chemin 3) : à interpréter ensemble.
- `stop_lost` : l'extension en aval n'est pas calculée (le GTF ne contient pas l'UTR 3').
- Départs non AUG : un changement du codon d'initiation est toujours `start_lost`, même entre codons
  quasi cognats.
- Conservation (PhyloCSF) et contrainte gnomAD dans les sORF : non implémentées (étape suivante).
- Pas de témoin réel : le catalogue ne contient pas de variant pathogène dans un micropeptide nucléaire
  (MTLN, LINC00467 et C15orf48 n'ont pas de variant de maladie). Le contrôle « bruit de fond » (variants rares
  de gnomAD dans les sORF) reste à faire sur Curta, pour mesure seulement.
