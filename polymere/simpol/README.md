# Simulation polymère gros-grain d'un locus, deux voies comparées

Deux modèles indépendants du repliement d'un locus, confrontés à la même
capture Hi-C et l'un à l'autre.

**Voie A — extrusion de boucles.** Des cohésines chargent sur la fibre,
poussent leurs deux jambes en sens opposé et s'arrêtent sur les sites CTCF
selon leur orientation. L'unique entrée est la position et l'orientation
des sites CTCF, tirées d'ENCODE. Le modèle est prédictif : il ne voit
jamais la matrice mesurée.

**Voie B — liants.** Chaque bille porte ou non l'affinité pour un type de
liant ; deux billes partageant un type sont pontées. Les affectations sont
inférées de la matrice elle-même, sans donnée externe et sans mécanisme
supposé. Le modèle est inverse : il voit la matrice.

Cette asymétrie est voulue et doit être dite. Comparer les qualités
d'ajustement brutes est défavorable à la voie A par construction. La
question qui a un sens est double : la voie A atteint-elle, avec CTCF pour
seule entrée, un accord comparable ; et les deux voies prédisent-elles la
même chose là où les données ne tranchent pas.

## Comment la comparaison est faite

La simulation produit une matrice tout-contre-tout, la mesure est une
capture. L'efficacité de capture d'un fragment est un facteur
multiplicatif sur toute sa ligne : normaliser chaque ligne par sa propre
somme le fait disparaître exactement, sans avoir à l'estimer. La
comparaison principale porte donc sur des profils de 4C virtuel pris aux
ancres annotées, et non sur les matrices entières. C'est la pratique du
modèle HiP-HoP, qui superpose profil simulé et profil mesuré point de vue
par point de vue.

La comparaison matricielle n'est gardée qu'en contrôle, sous forme de
corrélation stratifiée par distance : la décroissance avec la distance
génomique domine toute corrélation brute.

## Ce qui distingue les deux voies

Les contacts deux à deux sont probablement reproduits par l'une comme par
l'autre. Là où elles diffèrent, c'est sur la simultanéité : la fraction
des conformations où toutes les ancres du hub sont ensemble au même
instant. C'est une prédiction, et elle se teste par Tri-C ou Pore-C.

## Installation

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda create -y -p /scratch/$USER/envs/polymere -c conda-forge \
  python=3.11 numpy pandas scipy matplotlib openmm
```

Un environnement OpenMM pèse environ 1 Go. Sur le `/scratch`, attention à
la purge à 90 jours, qui porte sur la date de modification.

## Marche à suivre

```bash
# 0. réglages : chemins, locus, profondeur
vi config.sh

# 1. sites CTCF orientés depuis ENCODE          (une fois, ~20 min)
sbatch slurm/01_ctcf.sbatch

# 2. plan de campagne
$PY plan.py --locus tyrp1 --replicats 20
# → plan.tsv, et la commande de soumission à recopier

# 3. les trajectoires, en job array             (une tâche = un cœur)
sbatch --array=1-1080 slurm/02_trajectoires.sbatch

# 4. comparaison, quand tout est terminé
sbatch slurm/03_comparer.sbatch

# 5. contrôles et rendu 3D
sbatch slurm/04_verifier.sbatch
```

Après une interruption, quelle qu'en soit la cause :

```bash
$PY reprendre.py --locus tyrp1 --nettoyer
# → plan_reprise.tsv, les durées mesurées, et la commande à recopier
```

## Dimensionnement

Mesuré sur un cœur Curta, locus TYRP1, billes de 5 kb, 790 billes :
1 696 pas d'intégration par seconde. Une trajectoire de 2 × 10⁷ pas prend
donc 3,3 h et tient dans n'importe quelle walltime courte, sans découpage
ni reprise.

La campagne complète fait 1 080 tâches — deux conditions, deux voies,
27 jeux de paramètres, 20 réplicats — soit environ 3 500 heures-cœur. Les
tâches sont indépendantes : le temps écoulé ne dépend que du nombre de
cœurs obtenus. Les conformations occupent environ 10 Go au total.

Le GPU n'apporte pas grand-chose ici : à 790 billes le système
sous-utilise la carte, et sept GPU partagés sur tout Curta feront passer
plus de temps en file qu'en calcul.

## Fichiers

| | |
|---|---|
| `config.py` | locus, ancres, grilles de paramètres |
| `config.sh` | chemins propres à l'installation |
| `ctcf_encode.py` | pics ENCODE, balayage du motif MA0139, orientation |
| `extrusion.py` | voie A, simulation 1D des cohésines |
| `liants.py` | voie B, inférence des liants par recuit |
| `polymere.py` | moteur OpenMM, commun aux deux voies |
| `une_trajectoire.py` | unité de travail du job array |
| `plan.py` | énumère les tâches |
| `reprendre.py` | inventaire après interruption, plan de reprise |
| `comparer.py` | 4C virtuels, métriques, simultanéité |
| `bilan.py` | voie A contre voie B |
| `verifier.py` | éléments annotés, ancres absolues, différentiel et son nul |
| `rendu3d.py` | conformations représentatives, export PyMOL et ChimeraX |

## Ce que les simulations reçoivent, et ce qu'elles ne reçoivent pas

Vos éléments annotés **ne sont l'entrée d'aucune des deux voies**.
`une_trajectoire.py` ne lit que deux choses : le BED CTCF pour la voie A,
la matrice mesurée pour la voie B. Les ancres n'apparaissent qu'en aval,
dans `comparer.py` et `verifier.py`, comme points de vue 4C et comme
définition du hub. Les simulations peuvent donc être interrogées sur ces
positions sans circularité d'entrée.

La portée de la réponse diffère selon la voie. La voie A ne voit jamais la
matrice : si elle met en contact deux positions que vous avez annotées,
c'est une prédiction indépendante. La voie B a vu la matrice, donc le
signal qui a produit l'annotation ; elle ne peut pas confirmer l'existence
d'un contact, seulement établir qu'un liant est nécessaire à cet endroit
pour rendre compte des données. `verifier.py` étiquette chaque ligne en
conséquence.

## Ce que les scripts ne font pas

Les sites CTCF viennent d'ENCODE et non de vos cellules : c'est un
substitut, à écrire comme tel.

L'inférence des liants travaille sur un substitut analytique de la
probabilité de contact, et non sur une simulation à l'intérieur de la
boucle d'optimisation, ce qui serait hors de portée. C'est donc une
simplification de PRISMR ; la vérification se fait ensuite par dynamique
réelle. Un article doit dire : liants inférés par factorisation binaire de
la matrice observé sur attendu, puis validés par simulation.

La fraction de conformations à hub complet dépend du seuil de contact
retenu, `--seuil-hub`. C'est une prédiction relative entre voies et entre
conditions, pas une valeur absolue.

Une corrélation de différentiel ne se lit **jamais seule**. Mesuré sur ce
paquet : trois comparaisons structurellement nulles, où les deux
conditions avaient reçu des entrées identiques et ne différaient que par
la graine, ont donné des corrélations de +0,490, +0,215 et −0,325.
L'étendue du nul couvrait donc 0,8 à elle seule. `verifier.py` accompagne
toujours la valeur observée d'un nul par permutation des étiquettes, et
signale le cas où les deux conditions ont reçu le même BED CTCF — auquel
cas le différentiel de la voie A est un nul par construction.

Une image 3D est un tirage parmi des millions. Elle ne représente aucune
cellule et n'a pas valeur de preuve ; seules les moyennes d'ensemble en
ont. `rendu3d.py` écrit sous chaque conformation la fraction de l'ensemble
qu'elle représente, pour que la légende de figure le dise aussi.

L'échelle des fichiers PDB est arbitraire : la distance entre billes
voisines vaut 3,8 Å, seule valeur pour laquelle PyMOL et ChimeraX tracent
un ruban continu. Le facteur vers les nanomètres figure dans les REMARK.
