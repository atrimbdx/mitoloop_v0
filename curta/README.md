# Curta : inventaire et archivage du travail MITODARK / mitoloop_v0

Deux scripts, à lancer dans cet ordre, depuis n'importe quel dossier de Curta.

1. `inventaire_curta.sh` : inventaire en **lecture seule** des dossiers
   MITODARK / mitoloop (ou de ceux donnés en argument). Écrit
   `inventaire_<date>/resume.txt` (tailles, types, gros fichiers, dépôts git,
   liens symboliques) et `cles.tsv` (fichiers de résultats : lignes, gènes,
   promoteurs, termes), qui permettent de retrouver quels fichiers donnent les
   chiffres du texte ANR (2 193 boucles à 907 promoteurs, 3 135 éléments rE2G,
   208 / 1 057 boucles du muscle avec un élément).
2. `archiver_curta.py` : garde (1) la version derrière le texte ANR, copiée en
   lecture seule avec le code (git bundle, tag `anr-aapg2027`) et le PDF, et
   (2) la version en développement, avec ses listes finales copiées ; archive
   tout le reste (un tar.gz par racine). Simulation par défaut ; `--appliquer`
   crée la référence et les archives ; `--supprimer` retire les originaux,
   seulement après relecture de l'archive. Les cibles des liens symboliques
   utilisés (par exemple `~/reannot_git/data`) sont gardées ; `~/reannot` et
   `~/reannot_git` ne sont jamais touchés. Liste à garder : `garder.exemple.tsv`.
