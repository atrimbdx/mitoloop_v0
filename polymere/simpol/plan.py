#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enumere les taches de la campagne, une par ligne.

Le job array lit la ligne correspondant a son indice. Le plan est ecrit
une fois et relu par chaque tache : les indices restent stables si la
campagne est relancee, et une tache perdue se rejoue seule.
"""
import argparse

import config as cfg


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--locus", default="tyrp1", choices=sorted(cfg.LOCI))
    p.add_argument("--conditions", nargs="+", default=["MEL", "FIBRO"])
    p.add_argument("--voies", nargs="+", default=["A", "B"])
    p.add_argument("--replicats", type=int, default=20)
    p.add_argument("--sortie", default="plan.tsv")
    a = p.parse_args()

    lignes = ["condition\tvoie\tjeu\treplicat"]
    for cond in a.conditions:
        for voie in a.voies:
            for jeu in range(len(cfg.jeux(voie if voie != "AB" else "A"))):
                for rep in range(a.replicats):
                    lignes.append(f"{cond}\t{voie}\t{jeu}\t{rep}")
    with open(a.sortie, "w") as f:
        f.write("\n".join(lignes) + "\n")
    print(f"  {len(lignes) - 1} tâches → {a.sortie}")
    print(f"  soumettre : sbatch --array=1-{len(lignes) - 1} "
          f"slurm/02_trajectoires.sbatch")


if __name__ == "__main__":
    main()
