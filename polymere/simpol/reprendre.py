#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Inventaire d'une campagne interrompue, et plan de reprise.

Une tache est comptee terminee quand son .npy ET son .json existent, et que
le .npy a la forme attendue. L'ordre d'ecriture de une_trajectoire.py est
.npy puis .json : un .npy sans .json est donc une tache tuee en cours
d'ecriture, et son fichier est tronque. Exiger les deux ecarte ces
fichiers-la, qui passeraient sinon pour des resultats.

Les durees mesurees sur les taches terminees servent a proposer une
walltime. Elles viennent de votre materiel, ce qui vaut mieux que n'importe
quelle estimation faite ailleurs.

USAGE
  python3 reprendre.py                         inventaire + plan_reprise.tsv
  python3 reprendre.py --nettoyer              supprime aussi les .npy tronqués
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

import config as cfg


def etiquette(locus, cond, voie, jeu, rep):
    return f"{locus}_{cond}_{voie}_jeu{int(jeu):03d}_rep{int(rep):02d}"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--locus", default="tyrp1", choices=sorted(cfg.LOCI))
    p.add_argument("--bin", type=int, dest="taille_bin_cli")
    p.add_argument("--marge", type=int)
    p.add_argument("--plan", default="plan.tsv")
    p.add_argument("--resultats", default="resultats")
    p.add_argument("--sortie", default="plan_reprise.tsv")
    p.add_argument("--nettoyer", action="store_true",
                   help="supprime les .npy sans .json (écritures tronquées)")
    p.add_argument("--marge-temps", type=float, default=1.6, dest="marge_temps",
                   help="facteur appliqué à la durée la plus longue observée")
    a = p.parse_args()
    cfg.surcharger(a.locus, a.taille_bin_cli, a.marge)
    n_attendu = cfg.geometrie(a.locus)[0]

    plan = pd.read_csv(a.plan, sep="\t")
    res = Path(a.resultats)
    faits, manquants, tronques, durees = [], [], [], []
    vues = set()

    for i, r in plan.iterrows():
        e = etiquette(a.locus, r.condition, r.voie, r.jeu, r.replicat)
        npy, js = res / f"{e}.npy", res / f"{e}.json"
        if js.exists() and npy.exists():
            try:
                forme = np.load(npy, mmap_mode="r").shape
                m = json.loads(js.read_text())
            except Exception:
                forme = m = None
            # La forme est comparee au nombre de billes inscrit dans le
            # .json a l'execution, et non a la geometrie courante : lance
            # avec un --bin different de celui de la campagne, le script
            # declarerait sinon toutes les taches manquantes.
            if forme and m and forme[1] == m.get("n_billes"):
                vues.add(m.get("n_billes"))
                durees.append(dict(voie=r.voie, jeu=int(r.jeu),
                                   condition=r.condition,
                                   heures=m.get("secondes", np.nan) / 3600,
                                   blocs=forme[0]))
                faits.append(i)
                continue
        if npy.exists() and not js.exists():
            tronques.append(npy)
        manquants.append(i)

    if vues and n_attendu not in vues:
        print(f"\n  ATTENTION la campagne a été lancée avec "
              f"{sorted(vues)} billes, la géométrie courante en donne "
              f"{n_attendu}.\n  Vérifier --bin et --marge avant de relancer : "
              f"les deux lots ne seraient pas comparables.")
    print(f"\n  plan      : {len(plan):,} tâches  ({a.plan})")
    print(f"  terminées : {len(faits):,}")
    print(f"  à relancer: {len(manquants):,}")
    if tronques:
        print(f"  tronquées : {len(tronques):,} .npy sans .json "
              f"(tuées pendant l'écriture)")
        if a.nettoyer:
            for f in tronques:
                f.unlink()
            print(f"              supprimés")
        else:
            print(f"              relancer avec --nettoyer pour les supprimer")

    if durees:
        D = pd.DataFrame(durees)
        print(f"\n  durées mesurées sur les tâches terminées (heures) :")
        g = D.groupby("voie").heures
        for voie, s in g:
            print(f"    voie {voie} : médiane {s.median():.2f}  "
                  f"p95 {s.quantile(0.95):.2f}  max {s.max():.2f}  "
                  f"({len(s)} tâches)")
        pire = D.heures.max()
        reco = max(2, int(math.ceil(pire * a.marge_temps)))
        print(f"\n  la plus longue a pris {pire:.2f} h "
              f"→ walltime conseillée : {reco} h")
        # Les taches terminees l'ont ete sous l'ancienne walltime : la vraie
        # duree des taches tuees est inconnue et superieure. La marge porte
        # donc sur un maximum sous-estime, d'ou un facteur genereux.
        print(f"  (les tâches tuées étaient plus longues que ce maximum, "
              f"qui est donc minoré)")
    else:
        reco = 12
        print("\n  aucune tâche terminée : walltime conseillée 12 h par défaut")

    if not manquants:
        print("\n  rien à relancer. Suite :")
        print("    sbatch slurm/03_comparer.sbatch")
        return

    sous = plan.loc[manquants].reset_index(drop=True)
    sous.to_csv(a.sortie, sep="\t", index=False)
    print(f"\n  {len(sous):,} tâches écrites dans {a.sortie}")
    print(f"\n  relancer :")
    print(f"    PLAN={a.sortie} sbatch --time={reco}:00:00 "
          f"--array=1-{len(sous)}%200 slurm/02_trajectoires.sbatch")
    print(f"\n  puis, quand c'est fini :")
    print(f"    sbatch slurm/03_comparer.sbatch")
    print(f"    sbatch slurm/04_verifier.sbatch")

    reste = (sous.groupby(["condition", "voie"]).size()
             .reset_index(name="tâches"))
    print(f"\n  répartition de ce qui reste :")
    for _, r in reste.iterrows():
        print(f"    {r.condition:6s} voie {r.voie} : {r['tâches']:4d}")


if __name__ == "__main__":
    main()
