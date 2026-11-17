#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Une trajectoire : une voie, un jeu de parametres, un replicat.

C'est l'unite de travail du job array. Chaque appel occupe un coeur, ne
communique avec aucun autre, et ecrit un fichier de conformations.

USAGE
  python3 une_trajectoire.py --locus tyrp1 --condition MEL --voie A \\
      --jeu 7 --replicat 3 --sortie resultats/
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np


def graine_stable(*parties):
    """Graine reproductible d'un processus a l'autre.

    hash() de Python est sale par PYTHONHASHSEED des qu'une chaine entre
    dans le tuple : la meme ligne de commande donnait une graine differente
    a chaque appel, et la campagne n'etait donc pas reproductible. blake2b
    ne depend pas du processus.
    """
    cle = "|".join(str(p) for p in parties).encode()
    return int.from_bytes(hashlib.blake2b(cle, digest_size=4).digest(),
                          "big") % (2 ** 31)

import config as cfg
from extrusion import Extrusion, barrieres_depuis_bed
from polymere import simuler


def charger_matrice(chemin, chrom, debut, fin, taille_bin):
    """Matrice de contacts mesuree, projetee sur la grille simulee."""
    import pandas as pd
    d = pd.read_csv(chemin, sep="\t", comment="#")
    n = int(np.ceil((fin - debut) / taille_bin))
    M = np.zeros((n, n))
    vu = np.zeros((n, n), dtype=bool)

    def plage(s, e):
        return np.arange(max(0, int((s - debut) // taille_bin)),
                         min(n - 1, int((e - 1 - debut) // taille_bin)) + 1)

    for s1, e1, s2, e2, N in d[["startA", "endA", "startB", "endB",
                                "N"]].to_numpy():
        i, j = plage(s1, e1), plage(s2, e2)
        if i.size == 0 or j.size == 0:
            continue
        I, J = np.meshgrid(i, j, indexing="ij")
        np.add.at(M, (I, J), N / (i.size * j.size))
        np.add.at(M, (J, I), N / (i.size * j.size))
        vu[I, J] = vu[J, I] = True
    return M, vu


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--locus", default="tyrp1", choices=sorted(cfg.LOCI))
    p.add_argument("--bin", type=int, dest="taille_bin_cli",
                   help="change la résolution du locus")
    p.add_argument("--marge", type=int,
                   help="change le flanc simulé de part et d'autre")
    p.add_argument("--condition", required=True, help="MEL ou FIBRO")
    p.add_argument("--voie", required=True, choices=["A", "B", "AB"])
    p.add_argument("--jeu", type=int, required=True)
    p.add_argument("--replicat", type=int, default=0)
    p.add_argument("--matrice", help="TSV de la matrice mesurée (voie B)")
    p.add_argument("--ctcf", help="BED de sites CTCF orientés (voie A)")
    p.add_argument("--blocs", type=int, default=2000)
    p.add_argument("--pas-par-bloc", type=int, default=10_000,
                   dest="pas_par_bloc")
    p.add_argument("--plateforme", default=None)
    p.add_argument("--sortie", default="resultats")
    a = p.parse_args()
    cfg.surcharger(a.locus, a.taille_bin_cli, a.marge)

    c = cfg.LOCI[a.locus]
    n, d0, f0 = cfg.geometrie(a.locus)
    out = Path(a.sortie)
    out.mkdir(parents=True, exist_ok=True)
    etiq = f"{a.locus}_{a.condition}_{a.voie}_jeu{a.jeu:03d}_rep{a.replicat:02d}"
    # le locus et la condition entrent dans la graine : sans eux, deux
    # conditions tirees sur les memes entrees partent du meme etat
    graine = graine_stable(a.locus, a.condition, a.voie, a.jeu, a.replicat)

    traj_lef = liants = poids = None
    par = {}

    if a.voie in ("A", "AB"):
        jx = cfg.jeux("A")[a.jeu % len(cfg.jeux("A"))]
        par.update(jx)
        if a.ctcf:
            bp, bm = barrieres_depuis_bed(a.ctcf, c["chrom"], d0, f0,
                                          c["taille_bin"], n,
                                          force_max=jx["force_barriere"])
            # Empreinte du BED, relue par verifier.py : si les deux
            # conditions ont recu le meme fichier, la voie A n'a aucune
            # entree qui les distingue et leur differentiel est un nul.
            par["ctcf_empreinte"] = hashlib.blake2b(
                Path(a.ctcf).read_bytes(), digest_size=8).hexdigest()
            par["ctcf_fichier"] = Path(a.ctcf).name
        else:
            print("  aucun BED CTCF : extrusion sans barrière")
            bp = bm = np.zeros(n)
            par["ctcf_empreinte"] = "aucun"
        e = Extrusion(n, bp, bm,
                      separation_billes=jx["separation_kb"] * 1000
                      // c["taille_bin"],
                      processivite_billes=jx["processivite_kb"] * 1000
                      // c["taille_bin"],
                      graine=graine)
        # une position 1D par bloc tridimensionnel
        traj_lef = e.trajectoire(a.blocs * 4, 4)
        print(f"  voie A : {traj_lef.shape[1]} facteurs, "
              f"{len(traj_lef)} instants")

    if a.voie in ("B", "AB"):
        if not a.matrice:
            raise SystemExit("la voie B demande --matrice")
        jx = cfg.jeux("B")[a.jeu % len(cfg.jeux("B"))]
        par.update(jx)
        from liants import observe_sur_attendu, inferer
        M, vu = charger_matrice(a.matrice, c["chrom"], d0, f0, c["taille_bin"])
        L = observe_sur_attendu(M, vu)
        liants, poids, cout, _ = inferer(L, n_types=jx["n_types"],
                                         graine=graine, verbeux=False)
        print(f"  voie B : {liants.shape[1]} types retenus, coût {cout:.0f}")

    t0 = time.time()
    conf, plat = simuler(n, a.blocs, a.pas_par_bloc, traj_lef=traj_lef,
                         liants=liants, poids=poids,
                         eps_liant=par.get("eps_liant", 0.9),
                         graine=graine, plateforme=a.plateforme)
    dt = time.time() - t0

    np.save(out / f"{etiq}.npy", conf.astype(np.float32))
    meta = dict(locus=a.locus, condition=a.condition, voie=a.voie, jeu=a.jeu,
                replicat=a.replicat, n_billes=n, debut=d0, fin=f0,
                taille_bin=c["taille_bin"], parametres=par,
                blocs=a.blocs, pas_par_bloc=a.pas_par_bloc,
                plateforme=plat, secondes=round(dt, 1))
    (out / f"{etiq}.json").write_text(json.dumps(meta, indent=2))
    print(f"  {etiq} : {conf.shape} en {dt:.0f} s sur {plat}")


if __name__ == "__main__":
    main()
