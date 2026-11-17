#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Voie A — simulation unidimensionnelle de l'extrusion de boucles.

Des facteurs d'extrusion, assimiles aux complexes cohesine, se chargent sur
la fibre, poussent leurs deux jambes en sens oppose, s'arretent sur les
sites CTCF selon leur orientation, et se dechargent avec une duree de vie
finie. La sortie est une trajectoire de positions de jambes, qui fournit au
modele tridimensionnel la liste des liens transitoires a chaque instant.

CONVENTION D'ORIENTATION
  Une jambe qui progresse vers les coordonnees croissantes est arretee par
  un site d'orientation moins ; une jambe qui progresse vers les
  coordonnees decroissantes est arretee par un site d'orientation plus. Une
  boucle se forme donc entre deux sites convergents, un plus a gauche et un
  moins a droite, ce qui est la geometrie observee.

Aucune jambe ne traverse une autre jambe : deux facteurs qui se rencontrent
se bloquent mutuellement.
"""
import numpy as np


class Extrusion:
    def __init__(self, n_billes, barrieres_plus, barrieres_moins,
                 separation_billes=40, processivite_billes=60,
                 vitesse=1.0, graine=0, decharge_bloquee=1.0):
        """
        barrieres_plus / barrieres_moins : tableaux de taille n_billes
            donnant la probabilite d'arret d'une jambe arrivant sur la
            bille, pour chacune des deux orientations. Zero partout revient
            a une extrusion sans barriere.
        separation_billes : une cohesine toutes les N billes.
        processivite_billes : taille moyenne de boucle sans obstacle.
        decharge_bloquee : facteur multiplicatif sur le taux de decharge
            quand les deux jambes sont arretees.
        """
        self.n = int(n_billes)
        self.bp = np.asarray(barrieres_plus, dtype=float)
        self.bm = np.asarray(barrieres_moins, dtype=float)
        self.rng = np.random.default_rng(graine)
        self.n_lef = max(int(round(self.n / separation_billes)), 1)
        # une boucle croit de deux billes par pas : la duree de vie se
        # deduit de la processivite voulue
        self.duree = max(processivite_billes / (2.0 * vitesse), 1.0)
        self.vitesse = float(vitesse)
        self.decharge_bloquee = float(decharge_bloquee)
        self.pos = np.zeros((self.n_lef, 2), dtype=np.int64)
        self.occupe = np.zeros(self.n, dtype=bool)
        for k in range(self.n_lef):
            self._charger(k)

    # ------------------------------------------------------------------
    def _charger(self, k):
        for _ in range(200):
            i = int(self.rng.integers(0, self.n - 1))
            if not self.occupe[i] and not self.occupe[i + 1]:
                self.pos[k] = (i, i + 1)
                self.occupe[i] = self.occupe[i + 1] = True
                return
        self.pos[k] = (0, 1)          # fibre saturee : on insiste pas

    def _decharger(self, k):
        g, d = self.pos[k]
        self.occupe[g] = self.occupe[d] = False

    # ------------------------------------------------------------------
    def pas(self):
        """Un pas de temps : deplacement des jambes, puis decharge."""
        ordre = self.rng.permutation(self.n_lef)
        bloque = np.zeros(self.n_lef, dtype=bool)

        for k in ordre:
            g, d = self.pos[k]

            # jambe gauche, vers les coordonnees decroissantes
            bg = True
            if g > 0 and self.rng.random() < self.vitesse:
                if not self.occupe[g - 1] and self.rng.random() >= self.bp[g]:
                    self.occupe[g] = False
                    g -= 1
                    self.occupe[g] = True
                    bg = False

            # jambe droite, vers les coordonnees croissantes
            bd = True
            if d < self.n - 1 and self.rng.random() < self.vitesse:
                if not self.occupe[d + 1] and self.rng.random() >= self.bm[d]:
                    self.occupe[d] = False
                    d += 1
                    self.occupe[d] = True
                    bd = False

            self.pos[k] = (g, d)
            bloque[k] = bg and bd

        # decharge, acceleree quand les deux jambes sont arretees
        for k in range(self.n_lef):
            taux = 1.0 / self.duree
            if bloque[k]:
                taux *= self.decharge_bloquee
            if self.rng.random() < taux:
                self._decharger(k)
                self._charger(k)

    # ------------------------------------------------------------------
    def trajectoire(self, n_pas, intervalle, chauffe=None):
        """Renvoie un tableau (instants, facteurs, 2) de positions."""
        if chauffe is None:
            chauffe = int(10 * self.duree)
        for _ in range(chauffe):
            self.pas()
        sorties = []
        for t in range(n_pas):
            self.pas()
            if t % intervalle == 0:
                sorties.append(self.pos.copy())
        return np.array(sorties)


# ======================================================================
def barrieres_depuis_bed(chemin, chrom, debut, fin, taille_bin, n_billes,
                         force_max=0.95, colonne_signal=6):
    """Lit un BED de sites CTCF orientes et le projette sur les billes.

    Colonnes attendues : chrom, start, end, nom, score, brin, [signal].
    La probabilite d'arret croit avec le signal, bornee par force_max.
    """
    import pandas as pd
    d = pd.read_csv(chemin, sep="\t", header=None, comment="#")
    d = d[(d[0] == chrom) & (d[2] > debut) & (d[1] < fin)]
    bp = np.zeros(n_billes)
    bm = np.zeros(n_billes)
    if d.empty:
        return bp, bm

    sig = (d[colonne_signal].to_numpy(dtype=float)
           if d.shape[1] > colonne_signal else d[4].to_numpy(dtype=float))
    sig = sig / (np.percentile(sig, 95) + 1e-9)
    for (s, e, brin), v in zip(d[[1, 2, 5]].to_numpy(), sig):
        i = int(np.clip(((s + e) / 2 - debut) // taille_bin, 0, n_billes - 1))
        f = force_max * float(np.clip(v, 0.0, 1.0))
        if brin == "+":
            bp[i] = max(bp[i], f)
        elif brin == "-":
            bm[i] = max(bm[i], f)
        else:                       # orientation inconnue : barriere des deux cotes
            bp[i] = max(bp[i], f * 0.5)
            bm[i] = max(bm[i], f * 0.5)
    return bp, bm


def liens_par_instant(traj):
    """Convertit la trajectoire en listes de couples de billes."""
    return [[(int(a), int(b)) for a, b in inst] for inst in traj]
