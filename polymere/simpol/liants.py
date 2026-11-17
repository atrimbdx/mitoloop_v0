#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Voie B — inference des liants, dans l'esprit du modele a cordes et liants.

Chaque bille porte ou non l'affinite pour un type de liant donne. Deux
billes qui partagent un type sont ponte par un liant de ce type, ce qui
les rapproche. La matrice de contacts mesuree suffit a inferer ces
affectations : aucune donnee externe n'est requise, et aucun mecanisme
n'est suppose.

SUBSTITUT ANALYTIQUE
  L'inference travaille sur un substitut de la probabilite de contact :

      log2(observe / attendu)_ij  ~  somme_c  w_c * s_ic * s_jc

  ou s est binaire et w positif. Le recuit simule porte sur s, et w est
  reestime par moindres carres non negatifs a chaque palier. Ce substitut
  remplace la simulation a l'interieur de la boucle d'optimisation, qui
  serait sinon hors de portee ; la verification se fait ensuite par
  simulation tridimensionnelle reelle, dans polymere.py.

  C'est donc une simplification de PRISMR, et le texte d'un article doit
  le dire : les liants sont inferes par factorisation binaire de la
  matrice observe sur attendu, puis valides par dynamique.
"""
import numpy as np
from scipy.optimize import nnls


def observe_sur_attendu(M, masque_mesure=None):
    """log2 du rapport a la moyenne de la diagonale correspondante."""
    n = len(M)
    L = np.full((n, n), np.nan)
    for k in range(1, n):
        i = np.arange(n - k)
        v = M[i, i + k].astype(float)
        bon = np.isfinite(v) & (v > 0)
        if bon.sum() < 3:
            continue
        mu = v[bon].mean()
        r = np.where(bon, v / mu, np.nan)
        L[i, i + k] = L[i + k, i] = np.log2(np.maximum(r, 1e-3))
    if masque_mesure is not None:
        L[~masque_mesure] = np.nan
    return L


def _cout(L, M, P):
    d = np.where(M, L - P, 0.0)
    return float((d * d).sum())


def inferer(L, n_types=8, d_min=3, n_paliers=40, n_balayages=12,
            temp_debut=None, temp_fin=None, graine=0, verbeux=True):
    """Recuit simule sur l'affectation des billes aux types de liants.

    Renvoie s (n x n_types, binaire), w (n_types), et le cout final.
    """
    rng = np.random.default_rng(graine)
    n = len(L)
    k = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    M = np.isfinite(L) & (k >= d_min)
    Lz = np.where(M, L, 0.0)
    Mf = M.astype(float)

    s = (rng.random((n, n_types)) < 0.25).astype(float)
    ech0 = np.nanstd(L[np.isfinite(L)])
    plancher = 0.08 * ech0
    w = np.full(n_types, 0.5 * ech0)

    def refit(s, w):
        A = np.zeros((n_types, n_types))
        b = np.zeros(n_types)
        for c in range(n_types):
            b[c] = s[:, c] @ (Mf * Lz) @ s[:, c]
            for d in range(c, n_types):
                u = s[:, c] * s[:, d]
                A[c, d] = A[d, c] = u @ Mf @ u
        try:
            x, _ = nnls(A + 1e-6 * np.eye(n_types), b)
        except Exception:
            x = w
        # Plancher pendant le recuit : l'ajustement non negatif annule les
        # types faibles des le premier palier, et un type a poids nul ne
        # peut plus etre explore puisque ses basculements ne changent rien.
        # Ils ne sont elagues qu'a la fin.
        return np.maximum(x, plancher)

    w = refit(s, w)
    P = (s * w) @ s.T
    c0 = _cout(Lz, M, P)
    ech = np.sqrt(c0 / max(M.sum(), 1))
    T0 = temp_debut if temp_debut else ech * 2.0
    T1 = temp_fin if temp_fin else ech * 0.02

    for pal in range(n_paliers):
        T = T0 * (T1 / T0) ** (pal / max(n_paliers - 1, 1))
        for _ in range(n_balayages):
            for i in rng.permutation(n):
                c = int(rng.integers(0, n_types))
                signe = 1.0 - 2.0 * s[i, c]
                dP = signe * w[c] * s[:, c]          # variation de la ligne i
                mi = M[i]
                r = Lz[i] - P[i]
                # la cellule (i,i) n'est pas dans le masque, pas de
                # correction de double comptage a faire
                d = float((mi * (-2.0 * r * dP + dP * dP)).sum()) * 2.0
                if d <= 0 or rng.random() < np.exp(-d / max(T, 1e-12)):
                    s[i, c] += signe
                    P[i] += dP
                    P[:, i] += dP
                    P[i, i] += signe * w[c] * signe  # terme diagonal
        w = refit(s, w)
        P = (s * w) @ s.T
        if verbeux and pal % 8 == 0:
            print(f"    palier {pal:2d}  T={T:.3f}  coût={_cout(Lz, M, P):.1f}"
                  f"  types actifs={(w > 1e-3).sum()}")

    # elagage final : un type dont le poids reste au plancher n'a pas
    # trouve de structure a porter
    w = refit(s, w)
    garde = (w > 1.5 * plancher) & (s.sum(0) > 2)
    if not garde.any():
        garde = w == w.max()
    s, w = s[:, garde], w[garde]
    P = (s * w) @ s.T
    return s, w, _cout(Lz, M, P), M


def qualite(L, s, w, M):
    """Correlation entre le substitut ajuste et la matrice mesuree."""
    from scipy.stats import pearsonr, spearmanr
    P = (s * w) @ s.T
    a, b = L[M], P[M]
    return pearsonr(a, b).statistic, spearmanr(a, b).statistic
