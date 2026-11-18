#!/usr/bin/env python3
"""
evaluer_4c.py - agreement between simulated and measured contacts once the
distance decay is removed, for every (condition, parameter set) of path A.

simpol comparer.py correlates virtual 4C profiles normalised to their sum:
the decay with genomic distance dominates that correlation, so a featureless
polymer already scores high. Here both matrices are turned into
observed / expected (each diagonal divided by its mean over the analysed
window), and two things are reported:

  r_oe_<anchor>   Pearson correlation of the observed/expected profile from
                  each anchor, simulated against measured, over the analysed
                  window (simulated flanks left out), anchor bins +/- 2 left
                  out; r_oe_moyen is their mean
  oe_<a>__<b>     observed/expected at each pair of anchors, simulated; the
                  first row (condition MESURE) gives the measured values

Compare SAUVAGE with the ALEATOIRE<k> rows (CTCF sites at random positions):
the CTCF sites explain the locus only if SAUVAGE (or UNIFORME) beats the
random placements on r_oe and reproduces the anchor-pair values.

Usage (OpenMM environment or any with numpy, pandas):
  SIMPOL_LOCI=polymere/abat/loci.json python polymere/evaluer_4c.py \
      --dir polymere/abat --locus abat
Output: <dir>/evaluation_4c.tsv, sorted by r_oe_moyen.
"""

import argparse
import itertools
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "simpol"))
import config as cfg                                   # noqa: E402


def charger_matrice(chemin, debut, fin, taille_bin):
    """matrice mesurée sur la grille simulée (comme simpol une_trajectoire.py)"""
    import pandas as pd
    d = pd.read_csv(chemin, sep="\t", comment="#")
    n = int(np.ceil((fin - debut) / taille_bin))
    M = np.zeros((n, n))
    vu = np.zeros((n, n), dtype=bool)

    def plage(s, e):
        return np.arange(max(0, int((s - debut) // taille_bin)),
                         min(n - 1, int((e - 1 - debut) // taille_bin)) + 1)
    for s1, e1, s2, e2, N in d[["startA", "endA", "startB", "endB", "N"]].to_numpy():
        i, j = plage(s1, e1), plage(s2, e2)
        if i.size == 0 or j.size == 0:
            continue
        I, J = np.meshgrid(i, j, indexing="ij")
        np.add.at(M, (I, J), N / (i.size * j.size))
        np.add.at(M, (J, I), N / (i.size * j.size))
        vu[I, J] = vu[J, I] = True
    return M, vu


def matrice_simulee(conf, seuil=2.0, bloc=256):
    """fréquence de contact (simpol comparer.py, sans importer OpenMM)"""
    n = conf.shape[1]
    M = np.zeros((n, n))
    s2 = float(seuil) ** 2
    for k in range(0, len(conf), bloc):
        P = conf[k:k + bloc].astype(np.float64)
        q = (P * P).sum(2)
        d2 = q[:, :, None] + q[:, None, :] - 2.0 * (P @ P.transpose(0, 2, 1))
        M += (d2 < s2).sum(0)
    np.fill_diagonal(M, 0.0)
    return M / len(conf)


def oe(M, masque, lo, hi):
    """observé / attendu : chaque diagonale divisée par sa moyenne dans [lo, hi)"""
    n = len(M)
    out = np.full((n, n), np.nan)
    for k in range(1, hi - lo):
        i = np.arange(lo, hi - k)
        v = M[i, i + k]
        ok = masque[i, i + k] & np.isfinite(v)
        if ok.sum() < 3:
            continue
        mu = v[ok].mean()
        if mu > 0:
            out[i, i + k] = np.where(ok, v / mu, np.nan)
            out[i + k, i] = out[i, i + k]
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", required=True)
    p.add_argument("--locus", required=True)
    p.add_argument("--seuil", type=float, default=2.0)
    p.add_argument("--max-conformations", type=int, default=6000, dest="max_conf")
    p.add_argument("--conditions", default="", help="liste séparée par des virgules (défaut : toutes)")
    a = p.parse_args()
    d = Path(a.dir)
    c = cfg.LOCI[a.locus]
    n, d0, f0 = cfg.geometrie(a.locus)
    tb = c["taille_bin"]
    lo, hi = int((c["debut"] - d0) // tb), int(np.ceil((c["fin"] - d0) / tb))
    idx = cfg.indices_ancres(a.locus)
    noms = list(idx)
    Mes, vu = charger_matrice(d / "matrice.tsv", d0, f0, tb)
    OEm = oe(Mes, vu, lo, hi)
    paires = list(itertools.combinations(noms, 2))

    def profils(O):
        out = {}
        for nm, i in idx.items():
            v = O[i, lo:hi].copy()
            k = np.arange(lo, hi)
            v[np.abs(k - i) <= 2] = np.nan
            out[nm] = v
        return out
    pm = profils(OEm)
    groupes = defaultdict(list)
    voulu = set(a.conditions.split(",")) if a.conditions else None
    for f in sorted((d / "resultats").glob(f"{a.locus}_*_A_jeu*_rep*.npy")):
        m = re.match(rf"{re.escape(a.locus)}_(.+)_A_jeu(\d+)_rep\d+\.npy", f.name)
        if m and (voulu is None or m.group(1) in voulu):
            groupes[(m.group(1), int(m.group(2)))].append(f)
    if not groupes:
        sys.exit(f"aucune trajectoire dans {d / 'resultats'}")
    rows = [dict(condition="MESURE", jeu=".", replicats=".", r_oe_moyen=float("nan"),
                 **{f"r_oe_{nm}": float("nan") for nm in noms},
                 **{f"oe_{x}__{y}": OEm[idx[x], idx[y]] for x, y in paires})]
    for (cond, jeu), fs in sorted(groupes.items()):
        conf = np.concatenate([np.load(f) for f in fs])
        if a.max_conf and len(conf) > a.max_conf:
            conf = conf[np.linspace(0, len(conf) - 1, a.max_conf).astype(int)]
        OEs = oe(matrice_simulee(conf, a.seuil), np.ones((n, n), dtype=bool), lo, hi)
        ps = profils(OEs)
        r = {}
        for nm in noms:
            x, y = pm[nm], ps[nm]
            ok = np.isfinite(x) & np.isfinite(y)
            r[nm] = float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() >= 8 and x[ok].std() > 0 \
                and y[ok].std() > 0 else float("nan")
        rows.append(dict(condition=cond, jeu=jeu, replicats=len(fs),
                         r_oe_moyen=float(np.nanmean(list(r.values()))),
                         **{f"r_oe_{nm}": r[nm] for nm in noms},
                         **{f"oe_{x}__{y}": OEs[idx[x], idx[y]] for x, y in paires}))
    rows = rows[:1] + sorted(rows[1:], key=lambda r: -np.nan_to_num(r["r_oe_moyen"], nan=-9))
    cols = list(rows[0])
    with open(d / "evaluation_4c.tsv", "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            fo.write("\t".join(f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k]) for k in cols) + "\n")
    cle = [f"oe_{x}__{y}" for x, y in paires]
    print(f"{'condition':14s} {'jeu':>4s} {'rép':>4s} {'r_oe':>6s}  " + "  ".join(k[3:] for k in cle))
    for r in rows[:25]:
        print(f"{r['condition']:14s} {str(r['jeu']):>4s} {str(r['replicats']):>4s} "
              f"{r['r_oe_moyen']:6.3f}  " + "  ".join(f"{r[k]:.2f}" for k in cle))
    al = [r["r_oe_moyen"] for r in rows if str(r["condition"]).startswith("ALEATOIRE")]
    if al:
        print(f"\nALEATOIRE : r_oe de {min(al):.3f} à {max(al):.3f} ({len(al)} groupes)")
    print(f"-> {d / 'evaluation_4c.tsv'}")


if __name__ == "__main__":
    main()
