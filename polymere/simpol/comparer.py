#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Confrontation des deux voies a la mesure, et l'une a l'autre.

COMMENT ON COMPARE, ET POURQUOI AINSI
  La simulation produit une matrice tout-contre-tout, la mesure est une
  capture. L'efficacite de capture d'un fragment est un facteur
  multiplicatif sur toute sa ligne de la matrice : normaliser chaque ligne
  par sa propre somme le fait disparaitre exactement, sans avoir a
  l'estimer. C'est pourquoi la comparaison principale porte sur des profils
  de type 4C virtuel, pris aux ancres annotees, et non sur les matrices
  entieres. C'est la pratique du modele HiP-HoP, qui superpose profil
  simule et profil mesure point de vue par point de vue.

  La comparaison matricielle n'est gardee qu'en controle, sous forme de
  correlation stratifiee par distance : la decroissance avec la distance
  genomique domine toute correlation brute et la rendrait sans portee.

CE QUI DISTINGUE LES DEUX VOIES
  Les contacts deux a deux sont probablement reproduits par l'une comme
  par l'autre. La ou elles different, c'est sur la simultaneite : la
  fraction des conformations ou toutes les ancres du hub sont ensemble au
  meme instant. C'est une prediction, et elle se teste par Tri-C ou Pore-C.

USAGE
  python3 comparer.py --locus tyrp1 --condition MEL \\
      --matrice TYRP1_matrice_MEL_UNSTIM.tsv --resultats resultats/
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config as cfg
from une_trajectoire import charger_matrice

RED, BLUE, GOLD, VERT = "#B2182B", "#2166AC", "#C78A20", "#3F6B57"


# ======================================================================
def matrice_simulee(conformations, seuil=2.0, bloc=256):
    """Frequence de contact, par seuil sur la distance.

    Les distances passent par |x|^2 + |y|^2 - 2 x.y, ce qui remplace une
    difference diffusee de taille (n, n, 3) par un produit matriciel : a
    790 billes, 2,8 fois plus rapide. Le calcul est fait en float64 ; en
    float32 l'annulation faisait basculer une conformation sur mille
    exactement au seuil, et le resultat dependait alors du chemin
    arithmetique.
    """
    n = conformations.shape[1]
    M = np.zeros((n, n))
    s2 = float(seuil) ** 2
    for k in range(0, len(conformations), bloc):
        P = conformations[k:k + bloc].astype(np.float64)
        q = (P * P).sum(2)
        d2 = q[:, :, None] + q[:, None, :] - 2.0 * (P @ P.transpose(0, 2, 1))
        M += (d2 < s2).sum(0)
    np.fill_diagonal(M, 0.0)
    return M / len(conformations)


def sous_echantillon(conf, maximum):
    """Garde au plus `maximum` conformations, regulierement espacees.

    Les conformations d'une meme trajectoire sont correlees entre blocs
    voisins : l'effectif utile est tres inferieur a l'effectif brut, et
    40 000 conformations ne font pas une frequence de contact plus juste
    que 6 000. Le pas regulier fait contribuer tous les replicats.
    """
    if maximum and len(conf) > maximum:
        return conf[np.linspace(0, len(conf) - 1, maximum).astype(int)]
    return conf


def profils(M, idx, d_min=3):
    """4C virtuel a chaque ancre, chaque profil normalise a somme un.

    La normalisation par ligne est ce qui annule l'efficacite de capture
    propre a chaque fragment.
    """
    n = len(M)
    k = np.abs(np.arange(n)[None, :] - np.arange(n)[:, None])
    out = {}
    for nom, i in idx.items():
        v = np.where(k[i] >= d_min, M[i], np.nan)
        s = np.nansum(v)
        out[nom] = v / s if s > 0 else v
    return out


def correlation_stratifiee(A, B, masque=None, k_max=None):
    """Correlation moyenne par strate de distance, ponderee a la maniere
    de HiCRep : chaque strate pese son effectif et sa variabilite."""
    n = len(A)
    k_max = k_max or n - 1
    num = den = 0.0
    for k in range(1, k_max):
        i = np.arange(n - k)
        a, b = A[i, i + k], B[i, i + k]
        if masque is not None:
            m = masque[i, i + k]
            a, b = a[m], b[m]
        bon = np.isfinite(a) & np.isfinite(b)
        a, b = a[bon], b[bon]
        if len(a) < 5 or a.std() == 0 or b.std() == 0:
            continue
        w = len(a) * np.sqrt(a.std() * b.std())
        num += w * pearsonr(a, b).statistic
        den += w
    return num / den if den > 0 else np.nan


def simultaneite(conformations, idx, noms, seuil=2.5):
    """Fraction des conformations ou toutes les ancres listees sont
    simultanement a moins de `seuil` les unes des autres, et fraction des
    paires prises separement. L'ecart entre les deux est la cooperativite."""
    if len(noms) < 2:
        return np.nan, np.nan
    j = [idx[n] for n in noms]
    P = conformations[:, j, :]
    d = np.sqrt(((P[:, :, None, :] - P[:, None, :, :]) ** 2).sum(3))
    m = ~np.eye(len(j), dtype=bool)
    ensemble = float((d[:, m].max(1) < seuil).mean())
    paires = float((d[:, m] < seuil).mean())
    return ensemble, paires


# ======================================================================
def charger_resultats(dossier, locus, condition):
    """Regroupe les replicats par voie et par jeu."""
    grp = {}
    for j in sorted(Path(dossier).glob(f"{locus}_{condition}_*.json")):
        meta = json.loads(j.read_text())
        npy = j.with_suffix(".npy")
        if not npy.exists():
            continue
        grp.setdefault((meta["voie"], meta["jeu"]), []).append((meta, npy))
    return grp


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--locus", default="tyrp1", choices=sorted(cfg.LOCI))
    p.add_argument("--bin", type=int, dest="taille_bin_cli",
                   help="change la résolution du locus")
    p.add_argument("--marge", type=int,
                   help="change le flanc simulé de part et d'autre")
    p.add_argument("--condition", required=True)
    p.add_argument("--matrice", required=True)
    p.add_argument("--resultats", default="resultats")
    p.add_argument("--seuil-contact", type=float, default=2.0,
                   dest="seuil_contact")
    p.add_argument("--seuil-hub", type=float, default=2.5, dest="seuil_hub")
    p.add_argument("--max-conformations", type=int, default=6000,
                   dest="max_conf",
                   help="sous-échantillonnage régulier ; 0 pour tout garder")
    p.add_argument("--sortie", default="comparaison")
    a = p.parse_args()
    cfg.surcharger(a.locus, a.taille_bin_cli, a.marge)

    c = cfg.LOCI[a.locus]
    n, d0, f0 = cfg.geometrie(a.locus)
    idx = cfg.indices_ancres(a.locus)

    Mes, vu = charger_matrice(a.matrice, c["chrom"], d0, f0, c["taille_bin"])
    prof_mes = profils(np.where(vu, Mes, np.nan), idx)
    print(f"  mesure : {int(vu.sum()/2):,} cases renseignées sur "
          f"{n*(n-1)//2:,}")

    grp = charger_resultats(a.resultats, a.locus, a.condition)
    if not grp:
        raise SystemExit(f"aucun résultat dans {a.resultats}")
    print(f"  {len(grp)} couples (voie, jeu) trouvés")

    lignes, garde = [], {}
    for (voie, jeu), items in sorted(grp.items()):
        conf = sous_echantillon(
            np.concatenate([np.load(f) for _, f in items]), a.max_conf)
        Msim = matrice_simulee(conf, a.seuil_contact)
        prof_sim = profils(Msim, idx)

        # Une ancre dont le profil n'a pas assez de cases renseignees est
        # ecartee ; si aucune ne passe, la comparaison n'a pas de sens et
        # on le dit au lieu de renvoyer un nombre.
        r = []
        for nom in idx:
            x, y = prof_mes[nom], prof_sim[nom]
            bon = np.isfinite(x) & np.isfinite(y)
            if bon.sum() >= 8 and np.nanstd(x[bon]) > 0:
                r.append(pearsonr(x[bon], y[bon]).statistic)
        if not r:
            print(f"    voie {voie} jeu {jeu} : aucune ancre exploitable "
                  f"(résolution trop grossière ou couverture trop faible)")
            r = [np.nan]
        scc = correlation_stratifiee(Mes, Msim, vu)
        ens, par = simultaneite(conf, idx, c["hub"], a.seuil_hub)
        ensf, parf = simultaneite(conf, idx, c["hub_fibro"], a.seuil_hub)

        lignes.append(dict(voie=voie, jeu=jeu, replicats=len(items),
                           conformations=len(conf),
                           profils_pearson_moyen=float(np.mean(r)),
                           profils_pearson_min=float(np.min(r)),
                           matrice_scc=float(scc),
                           hub_simultane=ens, hub_paires=par,
                           hub_fibro_simultane=ensf,
                           parametres=json.dumps(items[0][0]["parametres"])))
        garde[(voie, jeu)] = (Msim, prof_sim, conf)

    T = pd.DataFrame(lignes).sort_values("profils_pearson_moyen",
                                         ascending=False)
    T.to_csv(f"{a.sortie}_{a.condition}.tsv", sep="\t", index=False)
    print(f"\n  meilleurs jeux par voie :")
    for voie in sorted(T.voie.unique()):
        b = T[T.voie == voie].iloc[0]
        print(f"    voie {voie}  jeu {b.jeu:3d}  profils r={b.profils_pearson_moyen:.3f}"
              f"  (min {b.profils_pearson_min:.3f})  matrice SCC={b.matrice_scc:.3f}")
        print(f"              hub simultané {b.hub_simultane:.3f} contre "
              f"{b.hub_paires:.3f} pour les paires séparées")
        print(f"              {b.parametres}")

    # ------------------------------------------------------------------
    meilleurs = {v: T[T.voie == v].iloc[0] for v in sorted(T.voie.unique())}
    coul = {"A": GOLD, "B": VERT, "AB": "#7A4FA3"}
    noms = list(idx)
    fig, axes = plt.subplots(len(noms), 1, figsize=(7.4, 1.25 * len(noms)),
                             sharex=True)
    axes = np.atleast_1d(axes)
    x = (d0 + (np.arange(n) + 0.5) * c["taille_bin"]) / 1e6
    for ax, nom in zip(axes, noms):
        ax.plot(x, prof_mes[nom], color="#1A1A1A", lw=1.1, label="mesure")
        for v, b in meilleurs.items():
            _, ps, _ = garde[(v, int(b.jeu))]
            ax.plot(x, ps[nom], color=coul.get(v, "#888"), lw=1.0,
                    alpha=0.85, label=f"voie {v}")
        ax.axvline((d0 + idx[nom] * c["taille_bin"]) / 1e6, color="#BBB",
                   lw=0.7, ls=(0, (2, 2)))
        ax.set_ylabel(nom.replace("_", "\n"), fontsize=6.2, rotation=0,
                      ha="right", va="center")
        ax.set_yticks([])
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(labelsize=7)
    axes[0].legend(frameon=False, fontsize=6.8, ncol=4, loc="upper left")
    axes[-1].set_xlabel(f"{c['chrom']}  (Mb)", fontsize=8)
    fig.suptitle(f"{a.locus.upper()} {a.condition} — virtual 4C profiles, "
                 "each normalised to its own sum", fontsize=9, x=0.012,
                 ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(f"{a.sortie}_{a.condition}_profils.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print(f"\n  {a.sortie}_{a.condition}.tsv")
    print(f"  {a.sortie}_{a.condition}_profils.png")


if __name__ == "__main__":
    main()
