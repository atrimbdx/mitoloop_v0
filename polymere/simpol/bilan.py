#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bilan : les deux voies face a la mesure, et face l'une a l'autre.

LECTURE DU TABLEAU, ET UNE ASYMETRIE A NE PAS OUBLIER
  La voie B est une methode inverse : elle voit la matrice mesuree et
  ajuste ses liants dessus. La voie A est predictive : son unique entree
  est la position et l'orientation des sites CTCF, elle ne voit jamais la
  matrice. Comparer leurs qualites d'ajustement brutes est donc
  defavorable a la voie A par construction, et le dire fait partie du
  resultat.

  La question qui a un sens est ailleurs : la voie A atteint-elle, avec
  CTCF pour seule entree, un accord comparable ? Et surtout, les deux
  voies predisent-elles la meme chose la ou les donnees ne tranchent pas,
  c'est-a-dire sur la simultaneite des contacts du hub.

SORTIES
  bilan.tsv   une ligne par voie et par condition, meilleur jeu retenu
  bilan.png   qualite d'ajustement, et prediction de simultaneite
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config as cfg

COUL = {"A": "#C78A20", "B": "#3F6B57", "AB": "#7A4FA3"}
NOM = {"A": "A — extrusion, CTCF ENCODE", "B": "B — liants inférés",
       "AB": "A+B — les deux"}


def nuage(ax, T, voies, conds, colonne, titre, ylabel, zero=False):
    rng = np.random.default_rng(0)
    for k, v in enumerate(voies):
        for j, c in enumerate(conds):
            s = T[(T.voie == v) & (T.condition == c)][colonne].dropna()
            if s.empty:
                continue
            x = k + (j - (len(conds) - 1) / 2) * 0.3
            ax.scatter(x + rng.normal(0, 0.03, len(s)), s, s=9,
                       color=COUL.get(v, "#888"), alpha=0.3 + 0.45 * j,
                       linewidths=0)
            ax.plot([x - 0.1, x + 0.1], [s.max()] * 2,
                    color=COUL.get(v, "#888"), lw=1.6)
    if zero:
        ax.axhline(0, color="#BBB", lw=0.7)
    ax.set_xticks(range(len(voies)))
    ax.set_xticklabels([f"voie {v}" for v in voies], fontsize=7.5)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.set_title(titre, fontsize=8, loc="left")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--locus", default="tyrp1", choices=sorted(cfg.LOCI))
    p.add_argument("--conditions", nargs="+", default=["MEL", "FIBRO"])
    p.add_argument("--prefixe", default="comparaison")
    p.add_argument("--sortie", default="bilan")
    a = p.parse_args()

    tabs = []
    for c in a.conditions:
        f = Path(f"{a.prefixe}_{c}.tsv")
        if not f.exists():
            print(f"  absent : {f}")
            continue
        t = pd.read_csv(f, sep="\t")
        t["condition"] = c
        tabs.append(t)
    if not tabs:
        raise SystemExit("aucun tableau de comparaison")
    T = pd.concat(tabs, ignore_index=True)

    best = (T.sort_values("profils_pearson_moyen", ascending=False)
              .groupby(["condition", "voie"], as_index=False).first())
    best.to_csv(f"{a.sortie}.tsv", sep="\t", index=False)

    print("\n  meilleur jeu par voie et condition :")
    for _, r in best.iterrows():
        print(f"    {r.condition:6s} voie {r.voie}  profils r="
              f"{r.profils_pearson_moyen:.3f} (min {r.profils_pearson_min:.3f})"
              f"  SCC={r.matrice_scc:.3f}")
        print(f"            hub simultané {r.hub_simultane:.3f}  "
              f"paires séparées {r.hub_paires:.3f}  coopérativité "
              f"{r.hub_simultane / max(r.hub_paires, 1e-9):.2f}")
        print(f"            {r.parametres}")

    voies = sorted(best.voie.unique())
    conds = [c for c in a.conditions if c in set(best.condition)]
    fig, axes = plt.subplots(1, 3, figsize=(8.6, 2.9))

    nuage(axes[0], T, voies, conds, "profils_pearson_moyen",
          "accord avec la mesure",
          "profils 4C virtuels\nPearson moyen aux ancres")
    nuage(axes[1], T, voies, conds, "matrice_scc", "contrôle matriciel",
          "matrice\ncorrélation stratifiée", zero=True)

    ax = axes[2]
    lar = 0.34
    for j, c in enumerate(conds):
        for k, v in enumerate(voies):
            r = best[(best.voie == v) & (best.condition == c)]
            if r.empty:
                continue
            r = r.iloc[0]
            x = j + (k - (len(voies) - 1) / 2) * lar
            ax.bar(x, r.hub_simultane, lar * 0.9, color=COUL.get(v, "#888"),
                   label=NOM.get(v, v) if j == 0 else None)
            ax.plot([x - lar * 0.45, x + lar * 0.45], [r.hub_paires] * 2,
                    color="#1A1A1A", lw=1.2)
    ax.set_xticks(range(len(conds)))
    ax.set_xticklabels(conds, fontsize=7.5)
    ax.set_ylabel("fraction des conformations\nà hub complet", fontsize=7)
    ax.set_title("prédiction discriminante", fontsize=8, loc="left")
    ax.set_ylim(0, max(best.hub_paires.max(), best.hub_simultane.max()) * 1.6)
    ax.legend(frameon=False, fontsize=6.2, loc="upper center",
              bbox_to_anchor=(0.5, 1.0), handlelength=1.1)
    ax.annotate("trait noir : paires prises séparément", xy=(0.5, 0.72),
                xycoords="axes fraction", ha="center", fontsize=6.0,
                color="#4A4A4A")

    for ax in axes:
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(labelsize=7)
    fig.suptitle(f"{a.locus.upper()} — voie A contre voie B", fontsize=9,
                 x=0.012, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(f"{a.sortie}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  {a.sortie}.tsv\n  {a.sortie}.png")


if __name__ == "__main__":
    main()
