#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trois controles que la comparaison par profils ne fait pas.

1. EXISTENCE DES ELEMENTS ANNOTES
   Les elements ont ete annotes a la main. Aucune des deux voies ne les
   recoit en entree : une_trajectoire.py ne lit jamais `ancres`. La voie A
   ne voit que les sites CTCF, la voie B ne voit que la matrice. On peut
   donc demander au modele s'il met ces positions en contact plus souvent
   qu'attendu a la meme distance genomique.

   La portee de la reponse n'est pas la meme selon la voie. La voie A est
   aveugle a la matrice : si elle place un contact la ou vous avez annote
   un element, c'est une prediction independante. La voie B a vu la
   matrice, donc le signal qui a conduit a l'annotation : elle ne peut pas
   confirmer l'existence du contact, seulement dire qu'un type de liant est
   necessaire a cette position pour rendre compte des donnees.

2. ABSOLU CONTRE DIFFERENTIEL
   Les elements viennent d'un seuil sur le rapport melanocyte sur
   fibroblaste. Ils marquent donc la ou les deux etats different, et non la
   ou l'architecture est forte dans l'un ou l'autre. Les simulations, elles,
   sont menees par condition et modelisent l'architecture absolue de chaque
   etat. Evaluer cette architecture absolue aux seules positions choisies
   sur un critere differentiel introduit un biais de selection dans
   l'evaluation.

   Deux corrections. Un second jeu de points de vue est choisi sur la force
   absolue des contacts dans chaque condition, independamment du
   differentiel. Et le differentiel simule est confronte au differentiel
   mesure, ce qui teste directement ce dont votre annotation parle.

3. OU LE MODELE SE TROMPE
   Les cellules les plus mal reproduites sont listees, pour voir si l'ecart
   se concentre quelque part.

USAGE
  python3 verifier.py --locus tyrp1 \\
      --mel TYRP1_matrice_MEL_UNSTIM.tsv --fibro TYRP1_matrice_FIBRO.tsv \\
      --resultats resultats/
"""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config as cfg
from comparer import (matrice_simulee, charger_resultats, profils,
                      sous_echantillon)
from une_trajectoire import charger_matrice
from liants import observe_sur_attendu

COUL = {"A": "#C78A20", "B": "#3F6B57", "AB": "#7A4FA3"}


# ======================================================================
def enrichissement(P, i, j, exclure, tol=0.12, n_min=30):
    """Contact simule en (i,j) rapporte aux paires de meme portee.

    Renvoie le rapport au mediane des temoins et le rang centile de la
    paire dans la distribution des temoins. Les positions annotees sont
    exclues des temoins, sans quoi elles se compareraient a elles-memes.
    """
    n = len(P)
    k = abs(i - j)
    lo, hi = int(k * (1 - tol)), int(np.ceil(k * (1 + tol)))
    temoins = []
    for a in range(n):
        for b in range(a + lo, min(a + hi + 1, n)):
            if a in exclure or b in exclure:
                continue
            temoins.append(P[a, b])
    temoins = np.array(temoins, dtype=float)
    temoins = temoins[np.isfinite(temoins)]
    if len(temoins) < n_min:
        return np.nan, np.nan, len(temoins)
    med = np.median(temoins)
    rapport = P[i, j] / med if med > 0 else np.nan
    centile = float((temoins < P[i, j]).mean())
    return rapport, centile, len(temoins)


def ancres_absolues(M, vu, n_top=6, d_min=6, ecart_min=4):
    """Points de vue choisis sur la force ABSOLUE du contact, sans aucun
    recours au differentiel. On prend les billes dont le signal observe
    sur attendu est le plus fort, en les espacant."""
    L = observe_sur_attendu(M, vu)
    n = len(L)
    k = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    L = np.where(k >= d_min, L, np.nan)
    score = np.nansum(np.where(np.isfinite(L) & (L > 0), L, 0.0), axis=1)
    out = []
    for i in np.argsort(score)[::-1]:
        if all(abs(i - j) >= ecart_min for j in out):
            out.append(int(i))
        if len(out) >= n_top:
            break
    return {f"abs_{p}": i for p, i in enumerate(sorted(out), start=1)}


def nul_differentiel(conf_a, conf_b, Dmes, masque, seuil, n_tirages=24,
                     graine=0):
    """Distribution de la correlation sous permutation des etiquettes.

    POURQUOI CE NUL EST INDISPENSABLE
      Mesure faite sur ce paquet : trois comparaisons structurellement
      nulles, ou les deux conditions avaient des entrees identiques et ne
      differaient que par la graine, ont donne des correlations de +0,490,
      +0,215 et -0,325. L'etendue du nul couvrait donc 0,8 a elle seule.
      Une correlation de differentiel annoncee sans son nul ne veut rien
      dire : elle est du meme ordre que le bruit de tirage.

      Le nul melange les conformations des deux conditions et les reseparte
      au hasard en deux groupes de meme effectif. Les effectifs, la
      longueur des trajectoires et la couverture sont donc conserves ; seul
      le contraste entre conditions est detruit.
    """
    rng = np.random.default_rng(graine)
    tout = np.concatenate([conf_a, conf_b])
    na = len(conf_a)
    out = []
    for _ in range(n_tirages):
        o = rng.permutation(len(tout))
        D = differentiel(matrice_simulee(tout[o[:na]], seuil),
                         matrice_simulee(tout[o[na:]], seuil))
        m = masque & np.isfinite(D)
        if m.sum() > 50:
            out.append(spearmanr(D[m], Dmes[m]).statistic)
    return np.array(out)


def differentiel(A, B, plancher=0.35):
    """log2 du rapport entre deux cartes, chacune ramenee a l'attendu de sa
    propre diagonale. Une pseudo-observation proportionnelle au niveau
    attendu amortit les cases peu couvertes."""
    n = len(A)
    out = np.full((n, n), np.nan)
    for k in range(1, n):
        i = np.arange(n - k)
        a, b = A[i, i + k].astype(float), B[i, i + k].astype(float)
        bon = np.isfinite(a) & np.isfinite(b)
        if bon.sum() < 5:
            continue
        ma = a[bon].mean() or 1.0
        mb = b[bon].mean() or 1.0
        p = plancher
        v = np.log2((a / ma + p) / (b / mb + p))
        out[i, i + k] = out[i + k, i] = np.where(bon, v, np.nan)
    return out


# ======================================================================
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--locus", default="tyrp1", choices=sorted(cfg.LOCI))
    p.add_argument("--bin", type=int, dest="taille_bin_cli")
    p.add_argument("--marge", type=int)
    p.add_argument("--mel", required=True)
    p.add_argument("--fibro", required=True)
    p.add_argument("--resultats", default="resultats")
    p.add_argument("--seuil-contact", type=float, default=2.0,
                   dest="seuil_contact")
    p.add_argument("--max-conformations", type=int, default=6000,
                   dest="max_conf",
                   help="sous-échantillonnage régulier ; 0 pour tout garder")
    p.add_argument("--tirages-nul", type=int, default=24, dest="tirages_nul",
                   help="permutations pour le nul du différentiel")
    p.add_argument("--jeu-mel", type=int, default=None, dest="jeu_mel",
                   help="jeu retenu ; par défaut le meilleur de bilan.tsv")
    p.add_argument("--jeu-fibro", type=int, default=None, dest="jeu_fibro")
    p.add_argument("--bilan", default="bilan.tsv")
    p.add_argument("--sortie", default="verification")
    a = p.parse_args()

    cfg.surcharger(a.locus, a.taille_bin_cli, a.marge)
    c = cfg.LOCI[a.locus]
    n, d0, f0 = cfg.geometrie(a.locus)
    idx = cfg.indices_ancres(a.locus)
    noms = list(idx)

    mes = {}
    for cond, chemin in (("MEL", a.mel), ("FIBRO", a.fibro)):
        M, vu = charger_matrice(chemin, c["chrom"], d0, f0, c["taille_bin"])
        mes[cond] = (M, vu)
    Dmes = differentiel(np.where(mes["MEL"][1], mes["MEL"][0], np.nan),
                        np.where(mes["FIBRO"][1], mes["FIBRO"][0], np.nan))

    # jeu retenu par voie : celui du bilan, sinon le premier trouve
    choix = {}
    try:
        b = pd.read_csv(a.bilan, sep="\t")
        for _, r in b.iterrows():
            choix[(r.condition, r.voie)] = int(r.jeu)
    except (FileNotFoundError, pd.errors.EmptyDataError, KeyError, OSError):
        print(f"  pas de {a.bilan} exploitable : tous les jeux sont gardés")

    sim = {}
    for cond in ("MEL", "FIBRO"):
        grp = charger_resultats(a.resultats, a.locus, cond)
        for (voie, jeu), items in grp.items():
            if choix and choix.get((cond, voie), jeu) != jeu:
                continue
            conf = sous_echantillon(
                np.concatenate([np.load(f) for _, f in items]), a.max_conf)
            sim[(cond, voie)] = (matrice_simulee(conf, a.seuil_contact), conf)
    if not sim:
        raise SystemExit(f"aucun résultat dans {a.resultats}")
    voies = sorted({v for _, v in sim})
    print(f"  voies trouvées : {voies}")

    # ------------------------------------------------------------------
    # 1. existence des elements annotes
    exclure = set(idx.values())
    lignes = []
    for cond in ("MEL", "FIBRO"):
        for voie in voies:
            if (cond, voie) not in sim:
                continue
            P, _ = sim[(cond, voie)]
            for na, nb in itertools.combinations(noms, 2):
                i, j = idx[na], idx[nb]
                if abs(i - j) < 4:
                    continue
                rap, cent, nt = enrichissement(P, i, j, exclure)
                lignes.append(dict(controle="éléments annotés", condition=cond,
                                   voie=voie, paire=f"{na} – {nb}",
                                   portee_billes=abs(i - j),
                                   rapport_aux_temoins=rap,
                                   centile=cent, n_temoins=nt,
                                   aveugle_a_la_matrice=(voie == "A")))
    E = pd.DataFrame(lignes)

    # ------------------------------------------------------------------
    # 2. points de vue choisis sur la force absolue, par condition
    abs_lignes = []
    for cond in ("MEL", "FIBRO"):
        M, vu = mes[cond]
        ia = ancres_absolues(M, vu)
        pm = profils(np.where(vu, M, np.nan), ia)
        for voie in voies:
            if (cond, voie) not in sim:
                continue
            P, _ = sim[(cond, voie)]
            ps = profils(P, ia)
            r = []
            for nom in ia:
                x, y = pm[nom], ps[nom]
                bon = np.isfinite(x) & np.isfinite(y)
                if bon.sum() >= 8 and np.nanstd(x[bon]) > 0:
                    r.append(pearsonr(x[bon], y[bon]).statistic)
            abs_lignes.append(dict(condition=cond, voie=voie,
                                   n_ancres=len(ia),
                                   pearson_moyen=float(np.mean(r)) if r else np.nan,
                                   positions=json.dumps(sorted(ia.values()))))
    A = pd.DataFrame(abs_lignes)

    # ------------------------------------------------------------------
    # 3. differentiel simule contre differentiel mesure
    diff_lignes = []
    Dsim, nuls = {}, {}
    k = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    for voie in voies:
        if ("MEL", voie) not in sim or ("FIBRO", voie) not in sim:
            continue
        Pm, cm = sim[("MEL", voie)]
        Pf, cf = sim[("FIBRO", voie)]
        D = differentiel(Pm, Pf)
        Dsim[voie] = D
        m = np.isfinite(D) & np.isfinite(Dmes) & (k >= 4)
        rho = spearmanr(D[m], Dmes[m]).statistic if m.sum() > 50 else np.nan

        base = np.isfinite(Dmes) & (k >= 4)
        z = nul_differentiel(cm, cf, Dmes, base, a.seuil_contact,
                             a.tirages_nul)
        nuls[voie] = z
        p = (float((np.abs(z) >= abs(rho)).mean()) if len(z) else np.nan)

        # La voie A ne distingue les conditions que si elle a recu des
        # sites CTCF differents. Meme fichier des deux cotes : le
        # differentiel simule est un nul par construction, et le dire est
        # plus utile que de publier le nombre.
        emps = []
        for c_ in ("MEL", "FIBRO"):
            g = charger_resultats(a.resultats, a.locus, c_)
            for (v_, j_), it in g.items():
                if v_ == voie:
                    emps.append(it[0][0]["parametres"].get("ctcf_empreinte"))
                    break
        memes = (voie in ("A", "AB") and len(set(emps)) == 1
                 and emps[0] not in (None, "aucun"))

        ligne = dict(voie=voie, spearman_differentiel=float(rho),
                     nul_median=float(np.median(z)) if len(z) else np.nan,
                     nul_p05=float(np.quantile(z, 0.05)) if len(z) else np.nan,
                     nul_p95=float(np.quantile(z, 0.95)) if len(z) else np.nan,
                     p_permutation=p, tirages_nul=len(z),
                     entrees_identiques_entre_conditions=bool(memes),
                     cases=int(m.sum()))
        for na, nb in itertools.combinations(noms, 2):
            i, j = idx[na], idx[nb]
            if abs(i - j) < 4:
                continue
            ligne[f"mes_{na[:6]}_{nb[:6]}"] = round(float(Dmes[i, j]), 2)
            ligne[f"sim_{na[:6]}_{nb[:6]}"] = round(float(D[i, j]), 2)
        diff_lignes.append(ligne)
    D_ = pd.DataFrame(diff_lignes)

    for nom, T in (("elements", E), ("ancres_absolues", A),
                   ("differentiel", D_)):
        T.to_csv(f"{a.sortie}_{nom}.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------
    print("\n  1. éléments annotés, contact simulé rapporté aux témoins "
          "de même portée")
    for voie in voies:
        s = E[(E.voie == voie) & (E.condition == "MEL")]
        if s.empty:
            continue
        au = "aveugle à la matrice" if voie == "A" else "a vu la matrice"
        print(f"     voie {voie} ({au}) : rapport médian "
              f"{s.rapport_aux_temoins.median():.2f}, "
              f"{(s.centile > 0.9).sum()}/{len(s)} paires au-dessus du "
              f"90e centile des témoins")
    print("\n  2. points de vue choisis sur la force absolue, hors "
          "différentiel")
    for _, r in A.iterrows():
        print(f"     {r.condition:6s} voie {r.voie} : Pearson moyen "
              f"{r.pearson_moyen:.3f} sur {r.n_ancres} ancres")
    print("\n  3. différentiel simulé contre différentiel mesuré")
    for _, r in D_.iterrows():
        print(f"     voie {r.voie} : Spearman {r.spearman_differentiel:+.3f} "
              f"sur {r.cases:,} cases")
        print(f"             nul par permutation [{r.nul_p05:+.3f}, "
              f"{r.nul_p95:+.3f}], médiane {r.nul_median:+.3f}, "
              f"p = {r.p_permutation:.3f} sur {r.tirages_nul} tirages")
        if r.p_permutation > 0.05:
            print(f"             dans le bruit de tirage : ne pas rapporter "
                  f"ce nombre seul")
        if r.entrees_identiques_entre_conditions:
            print(f"             ATTENTION les deux conditions ont reçu le "
                  f"même BED CTCF : la voie {r.voie} n'a aucune entrée qui "
                  f"les distingue,")
            print(f"             son différentiel est un nul par "
                  f"construction. Fournir un BED par condition.")

    # ------------------------------------------------------------------
    if Dsim:
        fig, axes = plt.subplots(1, len(Dsim) + 2,
                                 figsize=(3.1 * (len(Dsim) + 2), 3.0))
        axes = np.atleast_1d(axes)
        base = np.isfinite(Dmes) & (k >= 4)
        for ax, voie in zip(axes, sorted(Dsim)):
            m = base & np.isfinite(Dsim[voie])
            ax.scatter(Dmes[m], Dsim[voie][m], s=3, alpha=0.2,
                       color=COUL.get(voie, "#888"), linewidths=0)
            for na, nb in itertools.combinations(noms, 2):
                i, j = idx[na], idx[nb]
                if abs(i - j) >= 4 and np.isfinite(Dsim[voie][i, j]):
                    ax.scatter(Dmes[i, j], Dsim[voie][i, j], s=26,
                               facecolor="none", edgecolor="#1A1A1A", lw=0.9)
            ax.axhline(0, color="#CCC", lw=0.6)
            ax.axvline(0, color="#CCC", lw=0.6)
            ax.set_xlabel("différentiel mesuré\nlog2 MEL / FIBRO", fontsize=7)
            ax.set_ylabel("différentiel simulé", fontsize=7)
            ax.set_title(f"voie {voie}", fontsize=8, loc="left")
        ax = axes[-2]
        for voie in sorted(Dsim):
            s = E[(E.voie == voie) & (E.condition == "MEL")]
            ax.scatter(s.portee_billes, s.rapport_aux_temoins, s=22,
                       color=COUL.get(voie, "#888"), label=f"voie {voie}")
        ax.axhline(1, color="#1A1A1A", lw=0.8, ls=(0, (3, 3)))
        ax.set_xlabel("portée de la paire (billes)", fontsize=7)
        ax.set_ylabel("contact simulé / témoins\nde même portée", fontsize=7)
        ax.set_title("éléments annotés", fontsize=8, loc="left")
        ax.legend(frameon=False, fontsize=6.6)

        # Le differentiel observe, place dans son nul. C'est la seule
        # lecture honnete de ce nombre.
        ax = axes[-1]
        for q, voie in enumerate(sorted(nuls)):
            z = nuls[voie]
            if not len(z):
                continue
            ax.scatter(np.full(len(z), q) + np.random.default_rng(1)
                       .normal(0, 0.035, len(z)), z, s=11, alpha=0.45,
                       color="#9A9A9A", linewidths=0,
                       label="nul par permutation" if q == 0 else None)
            ax.plot([q - 0.2, q + 0.2], [np.median(z)] * 2, color="#5A5A5A",
                    lw=1.2)
            r = D_[D_.voie == voie].iloc[0]
            ax.scatter([q], [r.spearman_differentiel], s=58, marker="D",
                       color=COUL.get(voie, "#888"), zorder=4,
                       edgecolor="white", linewidths=0.8)
        ax.axhline(0, color="#CCC", lw=0.6)
        ax.set_xticks(range(len(nuls)))
        ax.set_xticklabels([f"voie {v}" for v in sorted(nuls)], fontsize=7.5)
        ax.set_ylabel("Spearman, différentiel simulé\ncontre mesuré",
                      fontsize=7)
        ax.set_title("différentiel, placé dans son nul", fontsize=8,
                     loc="left")
        # mandataire neutre : un losange colore par voie ferait croire que
        # la legende designe une voie plutot que la valeur observee
        from matplotlib.lines import Line2D
        mand = Line2D([], [], marker="D", ls="none", color="#4A4A4A",
                      markersize=5.4, label="observé")
        h, l = ax.get_legend_handles_labels()
        ax.legend(h + [mand], l + ["observé"], frameon=False, fontsize=6.4,
                  loc="lower right")
        for ax in axes:
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            ax.tick_params(labelsize=7)
        fig.tight_layout()
        fig.savefig(f"{a.sortie}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"\n  {a.sortie}.png")
    print(f"  {a.sortie}_elements.tsv  {a.sortie}_ancres_absolues.tsv  "
          f"{a.sortie}_differentiel.tsv")


if __name__ == "__main__":
    main()
