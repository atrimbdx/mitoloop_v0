#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rendu 3D des conformations simulees.

CE QU'UNE IMAGE 3D MONTRE, ET CE QU'ELLE NE MONTRE PAS
  Une conformation est un tirage parmi des millions. Elle ne represente
  aucune cellule en particulier et n'a pas de valeur de preuve : seules les
  moyennes d'ensemble en ont. Une image sert a montrer a quoi ressemble une
  conformation ou le hub est forme, pas a affirmer que le locus a cette
  forme.

  Pour que la figure reste honnete, trois conformations sont rendues et non
  une : la plus proche du comportement moyen de l'ensemble, une ou le hub
  est forme, une ou il ne l'est pas. La fraction de l'ensemble que chacune
  represente est inscrite sous l'image.

SORTIES
  <prefixe>_<condition>_<voie>.png        planche des trois conformations
  <prefixe>_<condition>_<voie>_<type>.xyz coordonnees en unites de bille
  <prefixe>_<condition>_<voie>_<type>.pdb pour PyMOL ou ChimeraX
  <prefixe>_ensemble.tsv                  ce que chaque image represente

PYMOL
  load conformation.pdb
  hide everything ; show cartoon ; cartoon tube ; set cartoon_tube_radius, 1.2
  spectrum b, rainbow            # colore selon la position genomique
  show spheres, name ZN ; color red, name ZN    # les elements annotes

CHIMERAX
  open conformation.pdb
  cartoon ; cartoon style modeh tube rad 1.2 ; color bfactor palette rainbow

USAGE
  python3 rendu3d.py --locus tyrp1 --resultats resultats/ \\
      --conditions MEL FIBRO --sortie figures/conf
"""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection

import config as cfg
from comparer import charger_resultats

# Echelle d'affichage du PDB. Les coordonnees simulees sont en diametres de
# bille ; on les met a l'echelle pour que la distance entre billes voisines
# vaille 3,8 A, seule valeur pour laquelle PyMOL et ChimeraX acceptent de
# tracer un ruban continu sur une chaine de carbones alpha. L'echelle n'a
# donc aucun sens physique : le facteur vers les nanometres est inscrit dans
# les REMARK du fichier.
ECHELLE_PDB = 3.8
NM_REF, PB_REF = 30.0, 5_000     # ordre de grandeur : 5 kb tiennent dans ~30 nm


def sans_accent(s):
    """Le format PDB est ASCII : un REMARK accentue casse les analyseurs
    stricts, et un nom de fichier accentue voyage mal."""
    import unicodedata
    return (unicodedata.normalize("NFKD", str(s))
            .encode("ascii", "ignore").decode())


def nm_par_bille(taille_bin):
    """Diametre d'une bille, a volume de chromatine proportionnel."""
    return NM_REF * (taille_bin / PB_REF) ** (1 / 3)


# ======================================================================
def aligner(p):
    """Repere propre de la conformation : plus grande extension en x."""
    p = p - p.mean(0)
    _, _, V = np.linalg.svd(p, full_matrices=False)
    q = p @ V.T
    # orientation stable d'une image a l'autre
    for k in range(3):
        if q[:, k].max() + q[:, k].min() < 0:
            q[:, k] *= -1
    return q


def lisser(p, facteur=6):
    """Spline cubique passant par les billes, pour un tube continu."""
    from scipy.interpolate import make_interp_spline
    n = len(p)
    t = np.arange(n)
    tt = np.linspace(0, n - 1, (n - 1) * facteur + 1)
    k = 3 if n > 3 else 1
    return np.column_stack([make_interp_spline(t, p[:, d], k=k)(tt)
                            for d in range(3)]), tt


def distances_ancres(conf, j):
    """Distances entre toutes les paires d'ancres, par conformation."""
    P = conf[:, j, :]
    return np.sqrt(((P[:, :, None, :] - P[:, None, :, :]) ** 2).sum(3))


def choisir(conf, j, seuil):
    """Trois conformations : mediane de l'ensemble, hub forme, hub ouvert.

    La mediane est celle dont le vecteur de distances entre ancres est le
    plus proche de la mediane de l'ensemble, et non une moyenne de
    coordonnees : moyenner des coordonnees produirait une forme qui
    n'existe dans aucune conformation.
    """
    if len(j) < 2:
        return {"médiane": 0}, {"médiane": 1.0}
    d = distances_ancres(conf, j)
    m = ~np.eye(len(j), dtype=bool)
    v = d[:, m]
    ferme = v.max(1) < seuil
    med = np.median(v, axis=0)
    ecart = np.abs(v - med).mean(1)

    out, frac = {}, {}
    out["médiane"] = int(np.argmin(ecart))
    frac["médiane"] = 1.0
    for nom, masque in (("hub formé", ferme), ("hub ouvert", ~ferme)):
        if not masque.any():
            continue
        k = np.where(masque)[0]
        o = k[np.argsort(ecart[k])]
        # La conformation mediane appartient forcement a l'un des deux
        # groupes : sans ce saut, elle etait affichee deux fois, une fois
        # comme mediane et une fois comme representante de son groupe.
        choix = next((int(i) for i in o if i not in out.values()), int(o[0]))
        out[nom] = choix
        frac[nom] = float(masque.mean())
    return out, frac


# ======================================================================
def ecrire_xyz(chemin, p, idx, meta):
    noms = {i: n for n, i in idx.items()}
    with open(chemin, "w") as f:
        f.write(f"{len(p)}\n")
        f.write(f"# {sans_accent(meta)}  unite = 1 diametre de bille\n")
        for i, (x, y, z) in enumerate(p):
            f.write(f"{'Fe' if i in noms else 'C':2s} {x:9.3f} {y:9.3f} "
                    f"{z:9.3f}\n")


def ecrire_pdb(chemin, p, idx, c, d0, titre=""):
    """Une bille par carbone alpha, pour que le ruban se trace. Le facteur
    B porte la position genomique : `spectrum b` colore alors selon la
    coordonnee le long du chromosome."""
    q = aligner(p)
    d = np.linalg.norm(np.diff(q, axis=0), axis=1).mean() or 1.0
    q = q * (ECHELLE_PDB / d)
    noms = {i: n for n, i in idx.items()}
    pos_mb = (d0 + np.arange(len(q)) * c["taille_bin"]) / 1e6
    b = 100 * (pos_mb - pos_mb.min()) / max(np.ptp(pos_mb), 1e-9)

    nm = nm_par_bille(c["taille_bin"])
    with open(chemin, "w") as f:
        f.write(f"REMARK   1 {sans_accent(titre)}\n")
        f.write(f"REMARK   2 {c['chrom']}:{d0:,}-{d0 + len(q)*c['taille_bin']:,}"
                f"  {c['taille_bin']:,} pb par bille\n")
        f.write("REMARK   3 ECHELLE ARBITRAIRE : la distance entre billes "
                "voisines vaut 3,8 A pour\n")
        f.write("REMARK   3 que le ruban se trace. 1 bille simulee vaut de "
                f"l'ordre de {nm:.0f} nm,\n")
        f.write(f"REMARK   3 soit un facteur {nm*10/ECHELLE_PDB:.0f} "
                "vers les angstroms reels.\n")
        f.write("REMARK   4 FACTEUR B = position genomique, de 0 a 100. "
                "PyMOL : spectrum b\n")
        f.write("REMARK   5 Les elements annotes sont des atomes ZN, "
                "residus nommes ELE.\n")
        for i, (x, y, z) in enumerate(q, start=1):
            k = i - 1
            anc = k in noms
            at = "ZN" if anc else "CA"
            res = "ELE" if anc else "ALA"
            f.write(f"ATOM  {i:5d}  {at:<3s} {res} A{i:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{b[k]:6.2f}"
                    f"          {'ZN' if anc else 'C':>2s}\n")
        for i in range(1, len(q)):
            f.write(f"CONECT{i:5d}{i+1:5d}\n")
        f.write("END\n")


# ======================================================================
def tracer(ax, p, idx, c, d0, titre, sous_titre):
    q = aligner(p)
    s, tt = lisser(q)
    pos = (d0 + tt * c["taille_bin"]) / 1e6
    seg = np.stack([s[:-1], s[1:]], axis=1)
    lc = Line3DCollection(seg, cmap="viridis", linewidths=2.2)
    lc.set_array((pos[:-1] + pos[1:]) / 2)
    ax.add_collection3d(lc)

    for nom, i in idx.items():
        ax.scatter(*q[i], s=64, color="#B2182B", edgecolor="white",
                   linewidths=0.7, depthshade=False, zorder=5)

    # Etiquettes ecartees le long de l'axe le plus court : posees sur la
    # bille, elles se recouvraient des que deux ancres etaient proches.
    ordre = sorted(idx.items(), key=lambda kv: q[kv[1], 1])
    etendue = np.ptp(q, axis=0)
    pas = etendue[1] * 0.075
    y0 = q[[i for _, i in ordre], 1].mean() - pas * (len(ordre) - 1) / 2
    for r_, (nom, i) in enumerate(ordre):
        cible = np.array([q[i, 0] + etendue[0] * 0.10, y0 + r_ * pas, q[i, 2]])
        ax.plot(*np.stack([q[i], cible], axis=1), color="#9A9A9A", lw=0.5)
        ax.text(*cible, " " + nom.replace("_", " "), fontsize=5.6,
                color="#1A1A1A", va="center")

    # Cadre ajuste a la conformation, meme echelle sur les trois axes :
    # un cube centre sur zero laissait la fibre occuper un dixieme du cadre.
    c0 = (q.max(0) + q.min(0)) / 2
    demi = etendue / 2 * 1.18
    demi = np.maximum(demi, etendue.max() * 0.06)
    for k, f in enumerate((ax.set_xlim, ax.set_ylim, ax.set_zlim)):
        f(c0[k] - demi[k], c0[k] + demi[k])
    # zoom : sans lui, les axes 3D reservent une region cubique et la
    # conformation n'occupait qu'une fraction du cadre
    ax.set_box_aspect(tuple(demi / demi.max()), zoom=1.45)
    ax.set_axis_off()
    ax.set_title(titre, fontsize=8, loc="left", pad=0)
    ax.text2D(0.0, -0.04, sous_titre, transform=ax.transAxes, fontsize=6.2,
              color="#4A4A4A", va="top")
    return lc


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--locus", default="tyrp1", choices=sorted(cfg.LOCI))
    p.add_argument("--bin", type=int, dest="taille_bin_cli")
    p.add_argument("--marge", type=int)
    p.add_argument("--resultats", default="resultats")
    p.add_argument("--conditions", nargs="+", default=["MEL", "FIBRO"])
    p.add_argument("--voies", nargs="+", default=None)
    p.add_argument("--bilan", default="bilan.tsv")
    p.add_argument("--seuil-hub", type=float, default=2.5, dest="seuil_hub")
    p.add_argument("--vue", type=float, nargs=2, default=(18.0, -62.0),
                   help="élévation et azimut")
    p.add_argument("--sortie", default="conformation")
    a = p.parse_args()
    cfg.surcharger(a.locus, a.taille_bin_cli, a.marge)

    c = cfg.LOCI[a.locus]
    n, d0, _ = cfg.geometrie(a.locus)
    idx = cfg.indices_ancres(a.locus)

    choix = {}
    if Path(a.bilan).exists():
        import pandas as pd
        for _, r in pd.read_csv(a.bilan, sep="\t").iterrows():
            choix[(r.condition, r.voie)] = int(r.jeu)

    Path(a.sortie).parent.mkdir(parents=True, exist_ok=True)
    lignes = []
    for cond in a.conditions:
        grp = charger_resultats(a.resultats, a.locus, cond)
        if not grp:
            print(f"  aucun résultat pour {cond}")
            continue
        voies = a.voies or sorted({v for v, _ in grp})
        for voie in voies:
            jeux = [j for (v, j) in grp if v == voie]
            if not jeux:
                continue
            jeu = choix.get((cond, voie), jeux[0])
            if jeu not in jeux:
                jeu = jeux[0]
            conf = np.concatenate([np.load(f) for _, f in grp[(voie, jeu)]])

            noms_hub = c["hub"] if cond == "MEL" else (c["hub_fibro"]
                                                       or c["hub"])
            j = [idx[x] for x in noms_hub if x in idx]
            sel, frac = choisir(conf, j, a.seuil_hub)

            fig = plt.figure(figsize=(3.2 * len(sel), 2.9))
            lc = None
            for k, (nom, i) in enumerate(sel.items(), start=1):
                ax = fig.add_subplot(1, len(sel), k, projection="3d")
                ax.view_init(*a.vue)
                st = ("conformation la plus proche du comportement moyen"
                      if nom == "médiane"
                      else f"{frac[nom]*100:.0f} % de l'ensemble")
                lc = tracer(ax, conf[i], idx, c, d0, nom, st)
                base = (f"{a.sortie}_{cond}_{voie}_"
                        f"{sans_accent(nom).replace(' ', '-')}")
                ecrire_xyz(f"{base}.xyz", conf[i], idx,
                           f"{a.locus} {cond} voie {voie} jeu {jeu} {nom}")
                ecrire_pdb(f"{base}.pdb", conf[i], idx, c, d0,
                           f"{a.locus.upper()} {cond} voie {voie} "
                           f"jeu {jeu} - {nom}")
                lignes.append(dict(locus=a.locus, condition=cond, voie=voie,
                                   jeu=jeu, type=nom,
                                   fraction_ensemble=round(frac[nom], 4),
                                   conformations=len(conf),
                                   indice=int(i), fichier=f"{base}.pdb"))
            if lc is not None:
                cb = fig.colorbar(lc, ax=fig.axes, fraction=0.016, pad=0.02)
                cb.set_label(f"{c['chrom']}  (Mb)", fontsize=6.5)
                cb.ax.tick_params(labelsize=6)
            fig.suptitle(f"{a.locus.upper()} {cond} — voie {voie}, jeu {jeu}, "
                         f"{len(conf):,} conformations.  Un tirage ne vaut pas "
                         "une mesure : seules les moyennes d'ensemble le font.",
                         fontsize=7.6, x=0.012, ha="left")
            f = f"{a.sortie}_{cond}_{voie}.png"
            fig.savefig(f, dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"  {f}")

    if lignes:
        import pandas as pd
        t = pd.DataFrame(lignes)
        t.to_csv(f"{a.sortie}_ensemble.tsv", sep="\t", index=False)
        print(f"  {a.sortie}_ensemble.tsv")
        for _, r in t[t.type == "hub formé"].iterrows():
            print(f"    {r.condition:6s} voie {r.voie} : hub complet dans "
                  f"{r.fraction_ensemble*100:.1f} % des conformations")


if __name__ == "__main__":
    main()
