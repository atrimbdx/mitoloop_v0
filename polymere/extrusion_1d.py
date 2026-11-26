#!/usr/bin/env python3
"""
extrusion_1d.py - does the 1D stage of simpol path A (loop extrusion stopped by
oriented CTCF sites) rank contacts as the 3D polymer does? Test on a locus
whose 3D trajectories are already computed (ABAT), before using the 1D stage
alone on every MitoLoop loop.

For every 3D trajectory of path A in <dir>/resultats (<locus>_<cond>_A_jeuNNN_repRR),
the 1D extrusion trajectory that drove it is regenerated exactly: same CTCF
bed (blake2b fingerprint checked against the .json of the run), same
parameter set, same seed (graine_stable(locus, condition, "A", jeu, replicat),
as simpol une_trajectoire.py), same number of steps. Nothing is fitted.

Two 1D scores per (condition, parameter set), replicates pooled:
  pont_<a>__<b>   direct bridging: fraction of instants where one cohesin has
                  one leg within --fenetre beads of anchor a and the other
                  within --fenetre beads of anchor b
  oe1d_<a>__<b>   observed/expected of a 1D contact proxy that also counts
                  indirect contacts: at each instant, genomic distance along
                  the shortest path of the chain plus the cohesin bridges
                  (s_eff), contact = s_eff^-alpha (alpha 1.5, ideal chain, set
                  beforehand); then each diagonal divided by its mean over the
                  analysed window, as evaluer_4c.py does for the 3D matrix
and r_oe_1d, the agreement of the 1D observed/expected profiles from each
anchor with the measured ones (same definition as r_oe in evaluer_4c.py).

The 3D values are read from <dir>/evaluation_4c.tsv (evaluer_4c.py). Across
the (condition, parameter set) groups, Spearman correlations 1D against 3D are
reported for r_oe and for each anchor pair. Rule fixed beforehand: the 1D
stage is taken as sufficient if rho >= 0.8 for r_oe and for the median anchor
pair, and if it orders SAUVAGE / UNIFORME against ALEATOIRE as 3D does.

Usage (any environment with numpy, pandas, scipy; OpenMM not needed):
  SIMPOL_LOCI=polymere/abat/loci.json python polymere/extrusion_1d.py --dir polymere/abat --locus abat
Output: <dir>/extrusion_1d.tsv (one row per group, 1D and 3D side by side),
<dir>/extrusion_1d_resume.tsv (correlations, verdict).
"""

import argparse
import csv
import hashlib
import itertools
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "simpol"))
sys.path.insert(0, str(HERE))
import config as cfg                                    # noqa: E402
from extrusion import Extrusion, barrieres_depuis_bed    # noqa: E402
from evaluer_4c import charger_matrice, oe               # noqa: E402

SEUIL_RHO = 0.8


def graine_stable(*parties):
    """copie de simpol une_trajectoire.py (qui importe OpenMM au chargement)"""
    cle = "|".join(str(p) for p in parties).encode()
    return int.from_bytes(hashlib.blake2b(cle, digest_size=4).digest(), "big") % (2 ** 31)


def empreinte(path):
    return hashlib.blake2b(Path(path).read_bytes(), digest_size=8).hexdigest()


def rangs(x):
    x = np.asarray(x, dtype=float)
    o = np.argsort(x, kind="mergesort")
    r = np.empty(len(x))
    r[o] = np.arange(len(x))
    for v in np.unique(x):                      # ex aequo : rang moyen
        m = x == v
        r[m] = r[m].mean()
    return r


def spearman(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return float("nan"), int(ok.sum())
    return float(np.corrcoef(rangs(x[ok]), rangs(y[ok]))[0, 1]), int(ok.sum())


def trajectoire_1d(d, locus, cond, jeu, rep, meta, blocs_defaut, geo):
    """positions des jambes (instants, facteurs, 2) de la trajectoire 1D qui a piloté le run 3D"""
    n, d0, f0, tb, chrom = geo
    par = meta.get("parametres") or {}
    jx = cfg.jeux("A")[jeu % len(cfg.jeux("A"))]
    for k in ("separation_kb", "processivite_kb", "force_barriere"):
        par.setdefault(k, jx[k])
    bed = d / (par.get("ctcf_fichier") or f"ctcf_{cond}.bed")
    if not bed.is_file():
        return None, f"{bed.name} absent"
    note = ""
    if par.get("ctcf_empreinte") and par["ctcf_empreinte"] != empreinte(bed):
        note = f"{bed.name} modifié depuis le run 3D"
    bp, bm = barrieres_depuis_bed(str(bed), chrom, d0, f0, tb, n, force_max=par["force_barriere"])
    e = Extrusion(n, bp, bm, separation_billes=par["separation_kb"] * 1000 // tb,
                  processivite_billes=par["processivite_kb"] * 1000 // tb,
                  graine=graine_stable(locus, cond, "A", jeu, rep))
    blocs = int(meta.get("blocs") or blocs_defaut)
    return e.trajectoire(blocs * 4, 4), note


def contacts_1d(traj, n, instants, alpha):
    """moyenne sur les instants de s_eff^-alpha (chemin le plus court : fibre + ponts de cohésine)"""
    M = np.zeros((n, n))
    sel = np.linspace(0, len(traj) - 1, min(instants, len(traj))).astype(int)
    i0 = np.arange(n - 1)
    for t in sel:
        g, dd = traj[t][:, 0], traj[t][:, 1]
        ok = dd > g + 1
        r = np.r_[i0, g[ok]]
        c = np.r_[i0 + 1, dd[ok]]
        G = coo_matrix((np.ones(len(r)), (r, c)), shape=(n, n)).tocsr()
        S = shortest_path(G, method="D", directed=False, unweighted=True)
        with np.errstate(divide="ignore"):
            P = np.where(S > 0, S ** (-alpha), 0.0)
        M += P
    return M / len(sel)


def pontage(traj, ia, ib, w):
    lo, hi = min(ia, ib), max(ia, ib)
    g, dd = traj[:, :, 0], traj[:, :, 1]
    pont = (np.abs(g - lo) <= w) & (np.abs(dd - hi) <= w)
    return float(pont.any(axis=1).mean())


def lire_3d(path):
    if not path.is_file():
        return {}
    with open(path) as fh:
        return {(r["condition"], r["jeu"]): r for r in csv.DictReader(fh, delimiter="\t")}


def fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def famille(cond):
    for f in ("ALEATOIRE", "UNIFORME", "SAUVAGE"):
        if cond.startswith(f):
            return f
    return "autre"


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", required=True)
    p.add_argument("--locus", required=True)
    p.add_argument("--instants", type=int, default=200, help="instants 1D par réplicat pour s_eff")
    p.add_argument("--alpha", type=float, default=1.5, help="exposant de la chaîne idéale (fixé a priori)")
    p.add_argument("--fenetre", type=int, default=1, help="tolérance (billes) du pontage direct")
    p.add_argument("--blocs", type=int, default=2000, help="blocs 3D si le .json du run manque")
    p.add_argument("--conditions", default="", help="liste séparée par des virgules (défaut : toutes)")
    a = p.parse_args()
    d = Path(a.dir)
    c = cfg.LOCI[a.locus]
    n, d0, f0 = cfg.geometrie(a.locus)
    tb = c["taille_bin"]
    geo = (n, d0, f0, tb, c["chrom"])
    lo, hi = int((c["debut"] - d0) // tb), int(np.ceil((c["fin"] - d0) / tb))
    idx = cfg.indices_ancres(a.locus)
    noms = list(idx)
    paires = list(itertools.combinations(noms, 2))
    Mes, vu = charger_matrice(d / "matrice.tsv", d0, f0, tb) if (d / "matrice.tsv").is_file() else (None, None)
    OEm = oe(Mes, vu, lo, hi) if Mes is not None else None

    def profils(O):
        out = {}
        for nm, i in idx.items():
            v = O[i, lo:hi].copy()
            v[np.abs(np.arange(lo, hi) - i) <= 2] = np.nan
            out[nm] = v
        return out
    pm = profils(OEm) if OEm is not None else None

    voulu = set(a.conditions.split(",")) if a.conditions else None
    groupes = defaultdict(list)
    for f in sorted((d / "resultats").glob(f"{a.locus}_*_A_jeu*_rep*.npy")):
        m = re.match(rf"{re.escape(a.locus)}_(.+)_A_jeu(\d+)_rep(\d+)\.npy", f.name)
        if m and (voulu is None or m.group(1) in voulu):
            groupes[(m.group(1), int(m.group(2)))].append(int(m.group(3)))
    if not groupes:
        sys.exit(f"aucune trajectoire de la voie A dans {d / 'resultats'}")
    trois_d = lire_3d(d / "evaluation_4c.tsv")
    if not trois_d:
        print(f"attention : {d / 'evaluation_4c.tsv'} absent, lancer evaluer_4c.py d'abord "
              "(sans lui, pas de comparaison 1D/3D)", file=sys.stderr)

    rows, notes = [], []
    for k, ((cond, jeu), reps) in enumerate(sorted(groupes.items()), 1):
        M, trajs = np.zeros((n, n)), []
        for rep in sorted(reps):
            jf = d / "resultats" / f"{a.locus}_{cond}_A_jeu{jeu:03d}_rep{rep:02d}.json"
            meta = json.loads(jf.read_text()) if jf.is_file() else {}
            if meta.get("n_billes") and int(meta["n_billes"]) != n:
                notes.append(f"{cond} jeu {jeu} rep {rep} : {meta['n_billes']} billes au run, {n} ici (--bin/--marge ?)")
                continue
            traj, note = trajectoire_1d(d, a.locus, cond, jeu, rep, meta, a.blocs, geo)
            if traj is None:
                notes.append(f"{cond} jeu {jeu} rep {rep} : {note}")
                continue
            if note:
                notes.append(f"{cond} jeu {jeu} rep {rep} : {note}")
            trajs.append(traj)
            M += contacts_1d(traj, n, a.instants, a.alpha)
        if not trajs:
            continue
        M /= len(trajs)
        T = np.concatenate(trajs)
        O1 = oe(M, np.ones((n, n), dtype=bool), lo, hi)
        r1 = {}
        if pm is not None:
            ps = profils(O1)
            for nm in noms:
                x, y = pm[nm], ps[nm]
                ok = np.isfinite(x) & np.isfinite(y)
                r1[nm] = float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() >= 8 and x[ok].std() > 0 \
                    and y[ok].std() > 0 else float("nan")
        t3 = trois_d.get((cond, str(jeu)), {})
        row = dict(condition=cond, jeu=jeu, replicats=len(trajs),
                   r_oe_1d=float(np.nanmean(list(r1.values()))) if r1 else float("nan"),
                   r_oe_3d=fl(t3.get("r_oe_moyen")))
        for x, y in paires:
            row[f"oe1d_{x}__{y}"] = O1[idx[x], idx[y]]
            row[f"oe3d_{x}__{y}"] = fl(t3.get(f"oe_{x}__{y}"))
            row[f"pont_{x}__{y}"] = pontage(T, idx[x], idx[y], a.fenetre)
        rows.append(row)
        print(f"  {k}/{len(groupes)} {cond} jeu {jeu} : {len(trajs)} réplicats", file=sys.stderr, flush=True)
    if not rows:
        sys.exit("aucun groupe relu :\n  " + "\n  ".join(notes))

    cols = list(rows[0])
    with open(d / "extrusion_1d.tsv", "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            fo.write("\t".join(f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k]) for k in cols) + "\n")

    # comparaison 1D / 3D sur les groupes
    res = []
    rho, m = spearman([r["r_oe_1d"] for r in rows], [r["r_oe_3d"] for r in rows])
    res.append(("r_oe", "oe1d", rho, m))
    rp = []
    for x, y in paires:
        for kind in ("oe1d", "pont"):
            rr, m = spearman([r[f"{kind}_{x}__{y}"] for r in rows], [r[f"oe3d_{x}__{y}"] for r in rows])
            res.append((f"{x}__{y}", kind, rr, m))
            if kind == "oe1d":
                rp.append(rr)
    med = float(np.nanmedian(rp)) if np.isfinite(rp).any() else float("nan")

    def ordre(cle):
        """famille -> meilleure valeur ; True si SAUVAGE ou UNIFORME dépasse le meilleur ALEATOIRE"""
        best = defaultdict(lambda: -np.inf)
        for r in rows:
            v = r[cle]
            if np.isfinite(v):
                best[famille(r["condition"])] = max(best[famille(r["condition"])], v)
        if "ALEATOIRE" not in best:
            return None, dict(best)
        return max(best.get("SAUVAGE", -np.inf), best.get("UNIFORME", -np.inf)) > best["ALEATOIRE"], dict(best)
    o1, b1 = ordre("r_oe_1d")
    o3, b3 = ordre("r_oe_3d")
    meme_ordre = None if o1 is None or o3 is None else (o1 == o3)
    suffit = (np.isfinite(rho) and rho >= SEUIL_RHO and np.isfinite(med) and med >= SEUIL_RHO
              and meme_ordre is not False)
    with open(d / "extrusion_1d_resume.tsv", "w") as fo:
        fo.write("mesure\tscore_1d\tspearman_vs_3d\tgroupes\n")
        for nom, kind, rr, m in res:
            fo.write(f"{nom}\t{kind}\t{rr:.3f}\t{m}\n")
        fo.write(f"paires_mediane\toe1d\t{med:.3f}\t.\n")
        fo.write(f"ordre_temoins\tr_oe\t{meme_ordre}\t.\n")
        fo.write(f"verdict\t.\t{'1D suffisant' if suffit else '1D insuffisant ou non testable'}\t.\n")

    print(f"{'mesure':42s} {'score 1D':>8s} {'rho vs 3D':>9s} {'groupes':>7s}")
    for nom, kind, rr, m in res:
        print(f"{nom:42s} {kind:>8s} {rr:9.3f} {m:7d}")
    print(f"\nmédiane des paires d'ancres (oe1d) : rho = {med:.3f}")
    if o1 is not None:
        f2 = lambda b: ", ".join(f"{k} {v:.3f}" for k, v in sorted(b.items()))  # noqa: E731
        print(f"meilleur r_oe par famille, 1D : {f2(b1)}")
        print(f"meilleur r_oe par famille, 3D : {f2(b3)}")
        print(f"SAUVAGE/UNIFORME au-dessus d'ALEATOIRE : 1D {o1}, 3D {o3}")
    print(f"verdict (rho >= {SEUIL_RHO} pour r_oe et la médiane des paires, même ordre face aux témoins) : "
          + ("l'étape 1D suffit" if suffit else "l'étape 1D ne suffit pas, ou comparaison incomplète"))
    for nt in notes[:20]:
        print(f"note : {nt}", file=sys.stderr)
    print(f"→ {d / 'extrusion_1d.tsv'}, {d / 'extrusion_1d_resume.tsv'}")


if __name__ == "__main__":
    main()
