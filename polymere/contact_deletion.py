#!/usr/bin/env python3
"""
contact_deletion.py - contact frequency between the promoter and each other
anchor of a locus, in the wild-type simulation (SAUVAGE) and in each
deletion condition, path A of simpol, one parameter set.

For each pair (promoter, anchor) and each condition, the contact frequency
is the share of conformations where the two beads are closer than
--seuil (simpol comparer.py, matrice_simulee), computed replicate by
replicate. The change is log2(deletion / SAUVAGE) of the replicate means,
with
  ci_bas, ci_haut   a 95 % interval by resampling the replicates of both
                    conditions (bootstrap)
  p_nul             the share of random splits of the SAUVAGE replicates
                    into two groups of the same sizes whose |log2 ratio|
                    reaches the observed one: the spread that seeds alone
                    produce, the null the simpol README asks for
A condition that removed no CTCF site has the same input as SAUVAGE and
is reported as null by construction (conditions.tsv of preparer_locus.py).

Outputs (in --dir)
  contacts_deletion.tsv   one row per (condition, anchor)
  contacts_deletion.png/.svg

Usage
  SIMPOL_LOCI=polymere/abat/loci.json python3 polymere/contact_deletion.py \
      --dir polymere/abat --locus abat
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "simpol"))
import config as cfg                                   # noqa: E402


def sous_echantillon(conf, maximum):
    """au plus `maximum` conformations, régulièrement espacées (comme simpol
    comparer.py, sans importer OpenMM)"""
    if maximum and len(conf) > maximum:
        return conf[np.linspace(0, len(conf) - 1, maximum).astype(int)]
    return conf


def contacts(conf, i, js, seuil):
    """fraction des conformations où la bille i est à moins de seuil de chaque bille j"""
    d = np.sqrt(((conf[:, [i], :] - conf[:, js, :]) ** 2).sum(2))
    return (d < seuil).mean(0)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", required=True)
    p.add_argument("--locus", required=True)
    p.add_argument("--jeu", type=int, default=None, help="défaut : jeu_retenu.txt")
    p.add_argument("--seuil", type=float, default=2.0, help="distance de contact (comparer.py)")
    p.add_argument("--max-conformations", type=int, default=3000, dest="max_conf",
                   help="par réplicat, sous-échantillonnage régulier")
    p.add_argument("--tirages", type=int, default=2000)
    p.add_argument("--graine", type=int, default=1)
    a = p.parse_args()
    d = Path(a.dir)
    jeu = a.jeu if a.jeu is not None else int((d / "jeu_retenu.txt").read_text().split()[0])
    c = cfg.LOCI[a.locus]
    idx = cfg.indices_ancres(a.locus)
    promo = c["ancres"][0][0]
    autres = [nm for nm in idx if nm != promo]
    js = [idx[nm] for nm in autres]
    retires = {}
    if (d / "conditions.tsv").exists():
        for ln in (d / "conditions.tsv").read_text().splitlines()[1:]:
            k, v = ln.split("\t")
            retires[k] = int(v)
    res = d / "resultats"
    par_cond = {}
    for cond in ["SAUVAGE"] + [k for k in retires if k != "SAUVAGE"]:
        fs = sorted(res.glob(f"{a.locus}_{cond}_A_jeu{jeu:03d}_rep*.npy"))
        if not fs:
            continue
        par_cond[cond] = np.array([contacts(sous_echantillon(np.load(f), a.max_conf),
                                            idx[promo], js, a.seuil) for f in fs])
        print(f"{cond:12s} {len(fs)} réplicats (jeu {jeu})")
    if "SAUVAGE" not in par_cond:
        sys.exit(f"aucun réplicat SAUVAGE du jeu {jeu} dans {res}")
    rng = np.random.default_rng(a.graine)
    W = par_cond["SAUVAGE"]
    eps = 1e-4
    rows = []
    for cond, X in par_cond.items():
        if cond == "SAUVAGE":
            continue
        obs = np.log2((X.mean(0) + eps) / (W.mean(0) + eps))
        boot = np.array([np.log2((X[rng.integers(0, len(X), len(X))].mean(0) + eps)
                                 / (W[rng.integers(0, len(W), len(W))].mean(0) + eps))
                         for _ in range(a.tirages)])
        lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
        # nul : deux groupes de réplicats SAUVAGE, tailles de la comparaison
        k = min(len(X), len(W) // 2)
        nul = []
        for _ in range(a.tirages):
            perm = rng.permutation(len(W))
            g1, g2 = W[perm[:k]], W[perm[k:2 * k]]
            nul.append(np.log2((g1.mean(0) + eps) / (g2.mean(0) + eps)))
        nul = np.abs(np.array(nul))
        pn = ((nul >= np.abs(obs)[None, :]).sum(0) + 1) / (len(nul) + 1)
        for j, nm in enumerate(autres):
            rows.append(dict(condition=cond, ancre=nm, sauvage=W.mean(0)[j], deletion=X.mean(0)[j],
                             log2_ratio=obs[j], ci_bas=lo[j], ci_haut=hi[j], p_nul=pn[j],
                             replicats_sauvage=len(W), replicats_deletion=len(X),
                             sites_ctcf_retires=retires.get(cond, "."),
                             nul_par_construction="oui" if retires.get(cond) == 0 else "non"))
    for cond, n in retires.items():
        if cond != "SAUVAGE" and n == 0:
            for nm in autres:
                rows.append(dict(condition=cond, ancre=nm, sauvage=W.mean(0)[autres.index(nm)],
                                 deletion=float("nan"), log2_ratio=0.0, ci_bas=float("nan"),
                                 ci_haut=float("nan"), p_nul=float("nan"),
                                 replicats_sauvage=len(W), replicats_deletion=0,
                                 sites_ctcf_retires=0, nul_par_construction="oui"))
    cols = list(rows[0]) if rows else []
    with open(d / "contacts_deletion.tsv", "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            fo.write("\t".join(f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k])
                               for k in cols) + "\n")
    for r in rows:
        print(f"  {r['condition']:12s} {promo} <-> {r['ancre']:16s} "
              f"{r['sauvage']:.3f} -> {r['deletion']:.3f}  log2 {r['log2_ratio']:+.2f} "
              f"[{r['ci_bas']:+.2f} ; {r['ci_haut']:+.2f}]  p_nul {r['p_nul']:.3g}"
              + ("  (nul par construction)" if r["nul_par_construction"] == "oui" else ""))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"svg.fonttype": "none", "font.size": 7})
    conds = [k for k in par_cond if k != "SAUVAGE"]
    fig, ax = plt.subplots(figsize=(1.2 + 0.9 * len(autres), 2.4))
    w = 0.8 / max(1, len(conds))
    pal = ["#2166AC", "#1B9E77", "#7570B3", "#E7298A"]
    coul = {"SANSELEMENT": "#C78A20"}
    for k, cond in enumerate(c_ for c_ in par_cond if c_ not in ("SAUVAGE", "SANSELEMENT")):
        coul[cond] = pal[k % len(pal)]
    for k, cond in enumerate(conds):
        rr = [r for r in rows if r["condition"] == cond]
        x = np.arange(len(rr)) + (k - (len(conds) - 1) / 2) * w
        y = [r["log2_ratio"] for r in rr]
        err = [[r["log2_ratio"] - r["ci_bas"] for r in rr], [r["ci_haut"] - r["log2_ratio"] for r in rr]]
        ax.bar(x, y, width=w * 0.9, color=coul.get(cond, "#888"), label=cond.lower())
        ax.errorbar(x, y, yerr=err, fmt="none", ecolor="black", lw=0.7, capsize=1.5)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(range(len(autres)))
    ax.set_xticklabels([nm.replace("_", "\n") for nm in autres], fontsize=6)
    ax.set_ylabel(f"contact avec {promo.replace('_', ' ')}\nlog2 délétion / sauvage")
    ax.set_title(f"{a.locus.upper()} - extrusion CTCF seule (voie A), jeu {jeu}", fontsize=7, loc="left")
    if conds:
        ax.legend(frameon=False, fontsize=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(d / "contacts_deletion.svg", bbox_inches="tight")
    fig.savefig(d / "contacts_deletion.png", dpi=300, bbox_inches="tight")
    print(f"-> {d / 'contacts_deletion.tsv'}, .png, .svg")


if __name__ == "__main__":
    main()
