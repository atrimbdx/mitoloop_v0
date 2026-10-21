#!/usr/bin/env python3
"""
validation_eqtl.py — Curta. Les arguments de ajout_mitoloop.py désignent-ils
les variants qui modifient réellement l'expression des gènes mitochondriaux ?

Données extérieures : les eQTL affinés de GTEx v8 (DAP-G, Wen et al. 2017 ;
GTEx_v8_finemapping_DAPG.txt.gz). Pour un même signal d'eQTL (un gène, un
tissu, un groupe de variants en déséquilibre de liaison), DAP-G donne à chaque
variant sa probabilité d'être causal (PIP). On compare, à l'intérieur de
chaque signal, les variants causaux probables (PIP >= --causal) aux variants
du même signal qui ne font que les accompagner (PIP <= --background). La
comparaison à l'intérieur d'un signal neutralise ce qui est commun à ses
variants : gène, tissu, région, fréquence, puissance de détection.

Seuls les signaux des gènes du paquet MitoLoop, et leurs variants situés dans
les régions régulatrices du même gène, sont gardés. L1 est une propriété du
gène, identique dans un signal : il n'est pas testé ici (il l'est sur les
fibroblastes de patients, dans results/l1_raw/validation).

Statistique, par argument : la différence moyenne, sur les signaux, entre la
part des variants causaux et la part des variants d'accompagnement qui le
portent ; p du test stratifié sur les signaux (approximation normale de la
permutation des étiquettes à l'intérieur de chaque signal, comme le test de
Cochran-Mantel-Haenszel). Même chose pour le nombre d'arguments.

Usage
    python3 reannot/validation_eqtl.py --data transfert_mitoloop \\
        --dapg data/GTEx_v8_finemapping_DAPG.txt.gz --fasta data/GRCh38.fa \\
        --out results/reannot/validation_eqtl.tsv
"""

from __future__ import annotations

import argparse
import gzip
import math
import sys
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ajout_mitoloop import (W_CPG, Intervals, annotate, fetch_phylop,  # noqa: E402
                            fetch_seqs, load_genes, read_motifs)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--data", required=True, type=Path)
    p.add_argument("--dapg", required=True, type=Path)
    p.add_argument("--fasta", required=True)
    p.add_argument("--samtools", default="samtools")
    p.add_argument("--causal", type=float, default=0.5)
    p.add_argument("--background", type=float, default=0.05)
    p.add_argument("--out", required=True, type=Path)
    a = p.parse_args()
    opt = SimpleNamespace(min_biosamples=5, min_rel=0.8, min_ic=1.0, min_delta=1.0,
                          min_phylop=2.27, site_fort_ic=1.5, site_fort_delta=2.0)
    D = a.data
    genes = load_genes(D / "genes.tsv")
    by_id = {r["gene_id"].split(".")[0]: g for g, r in genes.items() if r.get("gene_id")}
    regions = Intervals(D / "regions.bed.gz")
    ann = {f.name[:-7]: Intervals(f) for f in sorted((D / "annotations").glob("*.bed.gz"))}
    chrom_ann = [k for k in ann if k.startswith("chromhmm_")]
    mats = read_motifs(D / "motifs.jaspar")
    tfs = sorted(mats)

    # 1. signaux DAP-G des gènes du paquet, variants dans les régions du gène
    clusters = defaultdict(list)          # (tissu, gène, signal) -> [(variant, pip)]
    n_rows = 0
    with gzip.open(a.dapg, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        c = {h: i for i, h in enumerate(head)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            g = by_id.get(f[c["gene_id"]].split(".")[0])
            if g is None:
                continue
            n_rows += 1
            v = f[c["variant_id"]].split("_")
            if len(v) < 4:
                continue
            chrom, pos, ref, alt = v[0], int(v[1]), v[2], v[3]
            hit = regions.query(chrom, pos - 1, pos - 1 + len(ref))
            if not any(lab[0] == g for _s, _e, lab in hit):
                continue
            pip = float(f[c["variant_pip"]])
            clusters[(f[c["tissue_id"]], g, f[c["cluster_id"]])].append(
                ((chrom, pos, ref, alt), pip))
    keep = {}
    for k, vs in clusters.items():
        ca = [v for v, p_ in vs if p_ >= a.causal]
        bg = [v for v, p_ in vs if p_ <= a.background]
        if ca and bg:
            keep[k] = (ca, bg)
    variants = {v for ca, bg in keep.values() for v in ca + bg}
    print(f"DAP-G : {n_rows} lignes des gènes du paquet, {len(clusters)} signaux "
          f"avec des variants dans les régions, {len(keep)} avec au moins un variant "
          f"causal (PIP >= {a.causal}) et un d'accompagnement (PIP <= {a.background}), "
          f"{len(variants)} variants")
    if not keep:
        a.out.write_text("argument\tsignals\tcausal\tbackground\tdifference\tp\n")
        return

    # 2. arguments de chaque variant (gène du signal)
    need = {f"{v[0]}:{max(1, v[1] - W_CPG)}-{v[1] + len(v[2]) + W_CPG}" for v in variants}
    seqs = fetch_seqs(a.fasta, a.samtools, sorted(need))
    phy = fetch_phylop(D / "annotations" / "phylop.bedgraph.gz", "tabix",
                       {(v[0], v[1]) for v in variants})
    args_of = {}
    for (t, g, s), (ca, bg) in keep.items():
        for v in ca + bg:
            if (v, g) in args_of:
                continue
            row = [v[0], str(v[1]), v[2], v[3], g]
            hit = regions.query(v[0], v[1] - 1, v[1] - 1 + len(v[2]))
            hit = [h for h in hit if h[2][0] == g]
            r = annotate(row, 0, 1, 2, 3, 4, hit, genes, ann, chrom_ann,
                         mats, tfs, seqs, opt, phy)
            al = [x for x in r["MITO_arguments"].split(",") if x not in (".", "")]
            # L1 est le même pour tout le signal
            args_of[(v, g)] = {x for x in al if x != "L1_transmis"}
    names = sorted({x for s_ in args_of.values() for x in s_})
    # site_<TF> et fixation_<TF> regroupés aussi par famille
    fam = {"site_TF": lambda s_: any(x.startswith("site_") for x in s_),
           "fixation_TF": lambda s_: any(x.startswith("fixation_") for x in s_),
           "chromatine_active": lambda s_: any(x.startswith("chromatine_active") for x in s_)}

    # stratified on the signals: under the null (labels exchangeable within a
    # signal), the causal sum of a signal has mean nc * mean(y) and variance
    # nc * (n - nc) / (n - 1) * var(y); Z over the signals (the permutation
    # test's normal approximation, as in Cochran-Mantel-Haenszel)
    tests = [(n, (lambda n: lambda s_: float(n in s_))(n)) for n in names] + \
        [(k, (lambda f_: lambda s_: float(f_(s_)))(f_)) for k, f_ in fam.items()] + \
        [("n_arguments", lambda s_: float(len(s_)))]
    rows = []
    for name, fn in tests:
        T = E = V = 0.0
        diffs, ca_all, bg_all = [], [], []
        for k, (ca, bg) in keep.items():
            yc = [fn(args_of[(v, k[1])]) for v in ca]
            yb = [fn(args_of[(v, k[1])]) for v in bg]
            y = yc + yb
            n, nc = len(y), len(yc)
            m = sum(y) / n
            var = sum((x - m) ** 2 for x in y) / n
            T += sum(yc)
            E += nc * m
            V += nc * (n - nc) / (n - 1) * var if n > 1 else 0
            diffs.append(sum(yc) / nc - sum(yb) / len(yb))
            ca_all += yc
            bg_all += yb
        z = (T - E) / math.sqrt(V) if V > 0 else 0.0
        pval = math.erfc(abs(z) / math.sqrt(2))
        rows.append((name, len(keep), sum(ca_all) / len(ca_all),
                     sum(bg_all) / len(bg_all), sum(diffs) / len(diffs), z, pval))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("argument\tsignals\tcausal\tbackground\tdifference\tz\tp\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]:.4g}\t{r[3]:.4g}\t{r[4]:.4g}\t{r[5]:.3f}\t{r[6]:.3g}\n")
    print(f"\n{'argument':28s}{'causal':>9s}{'accomp.':>9s}{'diff':>9s}{'z':>8s}{'p':>9s}")
    for r in rows:
        print(f"{r[0]:28s}{r[2]:9.3f}{r[3]:9.3f}{r[4]:9.3f}{r[5]:8.2f}{r[6]:9.3g}")
    print(f"Résultat : {a.out}")


if __name__ == "__main__":
    main()
