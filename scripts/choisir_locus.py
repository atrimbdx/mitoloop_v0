#!/usr/bin/env python3
"""
choisir_locus.py - ranks the promoters of one tissue by how well they
illustrate, on one locus figure, a Hi-C map, the virtual 4C loop call and the
ENCODE-rE2G elements of the same tissue.

A good locus to show has one clear loop, called in the tissue and found again
in each of its libraries, that ends on an element ENCODE-rE2G links to the
same gene in biosamples of the same tissue, far enough from the promoter to
stand apart on the map, near enough to stay inside a readable window, and few
other loops around to blur it. The script measures these points for every
loop of the tissue and ranks the genes by their best loop. It reads only the
outputs already written; it draws nothing.

Inputs
  --dir        results/loops_tissus: loops_<tissu>.tsv (loops_par_tissu.py
               separer) and <tissu>/peaks_mito.tsv (v4c_loops.R call +
               reciprocal, for summit_ratio)
  --tissu      muscle, cardiomyocyte, neurone or rein
  --elements   re2g_elements.tsv (re2g_panel.py), with its biosamples column
               (file=CellType;...); elements linked to the gene of the loop
  --interet    optional list of genes to keep (MITO list, L1-sensitive genes)

Output (--out, default <dir>/locus_<tissu>.tsv), one row per gene, best loop
first, columns:
  gene chrom vp_pos summit distance_kb n_fichiers n_presents summit_ratio
  elements_tissu   biosamples of the tissue linking the element on the anchor
  elements_tous    biosamples of any tissue linking it
  boucles_gene     loops of this gene in the tissue
  boucles_fenetre  loops of any gene with both ends in the window
  tissus_avec      tissues where this gene has a loop on a same-tissue element
  score, commande  (the locus_figure.py call that draws it)

Ranking (by gene, its best loop; stated, not tuned on anything)
  required   the loop anchor (summit +/- --tol) holds an element linked to the
             gene in >= 1 biosample of the tissue
  score      presence (n_presents / n_fichiers, 1 when not measured)
             + min(biosamples of the tissue, 3) / 3
             + min(log2(summit_ratio), 2) / 2
             + 1 if --dmin <= distance <= --dmax
             - 0.25 per extra loop of the gene beyond 2
             - 0.1 per loop of other genes in the window beyond 3

Usage
  python3 scripts/choisir_locus.py --dir results/loops_tissus --tissu muscle \\
      --elements ~/reannot_git/results/paquet/re2g_tous/re2g_elements.tsv \\
      --interet data/mito_genes.txt --top 15
"""

import argparse
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

TISSUS = ("muscle", "cardiomyocyte", "neurone", "rein")
# CellType of the rE2G biosamples of each tissue (as ag2_balayage.py CELLTYPES)
CELLTYPES = {
    "muscle": (r"muscle|gastrocnemius|psoas|myotube|myoblast|myocyte", r"cardiac|heart|smooth"),
    "cardiomyocyte": (r"cardiac|heart|ventric|atri(um|al)", r"smooth"),
    "neurone": (r"brain|cortex|gyrus|caudate|putamen|cerebell|hippocamp|neuron|astrocyte|"
                r"nucleus accumbens|substantia nigra|amygdala",
                r"adrenal|kidney|smooth|vascul|endotheli|pericyte"),
    "rein": (r"kidney|\brenal|nephr|tubul|HK-2|RPTEC", r"adrenal"),
}


def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").lstrip("#").split("\t")
        for line in fh:
            if line.strip():
                f = line.rstrip("\n").split("\t")
                yield dict(zip(head, f + [""] * (len(head) - len(f))))


def num(x, default=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def de_tissu(cts, tissu):
    inc, exc = CELLTYPES[tissu]
    ri, rx = re.compile(inc, re.I), re.compile(exc, re.I)
    return sum(1 for c in cts if ri.search(c) and not rx.search(c))


def lire_elements(path):
    """gene -> [(chrom, start, end, [CellType per biosample])]"""
    el = defaultdict(list)
    for r in read_tsv(path):
        if r.get("biosamples"):
            cts = [x.split("=", 1)[1] for x in r["biosamples"].split(";") if "=" in x]
        else:
            cts = [x for x in r.get("celltypes", "").split(";") if x]
        el[r["gene"]].append((r["chrom"], int(num(r["start"], 0)), int(num(r["end"], 0)), cts))
    return el


def sur_ancre(elements, chrom, summit, tol):
    """elements of the gene within tol of the summit"""
    return [e for e in elements if e[0] == chrom and e[1] - tol <= summit <= e[2] + tol]


def boucles(d, tissu, interet):
    lp = d / f"loops_{tissu}.tsv"
    if not lp.exists():
        return []
    return [r for r in read_tsv(lp) if not interet or r["gene"] in interet]


def ratios(d, tissu):
    """(gene, summit) -> summit_ratio from the tissue's peaks_mito.tsv"""
    pm, out = d / tissu / "peaks_mito.tsv", {}
    if pm.exists():
        for r in read_tsv(pm):
            out[(r["gene"], int(num(r.get("summit"), 0)))] = num(r.get("summit_ratio"))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", required=True)
    p.add_argument("--tissu", required=True, choices=TISSUS)
    p.add_argument("--elements", required=True)
    p.add_argument("--interet", help="genes to keep (whitespace-separated)")
    p.add_argument("--tol", type=int, default=10_000, help="bp around the summit (one bin)")
    p.add_argument("--dmin", type=int, default=50_000)
    p.add_argument("--dmax", type=int, default=500_000)
    p.add_argument("--fenetre", type=int, default=0,
                   help="half-window of the figure (default: distance + 150 kb, at least 300 kb)")
    p.add_argument("--top", type=int, default=20, help="rows printed (all are written)")
    p.add_argument("--gtf", default="data/gencode.gtf.gz", help="written into the commands")
    p.add_argument("--out")
    a = p.parse_args()

    d = Path(a.dir)
    interet = set(Path(a.interet).read_text().split()) if a.interet else None
    el = lire_elements(a.elements)
    loops = {t: boucles(d, t, interet) for t in TISSUS}
    if not loops[a.tissu]:
        sys.exit(f"aucune boucle dans {d}/loops_{a.tissu}.tsv")
    # tissues where the gene has a loop on an element of that same tissue
    avec = defaultdict(set)
    for t, rows in loops.items():
        for r in rows:
            s = int(num(r["summit"], 0))
            if any(de_tissu(e[3], t) for e in sur_ancre(el.get(r["gene"], []), r["chrom"], s, a.tol)):
                avec[r["gene"]].add(t)
    rat = ratios(d, a.tissu)
    rows = loops[a.tissu]
    par_gene = defaultdict(list)
    for r in rows:
        par_gene[r["gene"]].append(r)
    meilleurs = []
    for g, rs in par_gene.items():
        best = None
        for r in rs:
            vp, s = int(num(r["vp_pos"], 0)), int(num(r["summit"], 0))
            dist = abs(s - vp)
            on = sur_ancre(el.get(g, []), r["chrom"], s, a.tol)
            n_t = max((de_tissu(e[3], a.tissu) for e in on), default=0)
            if not n_t:
                continue
            n_all = max((len(e[3]) for e in on), default=0)
            nf, npr = num(r["n_fichiers"], 1), num(r["n_presents"])
            pres = npr / nf if nf and not math.isnan(npr) else 1.0
            sr = rat.get((g, s), float("nan"))
            fen = a.fenetre or max(300_000, dist + 150_000)
            dans = sum(1 for x in rows if x["chrom"] == r["chrom"]
                       and abs(int(num(x["vp_pos"], 0)) - vp) <= fen
                       and abs(int(num(x["summit"], 0)) - vp) <= fen)
            autres = dans - len(rs)
            score = (pres + min(n_t, 3) / 3
                     + (min(math.log2(sr), 2) / 2 if sr > 1 else 0)
                     + (1 if a.dmin <= dist <= a.dmax else 0)
                     - 0.25 * max(0, len(rs) - 2) - 0.1 * max(0, autres - 3))
            if best is None or score > best[0]:
                acc = (r.get("fichiers", "").split("|") or [""])[0]
                cmd = (f"python3 scripts/locus_figure.py --tissues {d}/{a.tissu}/libraries.tsv "
                       f"--term-column accession --term {acc} --gene {g} --gtf {a.gtf} "
                       f"--loops {d}/{a.tissu}/loops.tsv --arcs --genes-track "
                       f"--elements {a.elements} --elements-tissu {a.tissu} "
                       f"--window {fen} --max-distance {fen} --norm SCALE "
                       f"--out figures/locus_{a.tissu}_{g}")
                best = (score, [g, r["chrom"], vp, s, f"{dist / 1000:.0f}", r["n_fichiers"],
                                r["n_presents"], f"{sr:.2f}", n_t, n_all, len(rs), dans,
                                ",".join(sorted(avec[g])), f"{score:.2f}", cmd])
        if best:
            meilleurs.append(best)
    meilleurs.sort(key=lambda x: -x[0])
    head = ["gene", "chrom", "vp_pos", "summit", "distance_kb", "n_fichiers", "n_presents",
            "summit_ratio", "elements_tissu", "elements_tous", "boucles_gene",
            "boucles_fenetre", "tissus_avec", "score", "commande"]
    out = Path(a.out) if a.out else d / f"locus_{a.tissu}.tsv"
    with open(out, "w") as fo:
        fo.write("\t".join(head) + "\n")
        for _, x in meilleurs:
            fo.write("\t".join(map(str, x)) + "\n")
    print(f"{a.tissu} : {len(par_gene)} promoteurs avec boucle, {len(meilleurs)} dont une "
          f"boucle finit sur un élément rE2G du même tissu -> {out}", file=sys.stderr)
    print("gene\tdist_kb\tpresents\tratio\tbios_tissu\tboucles\tfenetre\ttissus\tscore", file=sys.stderr)
    for _, x in meilleurs[:a.top]:
        print("\t".join(map(str, [x[0], x[4], f"{x[6]}/{x[5]}", x[7], x[8], x[10], x[11],
                                  x[12], x[13]])), file=sys.stderr)


if __name__ == "__main__":
    main()
