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
  --interet    optional list of genes to keep (MITO list)
  --frac-tissu an element counts when it is linked to the gene in at least
               this share of the tissue's rE2G biosamples (default 0.8, the
               rule of the reannot paquet and of the AlphaGenome targets)
  --l1         optional l1.tsv of `l1_score.py score` (results/l1_raw/l1.tsv,
               the one of the reannot paquet): only the genes whose l1_class
               is in --l1-classes are kept (default transmitted: RNA changes
               reach the protein, so a regulatory variant at the locus can
               matter); --l1-classes tous keeps every gene and only adds the
               class. A l1.tsv without l1_class is classed by the rule of
               l1_score.py (95% interval above --class-high: transmitted,
               below --class-low: not_transmitted)

Output (--out, default <dir>/locus_<tissu>.tsv), one row per gene, best loop
first, columns:
  gene l1_class l1 chrom vp_pos summit distance_kb n_fichiers n_presents summit_ratio
  elements_tissu   biosamples of the tissue linking the element on the anchor
  elements_tous    biosamples of any tissue linking it
  boucles_gene     loops of this gene in the tissue
  boucles_fenetre  loops of any gene with both ends in the window
  tissus_avec      tissues where this gene has a loop on a same-tissue element
  score, commande  (the locus_figure.py call that draws it)

Ranking (by gene, its best loop; stated, not tuned on anything)
  required   the loop anchor (peakC region start-end widened by --tol) holds an element linked to the
             gene in >= --frac-tissu of the tissue's biosamples
  score      presence (n_presents / n_fichiers, 1 when not measured)
             + share of the tissue's biosamples linking the element
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


def du_tissu(ct, tissu):
    inc, exc = CELLTYPES[tissu]
    return bool(re.search(inc, ct, re.I)) and not re.search(exc, ct, re.I)


def de_tissu(bios, tissu):
    """biosamples (fichiers distincts) du tissu qui relient l'élément"""
    return len({f for f, ct in bios if du_tissu(ct, tissu)})


def lire_elements(path):
    """gene -> [(chrom, start, end, [(fichier, CellType)])] et, par tissu, le nombre de
    biosamples (fichiers distincts de la colonne biosamples sur toute la table)"""
    el, tous = defaultdict(list), set()
    for r in read_tsv(path):
        if r.get("biosamples"):
            bios = [tuple(x.split("=", 1)) for x in r["biosamples"].split(";") if "=" in x]
        else:                                   # ancien format : un CellType par biosample
            bios = [(f"{r['gene']}:{i}", x) for i, x in
                    enumerate(x for x in r.get("celltypes", "").split(";") if x)]
        tous.update(bios)
        el[r["gene"]].append((r["chrom"], int(num(r["start"], 0)), int(num(r["end"], 0)), bios))
    taille = {t: len({f for f, ct in tous if du_tissu(ct, t)}) for t in TISSUS}
    return el, taille


def sur_ancre(elements, chrom, summit, tol, start=None, end=None):
    """elements of the gene on the loop anchor: overlapping the peakC region of the loop
    (start-end, summit included) widened by tol (one bin, the resolution of the call) on
    each side; the rule of the reannot AG2 campaign (ag2_balayage.py ancre_region)"""
    lo = min(summit, start) if start is not None else summit
    hi = max(summit, end) if end is not None else summit
    return [e for e in elements if e[0] == chrom and e[1] < hi + tol and lo - tol < e[2]]


def region(r):
    """(start, end) of a loop row, None when absent"""
    s, e = num(r.get("start")), num(r.get("end"))
    if math.isnan(s) or math.isnan(e):
        return None, None
    return int(s), int(e)


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


def classe_l1(r, high, low):
    """the rule of l1_score.py, for a l1.tsv without the l1_class column"""
    if r.get("l1_class"):
        return r["l1_class"]
    if not str(r.get("status", "")).startswith("estimable"):
        return "not_estimable"
    lo, hi = num(r.get("ci_low")), num(r.get("ci_high"))
    if math.isnan(lo) or math.isnan(hi):
        return "uncertain"
    return "transmitted" if lo > high else "not_transmitted" if hi < low else "uncertain"


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", required=True)
    p.add_argument("--tissu", required=True, choices=TISSUS)
    p.add_argument("--elements", required=True)
    p.add_argument("--interet", help="genes to keep (whitespace-separated)")
    p.add_argument("--l1", help="l1.tsv (l1_score.py score)")
    p.add_argument("--l1-classes", default="transmitted",
                   help="classes kept, comma-separated, or 'tous'")
    p.add_argument("--class-high", type=float, default=0.5)
    p.add_argument("--class-low", type=float, default=0.3)
    p.add_argument("--frac-tissu", type=float, default=0.8)
    p.add_argument("--tol", type=int, default=10_000,
                   help="bp added on each side of the loop region (one bin, the call resolution)")
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
    l1 = {}
    if a.l1:
        for r in read_tsv(a.l1):
            l1[r["gene"]] = (classe_l1(r, a.class_high, a.class_low),
                             r.get("l1") or r.get("score", ""))
        if a.l1_classes != "tous":
            garder = {x.strip() for x in a.l1_classes.split(",") if x.strip()}
            sel = {g for g, (c, _) in l1.items() if c in garder}
            print(f"L1 : {len(sel)} gènes de classe {', '.join(sorted(garder))} "
                  f"sur {len(l1)}", file=sys.stderr)
            interet = sel if interet is None else interet & sel
    el, taille = lire_elements(a.elements)
    seuil = {t: max(1, math.ceil(a.frac_tissu * taille[t])) for t in TISSUS}
    print("biosamples rE2G par tissu (seuil) : "
          + ", ".join(f"{t} {taille[t]} ({seuil[t]})" for t in TISSUS), file=sys.stderr)
    loops = {t: boucles(d, t, interet) for t in TISSUS}
    if not loops[a.tissu]:
        sys.exit(f"aucune boucle dans {d}/loops_{a.tissu}.tsv")
    # tissues where the gene has a loop on an element of that same tissue
    avec = defaultdict(set)
    for t, rows in loops.items():
        for r in rows:
            s = int(num(r["summit"], 0))
            if any(de_tissu(e[3], t) >= seuil[t]
                   for e in sur_ancre(el.get(r["gene"], []), r["chrom"], s, a.tol, *region(r))):
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
            on = sur_ancre(el.get(g, []), r["chrom"], s, a.tol, *region(r))
            n_t = max((de_tissu(e[3], a.tissu) for e in on), default=0)
            if n_t < seuil[a.tissu]:
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
            score = (pres + n_t / max(1, taille[a.tissu])
                     + (min(math.log2(sr), 2) / 2 if sr > 1 else 0)
                     + (1 if a.dmin <= dist <= a.dmax else 0)
                     - 0.25 * max(0, len(rs) - 2) - 0.1 * max(0, autres - 3))
            if best is None or score > best[0]:
                acc = (r.get("fichiers", "").split("|") or [""])[0]
                cmd = (f"python3 scripts/locus_figure.py --tissues {d}/{a.tissu}/libraries.tsv "
                       f"--term-column accession --term {acc} --gene {g} --gtf {a.gtf} "
                       f"--loops {d}/{a.tissu}/loops.tsv --arcs --genes-track "
                       f"--elements {a.elements} --elements-tissu {a.tissu} --elements-frac {a.frac_tissu} "
                       f"--window {fen} --max-distance {fen} --norm SCALE "
                       f"--out figures/locus_{a.tissu}_{g}")
                c1, v1 = l1.get(g, ("." if not a.l1 else "absent", ""))
                try:
                    v1 = f"{float(v1):.2f}"
                except ValueError:
                    v1 = v1 or "."
                best = (score, [g, c1, v1, r["chrom"], vp, s, f"{dist / 1000:.0f}", r["n_fichiers"],
                                r["n_presents"], f"{sr:.2f}", n_t, n_all, len(rs), dans,
                                ",".join(sorted(avec[g])), f"{score:.2f}", cmd])
        if best:
            meilleurs.append(best)
    meilleurs.sort(key=lambda x: -x[0])
    head = ["gene", "l1_class", "l1", "chrom", "vp_pos", "summit", "distance_kb", "n_fichiers", "n_presents",
            "summit_ratio", "elements_tissu", "elements_tous", "boucles_gene",
            "boucles_fenetre", "tissus_avec", "score", "commande"]
    out = Path(a.out) if a.out else d / f"locus_{a.tissu}.tsv"
    with open(out, "w") as fo:
        fo.write("\t".join(head) + "\n")
        for _, x in meilleurs:
            fo.write("\t".join(map(str, x)) + "\n")
    print(f"{a.tissu} : {len(par_gene)} promoteurs avec boucle, {len(meilleurs)} dont une "
          f"boucle finit sur un élément rE2G du même tissu → {out}", file=sys.stderr)
    print("gene\tL1\tdist_kb\tpresents\tratio\tbios_tissu\tboucles\tfenetre\ttissus\tscore", file=sys.stderr)
    for _, x in meilleurs[:a.top]:
        print("\t".join(map(str, [x[0], x[1], x[6], f"{x[8]}/{x[7]}", x[9], x[10], x[12],
                                  x[13], x[14], x[15]])), file=sys.stderr)


if __name__ == "__main__":
    main()
