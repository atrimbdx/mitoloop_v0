#!/usr/bin/env python3
"""
loops_par_tissu.py - loops of one tissue, replicated across that tissue's own
Hi-C libraries.

Input: peaks_lib.tsv from `v4c_loops.R call_lib` (peakC on each library on its
own) and the library table. For every tissue group (ENCODE terms matched by a
regular expression), the peaks of its libraries are matched by summit (within
--tol bp, same viewpoint) and a loop is kept when it is found in

    1 library  : its single library (flagged n_fichiers = 1)
    2 libraries: both (2/2)
    3 libraries: all three (3/3)
    > 3        : at least 3

No pooling with other tissues and no use of the multi-tissue loop catalogue:
a tissue's loops come from that tissue's libraries only. Reciprocity (the
peak seen again from the partner anchor) is not required here.

Outputs (--out): loops_<tissue>.tsv (gene, chrom, vp_pos, summit, start, end,
n_fichiers, n_presents, fichiers, termes, loop_id) and resume_tissus.tsv
(libraries, threshold, loops and viewpoints per tissue).

Usage
  python3 scripts/loops_par_tissu.py --peaks results/loops/peaks_lib.tsv \\
      --libraries results/libraries.tsv --genes genes.txt --out results/loops_tissus
  (--groupe NAME=REGEX, repeatable, replaces the default groups)
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

# ENCODE Hi-C terms of each tissue group: (included, excluded)
GROUPES = {
    "muscle": (r"gastrocnemius|psoas|skeletal muscle|muscle of (the )?(leg|arm|trunk|back)|"
               r"limb muscle|myotube|myoblast|quadriceps|vastus", r"smooth|cardiac|heart"),
    "cardiomyocyte": (r"heart|cardiac|ventric|atri(um|al)", r"smooth"),
    "neurone": (r"brain|cortex|gyrus|caudate|putamen|cerebell|hippocamp|neuron|"
                r"nucleus accumbens|substantia nigra|amygdala",
                r"adrenal|kidney|smooth|vascul|endotheli|pericyte"),
    "rein": (r"kidney|\brenal|nephr|tubul", r"adrenal"),     # « adrenal » contient « renal »
}


HEAD = ("gene\tchrom\tvp_pos\tsummit\tstart\tend\tn_fichiers\tn_presents\t"
        "fichiers\ttermes\tloop_id\n")


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def requis(n):
    """libraries in which a loop must be found, among n libraries of the tissue"""
    return n if n <= 3 else 3


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--peaks", required=True, help="peaks_lib.tsv (v4c_loops.R call_lib)")
    p.add_argument("--libraries", required=True, help="libraries.tsv (column kept used)")
    p.add_argument("--genes", default="", help="viewpoints kept (default: all)")
    p.add_argument("--groupe", action="append", default=[],
                   help="NAME=REGEX on the ENCODE term, repeatable (replaces the defaults)")
    p.add_argument("--tol", type=float, default=20e3, help="bp between summits of one loop")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    groups = {k: (v, "") for k, v in (g.split("=", 1) for g in a.groupe)} if a.groupe else GROUPES
    libs = [r for r in read_tsv(a.libraries) if r.get("kept", "1") in ("1", "TRUE", "True")]
    term_of = {r["accession"]: r["term"] for r in libs}
    want = set(Path(a.genes).read_text().split()) if a.genes else None
    peaks = [r for r in read_tsv(a.peaks) if r["accession"] in term_of
             and (want is None or r["gene"] in want)]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    resume = []
    for name, (inc, exc) in groups.items():
        ri, re_ = re.compile(inc, re.I), re.compile(exc, re.I) if exc else None
        accs = sorted(x for x, t in term_of.items() if ri.search(t) and not (re_ and re_.search(t)))
        n = len(accs)
        if not n:
            print(f"{name:14s} aucune librairie Hi-C (termes : {inc})", file=sys.stderr)
            resume.append([name, 0, 0, 0, 0, ""])
            rows = []
            with open(out / f"loops_{name}.tsv", "w") as fo:
                fo.write(HEAD)
            continue
        need = requis(n)
        by_gene = defaultdict(list)
        for r in peaks:
            if r["accession"] in accs:
                by_gene[r["gene"]].append(r)
        rows = []
        for g, pk in sorted(by_gene.items()):
            pk.sort(key=lambda r: float(r["summit"]))
            clusters, cur = [], []
            for r in pk:          # chaining summits within tol, never wider than 2 x tol
                s = float(r["summit"])
                if cur and (s - float(cur[-1]["summit"]) > a.tol or
                            s - float(cur[0]["summit"]) > 2 * a.tol):
                    clusters.append(cur)
                    cur = []
                cur.append(r)
            if cur:
                clusters.append(cur)
            for cl in clusters:
                present = sorted({r["accession"] for r in cl})
                if len(present) < need:
                    continue
                best = max(cl, key=lambda r: float(r["summit_ratio"]))
                summits = sorted(float(r["summit"]) for r in cl)
                summit = int(summits[len(summits) // 2])
                rows.append([g, best["chrom"], best["vp_pos"], summit,
                             int(min(float(r["start"]) for r in cl)),
                             int(max(float(r["end"]) for r in cl)), n, len(present),
                             "|".join(present), "|".join(sorted({term_of[x] for x in present})),
                             f"{g}:{summit}"])
        with open(out / f"loops_{name}.tsv", "w") as fo:
            fo.write(HEAD)
            for r in rows:
                fo.write("\t".join(map(str, r)) + "\n")
        terms = sorted({term_of[x] for x in accs})
        resume.append([name, n, need, len(rows), len({r[0] for r in rows}), "|".join(terms)])
        print(f"{name:14s} {n} librairies ({', '.join(terms)}) ; seuil {need}/{n} ; "
              f"{len(rows)} boucles à {len({r[0] for r in rows})} promoteurs", file=sys.stderr)
    with open(out / "resume_tissus.tsv", "w") as fo:
        fo.write("tissu\tlibrairies\tseuil\tboucles\tpromoteurs\ttermes\n")
        for r in resume:
            fo.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
