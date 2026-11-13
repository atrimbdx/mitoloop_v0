#!/usr/bin/env python3
"""
loops_par_tissu.py - loops of one tissue, called on that tissue's own Hi-C
libraries pooled as replicates (peakC combined.analysis through
v4c_loops.R call / reciprocal), one loop set per tissue group.

The libraries of a tissue group (ENCODE terms matched by a regular
expression: muscle, heart, brain, kidney) are treated as replicates of one
condition: each group gets its own library table, with the group as its term,
and the usual `call` then `reciprocal` steps of v4c_loops.R are run on it in
the group's own folder (separately per tissue, so that loop coordinates never
merge across tissues). A group with
one library is called with single.analysis (flagged n_fichiers = 1). No
pooling with other tissues and no use of the multi-tissue loop catalogue.

Subcommands
  librairies   DIR/<tissu>/libraries.tsv: the kept libraries of each group,
               term replaced by the group name (original term in
               term_encode); summary per group on stderr
  partenaires  partner anchors of the per-tissue calls (DIR/*/partners.tsv of
               v4c_loops.R call) missing from an existing partner export:
               partenaires_manquants.tsv, to extract with v4c_pipeline.py
  fusionner    merges two export folders library by library (rows of the
               second added when their name is new), so that `reciprocal`
               reads the old and the new partner profiles together
  separer      DIR/loops_<tissu>.tsv from each DIR/<tissu>/loops.tsv, with the
               number of the group's libraries; with --peaks-lib (v4c_loops.R call_lib, each
               library alone) the number of those libraries in which the loop
               is also found on its own (n_presents, control only, no filter)

Usage (mitoloop_v0 ; T = results/loops_tissus ; slurm/loops_tissus.sbatch runs
the R steps for the four tissues)
  python3 scripts/loops_par_tissu.py librairies --libraries libraries.tsv --out $T
  Rscript scripts/v4c_loops.R call --libraries $T/<tissu>/libraries.tsv \
      --export <export/both> --wsize 3 --out $T/<tissu>          (each tissue)
  python3 scripts/loops_par_tissu.py partenaires --dir $T --export <export/partners> \
      --out $T/partenaires_manquants.tsv
  (v4c_pipeline.py extract + export of the missing partners, then fusionner)
  Rscript scripts/v4c_loops.R reciprocal --libraries $T/<tissu>/libraries.tsv \
      --export <merged partner export> --wsize 3 --out $T/<tissu>   (each tissue)
  python3 scripts/loops_par_tissu.py separer --dir $T --genes genes.txt [--peaks-lib peaks_lib.tsv]
"""

import argparse
import csv
import gzip
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


def groups_of(a):
    return {k: (v, "") for k, v in (g.split("=", 1) for g in a.groupe)} if a.groupe else GROUPES


def group_of(term, groups):
    for k, (inc, exc) in groups.items():
        if re.search(inc, term, re.I) and not (exc and re.search(exc, term, re.I)):
            return k
    return None


def cmd_librairies(a):
    groups = groups_of(a)
    rows = read_tsv(a.libraries)
    out = Path(a.out)
    cols = list(rows[0]) + (["term_encode"] if "term_encode" not in rows[0] else [])
    n = defaultdict(list)
    for r in rows:
        if r.get("kept", "1") not in ("1", "TRUE", "True"):
            continue
        g = group_of(r["term"], groups)
        if g is not None:
            n[g].append(dict(r, term_encode=r["term"], term=g))
    for g in groups:
        terms = sorted({r["term_encode"] for r in n[g]})
        print(f"{g:14s} {len(n[g])} librairies"
              + (f" ({', '.join(terms)})" if terms else " : aucune")
              + (" -> appel sur une seule librairie (single.analysis)" if len(n[g]) == 1 else ""),
              file=sys.stderr)
        if not n[g]:
            continue
        (out / g).mkdir(parents=True, exist_ok=True)
        with open(out / g / "libraries.tsv", "w") as fo:
            fo.write("\t".join(cols) + "\n")
            for r in n[g]:
                fo.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"-> {out}/<tissu>/libraries.tsv (un appel call + reciprocal par tissu)", file=sys.stderr)


def export_names(path):
    """noms des points de vue d'un fichier d'export (première colonne)"""
    with gzip.open(path, "rt") as fh:
        fh.readline()
        return {ln.split("\t", 1)[0] for ln in fh}


def cmd_partenaires(a):
    # partenaires des pics des gènes voulus (--genes : sans les gènes témoins de l'export)
    want = set(Path(a.genes).read_text().split()) if a.genes else None
    keep_pos = None
    if want is not None:
        keep_pos = set()
        for f in sorted(Path(a.dir).glob("*/peaks_mito.tsv")):
            for r in read_tsv(f):
                if r["gene"] in want:
                    keep_pos.add((r["chrom"], str(int(float(r["summit"])))))
    parts, names = [], set()
    for f in sorted(Path(a.dir).glob("*/partners.tsv")):      # partenaires de tous les tissus
        for p in read_tsv(f):
            if keep_pos is not None and (p["chrom"], str(int(float(p["pos"])))) not in keep_pos:
                continue
            if p["name"] not in names:
                names.add(p["name"])
                parts.append(p)
    files = sorted(Path(a.export).glob("*.tsv.gz"))
    have = export_names(files[0]) if files else set()
    miss = [p for p in parts if p["name"] not in have]
    with open(a.out, "w") as fo:
        fo.write("name\tchrom\tpos\tset\n")
        for p in miss:
            fo.write(f"{p['name']}\t{p['chrom']}\t{p['pos']}\t{p.get('set', 'partner')}\n")
    print(f"{len(parts)} partenaires (tous tissus{', gènes de --genes' if want else ''}), "
          f"{len(parts) - len(miss)} déjà dans {a.export}, "
          f"{len(miss)} à extraire -> {a.out}", file=sys.stderr)


def cmd_fusionner(a):
    d1, d2, out = Path(a.export[0]), Path(a.export[1]), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    for f1 in sorted(d1.glob("*.tsv.gz")):
        f2 = d2 / f1.name
        with gzip.open(f1, "rt") as fh:
            head = fh.readline()
            lines = fh.readlines()
        names = {ln.split("\t", 1)[0] for ln in lines}
        if f2.exists():
            with gzip.open(f2, "rt") as fh:
                h2 = fh.readline()
                if h2 != head:
                    sys.exit(f"{f2} : colonnes différentes de {f1}")
                lines += [ln for ln in fh if ln.split("\t", 1)[0] not in names]
        with gzip.open(out / f1.name, "wt") as fo:
            fo.write(head)
            fo.writelines(lines)
        n += 1
    for f in (d1 / "decay.tsv", d2 / "decay.tsv"):
        if f.exists() and not (out / "decay.tsv").exists():
            (out / "decay.tsv").write_text(f.read_text())
    print(f"{n} librairies fusionnées -> {out}", file=sys.stderr)


def cmd_separer(a):
    groups = groups_of(a)
    d = Path(a.dir)
    by_g, term_enc, loops_g = defaultdict(list), {}, {}
    want = set(Path(a.genes).read_text().split()) if a.genes else None
    for g in groups:
        if (d / g / "libraries.tsv").exists():
            for r in read_tsv(d / g / "libraries.tsv"):
                by_g[g].append(r["accession"])
                term_enc[r["accession"]] = r.get("term_encode", r["term"])
        lf = d / g / "loops.tsv"
        loops_g[g] = [r for r in read_tsv(lf) if want is None or r["gene"] in want] if lf.exists() else []
        if by_g[g] and not lf.exists():
            print(f"{g} : {lf} absent (call + reciprocal pas encore faits)", file=sys.stderr)
    seen = defaultdict(list)                     # (gene, groupe) -> [(start, end, accession)]
    if a.peaks_lib:
        for p in read_tsv(a.peaks_lib):
            g = group_of(term_enc.get(p["accession"], p.get("term", "")), groups)
            if g:
                seen[(p["gene"], g)].append((float(p["start"]), float(p["end"]), p["accession"]))
    out = d
    resume = []
    for g in groups:
        accs = by_g.get(g, [])
        rows = []
        for r in loops_g[g]:
            s, e = float(r["start"]), float(r["end"])
            pres = sorted({x for ps, pe, x in seen.get((r["gene"], g), [])
                           if ps <= e + a.tol and s <= pe + a.tol}) if a.peaks_lib else []
            rows.append([r["gene"], r["chrom"], r["vp_pos"], r["summit"], r["start"], r["end"],
                         len(accs), len(pres) if a.peaks_lib else ".",
                         "|".join(pres) if a.peaks_lib else "|".join(accs),
                         "|".join(sorted({term_enc[x] for x in accs})), r["loop_id"]])
        with open(out / f"loops_{g}.tsv", "w") as fo:
            fo.write(HEAD)
            for x in rows:
                fo.write("\t".join(map(str, x)) + "\n")
        genes = {x[0] for x in rows}
        resume.append([g, len(accs), len(rows), len(genes)])
        print(f"{g:14s} {len(accs)} librairies ; {len(rows)} boucles réciproques à "
              f"{len(genes)} promoteurs" + ("  (une seule librairie)" if len(accs) == 1 else ""),
              file=sys.stderr)
    with open(out / "resume_tissus.tsv", "w") as fo:
        fo.write("tissu\tlibrairies\tboucles\tpromoteurs\n")
        for x in resume:
            fo.write("\t".join(map(str, x)) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--groupe", action="append", default=[],
                   help="NAME=REGEX on the ENCODE term, repeatable (replaces the default groups)")
    sub = p.add_subparsers(dest="cmd")          # Python 3.6 (Curta) : pas de required=
    s = sub.add_parser("librairies")
    s.add_argument("--libraries", required=True)
    s.add_argument("--out", required=True)
    s = sub.add_parser("partenaires")
    s.add_argument("--dir", required=True, help="dossier des tissus (<tissu>/partners.tsv)")
    s.add_argument("--export", required=True, help="existing partner export folder")
    s.add_argument("--genes", default="", help="viewpoints kept (genes.txt : mitochondrial only)")
    s.add_argument("--out", required=True)
    s = sub.add_parser("fusionner")
    s.add_argument("--export", nargs=2, required=True, help="two export folders")
    s.add_argument("--out", required=True)
    s = sub.add_parser("separer")
    s.add_argument("--dir", required=True, help="dossier des tissus (<tissu>/loops.tsv)")
    s.add_argument("--genes", default="")
    s.add_argument("--peaks-lib", default="", help="peaks_lib.tsv of v4c_loops.R call_lib (control)")
    s.add_argument("--tol", type=float, default=20e3)
    a = p.parse_args()
    if not a.cmd:
        p.error("sous-commande attendue : librairies, partenaires, fusionner ou separer")
    {"librairies": cmd_librairies, "partenaires": cmd_partenaires, "fusionner": cmd_fusionner,
     "separer": cmd_separer}[a.cmd](a)


if __name__ == "__main__":
    main()
