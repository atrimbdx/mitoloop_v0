#!/usr/bin/env python3
"""
mitoloop_paquet.py — Curta. Prépare le paquet MitoLoop du pipeline reannot :
les données que ajout_mitoloop.py lit sur la VM, réduites aux régions
régulatrices des gènes mitochondriaux.

Entrées (produites par run_mcia.sh de mitoloop_v0)
  --genes        ensemble de gènes (genes.txt)
  --l1           l1.tsv de L1 sans ajustement du contenu mitochondrial
                 (results/l1_raw/l1.tsv) : l1_class, score, ci_low, ci_high,
                 status, complex, protein_coupling, rna_coupling
  --fibroblast   optionnel, patient_genes.tsv de la validation (gene,
                 fibroblast_r)
  --regions      BED de `vcf_annotate.py regions` (chrom, start, end, gene,
                 class, ...) ; les classes cds, exon_nc et gene_body sont
                 écartées
  --epd          EPDnew BED (TSS de -49 à +10, séparé par des espaces)
  --cage         pics CAGE FANTOM5 (BED)
  --remap TF=FICHIER    pics ReMap non redondants d'un facteur, répétable
  --chromhmm TISSU=FICHIER   segmentation ChromHMM (état en 4e colonne),
                 répétable
  --jaspar       matrices JASPAR (format jaspar) ; --tfs : facteurs gardés
  --margin       pb autour des régions pour les extraits (défaut 1000)

Sortie (--out) : genes.tsv, regions.bed.gz, annotations/{epd,cage,
remap_<TF>,chromhmm_<tissu>}.bed.gz, motifs.jaspar, MANIFEST.tsv (source et
nombre de lignes de chaque fichier), puis, copiés par 13_mitoloop_paquet.sh,
les scripts VM et SHA256SUMS.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path

DROP = {"cds", "exon_nc", "gene_body"}


def opener(p):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def fields(line):
    line = line.rstrip("\n")
    return line.split("\t") if "\t" in line else line.split()


def chrom_key(c):
    c = c[3:] if c.startswith("chr") else c
    return "MT" if c in ("M", "MT") else c


def read_table(path):
    with opener(path) as fh:
        head = fh.readline().rstrip("\n").lstrip("#").split("\t")
        return [dict(zip(head, ln.rstrip("\n").split("\t"))) for ln in fh if ln.strip()]


class Windows:
    """fenêtres (régions +/- marge) pour filtrer les annotations"""

    def __init__(self, regions, margin):
        d = defaultdict(list)
        for c, s, e in regions:
            d[chrom_key(c)].append((max(0, s - margin), e + margin))
        self.d = {}
        for c, v in d.items():
            v.sort()
            merged = []
            for s, e in v:
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            self.d[c] = ([m[0] for m in merged], [m[1] for m in merged])

    def hit(self, c, s, e):
        x = self.d.get(chrom_key(c))
        if not x:
            return False
        st, en = x
        i = bisect_right(st, e - 1) - 1
        return i >= 0 and en[i] > s


def write_bed(path, rows):
    rows.sort(key=lambda r: (chrom_key(r[0]), int(r[1])))
    with gzip.open(path, "wt") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return len(rows)


def subset(path, win, convert):
    out = []
    with opener(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = fields(line)
            if len(f) >= 4 and f[0].isdigit() and not f[1].isdigit():
                f = f[1:]                                   # colonne bin UCSC
            try:
                s, e = int(f[1]), int(f[2])
            except (IndexError, ValueError):
                continue
            if win.hit(f[0], s, e):
                r = convert(f, s, e)
                if r:
                    out.append(r)
    return out


def l1_class(r, high, low):
    """la règle de l1_score.py, pour un l1.tsv qui n'a pas la colonne"""
    if not str(r.get("status", "")).startswith("estimable"):
        return "not_estimable"
    try:
        lo, hi = float(r.get("ci_low")), float(r.get("ci_high"))
    except (TypeError, ValueError):
        return "uncertain"
    if lo != lo or hi != hi:
        return "uncertain"
    return "transmitted" if lo > high else "not_transmitted" if hi < low \
        else "uncertain"


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--genes", required=True, type=Path)
    p.add_argument("--l1", required=True, type=Path)
    p.add_argument("--fibroblast", type=Path)
    p.add_argument("--regions", required=True, type=Path)
    p.add_argument("--epd", type=Path)
    p.add_argument("--cage", type=Path)
    p.add_argument("--remap", action="append", default=[])
    p.add_argument("--chromhmm", action="append", default=[])
    p.add_argument("--jaspar", type=Path)
    p.add_argument("--tfs", default="NRF1,GABPA,YY1,ESRRA,SP1,NFYA,CREB1,CTCF,TBP")
    p.add_argument("--margin", type=int, default=1000)
    p.add_argument("--class-high", type=float, default=0.5)
    p.add_argument("--class-low", type=float, default=0.3)
    p.add_argument("--out", required=True, type=Path)
    a = p.parse_args()
    out = a.out
    (out / "annotations").mkdir(parents=True, exist_ok=True)
    manifest = []

    # gènes : ensemble complet, L1 quand il existe
    gs = [g for g in re.split(r"\s+", a.genes.read_text()) if g and not g.startswith("#")]
    l1 = {r["gene"]: r for r in read_table(a.l1)}
    fib = {r["gene"]: r for r in read_table(a.fibroblast)} \
        if a.fibroblast and a.fibroblast.exists() else {}
    regs = [r for r in read_table(a.regions) if r.get("class") not in DROP]
    in_regions = {r["gene"] for r in regs}
    cols = ["gene", "gene_id", "l1_class", "l1", "ci_low", "ci_high", "status",
            "complex", "protein_coupling", "rna_coupling", "fibroblast_r"]
    n_cls = defaultdict(int)
    with open(out / "genes.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for g in sorted(set(gs) | in_regions):
            r = l1.get(g, {})
            cls = r.get("l1_class") or l1_class(r, a.class_high, a.class_low)
            n_cls[cls] += 1
            row = [g, r.get("gene_id", ""), cls, r.get("score", ""),
                   r.get("ci_low", ""), r.get("ci_high", ""),
                   r.get("status", "absent_du_proteome"), r.get("complex", ""),
                   r.get("protein_coupling", ""), r.get("rna_coupling", ""),
                   fib.get(g, {}).get("fibroblast_r", "")]
            fh.write("\t".join("" if x in ("NA", "nan") else str(x) for x in row) + "\n")
    print("classes L1 : " + ", ".join(f"{k} {v}" for k, v in sorted(n_cls.items())))
    manifest.append(("genes.tsv", sum(n_cls.values()), str(a.l1)))

    # régions
    rb = [(r["chrom"], int(r["start"]), int(r["end"]), r["gene"], r["class"])
          for r in regs]
    n = write_bed(out / "regions.bed.gz", [list(x) for x in rb])
    print(f"régions : {n} ({', '.join(sorted({x[4] for x in rb}))})")
    manifest.append(("regions.bed.gz", n, str(a.regions)))
    win = Windows([(c, s, e) for c, s, e, _g, _k in rb], a.margin)

    def put(name, src, conv):
        if not src or not Path(src).exists():
            print(f"{name} : {src} absent, ignoré")
            return
        rows = subset(src, win, conv)
        n = write_bed(out / "annotations" / f"{name}.bed.gz", rows)
        print(f"{name} : {n} intervalles dans les régions")
        manifest.append((f"annotations/{name}.bed.gz", n, str(src)))

    def epd(f, s, e):
        strand = f[5] if len(f) > 5 else "+"
        tss = s + 49 if strand == "+" else e - 50          # EPDnew : -49..+10
        return [f[0], tss, tss + 1, strand]

    put("epd", a.epd, epd)
    put("cage", a.cage, lambda f, s, e: [f[0], s, e, "."])
    for spec in a.remap:
        tf, path = spec.split("=", 1)

        def remap(f, s, e):
            nm = f[3] if len(f) > 3 else ""
            bio = nm.split(":", 1)[1] if ":" in nm else nm
            return [f[0], s, e, len({b for b in bio.split(",") if b}) or 1]
        put(f"remap_{tf.upper()}", path, remap)
    for spec in a.chromhmm:
        tissue, path = spec.split("=", 1)
        put(f"chromhmm_{tissue}", path, lambda f, s, e: [f[0], s, e, f[3]])

    # motifs JASPAR des facteurs gardés
    want = {t.upper() for t in a.tfs.split(",") if t}
    kept = []
    if a.jaspar and a.jaspar.exists():
        block, name = [], None
        with opener(a.jaspar) as fh, open(out / "motifs.jaspar", "w") as fo:
            for line in list(fh) + [">END"]:
                if line.startswith(">"):
                    if name and name.upper() in want:
                        fo.write("".join(block))
                        kept.append(name.upper())
                    f = line[1:].split()
                    name = f[1] if len(f) > 1 else (f[0] if f else None)
                    block = [line]
                else:
                    block.append(line)
        print(f"motifs JASPAR : {', '.join(kept) or 'aucun'} "
              f"(absents : {', '.join(sorted(want - set(kept))) or 'aucun'})")
        manifest.append(("motifs.jaspar", len(kept), str(a.jaspar)))
    else:
        print("JASPAR absent : pas de colonne MITO_motifs")
    with open(out / "MANIFEST.tsv", "w") as fh:
        fh.write("fichier\tlignes\tsource\n")
        for r in manifest:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"Paquet MitoLoop : {out}")


if __name__ == "__main__":
    sys.exit(main())
