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
  --gtf          GENCODE GTF : identifiant Ensembl de chaque symbole
  --gnomad       contraintes par gène de gnomAD v4.1
                 (gnomad.v4.1.constraint_metrics.tsv) : transcrit MANE
                 Select, sinon canonique ; LoF observés et attendus, o/e,
                 LOEUF, pLI, o/e et z faux-sens
  --gismo        dossier (ou table) GISMO de Liao et al. 2026 : GISMO =
                 upper.ci et son décile publiés (Supplementary Table 2 :
                 lower.ci, median, upper.ci, gene, decile ; 1 = le moins de
                 pertes d'orthologue, le plus contraint ; virgule décimale
                 acceptée) ; GISMO-mis et son décile 'dec' de la table de
                 consensus des auteurs (Zenodo 20547722,
                 gismo-mis_consensus_genenames_*.txt, transcrit MANE Select),
                 ou d'une table publiée avec une colonne GISMO-mis
  --phylop       bigWig phyloP 241 mammifères (Zoonomia) ; ou
                 --phylop-bedgraph, un bedGraph déjà réduit : extrait sur les
                 régions, bgzip + tabix (annotations/phylop.bedgraph.gz)

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


GNOMAD_IN = {"gnomad_lof_obs": "lof.obs", "gnomad_lof_exp": "lof.exp",
             "gnomad_lof_oe": "lof.oe", "gnomad_loeuf": "lof.oe_ci.upper",
             "gnomad_pli": "lof.pLI", "gnomad_mis_oe": "mis.oe",
             "gnomad_mis_z": "mis.z_score"}
GNOMAD_OUT = list(GNOMAD_IN)
GISMO_OUT = ["gismo", "gismo_decile", "gismo_mis", "gismo_mis_decile"]


def gtf_ids(gtf):
    """symbole -> identifiant Ensembl (gènes)"""
    out = {}
    with opener(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "gene":
                continue
            m1 = re.search(r'gene_id "([^"]+)"', f[8])
            m2 = re.search(r'gene_name "([^"]+)"', f[8])
            if m1 and m2:
                out.setdefault(m2.group(1), m1.group(1).split(".")[0])
    return out


def gnomad_table(path):
    """identifiant Ensembl et symbole -> contraintes du transcrit MANE Select
    (sinon canonique) de gnomAD"""
    best = {}
    with opener(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        c = {h: i for i, h in enumerate(head)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            gid = f[c["gene_id"]]
            if not gid.startswith("ENSG"):
                continue
            rank = 2 * (f[c["mane_select"]] == "true") + (f[c["canonical"]] == "true")
            if gid in best and best[gid][0] >= rank:
                continue
            row = {k: f[c[v]] if v in c and f[c[v]] != "NA" else ""
                   for k, v in GNOMAD_IN.items()}
            for k in ("gnomad_lof_exp", "gnomad_lof_oe", "gnomad_loeuf",
                      "gnomad_pli", "gnomad_mis_oe", "gnomad_mis_z"):
                try:
                    row[k] = f"{float(row[k]):.3g}" if row[k] else ""
                except ValueError:
                    pass
            best[gid] = (rank, f[c["gene"]], row)
    out = {gid: r for gid, (_k, _g, r) in best.items()}
    out.update({g: r for _gid, (_k, g, r) in best.items()})
    return out


def gtf_mane(gtf):
    """transcrits MANE Select (ENST sans version)"""
    out = set()
    with opener(gtf) as fh:
        for line in fh:
            if "MANE_Select" not in line or line[0] == "#":
                continue
            f = line.split("\t", 9)
            m = re.search(r'transcript_id "(ENST\d+)', f[8]) if len(f) > 8 and f[2] == "transcript" else None
            if m:
                out.add(m.group(1))
    return out


def gismo_table(path, mane=frozenset()):
    """identifiant Ensembl et symbole -> GISMO, GISMO-mis et leurs déciles ;
    `mane` : transcrits préférés quand la table donne plusieurs transcrits par gène.
    GISMO et GISMO-mis viennent chacun de la meilleure table qui les porte :
    Supplementary Table 2 de Liao et al. (lower.ci, median, upper.ci, gene,
    decile ; GISMO = upper.ci, décile repris tel quel), une table avec des
    colonnes GISMO / GISMO-mis, ou la table de consensus des auteurs (GISMO-mis)"""
    import pandas as pd
    # gismo_genes*.tsv : recalcul de contrôle (gismo_scores.py), jamais pris pour la table publiée
    files = [path] if path.is_file() else sorted(
        f for f in path.rglob("*") if f.is_file() and re.search(
            r"\.(tsv|csv|txt|xlsx)(\.gz)?$", f.name, re.I)
        and not f.name.startswith("gismo_genes"))
    norm = lambda x: re.sub(r"[^a-z0-9]", "", str(x).lower())
    num = lambda s: pd.to_numeric(s.astype(str).str.strip().str.replace(",", ".", regex=False),
                                  errors="coerce")   # virgule décimale des tables collées d'Excel
    best = {}                                         # mesure -> (rang, fichier, onglet, table, valeur, décile, clé, symbole)
    for f in files:
        try:
            if f.name.lower().endswith(".xlsx"):
                tabs = pd.read_excel(f, sheet_name=None)
            else:
                sep = "," if ".csv" in f.name.lower() else "\t"
                tabs = {"": pd.read_csv(f, sep=sep, low_memory=False)}
        except Exception as e:  # tables of every kind
            print(f"GISMO : {f.name} illisible ({type(e).__name__})")
            continue
        for sheet, t in tabs.items():
            cols = {norm(c): c for c in t.columns}
            g = next((cols[k] for k in ("gismo", "gismoscore") if k in cols), None)
            m = next((cols[k] for k in ("gismomis", "gismomisscore") if k in cols), None)
            gd, md = cols.get("gismodecile"), cols.get("gismomisdecile")
            rank = 2
            if g is None and m is None and all(k in cols for k in ("upperci", "decile")):
                g, gd = cols["upperci"], cols["decile"]   # Supplementary Table 2
            elif g is None and m is None and all(k in cols for k in ("miscal", "dec", "genename")):
                m, md, rank = cols["miscal"], cols["dec"], 1   # consensus des auteurs (Zenodo)
                tx = next((c for c in t.columns if t[c].astype(str).head(200)
                           .str.match(r"^ENST\d{11}").mean() > 0.5), None)
                if tx is not None:   # transcrit MANE Select d'abord pour chaque gène
                    first = ~t[tx].astype(str).str.split(".").str[0].isin(mane)
                    t = t.assign(_o=first).sort_values("_o", kind="stable").drop(columns="_o")
            if g is None and m is None:
                continue
            key = next((c for c in t.columns if t[c].astype(str).head(200)
                        .str.match(r"^ENSG\d{11}").mean() > 0.5), None)
            sym = next((cols[k] for k in ("gene", "genename", "symbol", "genesymbol",
                                          "hgncsymbol") if k in cols and cols[k] != key), None)
            if key is None and sym is None:
                continue
            for name, col, dcol in (("gismo", g, gd), ("gismo_mis", m, md)):
                n = (rank, len(t))
                if col is not None and (name not in best or n > best[name][0]):
                    best[name] = (n, f, sheet, t, col, dcol, key, sym)
    if not best:
        print("GISMO : aucune table avec une colonne GISMO ou GISMO-mis")
        return {}
    res, by_sym = {}, {}
    for name in ("gismo", "gismo_mis"):
        if name not in best:
            print(f"GISMO : pas de table pour {name}")
            continue
        n, f, sheet, t, col, dcol, key, sym = best[name]
        print(f"GISMO : {name} <- {f.name}{' / ' + str(sheet) if sheet else ''}, {n[1]} lignes, "
              f"colonne {col!r}, décile {dcol or 'calculé'!r}, clé {key or sym!r}")
        x = num(t[col])
        if dcol is not None:      # déciles publiés : gardés tels quels
            d = num(t[dcol])
        else:
            d = pd.Series(pd.NA, index=t.index, dtype="Float64")
            ok = x.notna()
            d[ok] = pd.qcut(x[ok].rank(method="first"), 10, labels=False) + 1
        val = x.map(lambda v: "" if v != v else f"{v:.4g}")
        dec = d.map(lambda v: "" if pd.isna(v) else str(int(v)))
        for i in t.index:
            r = {name: val[i], name + "_decile": dec[i]}
            if key is not None and str(t.at[i, key]).startswith("ENSG"):
                res.setdefault(str(t.at[i, key]).split(".")[0], {}).update(r)
            if sym is not None and not pd.isna(t.at[i, sym]):
                by_sym.setdefault(str(t.at[i, sym]), {}).setdefault(name, r)
    for s_, parts in by_sym.items():
        rs = res.setdefault(s_, {})
        for r in parts.values():
            for k, v in r.items():
                rs.setdefault(k, v)
    return {k: {c: v.get(c, "") for c in GISMO_OUT} for k, v in res.items()}


def phylop_cut(bw, bedgraph, rb, path):
    """phyloP des régions (fusionnées) en bedGraph : chrom, start, end,
    valeur ; les bases consécutives de même valeur forment un intervalle"""
    win = Windows([(c, s, e) for c, s, e, _g, _k in rb], 0)
    rows = []
    if bw and bw.exists():
        import pyBigWig
        f = pyBigWig.open(str(bw))
        chroms = f.chroms()
        for c, (st, en) in win.d.items():
            name = c if c in chroms else ("chr" + c if "chr" + c in chroms else
                                          ("chrM" if c == "MT" else None))
            if name is None:
                continue
            for s, e in zip(st, en):
                e = min(e, chroms[name])
                for a_, b_, v in f.intervals(name, s, e) or []:
                    rows.append([name, max(a_, s), min(b_, e), f"{v:.3f}"])
        f.close()
    elif bedgraph and bedgraph.exists():
        rows = subset(bedgraph, win, lambda f, s, e: [f[0], s, e, f[3]])
    if not rows:
        return 0
    rows.sort(key=lambda r: (chrom_key(r[0]), int(r[1])))
    merged = []
    for r in rows:
        if merged and merged[-1][0] == r[0] and merged[-1][2] == r[1] \
                and merged[-1][3] == r[3]:
            merged[-1][2] = r[2]
        else:
            merged.append(list(r))
    with open(path, "w") as fh:
        for r in merged:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return len(merged)


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
    p.add_argument("--gtf", type=Path)
    p.add_argument("--gnomad", type=Path)
    p.add_argument("--gismo", type=Path)
    p.add_argument("--phylop", type=Path)
    p.add_argument("--phylop-bedgraph", type=Path)
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
    ens = gtf_ids(a.gtf) if a.gtf and a.gtf.exists() else {}
    gn = gnomad_table(a.gnomad) if a.gnomad and a.gnomad.exists() else {}
    gm = gismo_table(a.gismo, gtf_mane(a.gtf) if a.gtf and a.gtf.exists() else frozenset()) \
        if a.gismo and a.gismo.exists() else {}
    cols = ["gene", "gene_id", "l1_class", "l1", "ci_low", "ci_high", "status",
            "complex", "protein_coupling", "rna_coupling", "fibroblast_r"] + \
        GNOMAD_OUT + GISMO_OUT
    n_gn = n_gm = 0
    n_cls = defaultdict(int)
    with open(out / "genes.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for g in sorted(set(gs) | in_regions):
            r = l1.get(g, {})
            cls = r.get("l1_class") or l1_class(r, a.class_high, a.class_low)
            n_cls[cls] += 1
            gid = (r.get("gene_id") or ens.get(g, "")).split(".")[0]
            x = gn.get(gid) or gn.get(g) or {}
            y = {c: (gm.get(gid, {}).get(c) or gm.get(g, {}).get(c, "")) for c in GISMO_OUT} \
                if gid in gm or g in gm else {}   # GISMO par ENSG, GISMO-mis par symbole
            n_gn += bool(x)
            n_gm += bool(y)
            row = [g, gid, cls, r.get("score", ""),
                   r.get("ci_low", ""), r.get("ci_high", ""),
                   r.get("status", "absent_du_proteome"), r.get("complex", ""),
                   r.get("protein_coupling", ""), r.get("rna_coupling", ""),
                   fib.get(g, {}).get("fibroblast_r", "")] + \
                [x.get(c, "") for c in GNOMAD_OUT] + [y.get(c, "") for c in GISMO_OUT]
            fh.write("\t".join("" if x in ("NA", "nan") else str(x) for x in row) + "\n")
    print("classes L1 : " + ", ".join(f"{k} {v}" for k, v in sorted(n_cls.items())))
    print(f"gnomAD : {n_gn} gènes ; GISMO : {n_gm} gènes")
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

    # phyloP 241 mammifères, base par base sur les régions (sans marge)
    pb = out / "annotations" / "phylop.bedgraph"
    n_ph = phylop_cut(a.phylop, a.phylop_bedgraph, rb, pb)
    if n_ph:
        import pysam
        pysam.tabix_compress(str(pb), str(pb) + ".gz", force=True)
        pysam.tabix_index(str(pb) + ".gz", preset="bed", force=True)
        pb.unlink()
        print(f"phyloP : {n_ph} intervalles (bgzip + tabix)")
        manifest.append(("annotations/phylop.bedgraph.gz", n_ph,
                         str(a.phylop or a.phylop_bedgraph)))
    else:
        print("phyloP absent : pas de colonne MITO_phyloP")

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
