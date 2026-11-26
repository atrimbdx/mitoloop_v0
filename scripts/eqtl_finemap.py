#!/usr/bin/env python3
"""
eqtl_finemap.py - fine-mapped GTEx v8 eQTLs (SuSiE credible sets, eQTL Catalogue) of a few
genes, for the eQTL track of locus_figure.py (--eqtl).

The eQTL Catalogue (EBI) reprocessed GTEx v8 with SuSiE: for each tissue and gene, credible sets
of variants, each variant with its posterior inclusion probability (PIP) and its effect (beta).
The lead variant of an independent eQTL is not the causal one when it is in linkage
disequilibrium with others; the credible set and the PIP say how sure the causal variant is.

  liste    the GTEx gene-expression datasets of the catalogue (tabix/tabix_ftp_paths.tsv of
           eQTL-Catalogue-resources, which also gives the path of each credible-set file): dataset_id, sample_group, tissue_label, sample_size;
           to check the tissue of each campaign type before fetch
  fetch    one credible_sets.tsv.gz per type (--datasets type=QTD...;...) into --datadir
           (ftp_cs_path of the table, read over https from ftp.ebi.ac.uk)
  extraire the credible-set variants of --genes (symbols, matched to Ensembl ids by the GTF)
           from the fetched files: type, gene, dataset, cs_id, cs_size, variant, chrom, pos
           (1-based), pip, beta, se

Column names are read from the header (gene_id or molecular_trait_id; variant; pip; beta; se;
cs_id; cs_size); a file without them stops with its header printed.
Needs internet for liste and fetch (login node). Python 3 standard library.

Usage (from mitoloop_v0)
  python3 scripts/eqtl_finemap.py liste
  python3 scripts/eqtl_finemap.py fetch --datasets "muscle=QTD000xxx;cardiomyocyte=QTD000yyy" --datadir data/eqtl_finemap
  python3 scripts/eqtl_finemap.py extraire --datadir data/eqtl_finemap --genes BOLA3,COX6A2,COQ9
      --gtf data/gencode.gtf.gz --out results/eqtl_finemap/eqtl_finemap.tsv
"""

import argparse
import csv
import gzip
import io
import shutil
import sys
import urllib.request
from pathlib import Path

META = ("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/"
        "tabix/tabix_ftp_paths.tsv")       # study_id, dataset_id, ..., quant_method, ftp_cs_path
SUSIE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/{study}/{ds}/{ds}.credible_sets.tsv.gz"
csv.field_size_limit(10 ** 9)


def lire_meta(url=META):
    with urllib.request.urlopen(url, timeout=120) as r:
        return list(csv.DictReader(io.TextIOWrapper(r, encoding="utf-8"), delimiter="\t"))


def gtex_ge(meta):
    return [r for r in meta if r.get("study_label", "").upper().startswith("GTEX")
            and r.get("quant_method") == "ge"]


def cmd_liste(a):
    rows = gtex_ge(lire_meta(a.meta))
    print("dataset_id\tstudy_id\tsample_group\ttissue_label\tsample_size")
    for r in sorted(rows, key=lambda r: r.get("sample_group", "")):
        print(f"{r.get('dataset_id')}\t{r.get('study_id')}\t{r.get('sample_group')}\t"
              f"{r.get('tissue_label')}\t{r.get('sample_size')}")
    print(f"{len(rows)} jeux GTEx (expression)", file=sys.stderr)


def cmd_fetch(a):
    meta = {r["dataset_id"]: r for r in lire_meta(a.meta)}
    d = Path(a.datadir)
    d.mkdir(parents=True, exist_ok=True)
    man = []
    for x in a.datasets.split(";"):
        if "=" not in x:
            continue
        t, ds = (y.strip() for y in x.split("=", 1))
        if ds not in meta:
            sys.exit(f"{ds} absent du catalogue (voir liste)")
        url = (meta[ds].get("ftp_cs_path") or SUSIE.format(study=meta[ds]["study_id"], ds=ds)) \
            .replace("ftp://", "https://", 1)
        f = d / f"{ds}.credible_sets.tsv.gz"
        if f.is_file() and f.stat().st_size > 0:
            print(f"  {t} : {f.name} déjà là")
        else:
            print(f"  {t} : {url}", flush=True)
            tmp = f.with_suffix(".part")
            with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as fo:
                shutil.copyfileobj(r, fo, 1 << 20)
            tmp.rename(f)
        man.append((t, ds, meta[ds].get("sample_group", ""), f.name))
    with open(d / "manifeste.tsv", "w") as fo:
        fo.write("type\tdataset_id\tsample_group\tfichier\n")
        for m in man:
            fo.write("\t".join(m) + "\n")
    print(f"→ {d / 'manifeste.tsv'}")


def gene_ids(gtf, symbols):
    out = {}
    op = gzip.open if str(gtf).endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t", 9)
            if len(f) < 9 or f[2] != "gene":
                continue
            at = f[8]
            nm = at.split('gene_name "', 1)[1].split('"', 1)[0] if 'gene_name "' in at else ""
            if nm in symbols:
                out[at.split('gene_id "', 1)[1].split('"', 1)[0].split(".")[0]] = nm
    return out


def col(h, *names):
    return next((n for n in names if n in h), None)


def cmd_extraire(a):
    d = Path(a.datadir)
    man = list(csv.DictReader(open(d / "manifeste.tsv"), delimiter="\t"))
    symbols = {x.strip() for x in a.genes.split(",") if x.strip()}
    ids = gene_ids(a.gtf, symbols)
    manque = symbols - set(ids.values())
    if manque:
        print(f"  absents du GTF : {', '.join(sorted(manque))}", file=sys.stderr)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["type", "gene", "dataset", "cs_id", "cs_size", "variant", "chrom", "pos", "pip", "beta", "se"]
    n = 0
    with open(out, "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for m in man:
            with gzip.open(d / m["fichier"], "rt") as fh:
                rd = csv.DictReader(fh, delimiter="\t")
                h = rd.fieldnames or []
                cg, cv, cp = col(h, "gene_id", "molecular_trait_id"), col(h, "variant", "variant_id"), col(h, "pip")
                if not (cg and cv and cp):
                    sys.exit(f"{m['fichier']} : colonnes inattendues : {h}")
                cb, cs_, cc, cz = col(h, "beta"), col(h, "se"), col(h, "cs_id"), col(h, "cs_size")
                k = 0
                for r in rd:
                    g = ids.get((r[cg] or "").split(".")[0])
                    if not g:
                        continue
                    v = r[cv].replace("chr", "", 1).split("_") if r[cv] else []
                    if len(v) < 2:
                        continue
                    fo.write("\t".join([m["type"], g, m["dataset_id"], r.get(cc, "") or "", r.get(cz, "") or "",
                                        r[cv], "chr" + v[0], v[1], r[cp], r.get(cb, "") or "",
                                        r.get(cs_, "") or ""]) + "\n")
                    k += 1
                print(f"  {m['type']} ({m['dataset_id']}) : {k} variants de credible sets", file=sys.stderr)
                n += k
    print(f"{n} variants → {out}", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("liste")
    s.add_argument("--meta", default=META)
    s = sub.add_parser("fetch")
    s.add_argument("--meta", default=META)
    s.add_argument("--datasets", required=True, help="type=dataset_id;... (voir liste)")
    s.add_argument("--datadir", required=True)
    s = sub.add_parser("extraire")
    s.add_argument("--datadir", required=True)
    s.add_argument("--genes", required=True, help="symboles séparés par des virgules")
    s.add_argument("--gtf", required=True)
    s.add_argument("--out", required=True)
    a = p.parse_args()
    {"liste": cmd_liste, "fetch": cmd_fetch, "extraire": cmd_extraire}[a.cmd](a)


if __name__ == "__main__":
    main()
