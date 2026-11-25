#!/usr/bin/env python3
"""
fetch_hic_encode.py - library table (tissues.tsv format) of the ENCODE Hi-C .hic files of chosen
biosamples, for v4c_pipeline.py extract and loops_par_tissu.py (new tissue groups such as blood
or fibroblasts). Nothing is downloaded: the .hic files are read at a distance by v4c_pipeline.py.

Search: released File objects, file_format hic, GRCh38, whose assay title matches --assay (intact
Hi-C and in situ Hi-C by default) and whose biosample term matches --biosamples and not
--exclure. --liste prints the candidates grouped by biosample, assay and output type, and writes
nothing; without it the table is written with one row per file of --output-type (one output type
per library, so that a library is not counted twice).

Output columns: accession, term_name, classification, assay, output_type, lab, file_size, dataset,
s3_url, download_url, gtex_tissue, mtcn, kept. gtex_tissue and mtcn are carried along by
v4c_pipeline.py (select drops libraries without a finite mtcn): --gtex-tissu and --mtcn fill them
for the whole table (mtcn is not used to call loops).

Usage
  python3 scripts/fetch_hic_encode.py --biosamples "GM12878|B cell|T cell" --liste
  python3 scripts/fetch_hic_encode.py --biosamples "fibroblast|IMR-90" --output-type "mapping quality thresholded chromatin interactions"
      --gtex-tissu "Cells - Cultured fibroblasts" --mtcn 1 --out data/hic/tissues_fibroblaste.tsv
Python 3 standard library ; needs access to www.encodeproject.org (login node).
"""

import argparse
import json
import re
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

PORTAL = "https://www.encodeproject.org"
FIELDS = ["accession", "file_format", "output_type", "assembly", "status", "href", "file_size",
          "cloud_metadata.url", "biosample_ontology.term_name", "biosample_ontology.classification",
          "assay_title", "dataset", "lab.title", "lab"]


def get(path, **params):
    q = urllib.parse.urlencode(params, doseq=True)
    req = urllib.request.Request(f"{PORTAL}{path}?{q}", headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.load(r)


def labo(f):
    """lab of a File object: title when embedded, else the @id (/labs/erez-aiden/)"""
    lb = f.get("lab") or ""
    return (lb.get("title") or lb.get("@id") or "") if isinstance(lb, dict) else str(lb)


def filtrer(graph, biosamples, exclure, assay, output_type="", lab=""):
    """rows of the File objects kept (dict per file)"""
    rx, rx_ex, rx_a = re.compile(biosamples, re.I), re.compile(exclure, re.I) if exclure else None, \
        re.compile(assay, re.I)
    out = []
    for f in graph:
        bo = f.get("biosample_ontology") or {}
        bo = bo[0] if isinstance(bo, list) and bo else bo
        term = bo.get("term_name", "") if isinstance(bo, dict) else ""
        if (f.get("file_format") != "hic" or f.get("status") != "released"
                or f.get("assembly", "GRCh38") != "GRCh38" or not rx.search(term)
                or (rx_ex and rx_ex.search(term)) or not rx_a.search(f.get("assay_title", ""))
                or (output_type and f.get("output_type") != output_type)
                or (lab and not re.search(lab, labo(f), re.I))):
            continue
        out.append({"accession": f.get("accession", ""), "term_name": term,
                    "classification": bo.get("classification", "") if isinstance(bo, dict) else "",
                    "assay": f.get("assay_title", ""), "output_type": f.get("output_type", ""),
                    "file_size": str(f.get("file_size", "")), "dataset": f.get("dataset", ""),
                    "lab": labo(f),
                    "s3_url": (f.get("cloud_metadata") or {}).get("url", ""),
                    "download_url": PORTAL + f["href"] if f.get("href") else ""})
    out.sort(key=lambda r: (r["term_name"], r["assay"], r["accession"]))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--biosamples", required=True, help="regex on the biosample term name")
    p.add_argument("--exclure", default="", help="regex of biosample terms left out")
    p.add_argument("--assay", default=r"intact Hi-C|in situ Hi-C", help="regex on the assay title")
    p.add_argument("--output-type", default="", help="output type kept in the table (see --liste)")
    p.add_argument("--lab", default="", help="regex on the lab of the experiment that produced the data "
                   "(e.g. aiden: Erez Aiden, Baylor), not the lab that processed the file")
    p.add_argument("--gtex-tissu", default="", help="gtex_tissue column")
    p.add_argument("--mtcn", default="1", help="mtcn column (not used to call loops)")
    p.add_argument("--liste", action="store_true", help="print the candidates, write nothing")
    p.add_argument("--out", default="", help="table to write (without --liste)")
    p.add_argument("--json", default="", help="ENCODE search result saved before (tests, offline)")
    a = p.parse_args()
    if a.json:
        graph = json.load(open(a.json)).get("@graph", [])
    else:
        print("recherche des fichiers .hic sur ENCODE (une à quelques minutes)...", file=sys.stderr, flush=True)
        graph = get("/search/", type="File", file_format="hic", assembly="GRCh38", status="released",
                    limit="all", format="json", field=FIELDS).get("@graph", [])
    if a.json:
        labs_exp = {e.get("@id", ""): labo(e) for e in json.load(open(a.json)).get("experiments", [])}
    else:
        print("laboratoires des expériences (producteurs des données)...", file=sys.stderr, flush=True)
        labs_exp = {e.get("@id", ""): labo(e) for e in get(
            "/search/", type="Experiment", assay_title=["intact Hi-C", "in situ Hi-C"], limit="all",
            format="json", field=["@id", "accession", "lab.title", "lab"]).get("@graph", [])}
    for f in graph:                     # le laboratoire qui a produit les données est celui de l'expérience ;
        f["lab_fichier"] = labo(f)      # celui du fichier est souvent « ENCODE Processing Pipeline »
        f["lab"] = labs_exp.get(f.get("dataset", ""), f["lab_fichier"])
    rows = filtrer(graph, a.biosamples, a.exclure, a.assay, "" if a.liste else a.output_type, a.lab)
    if a.liste or not a.output_type:
        c = Counter((r["term_name"], r["classification"], r["assay"], r["output_type"], r["lab"]) for r in rows)
        print("term_name\tclassification\tassay\toutput_type\tlab\tfichiers")
        for (t, cl, asy, ot, lb), n in sorted(c.items()):
            print(f"{t}\t{cl}\t{asy}\t{ot}\t{lb}\t{n}")
        print(f"{len(rows)} fichiers .hic ; {len({r['term_name'] for r in rows})} biosamples"
              + ("" if a.liste else " ; --output-type pour écrire la table"), file=sys.stderr)
        return
    if not a.out:
        sys.exit("--out requis sans --liste")
    cols = ["accession", "term_name", "classification", "assay", "output_type", "lab", "file_size", "dataset",
            "s3_url", "download_url", "gtex_tissue", "mtcn", "kept"]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            r.update(gtex_tissue=a.gtex_tissu, mtcn=a.mtcn, kept="1")
            fo.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"{len(rows)} librairies ({len({r['term_name'] for r in rows})} biosamples) → {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
