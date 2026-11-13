#!/usr/bin/env python3
"""
fetch_re2g_encode.py - adds ENCODE-rE2G link files of chosen biosamples
(kidney, heart, brain...) to the folder that re2g_panel.py reads.

The file type to fetch is copied from the files already in the folder: the
ENCFF accessions found there are looked up on the ENCODE portal (output
type, file format, and the annotation they belong to), and only files of the
same output type, from annotations of the same kind, are listed. Nothing is
downloaded with --liste.

Steps
  1. reference: metadata of up to 3 ENCFF files already in --dest
  2. search: released annotations of type "element gene regulatory
     interaction predictions", GRCh38, whose biosample term matches
     --biosamples and whose description mentions rE2G (or which share the
     annotation subtype of the reference files)
  3. files of these annotations with the reference output type and format,
     not already in --dest: listed (manifeste_<date>.tsv), then downloaded
     (*.part, renamed when complete) unless --liste

Usage
  python3 scripts/fetch_re2g_encode.py --dest data/re2g --liste
  python3 scripts/fetch_re2g_encode.py --dest data/re2g            # download
  (--biosamples REGEX replaces the default: kidney, renal, heart, cardiac,
   brain regions)
Python 3 standard library ; needs access to www.encodeproject.org.
"""

import argparse
import datetime
import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

PORTAL = "https://www.encodeproject.org"
DEFAULT = (r"kidney|renal|nephr|tubul|HK-2|RPTEC|heart|cardiac|ventric|atri(um|al)|"
           r"brain|cortex|gyrus|caudate|putamen|cerebell|hippocamp|neuron|astrocyte")


def get(path, **params):
    q = urllib.parse.urlencode(params, doseq=True)
    url = f"{PORTAL}{path}" + (f"?{q}" if q else "")
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--dest", required=True, help="folder of link files (re2g_panel.py --links-dir)")
    p.add_argument("--biosamples", default=DEFAULT, help="regex on the biosample term name")
    p.add_argument("--liste", action="store_true", help="list only, download nothing")
    p.add_argument("--max", type=int, default=0, help="at most N files (0: all)")
    a = p.parse_args()
    dest = Path(a.dest)
    have = sorted({m.group(0) for f in dest.glob("*") for m in [re.search(r"ENCFF[0-9A-Z]+", f.name)] if m})
    if not have:
        sys.exit(f"aucun fichier ENCFF dans {dest} : impossible de savoir quel type de fichier prendre")
    ref_types, ref_formats, ref_subtypes = set(), set(), set()
    for acc in have[:3]:
        f = get(f"/files/{acc}/", format="json", frame="object")
        ref_types.add(f.get("output_type", ""))
        ref_formats.add(f.get("file_format", ""))
        ds = f.get("dataset", "")
        if ds:
            d = get(ds, format="json", frame="object")
            ref_subtypes.add(d.get("annotation_subtype", ""))
            print(f"référence {acc} : {f.get('output_type')} / {f.get('file_format')} ; "
                  f"{d.get('accession')} {d.get('description', '')[:80]}", file=sys.stderr)
    res = get("/search/", type="Annotation",
              annotation_type="element gene regulatory interaction predictions",
              status="released", assembly="GRCh38", limit="all", format="json",
              field=["accession", "description", "annotation_subtype",
                     "biosample_ontology.term_name", "biosample_ontology.classification",
                     "files.accession", "files.output_type", "files.file_format",
                     "files.href", "files.status", "files.assembly"])
    rx = re.compile(a.biosamples, re.I)
    rows = []
    for an in res.get("@graph", []):
        bo = an.get("biosample_ontology") or {}
        bo = bo[0] if isinstance(bo, list) and bo else bo
        term = bo.get("term_name", "") if isinstance(bo, dict) else ""
        if not rx.search(term):
            continue
        desc = an.get("description", "") or ""
        if not (re.search(r"re2g", desc, re.I) or an.get("annotation_subtype") in ref_subtypes - {""}):
            continue
        for f in an.get("files", []):
            if (f.get("output_type") in ref_types and f.get("file_format") in ref_formats
                    and f.get("status") == "released" and f.get("assembly", "GRCh38") == "GRCh38"
                    and f.get("accession") not in have):
                rows.append((f["accession"], an.get("accession", ""), term,
                             bo.get("classification", "") if isinstance(bo, dict) else "",
                             f.get("output_type", ""), f.get("href", "")))
    rows.sort(key=lambda r: (r[2], r[0]))
    if a.max:
        rows = rows[:a.max]
    man = dest / f"manifeste_{datetime.date.today():%Y%m%d}.tsv"
    with open(man, "w") as fo:
        fo.write("fichier\tannotation\tbiosample\tclassification\toutput_type\thref\n")
        for r in rows:
            fo.write("\t".join(r) + "\n")
    terms = {}
    for r in rows:
        terms[r[2]] = terms.get(r[2], 0) + 1
    print(f"{len(rows)} fichiers à ajouter ({len(have)} déjà présents) :", file=sys.stderr)
    for t, n in sorted(terms.items()):
        print(f"   {n:3d}  {t}", file=sys.stderr)
    print(f"liste : {man}", file=sys.stderr)
    if a.liste:
        return
    for i, (acc, _an, term, _c, _o, href) in enumerate(rows, 1):
        tgt = dest / Path(href).name
        part = tgt.with_name(tgt.name + ".part")
        print(f"  [{i}/{len(rows)}] {acc} ({term})", file=sys.stderr)
        with urllib.request.urlopen(PORTAL + href, timeout=600) as r, open(part, "wb") as fo:
            while True:
                b = r.read(1 << 20)
                if not b:
                    break
                fo.write(b)
        part.rename(tgt)
    print("terminé : relancer re2g_panel.py (étape paquet de reannot) pour les prendre en compte",
          file=sys.stderr)


if __name__ == "__main__":
    main()
