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
     not already in --dest, biosamples outside the simulated cell types left
     out (--exclure), at most --max-par-biosample per biosample term, and
     only those whose model columns (Feature) are those of the files already
     present (header read from the first block, no full download):
     listed (manifeste_<date>.tsv, with sizes), then downloaded (*.part,
     renamed when complete) unless --liste

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
import zlib
from pathlib import Path

PORTAL = "https://www.encodeproject.org"
# biosamples hors des types simulés (muscle, cardiomyocyte, neurone, tubule rénal)
EXCLURE = (r"astrocyte|endotheli|pericyte|fibroblast|interstitium|pelvis|glomerul|progenitor|"
           r"stem cell|smooth|capillary")
DEFAULT = (r"kidney|\brenal|nephr|tubul|HK-2|RPTEC|heart|cardiac|ventric|atri(um|al)|"
           r"brain|cortex|gyrus|caudate|putamen|cerebell|hippocamp|neuron|astrocyte")


def get(path, **params):
    q = urllib.parse.urlencode(params, doseq=True)
    url = f"{PORTAL}{path}" + (f"?{q}" if q else "")
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def features_local(path):
    """colonnes de modèle (Feature) d'un fichier rE2G local"""
    import gzip
    with gzip.open(path, "rt") as fh:
        return sorted(c for c in fh.readline().rstrip("\n").split("\t") if "Feature" in c)


def features_distant(href, octets=1 << 16):
    """mêmes colonnes, lues dans le premier bloc du fichier sur ENCODE (pas de téléchargement
    complet)"""
    req = urllib.request.Request(PORTAL + href, headers={"Range": f"bytes=0-{octets - 1}"})
    with urllib.request.urlopen(req, timeout=120) as r:
        raw = r.read(octets)
    txt = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(raw).decode("utf-8", "replace")
    return sorted(c for c in txt.split("\n", 1)[0].split("\t") if "Feature" in c)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--dest", required=True, help="folder of link files (re2g_panel.py --links-dir)")
    p.add_argument("--biosamples", default=DEFAULT, help="regex on the biosample term name")
    p.add_argument("--liste", action="store_true", help="list only, download nothing")
    p.add_argument("--exclure", default=EXCLURE, help="regex of biosamples left out")
    p.add_argument("--max-par-biosample", type=int, default=6,
                   help="at most N files per biosample term (0: all ; enough for the >= 3 rule)")
    p.add_argument("--sans-controle-modele", action="store_true",
                   help="do not check that the model columns match the files already present")
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
                     "files.href", "files.status", "files.assembly", "files.file_size"])
    rx = re.compile(a.biosamples, re.I)
    rx_ex = re.compile(a.exclure, re.I) if a.exclure else None
    rows = []
    for an in res.get("@graph", []):
        bo = an.get("biosample_ontology") or {}
        bo = bo[0] if isinstance(bo, list) and bo else bo
        term = bo.get("term_name", "") if isinstance(bo, dict) else ""
        if not rx.search(term) or (rx_ex and rx_ex.search(term)):
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
                             f.get("output_type", ""), f.get("href", ""),
                             str(f.get("file_size", 0))))
    rows.sort(key=lambda r: (r[2], r[0]))
    # même modèle rE2G que les fichiers présents (re2g_panel.py applique le seuil d'un seul modèle)
    if not a.sans_controle_modele:
        ref = features_local(next(f for f in sorted(dest.glob("*.bed.gz"))))
        keep, n_out, n_err = [], 0, 0
        per_term = {}
        for r in rows:
            if a.max_par_biosample and per_term.get(r[2], 0) >= a.max_par_biosample:
                continue
            try:
                ok = features_distant(r[5]) == ref
            except Exception:                    # en-tête illisible : écarté, compté
                ok, n_err = False, n_err + 1
            if ok:
                keep.append(r)
                per_term[r[2]] = per_term.get(r[2], 0) + 1
            else:
                n_out += 1
        print(f"modèle : {len(ref)} colonnes Feature de référence ; {n_out} fichiers d'un autre modèle "
              f"écartés ({n_err} en-têtes illisibles)", file=sys.stderr)
        rows = keep
    elif a.max_par_biosample:
        per_term, keep = {}, []
        for r in rows:
            if per_term.get(r[2], 0) < a.max_par_biosample:
                keep.append(r)
                per_term[r[2]] = per_term.get(r[2], 0) + 1
        rows = keep
    if a.max:
        rows = rows[:a.max]
    man = dest / f"manifeste_{datetime.date.today():%Y%m%d}.tsv"
    with open(man, "w") as fo:
        fo.write("fichier\tannotation\tbiosample\tclassification\toutput_type\thref\ttaille\n")
        for r in rows:
            fo.write("\t".join(r) + "\n")
    terms = {}
    for r in rows:
        terms[r[2]] = terms.get(r[2], 0) + 1
    tot = sum(int(r[6] or 0) for r in rows)
    print(f"{len(rows)} fichiers à ajouter ({len(have)} déjà présents), {tot / 1e9:.1f} Go :",
          file=sys.stderr)
    for t, n in sorted(terms.items()):
        print(f"   {n:3d}  {t}", file=sys.stderr)
    print(f"liste : {man}", file=sys.stderr)
    if a.liste:
        return
    for i, (acc, _an, term, _c, _o, href, _t) in enumerate(rows, 1):
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
