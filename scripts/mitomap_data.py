#!/usr/bin/env python3
"""
mitomap_data.py - build the data file of the MitoMap viewer (app/mitomap.html)

Purpose
    Gather, for the genes of the set (MitoCarta3.0 and the PanelApp
    "Mitochondrial disorders" panel), where each protein sits in the cell,
    which MitoPathways it belongs to, and its functional links, into one
    JSON file read by the viewer.

Inputs
    --mitocarta      Human.MitoCarta3.0.xls, or its sheet A as TSV (Symbol,
                     Synonyms, MitoCarta3.0_SubMitoLocalization,
                     MitoCarta3.0_MitoPathways); --mitopathways sheet C as
                     TSV when sheet A gives leaf names only.
    --genes          optional extra genes (symbols, or a table with a gene
                     column); PanelApp genes outside MitoCarta go here.
    --pathways       optional pathways.tsv of regulome.py pathways: the
                     MitoPathways it inferred from GO for genes outside
                     MitoCarta are used and flagged as inferred.
    --hpa            optional Human Protein Atlas subcellular_location.tsv
                     (Gene name, Reliability, Main location, Additional
                     location): locations outside mitochondria, and the
                     location of genes outside MitoCarta.
    --gaf            optional GO annotations (goa_human.gaf): cellular
                     component terms for genes HPA does not place.
    --extra-locations optional TSV gene, compartment[, source]: curated
                     locations added to the others.
    --string-links   optional STRING protein.links(.detailed) file (space-
                     separated, gzipped or not); --string-info the matching
                     protein.info file (string_protein_id, preferred_name).

Output
    --out            JSON: meta; compartments (id, label); genes (s symbol,
                     e Ensembl id when known, syn synonyms, loc compartments
                     with the primary one first, ls the source of each
                     location, sec MitoPathways sections, pw pathways, src
                     mitocarta / panelapp, inf 1 when pathways are inferred);
                     links ([i, j, combined, experimental, database,
                     coexpression, textmining], scores 0-1000, the channels
                     absent from a plain links file are 0).

Method
    Compartments: matrix, nucleoid, imm, ims, omm, mito_membrane (MitoCarta
    "Membrane"), mito_unknown, cytosol, nucleus, er, peroxisome,
    plasma_membrane, other, unlocalized. A MitoCarta gene takes its
    sub-mitochondrial localisation; a matrix gene of a pathway matching
    --nucleoid-regex goes to the nucleoid. Locations outside mitochondria
    from HPA (reliability in --hpa-reliability) are added as further
    locations, and the viewer joins them by an arrow. A gene outside
    MitoCarta is placed by HPA, then GO cellular component, else
    unlocalized. At most --max-locations locations per gene are kept.

Usage
    python3 mitomap_data.py --mitocarta data/Human.MitoCarta3.0.xls \\
        --genes genes.txt --pathways results/regulome/pathways.tsv \\
        --hpa data/subcellular_location.tsv.zip \\
        --gaf data/goa_human.gaf.gz \\
        --string-links data/9606.protein.links.detailed.v12.0.txt.gz \\
        --string-info data/9606.protein.info.v12.0.txt.gz \\
        --out results/mitomap.json

References
    Rath S et al. 2021. MitoCarta3.0. Nucleic Acids Res 49:D1541-D1547.
    Thul PJ et al. 2017. A subcellular map of the human proteome. Science
        356:eaal3321.
    Szklarczyk D et al. 2023. The STRING database in 2023. Nucleic Acids
        Res 51:D638-D646.
    Martin AR et al. 2019. PanelApp. Nat Genet 51:1560-1565.
"""

from __future__ import annotations

import argparse
import datetime
import io
import json
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from re2g_panel import ALIASES, opener  # noqa: E402
from regulome import read_gene_list, read_mitocarta, read_table  # noqa: E402

COMPARTMENTS = [
    ("matrix", "Matrice"), ("nucleoid", "Nucléoïde (ADNmt)"),
    ("imm", "Membrane interne"), ("ims", "Espace intermembranaire"),
    ("omm", "Membrane externe"),
    ("mito_membrane", "Membrane mitochondriale, non précisée"),
    ("mito_unknown", "Mitochondrie, sous-compartiment inconnu"),
    ("cytosol", "Cytosol"), ("nucleus", "Noyau"),
    ("er", "Réticulum endoplasmique"), ("peroxisome", "Peroxysome"),
    ("plasma_membrane", "Membrane plasmique"), ("other", "Autre"),
    ("unlocalized", "Non localisé"),
]
MITO = {"matrix", "nucleoid", "imm", "ims", "omm", "mito_membrane",
        "mito_unknown"}
MC_LOC = {"matrix": "matrix", "mim": "imm", "ims": "ims", "mom": "omm",
          "membrane": "mito_membrane", "unknown": "mito_unknown",
          "inner membrane": "imm", "outer membrane": "omm",
          "intermembrane space": "ims"}
HPA_LOC = [
    (r"^mitochondria", "mito_unknown"),
    (r"^nucle|^nuclear|kinetochore|mitotic chromosome", "nucleus"),
    (r"^cytosol|^cytoplasmic bodies|^rods & rings|aggresome", "cytosol"),
    (r"^endoplasmic reticulum", "er"),
    (r"^peroxisomes", "peroxisome"),
    (r"^plasma membrane|^cell junctions|^focal adhesion", "plasma_membrane"),
]
GO_LOC = {
    "GO:0005759": "matrix", "GO:0042645": "nucleoid",
    "GO:0005743": "imm", "GO:0005758": "ims", "GO:0005741": "omm",
    "GO:0031966": "mito_membrane", "GO:0005739": "mito_unknown",
    "GO:0005829": "cytosol", "GO:0005634": "nucleus",
    "GO:0005654": "nucleus", "GO:0005783": "er", "GO:0005789": "er",
    "GO:0005777": "peroxisome", "GO:0005886": "plasma_membrane",
}


def hpa_compartment(name: str):
    n = name.strip().lower()
    for rx, c in HPA_LOC:
        if re.search(rx, n):
            return c
    return "other" if n else None


def read_hpa(path: Path, reliability: set) -> dict:
    """gene -> [(compartment, 'main'|'additional')]"""
    if str(path).endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            name = next(n for n in z.namelist() if n.endswith(".tsv"))
            text = io.TextIOWrapper(z.open(name), encoding="utf-8")
            lines = text.read().splitlines()
    else:
        with opener(path) as fh:
            lines = fh.read().splitlines()
    head = lines[0].split("\t")
    col = {h.strip().lower(): i for i, h in enumerate(head)}
    gi = col.get("gene name")
    ri = col.get("reliability")
    mi = col.get("main location")
    ai = col.get("additional location")
    if gi is None or mi is None:
        sys.exit(f"{path}: expected the columns Gene name, Main location")
    out = defaultdict(list)
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) <= max(gi, mi):
            continue
        if ri is not None and reliability and f[ri] not in reliability:
            continue
        for idx, kind in ((mi, "main"), (ai, "additional")):
            if idx is None or idx >= len(f):
                continue
            for loc in f[idx].split(";"):
                c = hpa_compartment(loc)
                if c and (c, kind) not in out[f[gi]]:
                    out[f[gi]].append((c, kind))
    return out


def read_gaf_locations(path: Path) -> dict:
    out = defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if line[0] == "!":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[8] != "C" or "NOT" in f[3] or f[6] == "ND":
                continue
            c = GO_LOC.get(f[4])
            if c and c not in out[f[2]]:
                out[f[2]].append(c)
    return out


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mitocarta", type=Path)
    p.add_argument("--mitopathways", type=Path)
    p.add_argument("--genes", type=Path)
    p.add_argument("--pathways", type=Path)
    p.add_argument("--hpa", type=Path)
    p.add_argument("--hpa-reliability", default="Enhanced,Supported,Approved",
                   help="HPA reliability levels kept (default %(default)s)")
    p.add_argument("--gaf", type=Path)
    p.add_argument("--extra-locations", type=Path)
    p.add_argument("--nucleoid-regex",
                   default=r"mtDNA nucleoid|mtDNA replication|"
                           r"mtDNA maintenance",
                   help="pathways whose matrix genes go to the nucleoid")
    p.add_argument("--string-links", type=Path)
    p.add_argument("--string-info", type=Path)
    p.add_argument("--min-string", type=int, default=400,
                   help="lowest STRING combined score kept, 0-1000 "
                        "(default %(default)s)")
    p.add_argument("--max-locations", type=int, default=3)
    p.add_argument("--out", required=True, type=Path)
    a = p.parse_args()

    genes = {}                       # symbol -> record
    order = []
    nucleoid = re.compile(a.nucleoid_regex, re.I)
    if a.mitocarta:
        mc, _h = read_mitocarta(a.mitocarta, a.mitopathways)
        for g, v in mc.items():
            locs = []
            for x in re.split(r"[|,;]", v.get("loc", "")):
                c = MC_LOC.get(x.strip().lower())
                if c and c not in locs:
                    locs.append(c)
            paths = [" > ".join(pth) for pth in v["paths"]]
            if "matrix" in locs and any(nucleoid.search(pth)
                                        for pth in paths):
                locs = ["nucleoid" if c == "matrix" else c for c in locs]
            genes[g] = {"s": g, "syn": v["synonyms"][:6],
                        "loc": locs or ["mito_unknown"],
                        "ls": ["MitoCarta3.0"] * max(1, len(locs)),
                        "sec": list(dict.fromkeys(pth[0] for pth in
                                                  v["paths"])),
                        "pw": paths, "src": "mitocarta", "inf": 0}
            order.append(g)
        print(f"MitoCarta: {len(order)} genes")

    syn_of = {}
    for g, r in genes.items():
        for s in r["syn"]:
            syn_of.setdefault(s, g)
    for old, new in ALIASES.items():
        for c in new:
            if c in genes:
                syn_of.setdefault(old, c)
    n_extra = 0
    if a.genes:
        for g in read_gene_list(a.genes):
            if g in genes or g in syn_of:
                continue
            if re.search(r"_(CTG|GAA|CGG|GCN)|-Loss|-Gain|ISCA-", g):
                continue
            genes[g] = {"s": g, "syn": [], "loc": [], "ls": [], "sec": [],
                        "pw": [], "src": "panelapp", "inf": 0}
            order.append(g)
            n_extra += 1
        print(f"{n_extra} genes outside MitoCarta added")

    if a.pathways:
        for r in read_table(a.pathways):
            g = r.get("gene", "")
            if g in genes and not genes[g]["sec"] and r.get("sections"):
                genes[g]["sec"] = r["sections"].split("|")
                genes[g]["pw"] = [x for x in r.get("pathways", "").split("|")
                                  if x]
                genes[g]["inf"] = int(r.get("source") == "go_inferred")

    hpa = read_hpa(a.hpa, set(a.hpa_reliability.split(","))) if a.hpa \
        else {}
    gaf = read_gaf_locations(a.gaf) if a.gaf else {}
    for g in order:
        r = genes[g]
        h = hpa.get(g, [])
        if r["src"] == "mitocarta":
            for c, _kind in h:
                if c not in MITO and c != "other" and c not in r["loc"]:
                    r["loc"].append(c)
                    r["ls"].append("HPA")
        else:
            for c, kind in sorted(h, key=lambda x: x[1] != "main"):
                if c not in r["loc"]:
                    r["loc"].append(c)
                    r["ls"].append("HPA")
            if not r["loc"]:
                for c in gaf.get(g, []):
                    if c not in r["loc"]:
                        r["loc"].append(c)
                        r["ls"].append("GO")
            # a specific mitochondrial compartment replaces 'mito_unknown'
            spec = [c for c in r["loc"] if c in MITO and c != "mito_unknown"]
            if spec and "mito_unknown" in r["loc"]:
                k = r["loc"].index("mito_unknown")
                del r["loc"][k], r["ls"][k]
    if a.extra_locations:
        for row in read_table(a.extra_locations):
            v = list(row.values())
            g, c = v[0], v[1] if len(v) > 1 else ""
            g = g if g in genes else syn_of.get(g, g)
            if g in genes and c in dict(COMPARTMENTS) \
                    and c not in genes[g]["loc"]:
                genes[g]["loc"].append(c)
                genes[g]["ls"].append(v[2] if len(v) > 2 and v[2]
                                      else "curated")
    for g in order:
        r = genes[g]
        if not r["loc"]:
            r["loc"], r["ls"] = ["unlocalized"], [""]
        r["loc"] = r["loc"][:a.max_locations]
        r["ls"] = r["ls"][:a.max_locations]

    links = []
    if a.string_links and a.string_info:
        name = {}
        with opener(a.string_info) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    name[f[0]] = f[1]
        idx = {g: i for i, g in enumerate(order)}
        for g in order:
            for s in genes[g]["syn"]:
                idx.setdefault(s, idx[g])
        with opener(a.string_links) as fh:
            head = fh.readline().split()
            col = {h: i for i, h in enumerate(head)}
            cs = col.get("combined_score", len(head) - 1)
            chan = [col.get(c) for c in ("experimental", "database",
                                         "coexpression", "textmining")]
            seen = set()
            for line in fh:
                f = line.split()
                try:
                    sc = int(f[cs])
                except (IndexError, ValueError):
                    continue
                if sc < a.min_string:
                    continue
                i = idx.get(name.get(f[0], ""))
                j = idx.get(name.get(f[1], ""))
                if i is None or j is None or i == j:
                    continue
                key = (min(i, j), max(i, j))
                if key in seen:
                    continue
                seen.add(key)
                links.append([key[0], key[1], sc] +
                             [int(f[c]) if c is not None else 0
                              for c in chan])
        print(f"STRING: {len(links)} links with a combined score >= "
              f"{a.min_string}")

    tally = defaultdict(int)
    for g in order:
        for c in genes[g]["loc"]:
            tally[c] += 1
    print("genes per compartment: " + ", ".join(
        f"{c} {tally[c]}" for c, _ in COMPARTMENTS if tally[c]))
    print(f"{sum(len(genes[g]['loc']) > 1 for g in order)} genes with "
          "several locations")
    data = {"meta": {"built": datetime.date.today().isoformat(),
                     "genes": len(order), "links": len(links),
                     "string_min": a.min_string if links else None,
                     "sources": [x for x, y in (
                         ("MitoCarta3.0", a.mitocarta), ("HPA", a.hpa),
                         ("GO", a.gaf), ("STRING", a.string_links),
                         ("regulome.py pathways", a.pathways)) if y]},
            "compartments": [{"id": c, "label": lab}
                             for c, lab in COMPARTMENTS],
            "genes": [genes[g] for g in order], "links": links}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(data, fh, ensure_ascii=False, separators=(",", ":"))
    print(f"Written to {a.out}")


if __name__ == "__main__":
    main()
