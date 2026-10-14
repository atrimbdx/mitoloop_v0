#!/usr/bin/env python3
"""
mitoshow_data.py - build the data file of the MitoShow viewer
(app/mitoshow.html)

Purpose
    Gather, for the genes of the set (MitoCarta3.0 and the PanelApp
    "Mitochondrial disorders" panel), where each protein sits in the cell,
    which MitoPathways it belongs to, where each subunit sits in its complex
    (from structures), and its functional links, into one JSON file read by
    the viewer. The viewer derives the nesting itself: compartment, then the
    MitoPathways levels, the OXPHOS complexes on a crista, the nucleoid and
    the mtDNA genome, and one copy of a protein per incompatible place.

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
    --extra-locations optional TSV gene, compartment[, path[, source]]:
                     curated places added to the others; path is a
                     ' > '-separated list of nested regions inside the
                     compartment (for instance oxphos and "Supercomplexes >
                     Respirasome assembly").
    --structures     optional TSV complex, file[, projection[, flip_u[,
                     flip_v]]]: one mmCIF file per complex (gzipped or not),
                     complex named as in MitoPathways ("Complex I"). The
                     position of each subunit across the complex (u along
                     the first axis, v across it) is taken from its C-alpha
                     atoms: projection 'z' uses x and z of a membrane-
                     oriented file (OPM, membrane normal along z), 'pca'
                     (default) the two first principal axes. Chains in
                     contact (at least --contact-min C-alpha pairs within
                     --contact-distance) become links of the structure
                     channel. Chains are mapped to genes through
                     _struct_ref UniProt accessions and --id-map, or by
                     --chain-map.
    --id-map         TSV with columns symbol and uniprot_ids (HGNC complete
                     set), for --structures.
    --chain-map      optional TSV complex, chain, gene, overriding --id-map.
    --string-links   optional STRING protein.links(.detailed) file (space-
                     separated, gzipped or not); --string-info the matching
                     protein.info file (string_protein_id, preferred_name).

Output
    --out            JSON: meta; compartments (id, label); genes (s symbol,
                     e Ensembl id when known, syn synonyms, loc compartments
                     with the primary one first, ls the source of each
                     location, sec MitoPathways sections, pw pathways, src
                     mitocarta / panelapp, inf 1 when pathways are inferred,
                     extra curated places [{c, p, src}], xy {complex: [u, v]}
                     from structures); links ([i, j, combined, experimental,
                     database, coexpression, textmining, structure], scores
                     0-1000, channels absent from the inputs are 0).

Method
    Compartments: matrix, imm, ims, omm, mito_membrane (MitoCarta
    "Membrane"), mito_unknown, cytosol, nucleus, er, peroxisome,
    plasma_membrane, other, unlocalized; the viewer adds nucleoid, oxphos
    and mtdna. A MitoCarta gene takes its sub-mitochondrial localisation.
    The 37 genes of the mtDNA keep no location of their own when they have
    none, and the viewer draws them on the mtDNA genome. Locations outside mitochondria
    from HPA (reliability in --hpa-reliability) are added as further
    locations, and the viewer joins them by an arrow. A gene outside
    MitoCarta is placed by HPA, then GO cellular component, else
    unlocalized. At most --max-locations locations per gene are kept.

Usage
    python3 mitoshow_data.py --mitocarta data/Human.MitoCarta3.0.xls \\
        --genes genes.txt --pathways results/regulome/pathways.tsv \\
        --hpa data/subcellular_location.tsv.zip \\
        --gaf data/goa_human.gaf.gz \\
        --extra-locations config/mitoshow_extra.tsv \\
        --structures config/mitoshow_structures.tsv \\
        --id-map data/hgnc_complete_set.txt \\
        --string-links data/9606.protein.links.detailed.v12.0.txt.gz \\
        --string-info data/9606.protein.info.v12.0.txt.gz \\
        --out results/mitoshow.json

References
    Rath S et al. 2021. MitoCarta3.0. Nucleic Acids Res 49:D1541-D1547.
    Thul PJ et al. 2017. A subcellular map of the human proteome. Science
        356:eaal3321.
    Szklarczyk D et al. 2023. The STRING database in 2023. Nucleic Acids
        Res 51:D638-D646.
    Berman HM et al. 2000. The Protein Data Bank. Nucleic Acids Res
        28:235-242.
    Lomize MA et al. 2012. OPM database and PPM web server. Nucleic Acids
        Res 40:D370-D376.
    Andrews RM et al. 1999. Reanalysis and revision of the Cambridge
        reference sequence for human mitochondrial DNA. Nat Genet 23:147.
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

MTDNA_GENES = {
    "MT-TF", "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1", "MT-ND1", "MT-TI",
    "MT-TQ", "MT-TM", "MT-ND2", "MT-TW", "MT-TA", "MT-TN", "MT-TC", "MT-TY",
    "MT-CO1", "MT-TS1", "MT-TD", "MT-CO2", "MT-TK", "MT-ATP8", "MT-ATP6",
    "MT-CO3", "MT-TG", "MT-ND3", "MT-TR", "MT-ND4L", "MT-ND4", "MT-TH",
    "MT-TS2", "MT-TL2", "MT-ND5", "MT-ND6", "MT-TE", "MT-CYB", "MT-TT",
    "MT-TP"}
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


def cif_tokens(path: Path):
    """Tokens of a CIF file, multi-line ';' fields and quotes handled."""
    with opener(path) as fh:
        lines = fh.read().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(";"):
            buf = [line[1:]]
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                buf.append(lines[i])
                i += 1
            yield "\n".join(buf).strip()
            i += 1
            continue
        j, n = 0, len(line)
        while j < n:
            c = line[j]
            if c.isspace():
                j += 1
            elif c == "#":
                break
            elif c in "'\"":
                k = j + 1
                while k < n and not (line[k] == c and (k + 1 == n or
                                                       line[k + 1].isspace())):
                    k += 1
                yield line[j + 1:k]
                j = k + 1
            else:
                k = j
                while k < n and not line[k].isspace():
                    k += 1
                yield line[j:k]
                j = k
        i += 1


def read_cif(path: Path, wanted: set) -> dict:
    """category -> list of rows (dict), for the categories in wanted."""
    out = {}
    toks = list(cif_tokens(path))
    i = 0
    while i < len(toks):
        t = toks[i]
        if t == "loop_":
            tags = []
            i += 1
            while i < len(toks) and toks[i].startswith("_"):
                tags.append(toks[i])
                i += 1
            vals = []
            while i < len(toks) and not toks[i].startswith("_") and \
                    toks[i] != "loop_" and not toks[i].startswith("data_"):
                vals.append(toks[i])
                i += 1
            cat = tags[0].split(".")[0] if tags else ""
            if cat in wanted:
                names = [x.split(".", 1)[1] for x in tags]
                rows = [dict(zip(names, vals[k:k + len(names)]))
                        for k in range(0, len(vals) - len(names) + 1,
                                       len(names))]
                out.setdefault(cat, []).extend(rows)
            continue
        if t.startswith("_") and "." in t:
            cat, name = t.split(".", 1)
            if cat in wanted and i + 1 < len(toks):
                rows = out.setdefault(cat, [{}])
                rows[0][name] = toks[i + 1]
            i += 2
            continue
        i += 1
    return out


def read_structures(a, genes, syn_of) -> dict:
    """xy per gene and complex; returns chain contacts as gene pairs."""
    import numpy as np
    acc2sym = {}
    if a.id_map:
        for r in read_table(a.id_map):
            sym = r.get("symbol", "")
            for acc in r.get("uniprot_ids", "").replace("|", ",").split(","):
                if acc.strip() and sym:
                    acc2sym[acc.strip()] = sym
    chain_map = {}
    if a.chain_map:
        for r in read_table(a.chain_map):
            v = list(r.values())
            chain_map[(v[0], v[1])] = v[2]
    contacts = {}
    for row in read_table(a.structures):
        v = list(row.values()) + ["", "", "", ""]
        cx, path, proj = v[0], v[1], (v[2] or "pca").lower()
        flip_u, flip_v = v[3] in ("1", "yes", "true"), v[4] in ("1", "yes",
                                                               "true")
        cif = read_cif(Path(path), {"_atom_site", "_struct_ref",
                                    "_struct_ref_seq", "_entity"})
        ent2gene = {}
        for r in cif.get("_struct_ref", []):
            sym = acc2sym.get(r.get("pdbx_db_accession", ""))
            if sym:
                ent2gene[r.get("entity_id")] = sym
        chains = {}
        for r in cif.get("_atom_site", []):
            if r.get("label_atom_id") != "CA":
                continue
            ch = r.get("auth_asym_id") or r.get("label_asym_id")
            g = chain_map.get((cx, ch)) or ent2gene.get(r.get("label_entity_id"))
            if not g:
                continue
            g = g if g in genes else syn_of.get(g, g)
            if g not in genes:
                continue
            chains.setdefault(ch, [g, []])[1].append(
                [float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])])
        first = {}
        for ch, (g, xyz) in chains.items():
            if g not in first:
                first[g] = np.array(xyz)
        if not first:
            print(f"structures: {cx} ({path}): no chain mapped to a gene")
            continue
        allxyz = np.vstack(list(first.values()))
        centre = allxyz.mean(0)
        if proj == "z":
            axes = np.array([[1, 0, 0], [0, 0, 1]], float)
        else:
            _u, _s, vt = np.linalg.svd(allxyz - centre, full_matrices=False)
            axes = vt[:2]
        uv = {g: (c.mean(0) - centre) @ axes.T for g, c in first.items()}
        arr = np.array(list(uv.values()))
        lo, hi = arr.min(0), arr.max(0)
        span = np.where(hi - lo > 0, hi - lo, 1)
        for g, p_ in uv.items():
            u, w = (p_ - lo) / span * 0.9 + 0.05
            if flip_u:
                u = 1 - u
            if flip_v:
                w = 1 - w
            genes[g].setdefault("xy", {})[cx] = [round(float(u), 3),
                                                 round(float(w), 3)]
        names = list(first)
        d2 = a.contact_distance ** 2
        for x in range(len(names)):
            for y in range(x + 1, len(names)):
                A, B = first[names[x]], first[names[y]]
                n = int((((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)
                         < d2).sum())
                if n >= a.contact_min:
                    key = tuple(sorted((names[x], names[y])))
                    contacts[key] = max(contacts.get(key, 0),
                                        min(999, 700 + 3 * n))
        print(f"structures: {cx}: {len(first)} subunits placed")
    return contacts


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
    p.add_argument("--structures", type=Path)
    p.add_argument("--id-map", type=Path)
    p.add_argument("--chain-map", type=Path)
    p.add_argument("--contact-distance", type=float, default=8.0)
    p.add_argument("--contact-min", type=int, default=8)
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
    if a.mitocarta:
        mc, _h = read_mitocarta(a.mitocarta, a.mitopathways)
        for g, v in mc.items():
            locs = []
            for x in re.split(r"[|,;]", v.get("loc", "")):
                c = MC_LOC.get(x.strip().lower())
                if c and c not in locs:
                    locs.append(c)
            paths = [" > ".join(pth) for pth in v["paths"]]
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
    places = dict(COMPARTMENTS)
    places.update({"oxphos": "", "nucleoid": "", "mtdna": ""})
    if a.extra_locations:
        for row in read_table(a.extra_locations):
            v = list(row.values()) + ["", "", ""]
            g, c, path, src = v[0], v[1], v[2], v[3]
            g = g if g in genes else syn_of.get(g, g)
            if g not in genes or c not in places:
                continue
            p_ = [x.strip() for x in path.split(">") if x.strip()]
            genes[g].setdefault("extra", []).append(
                {"c": c, "p": p_, "src": src or "curated"})
    for g in order:
        r = genes[g]
        if not r["loc"] and g not in MTDNA_GENES and not r.get("extra"):
            r["loc"], r["ls"] = ["unlocalized"], [""]
        r["loc"] = r["loc"][:a.max_locations]
        r["ls"] = r["ls"][:a.max_locations]

    links = []
    contacts = {}
    if a.structures:
        contacts = read_structures(a, genes, syn_of)
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
                              for c in chan] + [0])
        print(f"STRING: {len(links)} links with a combined score >= "
              f"{a.min_string}")
    if contacts:
        idx = {g: i for i, g in enumerate(order)}
        pos = {(l[0], l[1]): k for k, l in enumerate(links)}
        n_new = 0
        for (x, y), sc in contacts.items():
            if x not in idx or y not in idx:
                continue
            key = (min(idx[x], idx[y]), max(idx[x], idx[y]))
            if key in pos:
                links[pos[key]][7] = sc
                links[pos[key]][2] = max(links[pos[key]][2], sc)
            else:
                pos[key] = len(links)
                links.append([key[0], key[1], sc, 0, 0, 0, 0, sc])
                n_new += 1
        print(f"structures: {len(contacts)} chain contacts, {n_new} new "
              "links")

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
