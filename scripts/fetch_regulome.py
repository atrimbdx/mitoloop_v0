#!/usr/bin/env python3
"""
fetch_regulome.py - public inputs of the mitochondrial regulome

Purpose
    Download, into one directory, the public files that regulome.py and
    promoter_features.py read (the GTEx and GENCODE files are downloaded by
    `gtex_eqtl.py fetch`), and write the gene set: MitoCarta3.0 and the
    PanelApp "Mitochondrial disorders" panel, every rating.

Inputs
    network access (run it on a login node); --reuse DIR: directories
    searched first for files of the same name (for instance the data/ of a
    MitoShow run), which are linked instead of downloaded.

Outputs (in --datadir)
    Human.MitoCarta3.0.xls, goa_human.gaf.gz, go-basic.obo,
    GRCh38.fa (+ .fai built by promoter_features.py), cpgIslandExt.txt.gz,
    panelapp_genes_all.tsv (gene, confidence_level, entity_type),
    and with --validation: kopajtich/ (fibroblast RNA-seq and proteomics,
    Zenodo 4501904), attenuation.csv (Goncalves et al. 2017);
    MANIFEST.tsv (file, bytes, source or error);
    --genes (default genes.txt): MitoCarta3.0 symbols and PanelApp genes.

Usage
    python3 scripts/fetch_regulome.py --datadir data --reuse ../mitoshow/data
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
UA = {"User-Agent": "mitoloop_v0/0.1 (fetch_regulome.py)"}
SOURCES = {
    "Human.MitoCarta3.0.xls": [
        "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls",
        "https://www.broadinstitute.org/files/shared/metabolism/mitocarta/human.mitocarta3.0.xls",
        "https://data.broadinstitute.org/mitocarta/Human.MitoCarta3.0.xls"],
    "goa_human.gaf.gz": [
        "https://current.geneontology.org/annotations/goa_human.gaf.gz",
        "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz"],
    "go-basic.obo": [
        "https://current.geneontology.org/ontology/go-basic.obo",
        "https://purl.obolibrary.org/obo/go/go-basic.obo"],
    "cpgIslandExt.txt.gz": [
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz"],
    "hg38.fa.gz": [
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"],
}
VALIDATION = {
    "attenuation.csv": [
        "https://raw.githubusercontent.com/saezlab/protein_attenuation/master/tables/protein_attenuation_table.csv"],
}
# patient fibroblasts, RNA-seq and proteomics (Kopajtich et al.): the files
# of the Zenodo record whose name holds one of these words; their exact names
# are read from the record
ZENODO_RECORD = "4501904"
ZENODO_WANTED = ("count", "proteomics", "annotation", "sample")
PANELAPP = "https://panelapp.genomicsengland.co.uk/api/v1"


def get(url: str, dest: Path, tries=3, timeout=600):
    """download url to dest; None, or the last error"""
    tmp = dest.with_name(dest.name + ".part")
    err = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as fh:
                shutil.copyfileobj(r, fh, 1 << 20)
            if tmp.stat().st_size == 0:
                raise OSError("empty file")
            tmp.replace(dest)
            return None
        except Exception as e:  # network errors of every kind
            err = f"{type(e).__name__}: {e}"
            time.sleep(2 ** (k + 1))
    tmp.unlink(missing_ok=True)
    return err


def get_json(url: str):
    req = urllib.request.Request(url, headers=dict(UA, Accept="application/json"))
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def fetch(name: str, urls: list, d: Path, reuse: list, manifest: list) -> bool:
    dest = d / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        manifest.append((name, dest.stat().st_size, "present"))
        return True
    for r in reuse:
        src = r / Path(name).name
        if src.exists() and src.stat().st_size > 0:
            dest.symlink_to(src.resolve())
            manifest.append((name, src.stat().st_size, f"linked from {src}"))
            print(f"{name}: linked from {src}")
            return True
    err = None
    for u in urls:
        print(f"{name}: {u}", flush=True)
        err = get(u, dest)
        if err is None:
            manifest.append((name, dest.stat().st_size, u))
            return True
        print(f"  failed: {err}")
    manifest.append((name, 0, f"ERROR {err}"))
    return False


def zenodo(d: Path, record: str, wanted: tuple, manifest: list) -> None:
    """the files of a Zenodo record whose name contains one of `wanted`"""
    try:
        rec = get_json(f"https://zenodo.org/api/records/{record}")
    except Exception as e:
        manifest.append((f"zenodo:{record}", 0, f"ERROR {type(e).__name__}: {e}"))
        print(f"Zenodo {record}: {e}")
        return
    files = rec.get("files", [])
    print(f"Zenodo {record}: {[f.get('key') for f in files]}")
    for f in files:
        key = f.get("key", "")
        if not any(w in key.lower() for w in wanted):
            continue
        url = (f.get("links") or {}).get("self") or \
            f"https://zenodo.org/records/{record}/files/{key}?download=1"
        fetch(f"kopajtich/{key}", [url], d, [], manifest)


def panelapp(d: Path, name: str, manifest: list) -> list:
    out = d / "panelapp_genes_all.tsv"
    if not out.exists():
        panels, nxt = [], f"{PANELAPP}/panels/?page=1"
        while nxt:
            page = get_json(nxt)
            panels += page.get("results", [])
            nxt = page.get("next")
        hits = [p for p in panels if p["name"].lower().strip() == name.lower()] or \
            [p for p in panels if name.lower() in p["name"].lower()]
        if len(hits) != 1:
            raise OSError(f"{len(hits)} PanelApp panels match {name!r}")
        pid = hits[0]["id"]
        p = get_json(f"{PANELAPP}/panels/{pid}/")
        with open(out, "w") as fh:
            fh.write(f"# PanelApp {p.get('name')} id {pid} version {p.get('version')}\n")
            fh.write("gene\tconfidence_level\tentity_type\n")
            for g in p.get("genes", []):
                fh.write(f"{g['gene_data']['gene_symbol']}\t{g['confidence_level']}\t"
                         f"{g.get('entity_type', 'gene')}\n")
    rows = [ln.split("\t") for ln in out.read_text().splitlines()[2:] if ln.strip()]
    manifest.append((out.name, out.stat().st_size, "PanelApp API"))
    return sorted({r[0] for r in rows if r and r[0]})


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[1].strip(),
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--datadir", type=Path, default=Path("data"))
    p.add_argument("--reuse", type=Path, action="append", default=[],
                   help="directory searched first for the files (repeatable)")
    p.add_argument("--genes", type=Path, default=Path("genes.txt"))
    p.add_argument("--panel-name", default="Mitochondrial disorders")
    p.add_argument("--no-genome", action="store_true",
                   help="skip the GRCh38 FASTA (no promoter sequence features)")
    p.add_argument("--validation", action="store_true",
                   help="also the L1 validation data (Kopajtich, Goncalves)")
    a = p.parse_args()
    d = a.datadir
    d.mkdir(parents=True, exist_ok=True)
    reuse = [r for r in a.reuse if r.is_dir()]
    manifest, ok = [], True
    for name, urls in SOURCES.items():
        if name == "hg38.fa.gz":
            if a.no_genome or (d / "GRCh38.fa").exists():
                continue
            if fetch(name, urls, d, reuse, manifest):
                # promoter_features.py reads an uncompressed FASTA (or a
                # bgzipped one with pysam): plain gunzip
                print("hg38.fa.gz: decompressing to GRCh38.fa", flush=True)
                with gzip.open(d / name, "rb") as fi, open(d / "GRCh38.fa.part", "wb") as fo:
                    shutil.copyfileobj(fi, fo, 1 << 22)
                (d / "GRCh38.fa.part").replace(d / "GRCh38.fa")
                (d / name).unlink()
            else:
                ok = False
            continue
        ok &= fetch(name, urls, d, reuse, manifest)
    if a.validation:                           # optional: never a failure
        for name, urls in VALIDATION.items():
            fetch(name, urls, d, reuse, manifest)
        zenodo(d, ZENODO_RECORD, ZENODO_WANTED, manifest)
    # the gene set
    try:
        panel = panelapp(d, a.panel_name, manifest)    # kept when present
    except Exception as e:
        print(f"PanelApp failed ({type(e).__name__}: {e}); MitoCarta genes only")
        manifest.append(("panelapp_genes_all.tsv", 0, f"ERROR {e}"))
        panel = []
    mc = []
    if (d / "Human.MitoCarta3.0.xls").exists():
        from regulome import read_mitocarta
        try:
            genes, _ = read_mitocarta(d / "Human.MitoCarta3.0.xls", None)
            mc = sorted(genes)
        except (Exception, SystemExit) as e:
            print(f"MitoCarta3.0 unreadable ({type(e).__name__}: {e})")
            ok = False
    if not mc:
        print("no MitoCarta3.0 gene: genes.txt not written")
        ok = False
    gs = sorted(set(mc) | set(panel))
    if mc:
        a.genes.write_text("\n".join(gs) + "\n")
    print(f"gene set: {len(mc)} MitoCarta3.0 + {len(set(panel) - set(mc))} PanelApp "
          f"outside MitoCarta = {len(gs)} genes -> {a.genes}")
    with open(d / "MANIFEST.tsv", "w") as fh:
        fh.write("file\tbytes\tsource\n")
        for r in manifest:
            fh.write("\t".join(map(str, r)) + "\n")
    bad = [r[0] for r in manifest if str(r[2]).startswith("ERROR")]
    if bad:
        print("missing (put them in the data directory by hand, same name, then rerun): "
              + ", ".join(bad))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
