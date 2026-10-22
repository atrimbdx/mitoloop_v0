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
    with --promoter-annotations: promoters_exp/ (gnomAD v4.1 gene
    constraint, Zoonomia phyloP 241 mammals, EPDnew, FANTOM5 CAGE peaks,
    SCREEN cCREs, Roadmap ChromHMM of skeletal muscle, brain, left
    ventricle and liver, ReMap ChIP-seq peaks of 11 factors; see
    config/promoter_annotations.tsv),
    and with --validation: kopajtich/ (fibroblast RNA-seq and proteomics,
    datasets.zip of Zenodo 4501904, extracted), attenuation.csv,
    GTEx_v8_finemapping_DAPG.txt.gz (fine-mapped eQTL, validation of the
    reannot arguments), gismo/ (GISMO, Zenodo 20547722) (Goncalves et al. 2017);
    MANIFEST.tsv (file, bytes, source or error);
    --genes (default genes.txt): MitoCarta3.0 symbols and PanelApp genes,
    under the symbol of the --gtf when one is given (PanelApp keeps old HGNC
    symbols such as ATP5A1: without this, a gene can be counted twice and
    lost by the steps that match on the GTF symbol); genes_renamed.tsv in
    --datadir lists the renames.

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
# are read from the record. The record holds datasets.zip (processed tables),
# raw_data.zip (count and protein tables, read by l1_score.py validate once
# converted) and FRASER_counts (splicing counts, hg19): the first two. Of the
# archives, the tables are extracted, and every member is listed in
# kopajtich/CONTENTS.tsv.
ZENODO_RECORD = "4501904"
ZENODO_WANTED = ("datasets", "raw_data")
ZENODO_MAX_GB = 5
TABLE = (".tsv", ".csv", ".txt", ".tsv.gz", ".csv.gz", ".txt.gz", ".rds", ".xlsx")
# GISMO and GISMO-mis (Liao et al. 2026, Nat Genet): gene loss and missense
# constraint across 462 mammals, every file of the record, in gismo/
GISMO_RECORD = "20547722"
PANELAPP = "https://panelapp.genomicsengland.co.uk/api/v1"
# experimental promoter annotations (config/promoter_annotations.tsv), in
# promoters_exp/; ReMap: one file per factor
ROADMAP = ("https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/"
           "ChmmModels/coreMarks/jointModel/final/")
REMAP = "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/TF/"
REMAP_TFS = ("NRF1", "GABPA", "YY1", "ESRRA", "SP1", "NFYA", "CREB1", "CTCF", "TBP",
             "MYC", "RELA")             # MYC, RELA : recrutent P-TEFb (sortie de pause)
PROMOTER_SOURCES = {
    "Hs_EPDnew_006_hg38.bed": [
        "https://epd.expasy.org/ftp/epdnew/H_sapiens/006/Hs_EPDnew_006_hg38.bed",
        "ftp://ccg.epfl.ch/epdnew/H_sapiens/006/Hs_EPDnew_006_hg38.bed"],
    "hg38_fair+new_CAGE_peaks_phase1and2.bed.gz": [
        "https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/"
        "CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz"],
    "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt": [
        "https://jaspar.elixir.no/download/data/2024/CORE/"
        "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
        "https://jaspar.genereg.net/download/data/2024/CORE/"
        "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
        "https://jaspar2024.genereg.net/download/data/2024/CORE/"
        "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"],
    # gene constraint of gnomAD v4.1 (LoF obs/exp, o/e, LOEUF, pLI, missense)
    "gnomad.v4.1.constraint_metrics.tsv": [
        "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/"
        "constraint/gnomad.v4.1.constraint_metrics.tsv"],
    # base-wise conservation over 241 mammals (Zoonomia, phyloP; several GB,
    # cut to the regions of the gene set when the reannot package is built)
    "cactus241way.phyloP.bw": [
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus241way/cactus241way.phyloP.bw",
        "https://hgdownload.gi.ucsc.edu/goldenPath/hg38/cactus241way/cactus241way.phyloP.bw"],
    # gnomAD-SV v4.1 (1.3 GB; reduced to PASS SV, type and frequency, ~15 MB,
    # when the reannot package is built)
    "gnomad.v4.1.sv.sites.bed.gz": [
        "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/"
        "gnomad.v4.1.sv.sites.bed.gz"],
    "GRCh38-cCREs.bed": [
        "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed",
        "https://downloads.wenglab.org/V3/GRCh38-cCREs.bed"],
    **{f"{e}_15_coreMarks_hg38lift_mnemonics.bed.gz":
       [f"{ROADMAP}{e}_15_coreMarks_hg38lift_mnemonics.bed.gz"]
       for e in ("E108", "E073", "E095", "E066")},   # muscle, cerveau, ventricule gauche, foie
    **{f"remap2022_{t}_nr_macs2_hg38_v1_0.bed.gz":
       [f"{REMAP}{t}/remap2022_{t}_nr_macs2_hg38_v1_0.bed.gz"] for t in REMAP_TFS},
}


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


def zenodo(d: Path, record: str, wanted: tuple, manifest: list,
           sub: str = "kopajtich") -> None:
    """the files of a Zenodo record whose name contains one of `wanted`
    (every file when `wanted` is empty), in d/sub"""
    try:
        rec = get_json(f"https://zenodo.org/api/records/{record}")
    except Exception as e:
        manifest.append((f"zenodo:{record}", 0, f"ERROR {type(e).__name__}: {e}"))
        print(f"Zenodo {record}: {e}")
        return
    files = rec.get("files", [])
    (d / sub).mkdir(parents=True, exist_ok=True)
    (d / sub / "CONTENTS.tsv").write_text("archive\tmember\tbytes\n")
    print(f"Zenodo {record}: " + ", ".join(
        f"{f.get('key')} ({f.get('size', 0) / 1e6:,.0f} MB)" for f in files))
    for f in files:
        key = f.get("key", "")
        if wanted and not any(w in key.lower() for w in wanted):
            continue
        if f.get("size", 0) > ZENODO_MAX_GB * 1e9:
            manifest.append((f"{sub}/{key}", f.get("size", 0),
                             f"not downloaded: over {ZENODO_MAX_GB} GB"))
            print(f"{key}: over {ZENODO_MAX_GB} GB, not downloaded")
            continue
        url = (f.get("links") or {}).get("self") or \
            f"https://zenodo.org/records/{record}/files/{key}?download=1"
        if fetch(f"{sub}/{key}", [url], d, [], manifest) and key.endswith(".zip"):
            import zipfile
            try:
                with zipfile.ZipFile(d / sub / key) as z:
                    info = [i for i in z.infolist() if not i.is_dir()
                            and "__MACOSX" not in i.filename
                            and not i.filename.endswith(".DS_Store")]
                    with open(d / sub / "CONTENTS.tsv", "a") as fh:
                        for i in info:
                            fh.write(f"{key}\t{i.filename}\t{i.file_size}\n")
                    # tables only, under 3 GB each (macOS files left out)
                    names = [i.filename for i in info
                             if i.filename.lower().endswith(TABLE)
                             and i.file_size < 3e9]
                    z.extractall(d / sub, names)
                print(f"{key}: {len(names)} files extracted to {d / sub}")
                for n in names:
                    manifest.append((f"{sub}/{n}", (d / sub / n).stat().st_size,
                                     f"from {key}"))
            except Exception as e:
                manifest.append((f"{sub}/{key}", 0, f"unzip failed: {e}"))


DAPG_TAR = ("https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/"
            "fine-mapping-cis-eqtl/GTEx_v8_finemapping_DAPG.tar")
DAPG_MEMBER = "GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.txt.gz"


def gtex_dapg(d: Path, manifest: list) -> None:
    """GTEx v8 fine-mapped eQTL (DAP-G), the table only: read from the
    760 MB release tar as a stream, nothing else kept"""
    import tarfile
    dest = d / Path(DAPG_MEMBER).name
    if dest.exists() and dest.stat().st_size > 0:
        manifest.append((dest.name, dest.stat().st_size, "present"))
        return
    try:
        req = urllib.request.Request(DAPG_TAR, headers=UA)
        with urllib.request.urlopen(req, timeout=600) as r, \
                tarfile.open(fileobj=r, mode="r|") as tf:
            for m in tf:
                if m.name == DAPG_MEMBER:
                    tmp = dest.with_name(dest.name + ".part")
                    with tf.extractfile(m) as fi, open(tmp, "wb") as fo:
                        shutil.copyfileobj(fi, fo, 1 << 20)
                    tmp.replace(dest)
                    break
        manifest.append((dest.name, dest.stat().st_size, DAPG_TAR))
        print(f"{dest.name}: {dest.stat().st_size / 1e6:.0f} MB")
    except Exception as e:  # network errors of every kind
        manifest.append((dest.name, 0, f"ERROR {type(e).__name__}: {e}"))
        print(f"GTEx DAP-G failed: {e}")


def canonical(symbols: list, gtf: Path, d: Path) -> dict:
    """symbol -> its name in the GTF (renamed HGNC symbols resolved with
    the alias table of re2g_panel), None for panel entries that are not
    genes; symbols absent from the GTF are left out (kept as they are)"""
    from re2g_panel import build_aliases, gtf_gene_names, resolve_symbols
    m, unresolved, not_genes = resolve_symbols(symbols, gtf_gene_names(gtf),
                                               build_aliases())
    with open(d / "genes_renamed.tsv", "w") as fh:
        fh.write("symbol\tgtf_symbol\n")
        for g in symbols:
            if m.get(g, g) != g:
                fh.write(f"{g}\t{m[g]}\n")
            elif g in not_genes:
                fh.write(f"{g}\tnot a gene\n")
            elif g in unresolved:
                fh.write(f"{g}\tnot in the GTF\n")
    ren = dict(m, **{g: None for g in not_genes})
    print(f"symbols: {sum(m[g] != g for g in m)} renamed to the GTF symbol, "
          f"{len(not_genes)} not genes dropped, {len(unresolved)} not in the GTF "
          f"(kept); see {d / 'genes_renamed.tsv'}")
    if unresolved:
        print("  not in the GTF: " + " ".join(unresolved[:40]))
    return ren


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
    p.add_argument("--gtf", type=Path,
                   help="GENCODE GTF: write the genes under its symbols")
    p.add_argument("--no-genome", action="store_true",
                   help="skip the GRCh38 FASTA (no promoter sequence features)")
    p.add_argument("--promoter-annotations", action="store_true",
                   help="also the experimental promoter annotations (EPDnew, "
                        "FANTOM5, SCREEN, Roadmap ChromHMM, ReMap) in promoters_exp/")
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
    if a.promoter_annotations:                 # optional: never a failure
        for name, urls in PROMOTER_SOURCES.items():
            fetch(f"promoters_exp/{name}", urls, d, reuse, manifest)
    if a.validation:                           # optional: never a failure
        for name, urls in VALIDATION.items():
            fetch(name, urls, d, reuse, manifest)
        zenodo(d, ZENODO_RECORD, ZENODO_WANTED, manifest)
        zenodo(d, GISMO_RECORD, (), manifest, sub="gismo")
        gtex_dapg(d, manifest)
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
    ren = {}
    if mc and a.gtf and a.gtf.exists():
        ren = canonical(sorted(set(mc) | set(panel)), a.gtf, d)
    mc = {ren.get(g, g) for g in mc} - {None}
    panel = {ren.get(g, g) for g in panel} - {None}
    gs = sorted(mc | panel)
    if mc:
        a.genes.write_text("\n".join(gs) + "\n")
    print(f"gene set: {len(mc)} MitoCarta3.0 + {len(panel - mc)} PanelApp "
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
