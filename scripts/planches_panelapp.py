#!/usr/bin/env python3
"""
planches_panelapp.py - one plate per gene of the PanelApp "Mitochondrial
disorders" panel: the locus figure (Hi-C map, loops, rE2G elements, AG2 tracks,
genes, virtual 4C; locus_figure.py) on the left and the AG2 zooms of the gene
(10-bp deletions and mutagenesis of the promoter and of the scanned elements;
reannot mcia/ag2_figure_locus.py --gene) on the right, in one PNG.

For each gene of the panel (--niveaux: PanelApp confidence levels, 3 = green
by default), in one tissue (--tissu, the name shared by loops_par_tissu.py and
the AG2 campaign: muscle, cardiomyocyte, neurone, rein):
  - Hi-C library: the deepest kept library of <dir>/<tissu>/libraries.tsv, its
    .hic read at a distance from the URL of --bibliotheque (the library table
    that holds s3_url / download_url; the per-tissue table has no URL); without
    --bibliotheque, the first such .tsv under the current folder that names
    the library;
  - window: as choisir_locus.py, max(300 kb, farthest loop of the gene in the
    tissue + 150 kb), 500 kb for a gene without loop; capped at 1 Mb;
  - a gene without AG2 zoom (not in the campaign, nothing scanned in the
    tissue) keeps its locus figure alone (statut macro_seul).
Only the plate is kept (PNG, --dpi 150 by default, a few hundred kB); the
intermediate SVG/PNG/TSV are written in a temporary folder and removed.
<out>/index.tsv: gene, PanelApp level, loop in the tissue, status, plate.
--reprendre skips the genes whose plate exists, so a stopped run resumes.
--liste prints the genes and the expected count, without drawing anything.
--tranche K/N draws one share of the genes (one task of a SLURM array,
slurm/planches_panelapp.sbatch), with its own index_K.tsv; --fusionner then
writes index.tsv from them.

PanelApp symbols are old HGNC symbols in places (ATP5A1): they are renamed with
--renommes (genes_renamed.tsv of fetch_regulome.py, old -> GTF symbol).

Needs network access for the .hic (run it on a login node). Dependencies:
matplotlib, numpy, hicstraw, Pillow (the mitodark environment). No AG2 call.

Usage (from the mitoloop_v0 folder)
    python3 scripts/planches_panelapp.py --tissu muscle --bibliotheque libraries.tsv
        --elements ~/reannot_git/results/paquet/re2g_tous/re2g_elements.tsv
        --ag2-dir ~/reannot_git/results/ag2_v2 --reannot ~/reannot_git
        --out results/planches_panelapp/muscle
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_tsv(path):
    with open(path) as fh:
        lines = [l for l in fh if not l.startswith("# ")]
    head = lines[0].rstrip("\n").lstrip("#").split("\t")
    for line in lines[1:]:
        if line.strip():
            f = line.rstrip("\n").split("\t")
            yield dict(zip(head, f + [""] * (len(head) - len(f))))


def genes_panel(path, niveaux, renommes):
    """[(gene under the GTF symbol, level)] of the panel, levels kept, panel order"""
    ren = {}
    if renommes and Path(renommes).exists():
        for r in read_tsv(renommes):
            vals = list(r.values())
            if len(vals) >= 2 and vals[0] and vals[1]:
                ren[vals[0]] = vals[1]
    out, vus = [], set()
    for r in read_tsv(path):
        if r.get("entity_type", "gene") not in ("gene", ""):
            continue
        if r.get("confidence_level", "") not in niveaux:
            continue
        g = ren.get(r["gene"], r["gene"])
        if g not in vus:
            vus.add(g)
            out.append((g, r["confidence_level"]))
    return out


def librairie(d, tissu):
    """accession of the deepest kept library of the tissue"""
    rows = [r for r in read_tsv(d / tissu / "libraries.tsv") if r.get("kept", "1") in ("1", "oui", "True")]
    if not rows:
        sys.exit(f"aucune librairie gardée dans {d / tissu / 'libraries.tsv'}")
    return max(rows, key=lambda r: float(r.get("depth") or 0))["accession"]


def fenetres(d, tissu, plafond):
    """gene -> window of the locus figure (as choisir_locus.py)"""
    far = {}
    for r in read_tsv(d / tissu / "loops.tsv"):
        try:
            dist = abs(int(r["summit"]) - int(r["vp_pos"]))
        except (KeyError, ValueError):
            continue
        far[r["gene"]] = max(far.get(r["gene"], 0), dist)
    return {g: min(plafond, max(300_000, x + 150_000)) for g, x in far.items()}


def trouver_bibliotheque(racine, acc):
    """first .tsv under racine (results/ and data/ first) whose header has s3_url or
    download_url and that names the library acc (the per-tissue tables have no URL)"""
    vus = []
    for sous in ("data", "results", "."):
        base = racine / sous
        if base.is_dir():
            vus += sorted(base.rglob("*.tsv")) if sous != "." else sorted(base.glob("*.tsv"))
    for f in vus:
        try:
            with open(f) as fh:
                tete = fh.readline()
                if "s3_url" not in tete and "download_url" not in tete:
                    continue
                if any(acc in line for line in fh):
                    return str(f)
        except (OSError, UnicodeDecodeError):
            continue
    return ""


def planche(macro, zoom, out):
    """macro on the left, zoom on the right, top-aligned, white background"""
    from PIL import Image
    a = Image.open(macro).convert("RGB")
    if zoom is None:
        a.save(out, optimize=True)
        return
    b = Image.open(zoom).convert("RGB")
    if b.height > a.height:                         # the zoom never taller than the locus figure
        b = b.resize((round(b.width * a.height / b.height), a.height), Image.LANCZOS)
    c = Image.new("RGB", (a.width + b.width + 20, max(a.height, b.height)), "white")
    c.paste(a, (0, 0))
    c.paste(b, (a.width + 20, 0))
    c.save(out, optimize=True)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--panelapp", default="data/panelapp_genes_all.tsv")
    p.add_argument("--renommes", default="data/genes_renamed.tsv")
    p.add_argument("--niveaux", default="3", help="PanelApp levels kept, comma-separated (3 = green)")
    p.add_argument("--genes", default="", help="restrict to these genes (comma-separated)")
    p.add_argument("--tissu", required=True)
    p.add_argument("--dir", default="results/loops_tissus")
    p.add_argument("--bibliotheque", default="", help="library table with the .hic URL (s3_url); found under the current folder when not given")
    p.add_argument("--gtf", default="data/gencode.gtf.gz")
    p.add_argument("--elements", default="")
    p.add_argument("--elements-frac", default="0.8")
    p.add_argument("--ag2-dir", default="")
    p.add_argument("--reannot", default=str(Path.home() / "reannot_git"),
                   help="reannot clone, for mcia/ag2_figure_locus.py")
    p.add_argument("--fenetre-max", type=int, default=1_000_000)
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--out", required=True)
    p.add_argument("--reprendre", action="store_true")
    p.add_argument("--liste", action="store_true", help="list the genes, draw nothing")
    p.add_argument("--tranche", default="", help="K/N : only the genes of rank K, K+N, K+2N... (K from 0), "
                   "for a SLURM array; index written to index_K.tsv (merged by --fusionner)")
    p.add_argument("--fusionner", action="store_true", help="merge the index_K.tsv of --out into index.tsv")
    a = p.parse_args()

    if a.fusionner:
        out = Path(a.out)
        parts = sorted(out.glob("index_*.tsv"), key=lambda f: int(f.stem.split("_")[1]))
        rows = [r for f in parts for r in read_tsv(f)]
        cols = ["gene", "niveau", "boucle", "statut", "raison", "planche"]
        with open(out / "index.tsv", "w") as fo:
            fo.write("\t".join(cols) + "\n")
            for x in sorted(rows, key=lambda r: r["gene"]):
                fo.write("\t".join(x.get(c, "") for c in cols) + "\n")
        n = {st: sum(1 for x in rows if x.get("statut") == st) for st in ("ok", "macro_seul", "echec")}
        print(f"{len(parts)} tranches, {len(rows)} gènes : {n['ok']} complètes, {n['macro_seul']} sans zoom, "
              f"{n['echec']} en échec → {out / 'index.tsv'}", file=sys.stderr)
        return
    d = Path(a.dir)
    genes = genes_panel(a.panelapp, set(a.niveaux.split(",")), a.renommes)
    if a.genes:
        garder = set(a.genes.split(","))
        genes = [x for x in genes if x[0] in garder]
    if a.tranche:
        k, n = (int(x) for x in a.tranche.split("/"))
        genes = genes[k::n]
    fen = fenetres(d, a.tissu, a.fenetre_max)
    avec = sum(1 for g, _ in genes if g in fen)
    print(f"{len(genes)} gènes PanelApp (niveaux {a.niveaux}) ; {avec} avec boucle dans {a.tissu}",
          file=sys.stderr)
    if a.liste:
        for g, n in genes:
            print(f"{g}\t{n}\t{'oui' if g in fen else 'non'}")
        return
    acc = librairie(d, a.tissu)
    if not a.bibliotheque:
        a.bibliotheque = trouver_bibliotheque(Path.cwd(), acc)
        if not a.bibliotheque:
            sys.exit(f"aucune table .tsv avec une colonne s3_url ou download_url et la librairie {acc} "
                     f"sous {Path.cwd()} : la donner avec --bibliotheque")
        print(f"table des librairies (URL du .hic) : {a.bibliotheque}", file=sys.stderr)
    elif not Path(a.bibliotheque).exists():
        sys.exit(f"--bibliotheque : {a.bibliotheque} n'existe pas")
    zoom_py = Path(a.reannot) / "mcia" / "ag2_figure_locus.py"
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    index = out / (f"index_{a.tranche.split('/')[0]}.tsv" if a.tranche else "index.tsv")
    deja = {}
    if a.reprendre and index.exists():
        deja = {r["gene"]: r for r in read_tsv(index)}
    rows = []
    t0 = time.time()
    for i, (g, niv) in enumerate(genes, 1):
        png = out / f"{g}_{a.tissu}.png"
        if a.reprendre and png.exists() and g in deja:
            rows.append(deja[g])
            continue
        tmp = Path(tempfile.mkdtemp(prefix=f"planche_{g}_"))
        statut, raison = "ok", ""
        try:
            cmd = [sys.executable, str(HERE / "locus_figure.py"), "--tissues", a.bibliotheque,
                   "--term-column", "accession", "--term", acc, "--gene", g, "--gtf", a.gtf,
                   "--loops", str(d / a.tissu / "loops.tsv"), "--arcs", "--genes-track",
                   "--window", str(fen.get(g, 500_000)), "--max-distance", str(fen.get(g, 500_000)),
                   "--norm", "SCALE", "--dpi", str(a.dpi), "--out", str(tmp / "macro")]
            if a.elements:
                cmd += ["--elements", a.elements, "--elements-tissu", a.tissu, "--elements-frac", a.elements_frac]
            if a.ag2_dir:
                cmd += ["--ag2-dir", a.ag2_dir, "--ag2-type", a.tissu]
            r = subprocess.run(cmd, capture_output=True, text=True)
            macro = tmp / "macro.png"
            if r.returncode or not macro.exists():
                statut = "echec"
                raison = (r.stderr.strip().splitlines() or ["locus_figure.py sans sortie"])[-1][:200]
            else:
                zoom = None
                if a.ag2_dir and zoom_py.exists():
                    rz = subprocess.run([sys.executable, str(zoom_py), "--dir", a.ag2_dir, "--type", a.tissu,
                                         "--gene", g, "--dpi", str(a.dpi), "--out", str(tmp / "zoom")],
                                        capture_output=True, text=True)
                    if rz.returncode == 0 and (tmp / "zoom.png").exists():
                        zoom = tmp / "zoom.png"
                    else:
                        statut = "macro_seul"
                        raison = (rz.stderr.strip().splitlines() or ["pas de zoom"])[-1][:200]
                else:
                    statut = "macro_seul"
                    raison = "pas de --ag2-dir ou ag2_figure_locus.py introuvable"
                planche(macro, zoom, png)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        rows.append({"gene": g, "niveau": niv, "boucle": "oui" if g in fen else "non", "statut": statut,
                     "raison": raison, "planche": png.name if statut != "echec" else ""})
        el = time.time() - t0
        print(f"[{i}/{len(genes)}] {g} : {statut}{' (' + raison + ')' if raison else ''} ; "
              f"{el / 60:.1f} min écoulées", file=sys.stderr, flush=True)
        with open(index, "w") as fo:                   # rewritten at each gene: a stopped run keeps its index
            cols = ["gene", "niveau", "boucle", "statut", "raison", "planche"]
            fo.write("\t".join(cols) + "\n")
            for x in rows:
                fo.write("\t".join(str(x.get(c, "")) for c in cols) + "\n")
    n = {s: sum(1 for x in rows if x.get("statut") == s) for s in ("ok", "macro_seul", "echec")}
    print(f"planches : {n['ok']} complètes, {n['macro_seul']} sans zoom, {n['echec']} en échec → {out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
