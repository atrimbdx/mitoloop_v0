#!/usr/bin/env python3
"""
locus_figure.py - Hi-C map, loop arcs, genes and virtual 4C at one viewpoint

Purpose
    Draws, for one viewpoint, the views that show how a virtual 4C profile
    and its called loops arise from a Hi-C matrix. The panels share one x
    axis, from top to bottom:

      the contact map, as a triangle rotated by 45 degrees so that the
      diagonal lies along the bottom edge and each contact sits above the
      midpoint of its two positions (observed / expected by default)
      an optional arc panel of the loops of loops.tsv in the window, for the
      requested gene and for every other viewpoint that falls in the window
      an optional gene track read from a GTF
      the virtual 4C profile of the viewpoint, with the called loops marked

    The figure illustrates the method at one locus; it is not a measurement.

Inputs
    --hic        .hic (local path or URL, read with hicstraw) or .cool /
                 .mcool (read with cooler). Optional: without it only the
                 virtual 4C (and arc and gene panels) are drawn
    --tissues    library table (accession, term_name and a location column
                 s3_url or download_url), with --term, to take the .hic URL
                 of one ENCODE term instead of --hic
    --npz        profiles of one library, in either of two key schemes:
                   as written by `v4c_pipeline.py extract`: counts, names,
                   chrom, pos, offsets, accession, term, resolution,
                   anchor_half (among others)
                   or: profiles, genes, gene_set, offsets, accession, tissue,
                   resolution; the viewpoint position is then taken from
                   loops.tsv (vp_pos) or from the GTF
    --gene/--gtf gene symbol and GTF (GENCODE style attributes gene_name,
                 gene_type; transcripts tagged MANE_Select are preferred for
                 the TSS); or --chrom and --viewpoint
    --loops      loops.tsv from `v4c_loops.R reciprocal`: gene, chrom,
                 vp_pos, summit, n_terms (called, start, end, loop_id unused)

Outputs (prefix --out)
    <out>.svg, <out>.png   the figure (SVG text kept as text)
    <out>.tsv              chrom, start, end, distance, virtual_4c: the
                           profile drawn, one row per bin; distance is the bin
                           centre minus the viewpoint

Method
    The window spans the viewpoint +/- --window, snapped to bins. The map is
    the dense matrix of the window; by default each diagonal is divided by
    its mean (observed / expected), which removes the distance decay and
    leaves the loops visible. In the rotated view a bin (i, j) is a square
    turned by 45 degrees, centred at ((i + j) / 2, (j - i) / 2), so both axes
    must carry the same number of bp per inch: the map is drawn with an equal
    aspect and its height follows from --max-distance and the figure width.
    The equal aspect shrinks the map's box while the other panels keep the
    full width; once the figure is laid out, every panel is given the map's
    horizontal extent, so that one genomic coordinate falls at the same x in
    all panels.
    Without --npz, the virtual 4C is the viewpoint row of the matrix, summed
    over the anchor bins (viewpoint +/- --anchor-half), as in the pipeline.
    With --npz, the stored profile is drawn, so the panel shows the profile
    that the loop calling used. The anchor bins are masked, since they are
    not contacts. A called loop is the pair (viewpoint, summit); it is marked
    above the profile at the summit and, on the rotated map, above the pixel
    at their midpoint and half their separation.

Usage
    python3 locus_figure.py --hic data/ENCFF000AAA.hic --gene GENE1 \\
        --gtf data/gencode.gtf.gz --loops results/loops/loops.tsv \\
        --arcs --genes-track --out figures/locus_GENE1
    python3 locus_figure.py --npz results/raw/ENCFF000AAA.npz --gene GENE1 \\
        --loops results/loops/loops.tsv --out figures/v4c_GENE1
    python3 locus_figure.py --tissues tissues.tsv --term "psoas muscle" \\
        --gene GENE1 --gtf data/gencode.gtf.gz --out figures/locus_GENE1

References
    Durand NC et al. 2016. Juicer provides a one-click system for analyzing
        loop-resolution Hi-C experiments. Cell Systems. (.hic format,
        hicstraw)
    Abdennur N, Mirny LA. 2020. Cooler: scalable storage for Hi-C data and
        other genomically labeled arrays. Bioinformatics.
    Frankish A et al. 2023. GENCODE: reference annotation for the human and
        mouse genomes in 2023. Nucleic Acids Research. doi:10.1093/nar/gkac1071
    Morales J et al. 2022. A joint NCBI and EMBL-EBI transcript set for
        clinical genomics and research. Nature. (MANE Select)
    Hunter JD. 2007. Matplotlib: a 2D graphics environment. Computing in
        Science and Engineering.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np

RED, BLUE, GOLD, GREY = "#B2182B", "#2166AC", "#C78A20", "#5A5A5A"
# colour ramps of the map: raw contacts, and observed / expected around 1
RAW_RAMP = ["#FFFFFF", "#FFF0E0", "#F7C5A0", "#E8794B", RED, "#5E0E1B"]
OE_RAMP = [BLUE, "#FFFFFF", RED]


# ---------------------------------------------------------------------------
# readers

class HicReader:
    """A .hic through hicstraw. Returns a dense square of the window."""

    def __init__(self, path, norm, resolution):
        import hicstraw
        self.h = hicstraw.HiCFile(path)
        self.res = resolution
        self.norm = norm
        res_available = set(self.h.getResolutions())
        if resolution not in res_available:
            sys.exit(f"{resolution} bp absent from the file; available: "
                     f"{sorted(res_available)}")
        self.chroms = {c.name: c.length for c in self.h.getChromosomes()}

    def name(self, chrom):
        if chrom in self.chroms:
            return chrom
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if alt in self.chroms:
            return alt
        sys.exit(f"{chrom} absent; the file names its chromosomes "
                 f"{sorted(self.chroms)[:4]}...")

    def square(self, chrom, start, end):
        c = self.name(chrom)
        for norm in (self.norm, "NONE"):
            try:
                mzd = self.h.getMatrixZoomData(c, c, "observed", norm,
                                               "BP", self.res)
                m = mzd.getRecordsAsMatrix(start, end - 1, start, end - 1)
            except Exception as exc:
                print(f"  {norm}: {type(exc).__name__}: {exc}")
                continue
            if m is not None and np.size(m):
                if norm != self.norm:
                    print(f"  {self.norm} unavailable, using NONE")
                return np.asarray(m, dtype=float)
        sys.exit("empty matrix; check the locus and the resolution")


class CoolReader:
    """A .cool or .mcool through cooler."""

    def __init__(self, path, norm, resolution):
        import cooler
        if "::" not in path and path.endswith(".mcool"):
            path = f"{path}::resolutions/{resolution}"
        self.c = cooler.Cooler(path)
        self.res = self.c.binsize
        if self.res != resolution:
            print(f"  cooler holds {self.res} bp bins, using those")
        self.balance = norm.upper() not in ("NONE", "")

    def name(self, chrom):
        names = set(self.c.chromnames)
        if chrom in names:
            return chrom
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if alt in names:
            return alt
        sys.exit(f"{chrom} absent from the cooler")

    def square(self, chrom, start, end):
        c = self.name(chrom)
        reg = f"{c}:{start}-{end}"
        try:
            m = self.c.matrix(balance=self.balance).fetch(reg)
        except Exception:
            m = self.c.matrix(balance=False).fetch(reg)
            print("  no balancing weights, raw counts used")
        return np.asarray(m, dtype=float)


def open_reader(path, norm, resolution):
    p = path.split("::")[0].lower()
    if p.endswith(".hic"):
        return HicReader(path, norm, resolution)
    if p.endswith((".cool", ".mcool")):
        return CoolReader(path, norm, resolution)
    sys.exit("--hic must end in .hic, .cool or .mcool")


# ---------------------------------------------------------------------------

def tss_of(gtf, gene, tag="MANE_Select"):
    """TSS of the transcript carrying `tag` where there is one, otherwise of
    the gene record. Positions are 0-based."""
    best = None
    op = gzip.open if str(gtf).endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] not in ("gene", "transcript"):
                continue
            i = f[8].find('gene_name "')
            if i < 0 or f[8][i + 11:f[8].find('"', i + 11)] != gene:
                continue
            rank = 2 if (f[2] == "transcript" and tag and tag in f[8]) else \
                   (1 if f[2] == "gene" else 0)
            if rank and (best is None or rank > best[2]):
                pos = int(f[3]) - 1 if f[6] == "+" else int(f[4])
                best = (f[0], pos, rank, f[6])
    if best is None:
        sys.exit(f"{gene} not found in {gtf}")
    return best[0], best[1], best[3]


def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.strip():
                f = line.rstrip("\n").split("\t")
                yield dict(zip(head, f + [""] * (len(head) - len(f))))


def same_chrom(c, ref):
    """Chromosome names compared with or without the chr prefix."""
    return c in (ref, ref.lstrip("chr"), "chr" + ref.lstrip("chr"))


def pick(npz, what, *keys):
    """The first of `keys` present in the archive, or a readable failure."""
    for k in keys:
        if k in npz.files:
            return npz[k]
    have = ", ".join(f"{k}{tuple(np.shape(npz[k]))}" for k in npz.files)
    sys.exit(f"no {what} in the npz: none of {list(keys)} is there.\n"
             f"The archive holds: {have}")


def scalar(npz, key, default):
    if key in npz.files:
        v = npz[key]
        return v.item() if np.ndim(v) == 0 else v
    return default


def genes_in(gtf, chrom, start, end, types):
    """Gene records overlapping the window, as (name, start, end, strand)."""
    keep, want = [], {x.strip() for x in types.split(",") if x.strip()}
    alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
    op = gzip.open if str(gtf).endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f_ = line.split("\t")
            if len(f_) < 9 or f_[2] != "gene" or f_[0] not in (chrom, alt):
                continue
            s, e = int(f_[3]) - 1, int(f_[4])
            if e < start or s > end:
                continue
            a_ = f_[8]
            i = a_.find('gene_name "')
            if i < 0:
                continue
            name = a_[i + 11:a_.find('"', i + 11)]
            if want:
                j = a_.find('gene_type "')
                kind = a_[j + 11:a_.find('"', j + 11)] if j >= 0 else ""
                if kind not in want:
                    continue
            keep.append((name, s, e, f_[6]))
    return sorted(keep, key=lambda g: g[1])


def pack_rows(genes, span, n_rows):
    """Greedy rows so that neither a gene nor its label overlaps the next."""
    ends, rows = [], []
    for g in genes:
        pad = max(span * 0.012 * (len(g[0]) + 2), span * 0.02)
        for r, last in enumerate(ends):
            if g[1] > last:
                ends[r] = g[2] + pad
                rows.append(r)
                break
        else:
            if len(ends) >= n_rows:
                rows.append(None)        # no room, left out
                continue
            ends.append(g[2] + pad)
            rows.append(len(ends) - 1)
    return rows


def hic_url(tissues, term, term_column="term_name",
            location_columns=("s3_url", "download_url", "path"),
            accession_column="accession"):
    """The location of the first library of one ENCODE term in the library
    table: an http(s) URL or an existing local file, from the first of
    `location_columns` that holds one."""
    rows = read_tsv(tissues)
    hits = [r for r in rows
            if term.lower() in (r.get(term_column, "") or "").lower()]
    if not hits:
        terms = sorted({r.get(term_column, "") for r in rows})
        sys.exit(f"no {term_column} matching {term!r}; the file holds {terms}")
    if len(hits) > 1:
        print(f"  {len(hits)} libraries match {term!r}, taking "
              f"{hits[0].get(accession_column)}")
    for key in location_columns:
        u = (hits[0].get(key, "") or "").strip()
        if u.startswith(("http://", "https://")) or (u and Path(u).exists()):
            return u
    sys.exit(f"{hits[0].get(accession_column)} has no URL or existing file "
             f"in the columns {', '.join(location_columns)}")


def observed_over_expected(m):
    """Divide each diagonal by its mean, which removes the distance decay
    and leaves the loops visible."""
    n = len(m)
    out = np.array(m, dtype=float)
    for d in range(n):
        idx = np.arange(n - d)
        v = out[idx, idx + d]
        mu = np.nanmean(v) if np.isfinite(v).any() else np.nan
        if mu and np.isfinite(mu):
            out[idx, idx + d] = v / mu
            out[idx + d, idx] = v / mu
    return out


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--hic",
                   help="a .hic (local path or URL), .cool or .mcool for the "
                        "contact map. Optional: without it only the virtual 4C "
                        "panel is drawn")
    p.add_argument("--npz",
                   help="npz of one library, as written by `v4c_pipeline.py "
                        "extract` or with keys profiles/genes/offsets. The "
                        "virtual 4C is read from it, so the panel shows the "
                        "profile the loop calling used")
    p.add_argument("--tissues",
                   help="tissues.tsv, to resolve the .hic URL of --term "
                        "instead of giving --hic a path")
    p.add_argument("--term", help="ENCODE term name, with --tissues")
    p.add_argument("--term-column", default="term_name",
                   help="column of --tissues holding the term (default "
                        "%(default)s)")
    p.add_argument("--location-columns", default="s3_url,download_url,path",
                   help="columns of --tissues tried in order for a URL or a "
                        "local file (default %(default)s)")
    p.add_argument("--gene", help="gene symbol of the viewpoint")
    p.add_argument("--gtf", help="GTF annotation (plain or gzipped)")
    p.add_argument("--tss-tag", default="MANE_Select",
                   help="transcript tag preferred for the TSS; empty uses the "
                        "gene record")
    p.add_argument("--chrom", help="chromosome, with --viewpoint")
    p.add_argument("--viewpoint", type=int,
                   help="position of the anchor, instead of --gene")
    p.add_argument("--window", type=int, default=500_000,
                   help="half-width drawn, in bp")
    p.add_argument("--resolution", type=int, default=10_000,
                   help="bin size in bp; taken from the npz when it holds one")
    p.add_argument("--norm", default="KR",
                   help="KR, SCALE, VC or NONE; the map reads better balanced, "
                        "while the pipeline calls peaks on raw counts")
    p.add_argument("--raw", action="store_true",
                   help="draw the contacts as they are. By default each "
                        "diagonal is divided by its mean, which removes the "
                        "distance decay and leaves the loops visible")
    p.add_argument("--anchor-half", type=int, default=10_000,
                   help="the viewpoint spans the TSS plus and minus this, in "
                        "bp, as in the pipeline; taken from the npz when it "
                        "holds one")
    p.add_argument("--loops", help="loops.tsv, to mark the called loops")
    p.add_argument("--max-distance", type=int, default=0,
                   help="largest separation shown by the triangle, in bp; 0 "
                        "uses --window (at least 8 bins)")
    p.add_argument("--square", action="store_true",
                   help="draw the usual square map instead of the triangle")
    p.add_argument("--width", type=float, default=6.0,
                   help="figure width in inches; the map height follows")
    p.add_argument("--map-fraction", type=float, default=0.76,
                   help="fraction of the figure width given to the map when "
                        "its height is computed; the rest holds labels and "
                        "colour bar")
    p.add_argument("--profile-height", type=float, default=1.15,
                   help="height of the virtual 4C panel in inches")
    p.add_argument("--arrow", type=float, default=20,
                   help="length of the arrows marking the loops and the "
                        "viewpoint on the map, in points")
    p.add_argument("--arcs", action="store_true",
                   help="draw an arc panel of the loops of the window. Every "
                        "gene of the panel whose viewpoint falls in the "
                        "window contributes its own loops, not only --gene")
    p.add_argument("--arcs-height", type=float, default=0.95,
                   help="height of the arc panel in inches")
    p.add_argument("--arcs-min-terms", type=int, default=0,
                   help="keep only the loops called in at least this many "
                        "ENCODE terms")
    p.add_argument("--arcs-others", default="",
                   help="restrict the other genes drawn to this comma "
                        "separated list; empty keeps them all")
    p.add_argument("--genes-track", action="store_true",
                   help="draw the genes of the window under the map; needs "
                        "--gtf")
    p.add_argument("--gene-types", default="protein_coding",
                   help="gene_type values kept in the track, comma separated; "
                        "empty keeps them all")
    p.add_argument("--gene-rows", type=int, default=3,
                   help="rows available to the gene track")
    p.add_argument("--gene-row-height", type=float, default=0.20,
                   help="height of one gene row in inches")
    p.add_argument("--arcs-max-terms", type=int, default=38,
                   help="number of ENCODE terms at which an arc reaches its "
                        "largest line width")
    p.add_argument("--vmax", type=float, default=0,
                   help="upper end of the colour scale. Lowering it saturates "
                        "the map sooner, so a loop covers more coloured "
                        "pixels. 0 uses the 99th percentile (99.5th for raw "
                        "contacts)")
    p.add_argument("--vmin", type=float, default=None,
                   help="lower end of the colour scale. By default 2 - vmax "
                        "for observed / expected when vmax < 2 (a scale "
                        "centred on 1), else 0; the 20th percentile for raw "
                        "contacts")
    p.add_argument("--marker-size", type=float, default=5.5,
                   help="size of the loop markers, in points")
    p.add_argument("--font-scale", type=float, default=1.0,
                   help="factor applied to every font size")
    p.add_argument("--dpi", type=int, default=300, help="resolution of the PNG")
    p.add_argument("--title", default="",
                   help="title; defaults to the gene and the window")
    p.add_argument("--out", required=True, help="prefix of the output files")
    a = p.parse_args()

    npz = prof_npz = None
    if a.npz:
        npz = np.load(a.npz, allow_pickle=True)
        if not a.gene:
            sys.exit("--npz needs --gene, to pick the viewpoint")
        names = [str(x) for x in pick(npz, "gene names",
                                      "names", "genes", "gene", "viewpoints",
                                      "viewpoint_names", "symbols")]
        if a.gene not in names:
            near = [g for g in names if g.upper().startswith(a.gene[:3].upper())]
            sys.exit(f"{a.gene} absent from {a.npz}; it holds {len(names)} "
                     f"viewpoints" + (f", among them {near[:6]}" if near else ""))
        i = names.index(a.gene)
        chrom = vp = None
        for k in ("chrom", "chroms", "chrom_name", "chromosome"):
            if k in npz.files:
                chrom = str(npz[k][i])
                break
        for k in ("pos", "tss", "viewpoint", "vp_pos", "position"):
            if k in npz.files:
                vp = int(npz[k][i])
                break
        if chrom is None or vp is None:
            # the npz holds the profiles alone; the position comes from
            # loops.tsv, which carries it, or from the annotation
            source = "loops.tsv"
            if a.loops:
                for r in read_tsv(a.loops):
                    if r.get("gene") == a.gene and r.get("chrom"):
                        chrom = chrom or r["chrom"]
                        try:
                            vp = vp if vp is not None else int(float(r["vp_pos"]))
                        except (KeyError, ValueError):
                            pass
                        break
            if (chrom is None or vp is None) and a.gtf:
                chrom, vp, _ = tss_of(a.gtf, a.gene, a.tss_tag)
                source = "the annotation"
            if chrom is None or vp is None:
                sys.exit("the npz holds no viewpoint position; give --loops "
                         "(its vp_pos column) or --gtf")
            print(f"  position taken from {source}")
        prof_npz = np.asarray(
            pick(npz, "profiles", "counts", "profiles", "mat", "matrix",
                 "X", "data")[i], dtype=float)
        npz_off = np.asarray(pick(npz, "offsets", "offsets", "offset",
                                  "distance", "distances"), dtype=float)
        if not np.isfinite(prof_npz).any():
            sys.exit(f"{a.gene} has no profile in {a.npz}; the library was "
                     f"below the contact threshold at that viewpoint")
        a.resolution = int(scalar(npz, "resolution", a.resolution))
        a.anchor_half = int(scalar(npz, "anchor_half", a.anchor_half))
        label = f"{a.gene} TSS"
        print(f"{a.npz}: {scalar(npz, 'accession', '?')}  "
              f"{scalar(npz, 'term', scalar(npz, 'tissue', '?'))}  "
              f"{a.resolution} bp bins")
    elif a.gene:
        if not a.gtf:
            sys.exit("--gene needs --gtf")
        chrom, vp, strand = tss_of(a.gtf, a.gene, a.tss_tag)
        label = f"{a.gene} TSS"
    elif a.chrom and a.viewpoint:
        chrom, vp = a.chrom, a.viewpoint
        label = "viewpoint"
    else:
        sys.exit("give --npz with --gene, --gene with --gtf, or --chrom with "
                 "--viewpoint")

    hic = a.hic
    if not hic and a.tissues and a.term:
        hic = hic_url(a.tissues, a.term, a.term_column,
                      tuple(a.location_columns.split(',')))
        print(f"  {a.term}: {hic}")
    elif a.term and not a.tissues:
        sys.exit("--term needs --tissues")

    res = a.resolution
    start = max(0, (vp - a.window) // res * res)
    end = (vp + a.window) // res * res + res
    print(f"{label}  {chrom}:{vp:,}  window {start:,}-{end:,}  {res} bp bins")

    m = None
    if hic:
        reader = open_reader(hic, a.norm, res)
        m = reader.square(chrom, start, end)
        res = getattr(reader, "res", res)
        n = len(m)
        edges = start + np.arange(n + 1) * res
    else:
        n = int((end - start) // res)
        edges = start + np.arange(n + 1) * res
        print("  no contact map given, drawing the virtual 4C alone")
    centres = edges[:-1] + res / 2

    if prof_npz is not None:
        # the stored profile is indexed by offset to the viewpoint
        prof = np.full(n, np.nan)
        for k, off in enumerate(npz_off):
            j = int((vp + off - start) // res)   # the bin holding that offset
            if 0 <= j < n:
                prof[j] = prof_npz[k]
    elif m is not None:
        # the viewpoint row, summed over the anchor bins, is the virtual 4C
        lo = max(0, int((vp - a.anchor_half - start) // res))
        hi = min(n, int((vp + a.anchor_half - start) // res) + 1)
        if hi <= lo:
            sys.exit("the viewpoint falls outside the window")
        prof = np.nansum(m[lo:hi, :], axis=0)
    else:
        sys.exit("give --npz or a contact map; there is nothing to draw")
    lo = max(0, int((vp - a.anchor_half - start) // res))
    hi = min(n, int((vp + a.anchor_half - start) // res) + 1)
    prof[lo:hi] = np.nan                   # the anchor itself is not a contact

    oe = not a.raw
    drawn = observed_over_expected(m) if (oe and m is not None) else m

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
    from matplotlib.transforms import offset_copy

    def fs(size):
        return size * a.font_scale

    # text kept as text in the SVG, so that svg_to_pptx.py rebuilds it as
    # editable text boxes
    plt.rcParams.update({"svg.fonttype": "none",
                         "font.size": fs(10.0)})
    fall = LinearSegmentedColormap.from_list("fall", RAW_RAMP)
    diverging = LinearSegmentedColormap.from_list("oe", OE_RAMP)

    # a bin is a square rotated by 45 degrees, so the two axes must carry the
    # same number of bp per inch; the panel heights follow from that.
    # `top` is the largest separation shown. A bin sits at half its own
    # separation on the rotated axis, so the panel is half as tall as that.
    top = a.max_distance or max(res * 8, a.window)
    half = top / 2
    fig_w = a.width
    ax_w = fig_w * a.map_fraction           # what is left beside labels and bar
    map_h = ax_w * half / (edges[-1] - edges[0])
    prof_h = a.profile_height

    track = []
    if a.genes_track:
        if not a.gtf:
            sys.exit("--genes-track needs --gtf")
        track = genes_in(a.gtf, chrom, start, end, a.gene_types)
        print(f"  {len(track)} genes in the window")
    gene_rows = pack_rows(track, end - start, a.gene_rows) if track else []
    n_rows = (max(r for r in gene_rows if r is not None) + 1) \
        if any(r is not None for r in gene_rows) else 0
    gene_h = a.gene_row_height * n_rows + 0.08 if n_rows else 0

    arcs = []
    if a.arcs:
        if not a.loops:
            sys.exit("--arcs needs --loops")
        only = {x.strip() for x in a.arcs_others.split(",") if x.strip()}
        for r in read_tsv(a.loops):
            try:
                v_, s_ = float(r["vp_pos"]), float(r["summit"])
                nt = int(float(r.get("n_terms", 1) or 1))
            except (KeyError, ValueError):
                continue
            if not same_chrom(r.get("chrom", ""), chrom):
                continue
            if not (start <= v_ <= end and start <= s_ <= end):
                continue
            if nt < a.arcs_min_terms:
                continue
            g_ = r.get("gene", "")
            if g_ != a.gene and only and g_ not in only:
                continue
            arcs.append((g_, v_, s_, nt))
        others = sorted({g for g, *_ in arcs} - {a.gene})
        print(f"  {len(arcs)} loops in the window, "
              f"{sum(1 for g, *_ in arcs if g == a.gene)} at {a.gene}"
              + (f", the rest at {', '.join(others)}" if others else ""))
    arcs_h = a.arcs_height if arcs else 0

    gx = None
    ax2 = None
    if m is None and not n_rows and not arcs:
        fig, bx = plt.subplots(1, 1, figsize=(fig_w, prof_h + 0.95))
        ax = None
    else:
        heights = ([map_h] if m is not None else []) + \
                  ([arcs_h] if arcs else []) + \
                  ([gene_h] if n_rows else []) + [prof_h]
        fig, axes = plt.subplots(
            len(heights), 1, figsize=(fig_w, sum(heights) + 0.95), sharex=True,
            gridspec_kw={"height_ratios": heights, "hspace": 0.14})
        axes = list(np.atleast_1d(axes))
        ax = axes.pop(0) if m is not None else None
        ax2 = axes.pop(0) if arcs else None
        gx = axes.pop(0) if n_rows else None
        bx = axes.pop(0)

    v = drawn[np.isfinite(drawn) & (drawn > 0)] if m is not None \
        else np.array([])
    if oe:
        hi_v = a.vmax or (np.nanpercentile(v, 99) if v.size else 2.0)
        lo_v = a.vmin if a.vmin is not None else \
            (2 - hi_v if hi_v < 2 else 0)
        norm = Normalize(vmin=lo_v, vmax=hi_v)
        cmap = diverging
    else:
        norm = LogNorm(vmin=a.vmin if a.vmin is not None
                       else max(np.nanpercentile(v, 20), 1e-3),
                       vmax=a.vmax or np.nanpercentile(v, 99.5)) \
            if v.size else None
        cmap = fall

    if ax is None:
        pass
    elif a.square:
        im = ax.imshow(drawn, origin="lower", cmap=cmap, norm=norm,
                       extent=[edges[0], edges[-1], edges[0], edges[-1]],
                       aspect="equal", interpolation="nearest")
        ax.set_ylabel("position")
    else:
        K, L = np.meshgrid(edges, edges, indexing="ij")
        X = (K + L) / 2
        Y = (L - K) / 2
        M = np.array(drawn, dtype=float)
        M[np.tril_indices(n, -1)] = np.nan     # one half is enough
        im = ax.pcolormesh(X, Y, np.ma.masked_invalid(M),
                           cmap=cmap, norm=norm, shading="flat",
                           rasterized=True)
        ax.set_ylim(0, half)
        ax.set_aspect("equal", adjustable="box")
        ax.set_ylabel("separation (kb)")
        ax.set_yticks([0, half / 2, half])
        ax.set_yticklabels(["0", f"{top/2000:.0f}", f"{top/1000:.0f}"],
                           fontsize=fs(7))
    head = f"{a.gene}   " if a.gene else ""
    title = a.title or f"{head}{chrom}:{start:,}-{end:,}"
    if ax is not None:
        ax.set_xlim(edges[0], edges[-1])
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.annotate("", xy=(vp, 0), xytext=(0, -a.arrow),
                    textcoords="offset points", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-|>", color="black", lw=1.2,
                                    shrinkA=0, shrinkB=0, mutation_scale=10))
        ax.set_title(title, fontsize=fs(9), loc="left")
        cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.025)
        cb.ax.tick_params(labelsize=fs(7))
        cb.set_label("observed / expected" if oe else
                     ("balanced contacts" if a.norm.upper() != "NONE"
                      else "contacts"), fontsize=fs(7))
    else:
        bx.set_title(title, fontsize=fs(9), loc="left")

    if ax2 is not None:
        from matplotlib.path import Path as MPath
        from matplotlib.patches import PathPatch
        widest = max(abs(s_ - v_) for _, v_, s_, _ in arcs) or 1.0
        ax2.set_xlim(edges[0], edges[-1])
        ax2.set_ylim(0, 1.06)
        ax2.set_yticks([])
        ax2.set_ylabel("loops", fontsize=fs(7), rotation=0, ha="right",
                       va="center", labelpad=2)
        for s in ("top", "right", "left", "bottom"):
            ax2.spines[s].set_visible(False)
        ax2.tick_params(axis="x", length=0)
        # the widest arcs first, so the short ones stay readable on top
        for g_, v_, s_, nt in sorted(arcs, key=lambda x: -abs(x[2] - x[1])):
            mine = (g_ == a.gene)
            h = (abs(s_ - v_) / widest) ** 0.55
            col = GOLD if mine else GREY
            path = MPath([(v_, 0), (v_, h * 1.34), (s_, h * 1.34), (s_, 0)],
                         [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4,
                          MPath.CURVE4])
            ax2.add_patch(PathPatch(
                path, fill=False, edgecolor=col,
                lw=0.6 + 1.9 * min(nt, a.arcs_max_terms) / a.arcs_max_terms,
                alpha=1.0 if mine else 0.55,
                zorder=3 if mine else 2))
            ax2.plot([v_, s_], [0, 0], ls="none", marker="o", ms=2.4,
                     color=col, zorder=4 if mine else 2,
                     alpha=1.0 if mine else 0.55)

    if gx is not None:
        gx.set_xlim(edges[0], edges[-1])
        gx.set_ylim(n_rows, 0)
        gx.set_yticks([])
        gx.set_ylabel("genes", fontsize=fs(7), rotation=0, ha="right",
                      va="center", labelpad=2)
        for s in ("top", "right", "left", "bottom"):
            gx.spines[s].set_visible(False)
        gx.tick_params(axis="x", length=0)
        gx.plot([vp, vp], [0, n_rows], color="black", lw=0.7,
                ls=(0, (1, 2)), zorder=0)
        for (name, s_, e_, strand), r in zip(track, gene_rows):
            if r is None:
                continue
            y = r + 0.5
            is_vp = (name == a.gene)
            col = "black" if is_vp else GREY
            gx.plot([max(s_, start), min(e_, end)], [y, y],
                    color=col, lw=2.4 if is_vp else 1.6,
                    solid_capstyle="butt")
            tip = min(e_, end) if strand == "+" else max(s_, start)
            gx.plot([tip], [y], marker=">" if strand == "+" else "<",
                    ms=3.6, color=col)
            gx.annotate(name, xy=(min(e_, end) if strand == "+"
                                  else max(s_, start), y),
                        xytext=(5 if strand == "+" else -5, 0),
                        textcoords="offset points", fontsize=fs(6),
                        style="italic", color=col, va="center",
                        ha="left" if strand == "+" else "right")

    bx.fill_between(centres, 0, np.nan_to_num(prof), step="mid",
                    color=RED, alpha=0.25, linewidth=0)
    bx.step(centres, prof, where="mid", color=RED, lw=1.0)
    bx.axvline(vp, color=BLUE, lw=0.9, ls=(0, (3, 2)))
    bx.set_ylabel("virtual 4C")
    bx.set_xlabel(f"{chrom}  (Mb)")
    bx.set_xlim(edges[0], edges[-1])
    bx.set_ylim(bottom=0)
    for s in ("top", "right"):
        bx.spines[s].set_visible(False)
    bx.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, _: f"{x/1e6:.2f}"))

    n_peaks = 0
    if a.loops:
        ymax = np.nanmax(prof) if np.isfinite(prof).any() else 1.0
        for r in read_tsv(a.loops):
            if a.gene and r.get("gene") != a.gene:
                continue
            if not same_chrom(r.get("chrom", chrom), chrom):
                continue
            try:
                s_ = float(r["summit"])
            except (KeyError, ValueError):
                continue
            if not (start <= s_ <= end):
                continue
            n_peaks += 1
            bx.plot([s_], [ymax * 1.06], marker="v", ms=a.marker_size,
                    color=GOLD, clip_on=False)
            if ax is not None and not a.square:
                # the loop is the pair (viewpoint, summit); on the rotated
                # matrix it sits at their midpoint, at half their separation.
                # The same marker as the panel below, lifted a few points so
                # that it points at the pixel without covering it.
                ax.plot([(vp + s_) / 2], [abs(s_ - vp) / 2], marker="v",
                        ms=a.marker_size, color=GOLD, clip_on=False,
                        transform=offset_copy(ax.transData, fig,
                                              y=8, units="points"))
            elif ax is not None:
                ax.plot([s_], [vp], marker="v", ms=a.marker_size,
                        color=GOLD)
        bx.set_ylim(0, ymax * 1.16)
        print(f"  {n_peaks} called loops inside the window")

    bx.annotate(label, xy=(vp, bx.get_ylim()[1]), xytext=(4, -2),
                textcoords="offset points", ha="left", va="top",
                fontsize=fs(7), color=BLUE)

    # `set_aspect("equal")` shrinks the map's box to satisfy the aspect, while
    # sharex only shares the data limits: the other panels keep the full
    # width. Giving each panel the map's horizontal extent once the layout is
    # computed puts one genomic coordinate at the same x in every panel.
    fig.canvas.draw()
    if ax is not None:
        box = ax.get_position()
        for other in (ax2, gx, bx):
            if other is not None:
                q = other.get_position()
                other.set_position([box.x0, q.y0, box.width, q.height])

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=a.dpi, bbox_inches="tight")
    with open(out.with_suffix(".tsv"), "w") as fh:
        fh.write("chrom\tstart\tend\tdistance\tvirtual_4c\n")
        for i, c in enumerate(centres):
            fh.write(f"{chrom}\t{edges[i]:.0f}\t{edges[i+1]:.0f}\t"
                     f"{c - vp:.0f}\t{prof[i]:.4f}\n")
    print(f"  written {out.with_suffix('.svg')}, .png and .tsv")


if __name__ == "__main__":
    main()
