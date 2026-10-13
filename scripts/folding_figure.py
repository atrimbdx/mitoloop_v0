#!/usr/bin/env python3
"""
folding_figure.py - schematic chromatin folding of a locus, from loops.tsv

Purpose
    Draws one chromatin conformation consistent with the loops called in a
    window: the anchors linked by a contact are brought together into a
    rosette, the chromatin between two consecutive anchors of the rosette
    forms a petal, and the fibre outside the rosette meanders. The drawing is
    a schematic of the called contacts; it is not a resolved structure.

Inputs
    --loops      loops.tsv from `v4c_loops.R reciprocal`: gene, chrom,
                 vp_pos, summit, n_terms, called (`|`-separated ENCODE terms)
    --elements   optional re2g_elements.tsv from `re2g_panel.py`: chrom,
                 start, end, gene, prevalence
    --gtf        optional GTF (gene records with gene_name, gene_type), to
                 place on the fibre the genes of the window that carry no loop

Outputs (prefix --out)
    <out>.svg, <out>.png   the figure (SVG text kept as text)
    <out>.tsv              one row per petal: petal, anchor_a, anchor_b,
                           position_a, position_b, length, kind_a, kind_b
                           (kind: prom, enh, other, gene)

Method
    Model. Anchors linked by a contact are placed side by side on a short
    segment and form a rosette. The chromatin between two consecutive anchors
    of the rosette forms a petal whose size grows with the square root of the
    genomic length of the segment. Petals open upwards as a fan and the two
    tails leave downwards, so that nothing crosses. Outside the petals the
    fibre meanders: the data say nothing about how it folds there, and the
    drawing assumes nothing. Petals are cubic Bezier curves, tails are
    smoothed noise, and the whole polyline receives a final Gaussian
    smoothing, so the trace has no corners. Each --seed gives a different
    conformation compatible with the same constraints.

    Anchors. They are read from loops.tsv, and all the loops of the window
    are taken, not only those of the requested gene: a neighbouring viewpoint
    with its own loops enters the rosette in the same way. Positions closer
    than --merge are merged into one anchor. An anchor carrying the viewpoint
    of a gene takes the name of that gene and is drawn as a promoter.
    Otherwise, with --elements, an anchor overlapping an rE2G element of
    prevalence >= --min-prevalence is drawn as a regulatory element. Other
    anchors are named by their position relative to the requested gene.

    Labels. Labels start at a set distance from their anchor and are then
    relaxed under three forces: a pull back to the starting point, repulsion
    between labels, and repulsion from the nearest points of the fibre. The
    anchors of the rosette touch one another, so their labels would overlap
    at their starting points.

    Scope. The virtual 4C looks at one row of the matrix per viewpoint. A
    loop between two positions of which neither is a viewpoint of the panel
    is not in loops.tsv and does not appear here.

Usage
    python3 folding_figure.py --loops results/loops/loops.tsv --gene GENE1 \\
        --window 1000000 --elements results/re2g/re2g_elements.tsv \\
        --called-in "psoas muscle|gastrocnemius" --out figures/folding_GENE1
    python3 folding_figure.py --loops results/loops/loops.tsv --gene GENE1 \\
        --gtf data/gencode.gtf.gz --lang fr --seed 5 --out figures/folding_GENE1

References
    Gschwind AR et al. 2023. An encyclopedia of enhancer-gene regulatory
        interactions in the human genome. bioRxiv. (ENCODE-rE2G elements)
    Frankish A et al. 2023. GENCODE: reference annotation for the human and
        mouse genomes in 2023. Nucleic Acids Research. doi:10.1093/nar/gkac1071
    Virtanen P et al. 2020. SciPy 1.0: fundamental algorithms for scientific
        computing in Python. Nature Methods.
    Hunter JD. 2007. Matplotlib: a 2D graphics environment. Computing in
        Science and Engineering.
"""
import argparse
import gzip
import re
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RED, BLUE, GOLD, GREY = "#B2182B", "#2166AC", "#C78A20", "#5A5A5A"
OLIVE = "#7A8B2A"
PROMOTER_COLOUR, GENE_COLOUR = "#1A1A1A", "#9A9690"
LEADER_COLOUR = "#BDB8B4"

# figure labels; the French set is a label set only
LABELS = {
    "fr": dict(
        title="{g} — repliement de la chromatine, une conformation "
              "compatible avec les boucles mesurées",
        subtitle="{c}:{a}-{b}   {n} boucles, {m} ancres",
        element="élément régulateur réplicable", promoter="promoteur",
        other="ancre de boucle", gene="autre gène du locus",
        contact="contact entre deux ancres consécutives"),
    "en": dict(
        title="{g} — chromatin folding: one conformation consistent with "
              "the called loops",
        subtitle="{c}:{a}-{b}   {n} loops, {m} anchors",
        element="replicable regulatory element", promoter="promoter",
        other="loop anchor", gene="other gene of the locus",
        contact="contact between two consecutive anchors"),
}


# ==========================================================================
# Rosette and petals
# ==========================================================================
DELTA = 0.62            # spacing between two anchors in contact
LEFT_TAIL_ANGLE = 188.0     # direction of the left tail, in degrees
RIGHT_TAIL_ANGLE = -12.0    # direction of the right tail
FAN_HIGH, FAN_LOW = 142.0, 52.0   # bounds of the fan of petals
MAX_TAIL_LENGTH = 6.0


def _unit(deg):
    a = np.deg2rad(deg)
    return np.array([np.cos(a), np.sin(a)])


def _wave(C, amp, seed, n_waves=2.2):
    """Gentle perpendicular undulation that leaves the end points in place."""
    if amp <= 0 or len(C) < 8:
        return C
    rng = np.random.default_rng(seed)
    d = np.gradient(C, axis=0)
    L = np.linalg.norm(d, axis=1, keepdims=True) + 1e-9
    perp = np.column_stack([-d[:, 1], d[:, 0]]) / L
    w = gaussian_filter1d(rng.normal(0, 1, len(C)),
                          sigma=len(C) / (6 * n_waves), mode="nearest")
    t = np.linspace(0, 1, len(C))
    w = w - w[0] - t * (w[-1] - w[0])
    return C + perp * (w / (np.abs(w).max() + 1e-9) * amp)[:, None]


def petal(a, b, angle, size, seed, n=280):
    """Chromatin loop between two anchors in contact. The control points are
    pushed far along `angle` and spread on either side, which gives a full
    loop with a rounded end."""
    u = _unit(angle)
    perp = np.array([-u[1], u[0]])
    spread = size * 0.45
    c1 = a + u * size + perp * spread
    c2 = b + u * size - perp * spread
    t = np.linspace(0, 1, n)[:, None]
    C = ((1 - t) ** 3 * a + 3 * (1 - t) ** 2 * t * c1
         + 3 * (1 - t) * t ** 2 * c2 + t ** 3 * b)
    return _wave(C, size * 0.030, seed)


def tail(origin, angle, length, seed, n=220):
    """Chromatin outside the loops: a meandering line, without folding."""
    t = np.linspace(0, 1, n)[:, None]
    return _wave(origin + t * _unit(angle) * length, length * 0.080, seed)


def _tail_geometry(tail_anchors, rosette_anchor, centres, spacing,
                   margin=0.40):
    """Length of the tail and fractional position of each anchor on it.

    Positions are fractions of the actual (capped) tail length, so that a
    cap on the length keeps the anchors spread along the tail.
    """
    if not tail_anchors:
        return spacing * 0.8, []
    chain = [rosette_anchor] + list(tail_anchors)
    steps = np.array([abs(centres[chain[k + 1]] - centres[chain[k]]) / 1e6
                      for k in range(len(chain) - 1)]) ** 0.5
    cum = np.cumsum(steps)
    total = cum[-1] + margin
    return min(spacing * total, MAX_TAIL_LENGTH), list(cum / total)


def build(names, centres, contacts, seed=3, tail_spacing=3.1,
          petal_base=1.30, petal_slope=3.80):
    """Assemble the trace.

    `contacts` lists the retained contacts as pairs of anchors. Their union,
    in genomic order, forms the rosette; a petal is drawn between two
    consecutive members. If `contacts` is empty, the chromatin is a single
    line carrying every anchor.

    Returns the polyline, the index of each anchor on it, the base and
    direction of each petal, and the members of the rosette.
    """
    rank = {n: k for k, n in enumerate(names)}

    if not contacts:
        d = np.array([abs(centres[names[k + 1]] - centres[names[k]]) / 1e6
                      for k in range(len(names) - 1)]) ** 0.5
        cum = np.concatenate([[0.0], np.cumsum(d)])
        total = cum[-1] + 0.9
        L = min(tail_spacing * total, 2 * MAX_TAIL_LENGTH)
        C = tail(np.array([-L / 2.0, 0.0]), 0.0, L, seed + 1, n=420)
        C = gaussian_filter1d(C, sigma=2.0, axis=0, mode="nearest")
        idx = [int(np.clip((0.45 + c) / total, 0.0, 1.0) * (len(C) - 1))
               for c in cum]
        return C, idx, [], []

    members = sorted({x for p in contacts for x in p}, key=lambda n: rank[n])
    m = len(members)
    pos = {n: np.array([(k - (m - 1) / 2.0) * DELTA, 0.0])
           for k, n in enumerate(members)}

    n_petals = m - 1
    angles = (np.linspace(FAN_HIGH, FAN_LOW, n_petals) if n_petals > 1
              else np.array([95.0]))

    before = [n for n in names if rank[n] < rank[members[0]]][::-1]
    after = [n for n in names if rank[n] > rank[members[-1]]]

    pieces, marks, bases = [], {}, []

    # --- left tail, traced from the rosette outwards then reversed
    Lg, fg = _tail_geometry(before, members[0], centres, tail_spacing)
    Qg = tail(pos[members[0]], LEFT_TAIL_ANGLE, Lg, seed + 1)[::-1]
    pieces.append(Qg[:-1])
    for n, f in zip(before, fg):
        marks[n] = int((1.0 - min(f, 1.0)) * (len(Qg) - 1))

    # --- petals ----------------------------------------------------------
    for k in range(n_petals):
        a, b = members[k], members[k + 1]
        L = max(abs(centres[b] - centres[a]), 1.0)
        size = petal_base + petal_slope * (L / 1e6) ** 0.5
        marks[a] = sum(len(x) for x in pieces)
        P = petal(pos[a], pos[b], angles[k], size, seed + 10 + k)
        bases.append(((pos[a] + pos[b]) / 2.0, angles[k]))
        base_i = sum(len(x) for x in pieces)
        pieces.append(P[:-1])
        for n in [x for x in names if rank[a] < rank[x] < rank[b]]:
            f = (centres[n] - centres[a]) / (centres[b] - centres[a])
            marks[n] = base_i + int(np.clip(f, 0.04, 0.96) * (len(P) - 1))

    # --- right tail ------------------------------------------------------
    marks[members[-1]] = sum(len(x) for x in pieces)
    Ld, fd = _tail_geometry(after, members[-1], centres, tail_spacing)
    base_i = sum(len(x) for x in pieces)
    Qd = tail(pos[members[-1]], RIGHT_TAIL_ANGLE, Ld, seed + 2)
    pieces.append(Qd)
    for n, f in zip(after, fd):
        marks[n] = base_i + int(min(f, 1.0) * (len(Qd) - 1))

    C = np.vstack(pieces)
    C = gaussian_filter1d(C, sigma=2.0, axis=0, mode="nearest")
    idx = [int(np.clip(marks.get(n, 0), 0, len(C) - 1)) for n in names]
    return C, idx, bases, members


# ==========================================================================
# Label placement
# ==========================================================================
def place_labels(anchors, C, radius=1.9, pull=0.10, n_iter=900,
                 step=0.06, direction=None):
    """Move labels away from one another and from the fibre.

    Start: at distance `radius` from the anchor, along `direction` when
    given, otherwise away from the centre of the drawing. Then relaxation
    under three forces: a pull back to that starting point, repulsion
    between labels, and repulsion from the nearest points of the fibre.
    """
    anchors = np.asarray(anchors, dtype=float)
    M = len(anchors)
    ctr = C.mean(0)
    u = (anchors - ctr if direction is None
         else np.asarray(direction, dtype=float).reshape(M, 2))
    nu = np.linalg.norm(u, axis=1, keepdims=True)
    u = np.where(nu > 1e-6, u / np.maximum(nu, 1e-9), np.array([[0.0, 1.0]]))

    radius = np.broadcast_to(np.asarray(radius, dtype=float).reshape(-1, 1),
                             (M, 1)).copy()
    pull = np.broadcast_to(np.asarray(pull, dtype=float).reshape(-1, 1),
                           (M, 1)).copy()
    start = anchors + u * radius

    # Without their own direction, the anchors of one contact, which sit at
    # the same place, would all send their label to the same point: each
    # cluster is spread as a fan instead.
    if direction is None:
        D = np.linalg.norm(anchors[:, None] - anchors[None, :], axis=2)
        left = set(range(M))
        while left:
            k = left.pop()
            cluster, stack = [k], [k]
            while stack:
                i = stack.pop()
                for j in list(left):
                    if D[i, j] < 0.9:
                        left.discard(j)
                        cluster.append(j)
                        stack.append(j)
            if len(cluster) < 2:
                continue
            cluster.sort()
            ctr_c = anchors[cluster].mean(0)
            base = np.arctan2(*(ctr_c - ctr)[::-1])
            for t, i in enumerate(cluster):
                r = float(radius[i, 0]) * (1.0 + 0.07 * len(cluster))
                ang = base + (t - (len(cluster) - 1) / 2.0) * (
                    2.0 * np.pi / max(len(cluster), 3))
                start[i] = ctr_c + r * np.array([np.cos(ang), np.sin(ang)])

    q = start.copy()
    R_LBL = 1.9          # minimal distance between two labels
    R_FIB = 0.70         # minimal distance to the fibre
    F = C[:: max(len(C) // 400, 1)]

    for _ in range(n_iter):
        g = pull * (q - start)

        D = q[:, None] - q[None, :]
        d = np.linalg.norm(D, axis=2) + 1e-9
        np.fill_diagonal(d, 1e9)
        m = d < R_LBL
        if m.any():
            g -= 1.0 * np.where(m[..., None],
                                D / d[..., None] * (R_LBL - d[..., None]),
                                0).sum(1)

        DF = q[:, None] - F[None, :]
        df = np.linalg.norm(DF, axis=2) + 1e-9
        mf = df < R_FIB
        if mf.any():
            g -= 0.5 * np.where(mf[..., None],
                                DF / df[..., None] * (R_FIB - df[..., None]),
                                0).sum(1)

        gn = np.linalg.norm(g, axis=1, keepdims=True)
        g = np.where(gn > 1.0, g / np.maximum(gn, 1e-9), g)
        q = q - step * g
    return q


# ==========================================================================
# Data
# ==========================================================================
def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.strip():
                f = line.rstrip("\n").split("\t")
                yield dict(zip(head, f + [""] * (len(head) - len(f))))


def same_chrom(c, ref):
    return c.lstrip("chr") == ref.lstrip("chr")


def load_loops(path, chrom, start, end, min_terms, called_in):
    """The loops whose two anchors lie in the window."""
    rx = re.compile(called_in, re.I) if called_in else None
    out = []
    for r in read_tsv(path):
        if not same_chrom(r.get("chrom", ""), chrom):
            continue
        try:
            v, s = float(r["vp_pos"]), float(r["summit"])
            nt = int(float(r.get("n_terms", 1) or 1))
        except (KeyError, ValueError):
            continue
        if not (start <= v <= end and start <= s <= end) or nt < min_terms:
            continue
        if rx:
            terms = [t for t in r.get("called", "").split("|") if t]
            if terms and not any(rx.search(t) for t in terms):
                continue
        out.append(dict(gene=r.get("gene", ""), vp=v, summit=s, n_terms=nt))
    return out


def merge_positions(positions, tol):
    """Group close positions. Returns, for each position, its group centre."""
    order = sorted(set(positions))
    groups, cur = [], [order[0]]
    for p in order[1:]:
        if p - cur[-1] <= tol:
            cur.append(p)
        else:
            groups.append(cur)
            cur = [p]
    groups.append(cur)
    return {p: float(np.mean(g)) for g in groups for p in g}


def load_elements(path, chrom):
    els = []
    for r in read_tsv(path):
        try:
            if not same_chrom(r.get("chrom", ""), chrom):
                continue
            els.append((float(r["start"]), float(r["end"]), r.get("gene", ""),
                        int(float(r.get("prevalence", 0)))))
        except (KeyError, ValueError):
            continue
    return els


def load_genes(gtf, chrom, start, end, types, keep):
    """TSS of the genes of the window. They carry no contact: they are put
    on the fibre as landmarks, at their genomic place."""
    want = {x.strip() for x in types.split(",") if x.strip()}
    only = {x.strip() for x in keep.split(",") if x.strip()}
    out = {}
    op = gzip.open if str(gtf).endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f_ = line.split("\t")
            if len(f_) < 9 or f_[2] != "gene" or not same_chrom(f_[0], chrom):
                continue
            s, e = int(f_[3]) - 1, int(f_[4])
            tss = s if f_[6] == "+" else e
            if not (start <= tss <= end):
                continue
            a_ = f_[8]
            i = a_.find('gene_name "')
            if i < 0:
                continue
            name = a_[i + 11:a_.find('"', i + 11)]
            if only and name not in only:
                continue
            if want:
                j = a_.find('gene_type "')
                if (a_[j + 11:a_.find('"', j + 11)] if j >= 0 else "") \
                        not in want:
                    continue
            out[name] = float(tss)
    return out


def build_anchors(loops, tol, gene_vp, elements, min_prev, landmarks=None):
    """Merge positions into anchors, name them and classify them."""
    centre = merge_positions([b["vp"] for b in loops] +
                             [b["summit"] for b in loops], tol)

    promoters = {}                       # centre -> gene
    for b in loops:
        promoters.setdefault(centre[b["vp"]], b["gene"])

    names, info = {}, {}
    for c in sorted(set(centre.values())):
        carried = [e for e in elements
                   if e[0] - tol <= c <= e[1] + tol and e[3] >= min_prev]
        if c in promoters:
            name, kind = promoters[c], "prom"
        else:
            name = f"{(c - gene_vp) / 1000.0:+.0f} kb"
            kind = "enh" if carried else "other"
        while name in names.values():    # two anchors never share a name
            name += "'"
        names[c] = name
        info[name] = dict(centre=c, kind=kind, elements=len(carried))

    # Genes of the locus that carry no loop are added as landmarks at their
    # genomic place, unless they fall on an anchor already placed, which
    # already carries their name or that of its element.
    placed = sorted(names)
    for g, pos in sorted((landmarks or {}).items(), key=lambda x: x[1]):
        if g in info:
            continue
        if any(abs(pos - c) <= tol for c in placed):
            continue
        info[g] = dict(centre=pos, kind="gene", elements=0)
    return centre, names, info


# ==========================================================================
# Drawing
# ==========================================================================
def folding_panel(ax, names, centres, kinds, contacts, seed,
                  fibre_colour=GREY, fibre_width=2.5, label_size=6.6,
                  number_size=6.2, geometry=None):
    C, idx, bases, members = build(names, centres, contacts, seed=seed,
                                   **(geometry or {}))
    ax.plot(C[:, 0], C[:, 1], color=fibre_colour, lw=fibre_width, zorder=2,
            solid_capstyle="round", solid_joinstyle="round")

    anchors = np.array([C[i] for i in idx])
    mids = (np.array([b + _unit(a) * 0.85 for b, a in bases])
            if bases else np.zeros((0, 2)))

    # Starting direction of each label. An anchor of the rosette points down,
    # fanned from left to right, since the petals occupy the top; an anchor
    # on a tail points perpendicular to the fibre, downwards and away from
    # the centre of the drawing.
    ctr_r = (np.mean([C[idx[names.index(n)]] for n in members], axis=0)
             if members else C.mean(0))
    ang_ros = (np.linspace(212.0, -32.0, len(members)) if len(members) > 1
               else np.array([270.0]))
    dirs = []
    for k, n in enumerate(names):
        if n in members:
            a = np.deg2rad(ang_ros[members.index(n)])
            dirs.append(np.array([np.cos(a), -abs(np.sin(a))]))
        else:
            i = idx[k]
            t = C[min(i + 6, len(C) - 1)] - C[max(i - 6, 0)]
            t = t / (np.linalg.norm(t) + 1e-9)
            perp = np.array([-t[1], t[0]])
            if perp[1] > 0:
                perp = -perp
            lat = np.sign(anchors[k, 0] - ctr_r[0]) or 1.0
            dirs.append(perp + np.array([lat * 0.95, 0.0]))
    dirs = np.array(dirs)

    # Petal numbers stay at the base of their petal; they are added to the
    # fibre as obstacles, so that labels avoid them.
    obstacle = np.vstack([C, mids]) if len(mids) else C
    lbl = place_labels(anchors, obstacle, radius=1.55, pull=0.42,
                       direction=dirs)

    for k, n in enumerate(names):
        i = idx[k]
        d = C[min(i + 5, len(C) - 1)] - C[max(i - 5, 0)]
        d = d / (np.linalg.norm(d) + 1e-9)
        pe = np.array([-d[1], d[0]])
        if kinds[n] in ("prom", "gene"):
            col = PROMOTER_COLOUR if kinds[n] == "prom" else GENE_COLOUR
            ax.annotate("", xy=C[i] + d * 0.55 + pe * 0.22,
                        xytext=C[i] + pe * 0.22,
                        arrowprops=dict(arrowstyle="-|>", color=col,
                                        lw=1.6 if kinds[n] == "prom" else 1.1,
                                        mutation_scale=10), zorder=9)
        else:
            ax.plot(*np.array([C[i] - pe * 0.19, C[i] + pe * 0.19]).T,
                    color=OLIVE if kinds[n] == "enh" else GREY,
                    lw=3.2, solid_capstyle="round", zorder=9)
        ax.plot([anchors[k, 0], lbl[k, 0]], [anchors[k, 1], lbl[k, 1]],
                color=LEADER_COLOUR, lw=0.5, zorder=4)
        ax.annotate(n, lbl[k], ha="center", va="center", fontsize=label_size,
                    color="#8A8680" if kinds[n] == "gene" else "#2A2A2A",
                    style="italic" if kinds[n] in ("prom", "gene") else "normal",
                    zorder=11, linespacing=1.15,
                    bbox=dict(boxstyle="round,pad=0.14", fc="white",
                              ec="none", alpha=0.92))

    for k in range(len(mids)):
        ax.annotate(str(k + 1), xy=mids[k], ha="center", va="center",
                    fontsize=number_size, color="white",
                    fontweight="bold", zorder=12,
                    bbox=dict(boxstyle="circle,pad=0.18", fc=GOLD,
                              ec="white", lw=0.8))

    ax.set_aspect("equal")
    ax.axis("off")
    return (np.vstack([C, lbl, mids]) if len(mids) else np.vstack([C, lbl]),
            members)


# ==========================================================================
def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--loops", required=True,
                   help="loops.tsv from `v4c_loops.R reciprocal`")
    p.add_argument("--gene", required=True,
                   help="gene of the locus; centres the window and is the "
                        "origin of relative positions")
    p.add_argument("--window", type=int, default=1_000_000,
                   help="half-width of the window around the gene viewpoint, "
                        "in bp")
    p.add_argument("--elements",
                   help="re2g_elements.tsv, to mark the anchors that carry a "
                        "replicable element")
    p.add_argument("--min-prevalence", type=int, default=51,
                   help="minimal prevalence (column `prevalence`, a number "
                        "of biosamples) of an element for its anchor to be "
                        "drawn as regulatory; set it relative to the number "
                        "of biosamples kept by re2g_panel.py --keep")
    p.add_argument("--min-terms", type=int, default=3,
                   help="keep the loops called in at least this many ENCODE "
                        "terms")
    p.add_argument("--called-in", default="",
                   help="regular expression on ENCODE terms; a loop is kept "
                        "when one of its `called` terms matches")
    p.add_argument("--merge", type=int, default=20_000,
                   help="two anchor positions closer than this, in bp, are "
                        "merged")
    p.add_argument("--only-gene", action="store_true",
                   help="keep only the loops of the requested gene")
    p.add_argument("--max-anchors", type=int, default=0,
                   help="when the rosette is crowded, keep only the most "
                        "reproduced loops up to this number of anchors; 0 "
                        "keeps all")
    p.add_argument("--gtf",
                   help="annotation, to place on the fibre the genes of the "
                        "locus that carry no loop")
    p.add_argument("--gene-types", default="protein_coding",
                   help="gene_type values kept, comma separated; empty keeps "
                        "them all")
    p.add_argument("--gene-list", default="",
                   help="comma separated genes, to place only these")
    p.add_argument("--seed", type=int, default=3,
                   help="seed of the meandering and of the petal undulation")
    p.add_argument("--petal-base", type=float, default=1.30,
                   help="size of a petal of zero genomic length, in drawing "
                        "units")
    p.add_argument("--petal-slope", type=float, default=3.80,
                   help="growth of the petal size per square root of Mb")
    p.add_argument("--tail-spacing", type=float, default=3.1,
                   help="tail length per square root of Mb between anchors, "
                        "in drawing units")
    p.add_argument("--fibre-width", type=float, default=2.5,
                   help="line width of the fibre, in points")
    p.add_argument("--width", type=float, default=7.0,
                   help="figure width in inches")
    p.add_argument("--height", type=float, default=4.6,
                   help="figure height in inches")
    p.add_argument("--font-size", type=float, default=8.0,
                   help="base font size in points; titles, labels and legend "
                        "scale with it")
    p.add_argument("--dpi", type=int, default=300, help="resolution of the PNG")
    p.add_argument("--lang", default="en", choices=["en", "fr"],
                   help="language of the figure labels")
    p.add_argument("--out", required=True, help="prefix of the output files")
    a = p.parse_args()

    scale = a.font_size / 8.0
    plt.rcParams.update({
        "font.size": a.font_size, "axes.linewidth": 0.7,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })
    T = LABELS[a.lang]

    vp = chrom = None
    for r in read_tsv(a.loops):
        if r.get("gene") == a.gene:
            try:
                vp, chrom = float(r["vp_pos"]), r["chrom"]
            except (KeyError, ValueError):
                continue
            break
    if vp is None:
        sys.exit(f"{a.gene} has no loop in {a.loops}")
    start, end = vp - a.window, vp + a.window

    loops = load_loops(a.loops, chrom, start, end, a.min_terms, a.called_in)
    if a.only_gene:
        loops = [b for b in loops if b["gene"] == a.gene]
    if not loops:
        sys.exit("no loop left in the window after filtering")
    if a.max_anchors:
        loops.sort(key=lambda b: -b["n_terms"])
        kept, seen = [], set()
        for b in loops:
            seen |= {b["vp"], b["summit"]}
            kept.append(b)
            if len(seen) >= a.max_anchors:
                break
        loops = kept
    genes = sorted({b["gene"] for b in loops})
    print(f"  {len(loops)} loops, viewpoints: {', '.join(genes)}")

    elements = load_elements(a.elements, chrom) if a.elements else []
    if a.elements:
        print(f"  {len(elements)} elements read on {chrom}")

    landmarks = {}
    if a.gtf:
        landmarks = load_genes(a.gtf, chrom, start, end, a.gene_types,
                               a.gene_list)
        print(f"  {len(landmarks)} genes of the locus in the window")

    centre, name_of_centre, info = build_anchors(
        loops, a.merge, vp, elements, a.min_prevalence, landmarks)

    names = sorted(info, key=lambda n: info[n]["centre"])
    centres = {n: info[n]["centre"] for n in names}
    kinds = {n: info[n]["kind"] for n in names}
    n_el = sum(1 for n in names if kinds[n] == "enh")
    print(f"  {len(names)} anchors after merging at {a.merge/1000:.0f} kb, "
          f"{n_el} of them carrying a replicable element")

    contacts, seen = [], set()
    for b in loops:
        pa = name_of_centre[centre[b["vp"]]]
        pb = name_of_centre[centre[b["summit"]]]
        if pa == pb:
            continue
        k = frozenset((pa, pb))
        if k in seen:
            continue
        seen.add(k)
        contacts.append((pa, pb))
    if not contacts:
        sys.exit("no contact between two distinct anchors")

    fig, ax = plt.subplots(1, 1, figsize=(a.width, a.height))
    P, members = folding_panel(
        ax, names, centres, kinds, contacts, a.seed,
        fibre_width=a.fibre_width, label_size=6.6 * scale,
        number_size=6.2 * scale,
        geometry=dict(tail_spacing=a.tail_spacing, petal_base=a.petal_base,
                      petal_slope=a.petal_slope))

    # Framing: the content box is extended to the actual aspect ratio of the
    # axes, so that a drawing taller than wide keeps its top and bottom.
    fig.canvas.draw()
    bb = ax.get_position()
    ratio = (bb.width * a.width) / (bb.height * a.height)
    c = (P.min(0) + P.max(0)) / 2.0
    half = (P.max(0) - P.min(0)) / 2.0 * 1.10 + 0.25
    if half[0] / half[1] < ratio:
        half[0] = half[1] * ratio
    else:
        half[1] = half[0] / ratio
    ax.set_xlim(c[0] - half[0], c[0] + half[0])
    ax.set_ylim(c[1] - half[1], c[1] + half[1])

    fig.text(0.012, 0.985, T["title"].format(g=a.gene), fontsize=8.5 * scale,
             ha="left", va="top")
    fig.text(0.012, 0.942,
             T["subtitle"].format(c=chrom, a=f"{start:,.0f}", b=f"{end:,.0f}",
                                  n=len(loops), m=len(names)),
             fontsize=7 * scale, ha="left", va="top", color="#5A5A5A")
    handles = [Line2D([], [], color=PROMOTER_COLOUR, lw=1.6, marker=">", ms=4,
                      label=T["promoter"])]
    if n_el:
        handles.append(Line2D([], [], color=OLIVE, lw=3.2, label=T["element"]))
    if any(g == "gene" for g in kinds.values()):
        handles.append(Line2D([], [], color=GENE_COLOUR, lw=1.1, marker=">",
                              ms=4, label=T["gene"]))
    handles += [Line2D([], [], color=GREY, lw=3.2, label=T["other"]),
                Line2D([], [], color=GOLD, lw=0, marker="o", mfc=GOLD,
                       mec="white", ms=7, label=T["contact"])]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.52, 0.0), ncol=len(handles), frameon=False,
               fontsize=6.5 * scale, handlelength=1.5, columnspacing=1.6,
               handletextpad=0.5, labelcolor="#3A3A3A")

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=a.dpi, bbox_inches="tight")

    with open(out.with_suffix(".tsv"), "w") as fh:
        fh.write("petal\tanchor_a\tanchor_b\tposition_a\tposition_b\t"
                 "length\tkind_a\tkind_b\n")
        for k in range(len(members) - 1):
            x, y = members[k], members[k + 1]
            fh.write(f"{k+1}\t{x}\t{y}\t{centres[x]:.0f}\t{centres[y]:.0f}\t"
                     f"{centres[y]-centres[x]:.0f}\t{kinds[x]}\t"
                     f"{kinds[y]}\n")
    print(f"  rosette: {' - '.join(members)}")
    print(f"  written {out.with_suffix('.svg')}, .png and .tsv")
    plt.close(fig)


if __name__ == "__main__":
    main()
