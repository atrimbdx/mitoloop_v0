#!/usr/bin/env python3
"""
peaks_vs_elements.py - overlap of loop anchors with replicable rE2G elements

Purpose
    Cross two sets of coordinates around the same promoters: the loop summits
    called by peakC on tissue Hi-C virtual 4C profiles (loops.tsv), optionally
    restricted with --called-in to the ENCODE terms that called them, and the
    elements ENCODE-rE2G links to the same gene in at least --min-prevalence
    biosamples (re2g_elements.tsv). Report how often they coincide, for the
    target (mito) and control promoters, and test the coincidence against
    two permutation nulls.

Inputs
    --elements   re2g_elements.tsv written by re2g_panel.py: gene, set, chrom,
                 start, end, distance (element midpoint minus TSS, bp),
                 prevalence (biosamples linking the element to the gene).
    --loops      loops.tsv written by v4c_loops.R reciprocal: gene, chrom,
                 vp_pos, summit, and called (a '|'-separated list of ENCODE
                 terms, needed by --called-in and --list-terms). A table with
                 a distance column (summit minus viewpoint, bp), such as
                 assoc_loops.tsv, is also accepted; without that column the
                 distance is computed as summit - vp_pos.

Outputs
    printed summary, and with --out, <out>/peaks_vs_elements.tsv, one row per
    gene set: set, genes, replicable_elements, elements_on_a_peak,
    share_elements_on_a_peak, peaks, peaks_on_an_element,
    share_peaks_on_an_element, expected_share_if_anywhere and, with
    --perms > 0, observed_peaks_on_an_element_by_distance,
    observed_elements_on_a_peak_by_distance, then for mode in (uniform,
    distance): null_<mode>_share_peaks, p_<mode>_peaks,
    null_<mode>_share_elements, p_<mode>_elements.

Method
    A loop summit and an element coincide when they are on the same
    chromosome and the summit lies within --tol of the element. The overlap
    is read in two directions, which answer different questions:

      elements on a peak   of the replicable elements of a gene, the share
                           that falls inside one of its loop anchors
      peaks on an element  of the loops of a gene, the share whose anchor
                           holds one of its replicable elements

    expected_share_if_anywhere is the anchor length over the window length
    (two sides of --window beyond --min-distance), averaged over genes.

    Two permutation nulls are drawn, --perms times each, gene by gene,
    keeping the loops fixed and the number of elements of each gene:

      uniform   element distances to the TSS drawn uniformly between
                --min-distance and --window, with a random sign; a test
                against a flat distribution.
      distance  element distances redrawn with replacement from the pooled
                distances of all replicable elements of the same gene set,
                sign included. The crowding of elements near the promoter is
                preserved and the pairing between a gene and its own
                elements is broken. This is the conservative null, since
                elements and loop anchors both concentrate at short
                distances.

    The empirical p-value is (1 + number of permutations with a count at
    least the observed one) / (perms + 1). An overlap that departs from the
    uniform null and not from the distance-matched one says only that both
    quantities crowd at the same distance from the promoter.

    The test works on distances to the TSS (element midpoint against summit,
    within --tol), which is equivalent to the coordinate test at this
    tolerance; its observed counts are printed beside the coordinate counts
    and the two can differ slightly at element edges. Elements closer than
    --min-distance to the TSS can be left out when loops are called only
    beyond a minimum distance, since they cannot fall on an anchor.

Usage
    python3 peaks_vs_elements.py --elements results/re2g/re2g_elements.tsv \\
        --loops results/loops/loops.tsv --list-terms

    python3 peaks_vs_elements.py --elements results/re2g/re2g_elements.tsv \\
        --loops results/loops/loops.tsv --min-prevalence 51 \\
        --min-distance 50000 --called-in "psoas muscle|gastrocnemius" \\
        --out results/overlap

References
    Gschwind AR et al. 2023. An encyclopedia of enhancer-gene regulatory
        interactions in the human genome. bioRxiv.
    ENCODE Project Consortium et al. 2020. Expanded encyclopaedias of DNA
        elements in the human and mouse genomes (Registry of candidate
        cis-regulatory elements). Nature.
    Frankish A et al. 2023. GENCODE: reference annotation for the human and
        mouse genomes in 2023. Nucleic Acids Res.
    Geeven G et al. 2018. peakC: a flexible, non-parametric peak calling
        package for 4C and Capture-C data. Nucleic Acids Res.
    Phipson B, Smyth GK. 2010. Permutation P-values should never be zero:
        calculating exact P-values when permutations are randomly drawn.
        Stat Appl Genet Mol Biol.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.strip():
                f = line.rstrip("\n").split("\t")
                yield dict(zip(head, f + [""] * (len(head) - len(f))))


def count_hits(el, pk, tol):
    """(peaks holding an element, elements on a peak) over a dict of arrays."""
    import numpy as np
    n_pk = n_el = 0
    for g, e in el.items():
        p_ = pk.get(g)
        if p_ is None or len(p_) == 0 or len(e) == 0:
            continue
        m = np.abs(e[:, None] - p_[None, :]) <= tol
        n_pk += int(m.any(axis=0).sum())
        n_el += int(m.any(axis=1).sum())
    return n_pk, n_el


def permute(el, pk, tol, perms, window, dmin, rng, mode, pooled):
    """Empirical null for both statistics. mode is 'uniform' or 'distance'."""
    import numpy as np
    out_pk = np.empty(perms)
    out_el = np.empty(perms)
    for i in range(perms):
        sim = {}
        for g, e in el.items():
            n = len(e)
            if mode == "uniform":
                d = rng.uniform(dmin, window, n) * rng.choice((-1.0, 1.0), n)
            else:
                d = rng.choice(pooled, size=n, replace=True)
            sim[g] = d
        out_pk[i], out_el[i] = count_hits(sim, pk, tol)
    return out_pk, out_el


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--elements", required=True,
                   help="re2g_elements.tsv from re2g_panel.py")
    p.add_argument("--loops", required=True,
                   help="loops.tsv (carries the 'called' column and so allows "
                        "--called-in), or a loop table with a distance column")
    p.add_argument("--called-in", default="",
                   help="regular expression on the ENCODE terms that called "
                        "the loop, case ignored; keeps only the loops called "
                        "in a matching tissue. Without it, loops of every "
                        "tissue are compared with elements of the biosamples "
                        "selected in re2g_panel.py")
    p.add_argument("--list-terms", action="store_true",
                   help="print the terms found in the 'called' column and stop")
    p.add_argument("--min-prevalence", type=int, default=51,
                   help="biosamples a link must appear in to count as "
                        "replicable; set it relative to the number of "
                        "biosamples kept by re2g_panel.py --keep (default "
                        "%(default)s)")
    p.add_argument("--tol", type=float, default=10_000,
                   help="half-width of a loop anchor, bp; one bin of the "
                        "contact map (default %(default)s)")
    p.add_argument("--window", type=float, default=1_000_000,
                   help="half-width of the window around the TSS, bp "
                        "(default %(default)s)")
    p.add_argument("--min-distance", type=float, default=0,
                   help="drop elements closer than this to the TSS, bp; when "
                        "loops are called only beyond a minimum distance, "
                        "closer elements cannot fall on an anchor and only "
                        "dilute the count (default %(default)s)")
    p.add_argument("--perms", type=int, default=1000,
                   help="permutations for each null; 0 skips the test "
                        "(default %(default)s)")
    p.add_argument("--seed", type=int, default=0,
                   help="seed of the permutations (default %(default)s)")
    p.add_argument("--out", default=None,
                   help="output directory for peaks_vs_elements.tsv; "
                        "optional")
    a = p.parse_args()

    import re
    rows = list(read_tsv(a.loops))
    has_called = bool(rows) and "called" in rows[0]
    if a.list_terms:
        if not has_called:
            sys.exit(f"{a.loops} has no 'called' column; use loops.tsv")
        tally = defaultdict(int)
        for r in rows:
            for t in r["called"].split("|"):
                if t:
                    tally[t] += 1
        print(f"\n{len(tally)} terms called at least one loop:")
        for t, n in sorted(tally.items(), key=lambda x: -x[1]):
            print(f"  {n:6d}  {t}")
        return
    if a.called_in:
        if not has_called:
            sys.exit(f"--called-in needs the 'called' column, absent from "
                     f"{a.loops}; pass loops.tsv")
        rx = re.compile(a.called_in, re.I)
        n0 = len(rows)
        rows = [r for r in rows
                if any(rx.search(t) for t in r["called"].split("|") if t)]
        print(f"--called-in {a.called_in!r}: {len(rows)} of {n0} loops kept")
        if not rows:
            sys.exit("no loop left; --list-terms shows the available terms")

    peaks = defaultdict(list)          # gene -> [(chrom, summit)]
    peak_d = defaultdict(list)         # gene -> [distance to TSS]
    gene_set = {}
    for r in rows:
        try:
            peaks[r["gene"]].append((r["chrom"], float(r["summit"])))
        except (KeyError, ValueError):
            continue
        # loops.tsv carries no distance column: the distance of a summit to
        # its viewpoint is then summit - vp_pos, the definition used by the
        # loop tables that carry one
        try:
            if r.get("distance", "") != "":
                peak_d[r["gene"]].append(float(r["distance"]))
            else:
                peak_d[r["gene"]].append(float(r["summit"])
                                         - float(r["vp_pos"]))
        except (KeyError, ValueError):
            pass

    els = defaultdict(list)            # gene -> [(chrom, start, end)]
    el_d = defaultdict(list)           # gene -> [distance to TSS]
    for r in read_tsv(a.elements):
        try:
            if int(float(r["prevalence"])) < a.min_prevalence:
                continue
            if abs(float(r["distance"])) < a.min_distance:
                continue
            els[r["gene"]].append((r["chrom"], float(r["start"]),
                                   float(r["end"])))
            el_d[r["gene"]].append(float(r["distance"]))
            gene_set[r["gene"]] = r.get("set", "mito")
        except (KeyError, ValueError):
            continue

    genes = sorted(set(peaks) & set(els))
    if not genes:
        sys.exit("no gene carries both peaks and replicable elements; check "
                 "that the two tables use the same gene symbols")

    out = []
    for label in ("mito", "control"):
        gs = [g for g in genes if gene_set.get(g) == label]
        if not gs:
            continue
        n_el = n_el_hit = n_pk = n_pk_hit = 0
        cover = 0.0
        for g in gs:
            P = peaks[g]
            E = els[g]
            n_el += len(E)
            n_pk += len(P)
            cover += len(P) * 2 * a.tol / (2 * (a.window - a.min_distance))
            on = set()
            for i, (pc, ps) in enumerate(P):
                for c, s, e in E:
                    if c == pc and s - a.tol <= ps <= e + a.tol:
                        on.add(i)
            n_pk_hit += len(on)
            for c, s, e in E:
                if any(pc == c and s - a.tol <= ps <= e + a.tol
                       for pc, ps in P):
                    n_el_hit += 1
        out.append({
            "set": label, "genes": len(gs),
            "replicable_elements": n_el, "elements_on_a_peak": n_el_hit,
            "share_elements_on_a_peak": n_el_hit / n_el if n_el else 0.0,
            "peaks": n_pk, "peaks_on_an_element": n_pk_hit,
            "share_peaks_on_an_element": n_pk_hit / n_pk if n_pk else 0.0,
            "expected_share_if_anywhere": cover / len(gs)})

    if a.perms:
        import numpy as np
        rng = np.random.default_rng(a.seed)
        for r in out:
            gs = [g for g in genes if gene_set.get(g) == r["set"]]
            el = {g: np.asarray(el_d[g], float) for g in gs if el_d.get(g)}
            pk = {g: np.asarray(peak_d[g], float) for g in gs if peak_d.get(g)}
            if not el or not pk:
                print(f"  {r['set']}: no distance to the TSS for the loops "
                      "or the elements, permutation test skipped")
                continue
            pooled = np.concatenate(list(el.values()))
            obs_pk, obs_el = count_hits(el, pk, a.tol)
            r["observed_peaks_on_an_element_by_distance"] = obs_pk
            r["observed_elements_on_a_peak_by_distance"] = obs_el
            n_pk = sum(len(v) for v in pk.values())
            n_el = sum(len(v) for v in el.values())
            for mode in ("uniform", "distance"):
                s_pk, s_el = permute(el, pk, a.tol, a.perms, a.window,
                                     a.min_distance, rng, mode, pooled)
                r[f"null_{mode}_share_peaks"] = float(s_pk.mean()) / n_pk
                r[f"p_{mode}_peaks"] = (
                    1 + int((s_pk >= obs_pk).sum())) / (a.perms + 1)
                r[f"null_{mode}_share_elements"] = float(s_el.mean()) / n_el
                r[f"p_{mode}_elements"] = (
                    1 + int((s_el >= obs_el).sum())) / (a.perms + 1)

    print(f"\n{len(genes)} genes with both loop summits and elements replicable "
          f"in >= {a.min_prevalence} biosamples, anchor +/-{a.tol/1000:.0f} kb\n")
    for r in out:
        print(f"  {r['set']}  ({r['genes']} genes)")
        print(f"    {r['replicable_elements']:,} replicable elements, "
              f"{r['elements_on_a_peak']:,} inside a peak "
              f"({100*r['share_elements_on_a_peak']:.1f}%)")
        print(f"    {r['peaks']:,} peaks, {r['peaks_on_an_element']:,} holding "
              f"one ({100*r['share_peaks_on_an_element']:.1f}%)")
        print(f"    if elements fell anywhere in the window: "
              f"{100*r['expected_share_if_anywhere']:.1f}%")
        if "p_uniform_peaks" in r:
            print(f"    on distances, {r['observed_peaks_on_an_element_by_distance']:,}"
                  f" peaks and "
                  f"{r['observed_elements_on_a_peak_by_distance']:,} elements")
            for mode in ("uniform", "distance"):
                print(f"    {mode:8s} null: peaks "
                      f"{100*r[f'null_{mode}_share_peaks']:.1f}%, "
                      f"p = {r[f'p_{mode}_peaks']:.4f}   elements "
                      f"{100*r[f'null_{mode}_share_elements']:.1f}%, "
                      f"p = {r[f'p_{mode}_elements']:.4f}")
        print()

    if a.out:
        d = Path(a.out)
        d.mkdir(parents=True, exist_ok=True)
        cols = list(out[0])
        for r in out:
            for c in r:
                if c not in cols:
                    cols.append(c)
        for r in out:
            for c in cols:
                r.setdefault(c, "")
        with open(d / "peaks_vs_elements.tsv", "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in out:
                fh.write("\t".join(
                    f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c])
                    for c in cols) + "\n")
        print(f"Written to {d}/peaks_vs_elements.tsv")


if __name__ == "__main__":
    main()
