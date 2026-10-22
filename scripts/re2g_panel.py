#!/usr/bin/env python3
"""
re2g_panel.py - number and replicability of ENCODE-rE2G links per gene

Purpose
    For every gene of a target set (nuclear-encoded mitochondrial genes) and,
    optionally, of a control set, count the candidate elements of the ENCODE
    cCRE registry within --window of the promoter and the elements that
    ENCODE-rE2G links to the gene, then measure how replicable these links
    are across biosamples. The quantities reported per gene set are:

      1. the number of registry elements within +/- --window of the TSS;
      2. the number of elements linked to the gene in one biosample (median
         over biosamples), over all biosamples (union) and in every
         biosample at once (intersection);
      3. the share of candidates never linked to the gene;
      4. the share of genes that receive no link in any biosample;
      5. among the genes with at least one link, the share whose
         intersection over the biosamples is empty.

    With --controls the same quantities are computed on a control gene set,
    which says whether the annotation behaves differently at the target
    promoters or the same way everywhere.

Inputs
    --links-dir   directory of ENCODE-rE2G element-gene link files, one per
                  biosample, tab-separated, gzipped or not. The first line is
                  a header (a leading '#' is ignored), so the column order
                  does not matter. Columns 1-3 are chrom, start, end of the
                  element. Required columns: TargetGene, Score. Used when
                  present: CellType (biosample, for --keep), isSelfPromoter,
                  TargetGeneEnsemblID, TargetGeneTSS, distanceToTSS.Feature
                  and the other '<name>.Feature' columns, from which the
                  model is recognised. --links-glob selects the files.
    --manifest    optional TSV restricting the files, one row per biosample,
                  with a column holding the file name or accession
                  (--manifest-column, guessed when absent) and optionally a
                  biosample_type column, used only for reporting.
    --ccre        ENCODE cCRE registry, BED(.gz) on the same assembly
                  (chrom, start, end). Optional; without it the candidate
                  counts are left out and elements keep their own
                  coordinates.
    --gtf         GENCODE annotation (GTF, gzipped or not), for the TSS and
                  the Ensembl gene id of each gene.
    --genes       target gene symbols, whitespace-separated.
    --controls    control gene symbols, whitespace-separated. Optional.
    --aliases     optional extra symbol table, two tab-separated columns
                  (old symbol, current symbol), tried before the built-in
                  table; --no-builtin-aliases disables the built-in one.

Outputs (in --out)
    re2g_elements.tsv       one row per gene and element: gene, set, chrom,
                            start, end, distance (element midpoint minus TSS,
                            bp), prevalence (biosamples linking it),
                            celltypes (their CellType values, ';')
    re2g_per_gene.tsv       one row per gene: gene, set, chrom, tss,
                            n_ccre_window, links_per_biosample_median,
                            links_per_biosample_max, n_union, n_intersection,
                            n_biosamples_with_link, max_prevalence,
                            frac_never_linked
    re2g_summary.tsv        the quantities 1-5 per gene set
    re2g_sentence.txt       a short text giving the quantities of the target
                            set, with the score threshold used
    re2g_core_curve.tsv     k, then per set the share of genes with an element
                            linked in at least k biosamples (observed_*) and
                            in every biosample of random k-subsets
                            (subsampled_*)
    re2g_core_curve.svg     the same as a curve
    re2g_elements_per_gene.svg
                            cumulative distributions of registry candidates
                            per gene and of linked elements per gene (one
                            biosample, union, intersection)
    biosamples_found.tsv    file, CellType; written when --keep is empty or
                            with --list-biosamples
    genes_without_tss.txt   symbols left out for lack of a TSS

Method
    Score threshold. The distributed files carry every candidate pair with
    its score. ENCODE-rE2G binarises its predictions at the score giving 70%
    recall on its CRISPR benchmark, and that operating point differs between
    models. The model is recognised from the '.Feature' columns of the first
    file (Jaccard similarity with the feature set of each model, at least
    --model-min-similarity) and its published threshold is applied. Some
    models share the same feature set (megamap and intact Hi-C variants);
    the first one listed in the table is then taken, and --model or
    --threshold sets the operating point explicitly. --threshold 0 keeps
    every row, which shows how much the choice matters.

    Distal links. A link counts as distal when it is not flagged
    isSelfPromoter and lies at least --min-distance from the TSS; only
    distal links enter the counts. Links beyond --window are ignored.

    Coordinate system. rE2G builds its candidates from the accessible peaks
    of each biosample, so one element carries slightly different coordinates
    from one biosample to the next. Each element is mapped onto the registry
    element it overlaps most, or failing an overlap the nearest one within
    --ccre-tolerance, so that biosamples share one coordinate system; without
    it the union is inflated and the intersection empty by construction.

    Gene symbols. Genes are matched in the link files by Ensembl gene id
    first, which is stable across renames, then by symbol. Symbols absent
    from the GTF are resolved through an alias table of HGNC renames; the
    first candidate present in the GTF is used. Entries of gene panels that
    are not genes (repeat expansions, copy-number regions) are dropped
    silently. Genes on the mitochondrial genome are left out, since rE2G
    makes no prediction there.

    Core curve. For k = 1..n biosamples: the exact share of genes whose most
    recurrent element is linked in at least k biosamples, and the share of
    genes with one element linked in every biosample of a random subset of
    size k (--curve-draws subsets per k). The second reading tells whether
    the value at k = n is a plateau or a point on a decline. Each element
    carries one bit per biosample, so "linked in every biosample of the
    subset" is one integer operation.

Usage
    python3 re2g_panel.py --links-dir re2g --list-biosamples --gtf \\
        data/gencode.gtf.gz --genes genes.txt --out results/re2g

    python3 re2g_panel.py --links-dir re2g --ccre data/cCREs.bed.gz \\
        --gtf data/gencode.gtf.gz --genes genes.txt --controls controls.txt \\
        --keep "skeletal muscle|psoas muscle|gastrocnemius" \\
        --out results/re2g --cores 8

References
    Gschwind AR et al. 2023. An encyclopedia of enhancer-gene regulatory
        interactions in the human genome. bioRxiv.
    EngreitzLab/ENCODE_rE2G, GitHub repository: each directory of models/
        holds a file named threshold_<value> giving the operating point.
    Fulco CP et al. 2019. Activity-by-contact model of enhancer-promoter
        regulation from thousands of CRISPR perturbations. Nat Genet.
    ENCODE Project Consortium et al. 2020. Expanded encyclopaedias of DNA
        elements in the human and mouse genomes (Registry of candidate
        cis-regulatory elements). Nature.
    Frankish A et al. 2023. GENCODE: reference annotation for the human and
        mouse genomes in 2023. Nucleic Acids Res.
    Morales J et al. 2022. A joint NCBI and EMBL-EBI transcript set for
        clinical genomics and research (MANE). Nature.
    Seal RL et al. 2023. Genenames.org: the HGNC resources in 2023. Nucleic
        Acids Res.
    Martin AR et al. 2019. PanelApp crowdsources expert knowledge to
        establish consensus diagnostic gene panels. Nat Genet.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import sys
from collections import defaultdict
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np

WINDOW = 1_000_000
BLUE, RED, GOLD = "#2166AC", "#B2182B", "#C78A20"


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


# ---------------------------------------------------------------------------
# annotation

def load_tss(gtf: Path, wanted: set, ens_out: dict | None = None) -> dict:
    """MANE Select TSS where available, otherwise the gene record.

    With ens_out, also fills symbol -> Ensembl gene id (no version), the key
    used to find the gene in the link files whatever symbol they carry.
    """
    best = {}
    with opener(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] not in ("gene", "transcript"):
                continue
            a = f[8]
            i = a.find('gene_name "')
            if i < 0:
                continue
            name = a[i + 11:a.find('"', i + 11)]
            if name not in wanted:
                continue
            rank = 2 if (f[2] == "transcript" and "MANE_Select" in a) else (
                1 if f[2] == "gene" else 0)
            if rank == 0:
                continue
            tss = int(f[3]) - 1 if f[6] == "+" else int(f[4])
            j = a.find('gene_id "')
            ensg = a[j + 9:a.find('"', j + 9)].split(".")[0] if j >= 0 else ""
            prev = best.get(name)
            if prev is None or rank > prev[2]:
                best[name] = (f[0], tss, rank, ensg)
    if ens_out is not None:
        for k, v in best.items():
            if v[3]:
                ens_out[k] = v[3]
    return {k: (v[0], v[1]) for k, v in best.items()}


def load_ccre(path: Path):
    """Registry per chromosome: sorted starts, ends, and sorted midpoints."""
    iv = defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                iv[f[0]].append((int(f[1]), int(f[2])))
            except ValueError:
                continue
    reg = {}
    for c, v in iv.items():
        v.sort()
        st = np.array([x[0] for x in v], dtype=np.int64)
        en = np.array([x[1] for x in v], dtype=np.int64)
        reg[c] = {"start": st, "end": en,
                  "mid_sorted": np.sort((st + en) // 2)}
    return reg


def _reg_for(reg, chrom):
    for c in (chrom, chrom.replace("chr", ""), "chr" + chrom):
        if c in reg:
            return c, reg[c]
    return None, None


def map_to_ccre(reg, chrom, s, e, tol=500):
    """The registry element an rE2G element stands for.

    rE2G defines its candidates from the DNase peaks of each biosample, so the
    same element carries slightly different coordinates from one biosample to
    the next. Mapping every element onto the registry gives one coordinate
    system for all biosamples, without which the union is inflated and the
    intersection empty by construction. The element goes to the registry entry
    it overlaps most; failing an overlap, to the nearest one within `tol`;
    failing that, it keeps its own coordinates and is reported as unmapped.
    """
    c, r = _reg_for(reg, chrom)
    if r is None:
        return f"{chrom}:{s}-{e}", False
    st, en = r["start"], r["end"]
    i = int(np.searchsorted(st, e, "left")) - 1
    best, best_ov = None, 0
    k = i
    while k >= 0 and k > i - 6:
        ov = min(e, en[k]) - max(s, st[k])
        if ov > best_ov:
            best, best_ov = k, ov
        if en[k] < s - 2000:
            break
        k -= 1
    if best is None:
        m = (s + e) // 2
        mids = (st + en) // 2
        j = int(np.searchsorted(st, m))
        cand = [x for x in (j - 1, j) if 0 <= x < len(st)]
        if cand:
            x = min(cand, key=lambda q: abs(mids[q] - m))
            if abs(mids[x] - m) <= tol:
                best = x
    if best is None:
        return f"{chrom}:{s}-{e}", False
    return f"{c}:{st[best]}-{en[best]}", True


def ccre_in_window(reg, chrom, tss, window):
    _c, r = _reg_for(reg, chrom)
    if r is None or len(r["mid_sorted"]) == 0:
        return np.nan
    a = r["mid_sorted"]
    lo = bisect.bisect_left(a, tss - window)
    hi = bisect.bisect_right(a, tss + window)
    return int(hi - lo)


# Built-in alias table: HGNC renames frequent in mitochondrial gene lists
# (genenames.org). Each old symbol maps to candidates tried in order; the
# first one present in the GTF is used, which copes with GENCODE releases on
# either side of a rename. A symbol is never mapped to a name the GTF does not
# hold. --aliases adds entries, --no-builtin-aliases disables this table.
ALIASES = {
    "ATP5A1": ["ATP5F1A"], "ATP5B": ["ATP5F1B"], "ATP5C1": ["ATP5F1C"],
    "ATP5D": ["ATP5F1D"], "ATP5E": ["ATP5F1E"], "ATP5O": ["ATP5PO"],
    "ATP5F1": ["ATP5PB"], "ATP5H": ["ATP5PD"], "ATP5J": ["ATP5PF"],
    "ATP5J2": ["ATP5MF"], "ATP5I": ["ATP5ME"], "ATP5L": ["ATP5MG"],
    "ATP5G1": ["ATP5MC1"], "ATP5G2": ["ATP5MC2"], "ATP5G3": ["ATP5MC3"],
    "ATP5S": ["DMAC2L"], "ATP5SL": ["DMAC2"], "ATPIF1": ["ATP5IF1"],
    "USMG5": ["ATP5MK", "ATP5MD"], "C14orf2": ["ATP5MJ", "ATP5MPL"],
    "ATP5MD": ["ATP5MK"], "ATP5MPL": ["ATP5MJ"],
    "APOPT1": ["COA8"], "C10orf2": ["TWNK"], "TAZ": ["TAFAZZIN"],
    "C12orf65": ["MTRFR"], "C6orf203": ["MTRES1"], "C7orf55": ["FMC1"],
    "C20orf7": ["NDUFAF5"], "C8orf38": ["NDUFAF6"], "C17orf89": ["NDUFAF8"],
    "KIAA0391": ["PRORP"], "MTERFD1": ["MTERF3"], "MTERFD2": ["MTERF4"],
    "MTERFD3": ["MTERF2"], "ADCK3": ["COQ8A"], "ADCK4": ["COQ8B"],
    "CCDC90A": ["MCUR1"], "CCDC109A": ["MCU"], "CCDC109B": ["MCUB"],
    "C2orf47": ["MAIP1"], "C1orf31": ["COA6"], "FAM36A": ["COX20"],
    "CCDC56": ["COA3"], "C2orf64": ["COA5"], "C11orf83": ["UQCC3"],
    "MNF1": ["UQCC2"], "UQCC": ["UQCC1"], "C16orf91": ["UQCC4"],
    "KIAA0564": ["VWA8"], "FAM73A": ["MIGA1"], "FAM73B": ["MIGA2"],
    "QIL1": ["MICOS13"], "C19orf70": ["MICOS13"], "APOO": ["MICOS26"],
    "APOOL": ["MICOS27"], "MINOS1": ["MICOS10"], "GLTSCR2": ["NOP53"],
    "C15orf48": ["NMES1"], "MTL5": ["TESMIN"],
    # aminoacyl-tRNA synthetases and other HGNC renames of 2020-2023
    "DARS": ["DARS1"], "GARS": ["GARS1"], "IARS": ["IARS1"], "KARS": ["KARS1"],
    "LARS": ["LARS1"], "QARS": ["QARS1"], "AARS": ["AARS1"], "CARS": ["CARS1"],
    "EPRS": ["EPRS1"], "HARS": ["HARS1"], "MARS": ["MARS1"], "NARS": ["NARS1"],
    "RARS": ["RARS1"], "SARS": ["SARS1"], "TARS": ["TARS1"], "VARS": ["VARS1"],
    "WARS": ["WARS1"], "YARS": ["YARS1"],
    "PHB": ["PHB1"], "G6PC": ["G6PC1"], "SPATA5": ["AFG2A"],
    "SPATA5L1": ["AFG2B"], "ATP5L2": ["ATP5MGL"], "MRPS36": ["KGD4"],
    "CCDC58": ["TMEM242"], "GATD3A": ["GATD3"], "ZADH2": ["PTGR3"],
    "ADPRHL2": ["ADPRS"], "C19orf12": ["MPAN"], "FAM210B": ["ATAD3C"],
}

# Gene-panel entries (PanelApp) that are not genes: repeat expansions and
# copy-number regions, recognised by a substring, and clone-based names,
# recognised by a prefix. They are dropped without a warning, which would
# otherwise suggest an annotation problem. Defaults of --not-gene-substrings
# and --not-gene-prefixes.
NOT_GENES = ("_CTG", "_GAA", "_CGG", "_GCN", "-Loss", "-Gain", "ISCA-")
NOT_GENE_PREFIXES = ("RP11_", "RP11-")


def build_aliases(extra=None, builtin=True):
    """old symbol -> candidate current symbols, extra entries tried first."""
    table = defaultdict(list)
    if extra:
        with opener(extra) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 2 or not f[0].strip() or not f[1].strip():
                    continue
                old, new = f[0].strip(), f[1].strip()
                if new not in table[old]:
                    table[old].append(new)
    if builtin:
        for old, cands in ALIASES.items():
            for c in cands:
                if c not in table[old]:
                    table[old].append(c)
        # and back, after the forward entries: a GTF older than the rename
        # (GENCODE v26 of GTEx v8) holds the old symbol only
        for old, cands in ALIASES.items():
            for c in cands:
                if old not in table[c]:
                    table[c].append(old)
    return dict(table)


def resolve_symbols(symbols, gtf_names, aliases=None, substrings=NOT_GENES,
                    prefixes=NOT_GENE_PREFIXES):
    """symbol -> name found in the GTF (itself, or a renamed form)."""
    aliases = ALIASES if aliases is None else aliases
    substrings, prefixes = tuple(substrings), tuple(prefixes)
    out, unresolved, not_genes = {}, [], []
    for g in symbols:
        if any(t in g for t in substrings) or (prefixes and
                                                g.startswith(prefixes)):
            not_genes.append(g)
            continue
        if g in gtf_names:
            out[g] = g
            continue
        hit = next((c for c in aliases.get(g, []) if c in gtf_names), None)
        if hit:
            out[g] = hit
        else:
            unresolved.append(g)
    return out, unresolved, not_genes


def gtf_gene_names(gtf: Path) -> set:
    names = set()
    with opener(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            i = line.find('gene_name "')
            if i >= 0:
                names.add(line[i + 11:line.find('"', i + 11)])
    return names


def first_data_row(path):
    """CellType of a link file, read from its first data line."""
    try:
        with opener(path) as fh:
            head = fh.readline().rstrip("\n").lstrip("#").split("\t")
            row = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(head)}
        i = col.get("CellType")
        return row[i] if i is not None and i < len(row) else ""
    except Exception:                                     # noqa: BLE001
        return ""


# ---------------------------------------------------------------------------
# which model produced a link file, and the threshold that goes with it
#
# ENCODE-rE2G binarises its predictions at the score giving 70% recall on the
# CRISPR benchmark (Gschwind et al. 2023). That operating point differs
# between models, and each model is recognisable by the feature columns its
# files carry. Values from the EngreitzLab/ENCODE_rE2G GitHub repository,
# where each directory models/<name>/ holds a file named threshold_<value>.
# --model picks a row of this table, --threshold overrides it.

MODEL_THRESHOLDS = {
    "dhs_megamap": (0.201, {"ABC.Score", "contactFrequency", "distanceToTSS",
                            "normalizedDNase_prom", "numCandidateEnhGene",
                            "numTSSEnhGene", "sumNearbyEnhancers",
                            "ubiquitousExpressedGene"}),
    "dhs_intact_hic": (0.243, {"ABC.Score", "contactFrequency", "distanceToTSS",
                               "normalizedDNase_prom", "numCandidateEnhGene",
                               "numTSSEnhGene", "sumNearbyEnhancers",
                               "ubiquitousExpressedGene"}),
    "atac_megamap": (0.179, {"ABC.Score", "contactFrequency", "distanceToTSS",
                             "normalizedATAC_prom", "numCandidateEnhGene",
                             "numTSSEnhGene", "sumNearbyEnhancers",
                             "ubiquitousExpressedGene"}),
    "atac_intact_hic": (0.232, {"ABC.Score", "contactFrequency", "distanceToTSS",
                                "normalizedATAC_prom", "numCandidateEnhGene",
                                "numTSSEnhGene", "sumNearbyEnhancers",
                                "ubiquitousExpressedGene"}),
    "atac_h3k27ac_megamap": (0.203, {"ABC.Score", "contactFrequency",
                                     "distanceToTSS", "normalizedH3K27ac_prom",
                                     "numCandidateEnhGene", "numTSSEnhGene",
                                     "sumNearbyEnhancers",
                                     "ubiquitousExpressedGene"}),
    "atac_h3k27ac_intact_hic": (0.259, {"ABC.Score", "contactFrequency",
                                        "distanceToTSS",
                                        "normalizedH3K27ac_prom",
                                        "numCandidateEnhGene", "numTSSEnhGene",
                                        "sumNearbyEnhancers",
                                        "ubiquitousExpressedGene"}),
    "dhs_h3k27ac_megamap": (0.217, {"ABC.Score", "activity_prom",
                                    "contactFrequency", "distanceToTSS",
                                    "numCandidateEnhGene", "numTSSEnhGene",
                                    "sumNearbyEnhancers",
                                    "ubiquitousExpressedGene"}),
    "dhs_h3k27ac_intact_hic": (0.298, {"ABC.Score", "activity_prom",
                                       "contactFrequency", "distanceToTSS",
                                       "numCandidateEnhGene", "numTSSEnhGene",
                                       "sumNearbyEnhancers",
                                       "ubiquitousExpressedGene"}),
    "extended": (0.336, {"averageCorrWeighted", "P2PromoterClass",
                         "normalizedEP300_enh", "inTAD", "enhCTCF"}),
}
# "3DContact" in a file is the "contactFrequency" of the feature table
SYNONYMS = {"3DContact": "contactFrequency", "ABCScore": "ABC.Score"}


def detect_model(path, min_similarity=0.6):
    """(model name, threshold) from the feature columns of a link file.

    Models with identical feature sets cannot be told apart; the first one
    listed in MODEL_THRESHOLDS is returned.
    """
    try:
        with opener(path) as fh:
            head = fh.readline().rstrip("\n").lstrip("#").split("\t")
    except Exception:                                     # noqa: BLE001
        return None, None
    have = {SYNONYMS.get(c[:-8], c[:-8]) for c in head if c.endswith(".Feature")}
    if not have:
        return None, None
    best, best_score = None, -1.0
    for name, (_thr, feats) in MODEL_THRESHOLDS.items():
        inter = len(have & feats)
        score = inter / max(len(have | feats), 1)
        if score > best_score:
            best, best_score = name, score
    if best_score < min_similarity:
        return None, None
    return best, MODEL_THRESHOLDS[best][0]


# ---------------------------------------------------------------------------
# the link files

def scan_file(path, wanted, threshold, min_dist, window):
    """{gene: (set of element ids, set of distal element ids)} for one file."""
    out = defaultdict(lambda: (set(), set()))
    n_rows = n_kept = 0
    with opener(path) as fh:
        head = fh.readline().rstrip("\n").lstrip("#").split("\t")
        col = {name: i for i, name in enumerate(head)}
        need = ("TargetGene", "Score")
        miss = [c for c in need if c not in col]
        if miss:
            return Path(path).name, {}, 0, 0, f"columns missing: {miss}"
        i_g, i_s = col["TargetGene"], col["Score"]
        i_self = col.get("isSelfPromoter")
        i_d = col.get("distanceToTSS.Feature")
        i_tss = col.get("TargetGeneTSS")
        i_ens = col.get("TargetGeneEnsemblID")
        # The files hold every candidate pair, and almost all of them score
        # far below any usable threshold. When Score is the last column, read
        # it from the end of the line and drop the row before splitting the
        # rest, which matters on files of several million rows.
        fast = (i_s == len(head) - 1) and threshold > 0
        for line in fh:
            n_rows += 1
            if fast:
                try:
                    if float(line[line.rfind("\t") + 1:]) < threshold:
                        continue
                except ValueError:
                    continue
            f = line.rstrip("\n").split("\t")
            if len(f) <= i_s:
                continue
            g = None
            if i_ens is not None and i_ens < len(f):
                g = wanted.get(f[i_ens].split(".")[0])
            if g is None:
                g = wanted.get(f[i_g])
            if g is None:
                continue
            if not fast:
                try:
                    if float(f[i_s]) < threshold:
                        continue
                except ValueError:
                    continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            eid = (f[0], start, end)
            if i_d is not None and f[i_d] not in ("", "NA"):
                try:
                    dist = abs(float(f[i_d]))
                except ValueError:
                    dist = np.nan
            elif i_tss is not None:
                try:
                    dist = abs((start + end) / 2 - float(f[i_tss]))
                except ValueError:
                    dist = np.nan
            else:
                dist = np.nan
            if np.isfinite(dist) and dist > window:
                continue
            n_kept += 1
            allset, distset = out[g]
            allset.add(eid)
            selfp = (i_self is not None and f[i_self].upper() in ("TRUE", "1"))
            if not selfp and (not np.isfinite(dist) or dist >= min_dist):
                distset.add(eid)
    return (Path(path).name, {g: (a, d) for g, (a, d) in out.items()},
            n_rows, n_kept, None)


# ---------------------------------------------------------------------------

def summarise(per_gene, label):
    """The headline numbers for one gene set."""
    g = [r for r in per_gene if r["set"] == label]
    if not g:
        return None
    arr = lambda k: np.array([r[k] for r in g], dtype=float)  # noqa: E731
    med = lambda v: float(np.nanmedian(v)) if np.isfinite(v).any() else np.nan  # noqa: E731
    ccre = arr("n_ccre_window")
    linked_any = arr("n_union")
    with_link = linked_any > 0
    never = arr("frac_never_linked")
    inter = arr("n_intersection")
    out = {
        "set": label,
        "n_genes": len(g),
        "ccre_per_gene_median": med(ccre),
        "links_per_biosample_median": med(arr("links_per_biosample_median")),
        "links_per_biosample_mean": float(np.nanmean(arr("links_per_biosample_median"))),
        "union_median": med(linked_any),
        "intersection_median": med(inter),
        "pct_candidates_never_linked": 100 * med(never),
        "pct_genes_without_any_link": 100 * float(np.mean(~with_link)),
        "pct_linked_genes_with_empty_intersection":
            100 * float(np.mean(inter[with_link] == 0)) if with_link.any() else np.nan,
        "pct_genes_with_link_in_every_biosample":
            100 * float(np.mean(inter[with_link] > 0)) if with_link.any() else np.nan,
    }
    return out


def window_text(window):
    """+/- window written in Mb or kb."""
    return (f"{window / 1e6:g} Mb" if window >= 1e6 else f"{window / 1e3:g} kb")


def sentence(s, n_biosamples, n_types, threshold, window=WINDOW):
    """The quantities of one gene set, as a short text."""
    if s is None:
        return ""
    def f(x, d=0):
        return "[X]" if not np.isfinite(x) else f"{x:,.{d}f}"
    return (
        f"Within +/-{window_text(window)} of the TSS, the ENCODE cCRE registry "
        f"lists a median of {f(s['ccre_per_gene_median'])} candidate elements "
        f"per gene ({s['n_genes']} genes of the '{s['set']}' set). Across "
        f"{n_biosamples} biosamples of {n_types} types, ENCODE-rE2G links a "
        f"median of {f(s['links_per_biosample_median'], 1)} distal elements per "
        "gene per biosample, and for the median gene "
        f"{f(s['pct_candidates_never_linked'], 1)}% of candidates are never "
        f"linked to it. {f(s['pct_genes_without_any_link'], 1)}% of genes "
        "receive no link in any biosample; among those that receive one, "
        f"{f(s['pct_linked_genes_with_empty_intersection'], 1)}% have no "
        f"element linked in all {n_biosamples} biosamples.\n\n"
        f"(score threshold {threshold})\n")


def core_curve(rows, masks, mask_gene, n_bio, draws, rng):
    """Share of genes with an element shared by k biosamples, k = 1..n_bio.

    Two readings of the same question.

    `observed` is exact: the share of genes whose most recurrent element is
    linked in at least k of the n_bio biosamples.

    `subsampled` answers what the analysis would have returned with only k
    biosamples: random subsets of size k are drawn, and a gene counts when one
    element is linked in every biosample of the subset. It is the quantity that
    decides whether the value at k = n_bio is a plateau or a point on a slope
    that would keep falling with more biosamples.

    Each element carries a bit per biosample, so "linked in every biosample of
    the subset" is one integer operation: mask & subset == subset. The bits
    are held in as many 64-bit words as the number of biosamples requires.
    """
    out = []
    sets = np.array([r["set"] for r in rows])
    prev = np.array([r["max_prevalence"] for r in rows])
    n_words = max(1, (n_bio + 63) // 64)
    full = (1 << 64) - 1
    m = np.array([[(w >> (64 * j)) & full for j in range(n_words)]
                  for w in masks], dtype=np.uint64).reshape(-1, n_words)
    gi = np.array(mask_gene, dtype=np.int64)
    n_genes = len(rows)
    order = np.argsort(gi, kind="stable")
    m, gi = m[order], gi[order]
    starts = np.searchsorted(gi, np.arange(n_genes))

    for k in range(1, n_bio + 1):
        row = {"k": k}
        for label in ("mito", "control"):
            sel = sets == label
            if not sel.any():
                continue
            row[f"observed_{label}"] = float((prev[sel] >= k).mean())
        if draws > 0:
            hits = np.zeros(n_genes, dtype=bool)
            tally = np.zeros(n_genes, dtype=np.int64)
            for _ in range(draws):
                idx = rng.choice(n_bio, size=k, replace=False)
                sub = np.zeros(n_words, dtype=np.uint64)
                for b in idx:
                    b = int(b)
                    sub[b // 64] |= np.uint64(1) << np.uint64(b % 64)
                ok = ((m & sub) == sub).all(axis=1)
                hits[:] = False
                if ok.any():
                    hits[np.unique(gi[ok])] = True
                tally += hits
            share = tally / draws
            for label in ("mito", "control"):
                sel = sets == label
                if sel.any():
                    row[f"subsampled_{label}"] = float(share[sel].mean())
        out.append(row)
    return out


def curve_figure(curve, path, n_bio):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    k = [c["k"] for c in curve]
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="0.88", lw=0.6)
    ax.set_axisbelow(True)
    drawn = False
    for label, colour in (("mito", RED), ("control", BLUE)):
        key = f"subsampled_{label}"
        if key in curve[0]:
            ax.plot(k, [c[key] for c in curve], color=colour, lw=2,
                    label=f"{label}, random subsets")
            drawn = True
        key = f"observed_{label}"
        if key in curve[0]:
            ax.plot(k, [c[key] for c in curve], color=colour, lw=1.2, ls="--",
                    label=f"{label}, at least k of {n_bio}")
            drawn = True
    if not drawn:
        plt.close(fig)
        return
    ax.set_xlabel("biosamples required, k", fontsize=9)
    ax.set_ylabel("share of genes with a shared element", fontsize=9)
    ax.set_ylim(0, 1)
    ax.set_xlim(1, n_bio)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(path, format="svg")
    plt.close(fig)

def figure(per_gene, path, n_biosamples, window=WINDOW,
           labels=("mitochondrial genes", "matched controls")):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib missing, figure skipped")
        return
    mito = [r for r in per_gene if r["set"] == "mito"]
    ctrl = [r for r in per_gene if r["set"] == "control"]
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.8))
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
        a.tick_params(length=3, labelsize=9)
        a.grid(axis="y", color="0.88", lw=0.6)
        a.set_axisbelow(True)

    def ecdf(a, v, colour, label, ls="-"):
        v = np.array([x for x in v if np.isfinite(x)], dtype=float)
        if v.size == 0:
            return None
        x = np.sort(v)
        y = np.arange(1, x.size + 1) / x.size
        a.step(x, y, where="post", color=colour, lw=2, ls=ls, label=label)
        return x, y

    ccre_m = [r["n_ccre_window"] for r in mito]
    if np.isfinite(ccre_m).any():
        ecdf(ax[0], ccre_m, RED, labels[0])
        if ctrl:
            ecdf(ax[0], [r["n_ccre_window"] for r in ctrl], BLUE,
                 labels[1])
        ax[0].set_xscale("log")
        ax[0].set_xlabel("candidate elements within $\\pm$"
                         f"{window_text(window)}", fontsize=9)
        ax[0].set_ylabel("share of genes", fontsize=9)
        ax[0].legend(frameon=True, framealpha=0.9, edgecolor="none",
                     fontsize=8, loc="best")
        med = np.nanmedian(ccre_m)
        if np.isfinite(med):
            ax[0].axvline(med, color=RED, lw=0.8, ls=":")
            ax[0].annotate(f"median {med:,.0f}", (med, 0.06),
                           xytext=(4, 0), textcoords="offset points",
                           fontsize=8, color=RED, va="bottom")
    else:
        ax[0].text(0.5, 0.5, "no cCRE registry given", ha="center",
                   va="center", fontsize=9, color="0.4")
        ax[0].set_xticks([]); ax[0].set_yticks([])

    series = [("elements linked, union over biosamples", "n_union", BLUE),
              ("in one biosample (median)", "links_per_biosample_median", GOLD),
              (f"in all {n_biosamples} biosamples", "n_intersection", RED)]
    xmax = 1
    for lab, key, colour in series:
        v = [r[key] for r in mito]
        r = ecdf(ax[1], v, colour, lab)
        if r is not None:
            xmax = max(xmax, float(np.nanpercentile(r[0], 99)) + 1)
    ax[1].set_xlim(-0.3, max(2.0, xmax))
    ax[1].set_xlabel("elements linked to the gene", fontsize=9)
    ax[1].set_ylabel(f"share of {labels[0]}", fontsize=9)
    # "best" so the box avoids the steps, whose position depends on the data
    ax[1].legend(frameon=True, framealpha=0.9, edgecolor="none",
                 fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(path, format="svg")
    plt.close(fig)
    print(f"  figure -> {path}")


# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--links-dir", required=True,
                   help="directory of ENCODE-rE2G link files, one per "
                        "biosample")
    p.add_argument("--links-glob", default="",
                   help="pattern of the link files inside --links-dir; "
                        "*.bed.gz, *.bed, *.tsv.gz, *.tsv, *.txt.gz, *.gz are "
                        "tried in turn when this is left out")
    p.add_argument("--manifest", default=None,
                   help="optional TSV restricting the files to its rows")
    p.add_argument("--manifest-column", default="",
                   help="column of --manifest holding the file name or "
                        "accession; guessed when absent")
    p.add_argument("--ccre", default=None,
                   help="ENCODE cCRE registry BED(.gz); optional")
    p.add_argument("--ccre-tolerance", type=int, default=500,
                   help="an element overlapping no registry element is mapped "
                        "to the nearest one whose midpoint lies within this "
                        "many bp (default %(default)s)")
    p.add_argument("--gtf", required=True, help="GENCODE GTF(.gz)")
    p.add_argument("--genes", required=True,
                   help="target gene symbols, whitespace-separated")
    p.add_argument("--controls", default=None,
                   help="control gene symbols, whitespace-separated; optional")
    p.add_argument("--aliases", default=None,
                   help="extra alias table, two tab-separated columns (old "
                        "symbol, current symbol), tried before the built-in "
                        "one; optional")
    p.add_argument("--no-builtin-aliases", action="store_true",
                   help="do not use the built-in table of HGNC renames")
    p.add_argument("--not-gene-substrings", default=",".join(NOT_GENES),
                   help="comma-separated substrings marking list entries that "
                        "are not genes (default %(default)s); empty disables")
    p.add_argument("--not-gene-prefixes", default=",".join(NOT_GENE_PREFIXES),
                   help="comma-separated prefixes marking list entries that "
                        "are not genes (default %(default)s); empty disables")
    p.add_argument("--threshold", type=float, default=None,
                   help="score threshold; overrides --model and the model "
                        "detected from the feature columns")
    p.add_argument("--model", default=None, choices=sorted(MODEL_THRESHOLDS),
                   help="take the published threshold of this rE2G model "
                        "instead of the detected one")
    p.add_argument("--model-min-similarity", type=float, default=0.6,
                   help="minimum Jaccard similarity between the feature "
                        "columns of a file and the feature set of a model for "
                        "the model to be recognised (default %(default)s)")
    p.add_argument("--min-distance", type=float, default=10_000,
                   help="a link closer than this to the TSS (bp) is not "
                        "distal (default %(default)s)")
    p.add_argument("--window", type=float, default=WINDOW,
                   help="half-width of the window around the TSS, bp "
                        "(default %(default)s)")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--cores", type=int, default=4,
                   help="link files read in parallel (default %(default)s)")
    p.add_argument("--curve-draws", type=int, default=40,
                   help="random subsets drawn per k for the core-element "
                        "curve; 0 keeps only the exact curve, -1 skips it "
                        "(default %(default)s)")
    p.add_argument("--mito-chroms", default="chrM,MT,M",
                   help="names of the mitochondrial chromosome; genes encoded "
                        "there are left out (default %(default)s)")
    p.add_argument("--set-labels",
                   default="mitochondrial genes,matched controls",
                   help="figure labels of the target and control sets, "
                        "comma separated (default %(default)s)")
    p.add_argument("--seed", type=int, default=0,
                   help="seed of the random subsets (default %(default)s)")
    p.add_argument("--keep", default="",
                   help="regular expression on the CellType of each file "
                        "(case ignored); only matching biosamples are used. "
                        "Empty uses every file. Example for skeletal muscle: "
                        "'skeletal muscle|psoas muscle|gastrocnemius'")
    p.add_argument("--list-biosamples", action="store_true",
                   help="print the CellType of every file and stop")
    a = p.parse_args()

    # every input is checked before the long part starts
    for flag, path in (("--gtf", a.gtf), ("--genes", a.genes),
                       ("--controls", a.controls), ("--ccre", a.ccre),
                       ("--manifest", a.manifest), ("--aliases", a.aliases)):
        if path and not Path(path).exists():
            sys.exit(f"{flag} {path}: no such file "
                     f"(working directory: {Path.cwd()})")

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    d = Path(a.links_dir)
    if not d.is_dir():
        sys.exit(f"--links-dir {d} is not a directory "
                 f"(working directory: {Path.cwd()})")
    patterns = ([a.links_glob] if a.links_glob else
                ["*.bed.gz", "*.bed", "*.tsv.gz", "*.tsv", "*.txt.gz", "*.gz"])
    files, used = [], None
    for pat in patterns:
        files = sorted(d.glob(pat))
        if a.ccre:                  # the registry may sit in the same folder
            ccre_path = Path(a.ccre).resolve()
            files = [f for f in files if f.resolve() != ccre_path]
        if files:
            used = pat
            break
    if not files:
        here = sorted(p.name for p in d.iterdir())
        sys.exit(
            f"nothing matched in {d.resolve()}\n"
            f"  patterns tried: {', '.join(patterns)}\n"
            f"  the directory holds {len(here)} entries"
            + (": " + ", ".join(here[:15]) + ("..." if len(here) > 15 else "")
               if here else " (it is empty)")
            + "\n  pass --links-glob with the pattern that matches, or point "
              "--links-dir at the right folder.")
    print(f"{len(files)} files matching {used!r} in {d.resolve()}")

    celltype = {f: first_data_row(f) for f in files}
    if a.list_biosamples or not a.keep:
        tally = defaultdict(list)
        for f, c in celltype.items():
            tally[c or "(no CellType column)"].append(f.name.split(".")[0])
        print(f"\n{len(tally)} distinct CellType values:")
        for c, fs in sorted(tally.items(), key=lambda x: (-len(x[1]), x[0])):
            print(f"  {len(fs):3d}  {c}")
        with open(out / "biosamples_found.tsv", "w") as fh:
            fh.write("file\tCellType\n")
            for f in files:
                fh.write(f"{f.name}\t{celltype[f]}\n")
        print(f"  (listed in {out / 'biosamples_found.tsv'})\n")
        if a.list_biosamples:
            return
        print("No --keep given: every file is used. If the folder mixes "
              "biosample types, rerun with --keep.\n")
    model, auto = detect_model(files[0], a.model_min_similarity)
    if a.threshold is None and a.model:
        a.threshold = MODEL_THRESHOLDS[a.model][0]
        print(f"--model {a.model}: threshold {a.threshold}"
              + (f" (feature columns point to {model})"
                 if model and model != a.model else ""))
    elif a.threshold is None:
        if auto is None:
            sys.exit("cannot tell which rE2G model produced these files; "
                     "give --threshold explicitly (ENCODE_rE2G models/*/ "
                     "carry theirs in a file named threshold_<value>)")
        a.threshold = auto
        print(f"model recognised from the feature columns: {model}, "
              f"threshold {auto}")
    elif auto is not None and abs(a.threshold - auto) > 1e-9:
        print(f"note: the feature columns point to the {model} model, whose "
              f"threshold is {auto}, but --threshold {a.threshold} was given")

    if a.keep:
        import re
        rx = re.compile(a.keep, re.I)
        files = [f for f in files if rx.search(celltype[f] or "")]
        kept_types = sorted({celltype[f] for f in files})
        print(f"--keep {a.keep!r}: {len(files)} files, "
              f"{len(kept_types)} CellType values")
        for c in kept_types:
            print(f"    {c}")
        if not files:
            sys.exit("--keep matched nothing; see --list-biosamples")
    if a.manifest:
        rows = []
        with open(a.manifest) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                if line.strip():
                    f = line.rstrip("\n").split("\t")
                    rows.append(dict(zip(head, f + [""] * (len(head) - len(f)))))
        key = a.manifest_column or next(
            (c for c in head if "file" in c.lower() or "accession" in c.lower()),
            head[0])
        wanted_files = {r[key].split("/")[-1].split(".")[0] for r in rows}
        files = [f for f in files if f.name.split(".")[0] in wanted_files]
        types = {r.get("biosample_type") or r.get("biosample") or ""
                 for r in rows}
        n_types = len([t for t in types if t]) or 0
        print(f"manifest {a.manifest}: column {key!r}, "
              f"{len(wanted_files)} biosamples, {len(files)} files found")
        if not files:
            sys.exit(f"the manifest filter removed every file. Column {key!r} "
                     f"holds e.g. {sorted(wanted_files)[:3]}, while the "
                     f"directory holds e.g. "
                     f"{[f.name for f in sorted(d.glob(used))][:3]}. Give "
                     "--manifest-column the column that carries the file name.")
    else:
        n_types = len({celltype[f] for f in files if celltype[f]})
    mito = {g for g in Path(a.genes).read_text().split()}
    ctrl = ({g for g in Path(a.controls).read_text().split()} - mito
            if a.controls else set())
    wanted = mito | ctrl
    print(f"{len(mito)} target genes"
          + (f", {len(ctrl)} controls" if ctrl else ""))

    names = gtf_gene_names(Path(a.gtf))
    aliases = build_aliases(a.aliases, builtin=not a.no_builtin_aliases)
    split = lambda x: [t for t in x.split(",") if t]  # noqa: E731
    mapping, unresolved, not_genes = resolve_symbols(
        sorted(wanted), names, aliases, split(a.not_gene_substrings),
        split(a.not_gene_prefixes))
    if not_genes:
        print(f"  {len(not_genes)} entries left out, not genes "
              f"({', '.join(not_genes[:4])})")
    renamed = {k: v for k, v in mapping.items() if k != v}
    if renamed:
        print(f"  {len(renamed)} symbols renamed to their current form "
              f"(e.g. {', '.join(f'{k}->{v}' for k, v in list(renamed.items())[:4])})")
    # the set membership follows the canonical name
    mito = {mapping[g] for g in mito if g in mapping}
    ctrl = {mapping[g] for g in ctrl if g in mapping} - mito
    ens = {}
    tss = load_tss(Path(a.gtf), mito | ctrl, ens)
    chrm = sorted(g for g, (c, _p) in tss.items() if c in set(a.mito_chroms.split(",")))
    if chrm:
        # rE2G predicts nothing for the mitochondrial genome, so these genes
        # would count as "never linked" for a reason unrelated to annotation
        for g in chrm:
            tss.pop(g)
        mito -= set(chrm); ctrl -= set(chrm)
        print(f"  {len(chrm)} genes on the mitochondrial genome left out "
              f"(e.g. {', '.join(chrm[:4])})")
    print(f"{len(tss)}/{len(wanted)} genes with a TSS in {Path(a.gtf).name}")
    missing = sorted(set(unresolved) | ({mapping[g] for g in mapping} - set(tss)))
    if missing:
        (out / "genes_without_tss.txt").write_text("\n".join(missing) + "\n")
        print(f"  without a TSS, dropped: {len(missing)} "
              f"({', '.join(missing[:6])}{'...' if len(missing) > 6 else ''}); "
              f"full list in {out / 'genes_without_tss.txt'}")
    # rE2G may use either the old or the new symbol: accept both
    wanted = {}
    for old, new in mapping.items():
        if new in tss:
            wanted[old] = new
            wanted[new] = new
    # the Ensembl id is stable across renames: the link files may have been
    # built on an older GENCODE and carry symbols this GTF does not use
    for name, eid in ens.items():
        if name in tss:
            wanted[eid] = name
    print(f"  genes looked up by Ensembl id ({len(ens)}) and by symbol")

    mids = {}
    if a.ccre:
        mids = load_ccre(Path(a.ccre))
        print(f"cCRE registry: {sum(len(v['start']) for v in mids.values()):,} elements "
              f"on {len(mids)} chromosomes")

    print(f"\n{len(files)} link files, score >= {a.threshold}")
    work = partial(scan_file, wanted=wanted, threshold=a.threshold,
                   min_dist=a.min_distance, window=a.window)
    if a.cores > 1:
        with Pool(a.cores) as pool:
            res = pool.map(work, files)
    else:
        res = [work(f) for f in files]

    per_bio = {}
    for name, d, n_rows, n_kept, err in res:
        if err:
            print(f"  [{name}] {err}")
            continue
        per_bio[name] = d
        print(f"  [{name}] {n_rows:,} rows, {n_kept:,} links kept "
              f"on {len(d)} genes", flush=True)
    if not per_bio:
        sys.exit("no usable link file")
    n_bio = len(per_bio)

    # every element onto the registry, once
    cache, n_map, n_tot = {}, 0, 0

    def as_ccre(el):
        nonlocal n_map, n_tot
        if el not in cache:
            if mids:
                cache[el] = map_to_ccre(mids, *el, tol=a.ccre_tolerance)
            else:
                cache[el] = (f"{el[0]}:{el[1]}-{el[2]}", False)
            n_tot += 1
            n_map += cache[el][1]
        return cache[el][0]

    bio_order = list(per_bio)
    # biosample label of each file: its CellType, else the file name
    ct_of = {f.name: (celltype.get(f) or f.name).replace(";", ",").replace("\t", " ")
             for f in files}
    ct_of.update({lab: lab for lab in bio_order if lab not in ct_of})
    elements = []
    masks, mask_gene = [], []          # one bit mask per element, see core_curve
    rows = []
    for gi, g in enumerate(sorted(tss)):
        chrom, pos = tss[g]
        counts, unions, inters = [], set(), None
        prevalence = defaultdict(int)  # element -> biosamples that link it
        bit = {}
        for bi, lab in enumerate(bio_order):
            s = {as_ccre(el) for el in per_bio[lab].get(g, (set(), set()))[1]}
            counts.append(len(s))
            unions |= s
            inters = set(s) if inters is None else (inters & s)
            for e in s:
                prevalence[e] += 1
                bit[e] = bit.get(e, 0) | (1 << bi)
        inters = inters or set()
        for e, w in bit.items():
            masks.append(w)
            mask_gene.append(gi)
        # sorted by position: set iteration order depends on string hashing,
        # which varies between runs
        parsed = []
        for e, n in prevalence.items():
            c_, se = e.split(":")
            s_, e_ = se.split("-")
            parsed.append((c_, int(s_), int(e_), n))
        for c_, s_, e_, n in sorted(parsed):
            w = bit[f"{c_}:{s_}-{e_}"]
            types = sorted({ct_of[lab] for bi, lab in enumerate(bio_order) if w >> bi & 1})
            elements.append((g, "mito" if g in mito else "control", c_,
                             s_, e_, (s_ + e_) // 2 - pos, n, ";".join(types)))
        n_ccre = (ccre_in_window(mids, chrom, pos, a.window) if mids
                  else np.nan)
        never = (1 - min(1.0, len(unions) / n_ccre)) if (
            mids and np.isfinite(n_ccre) and n_ccre > 0) else np.nan
        rows.append({
            "gene": g, "set": "mito" if g in mito else "control",
            "chrom": chrom, "tss": pos,
            "n_ccre_window": n_ccre,
            "links_per_biosample_median": float(np.median(counts)),
            "links_per_biosample_max": int(np.max(counts)),
            "n_union": len(unions), "n_intersection": len(inters),
            "n_biosamples_with_link": int(sum(c > 0 for c in counts)),
            "max_prevalence": max(prevalence.values()) if prevalence else 0,
            "frac_never_linked": never})

    if mids:
        print(f"\n{n_tot:,} distinct rE2G elements linked to these genes, "
              f"{n_map / max(1, n_tot):.1%} mapped onto the registry")
        if n_tot and n_map / n_tot < 0.9:
            print("  !! under 90% mapped: check that the registry and the link "
                  "files are on the same assembly")
    else:
        print("\nno registry: elements compared on their own coordinates, so "
              "the union is inflated and the intersection deflated")
    with open(out / "re2g_elements.tsv", "w") as fh:
        fh.write("gene\tset\tchrom\tstart\tend\tdistance\tprevalence\tcelltypes\n")
        for r in elements:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"{len(elements):,} gene-element pairs -> "
          f"{out / 're2g_elements.tsv'}")

    cols = list(rows[0])
    with open(out / "re2g_per_gene.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                "NA" if isinstance(r[c], float) and not np.isfinite(r[c])
                else (f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]))
                for c in cols) + "\n")

    summaries = [s for s in (summarise(rows, "mito"), summarise(rows, "control"))
                 if s is not None]
    scols = list(summaries[0])
    with open(out / "re2g_summary.tsv", "w") as fh:
        fh.write("\t".join(scols) + "\n")
        for s in summaries:
            fh.write("\t".join(
                "NA" if isinstance(s[c], float) and not np.isfinite(s[c])
                else (f"{s[c]:.4f}" if isinstance(s[c], float) else str(s[c]))
                for c in scols) + "\n")

    print(f"\n{n_bio} biosamples, {len(rows)} genes\n")
    for s in summaries:
        print(f"  {s['set']}")
        for k in ("ccre_per_gene_median", "links_per_biosample_median",
                  "union_median", "intersection_median",
                  "pct_candidates_never_linked", "pct_genes_without_any_link",
                  "pct_linked_genes_with_empty_intersection"):
            v = s[k]
            print(f"    {k:<44} {'NA' if not np.isfinite(v) else f'{v:,.2f}'}")
        print()

    if a.curve_draws >= 0:
        rng = np.random.default_rng(a.seed)
        curve = core_curve(rows, masks, mask_gene, n_bio, a.curve_draws, rng)
        cols_c = list(curve[0])
        with open(out / "re2g_core_curve.tsv", "w") as fh:
            fh.write("\t".join(cols_c) + "\n")
            for c in curve:
                fh.write("\t".join(
                    f"{c[x]:.4f}" if isinstance(c[x], float) else str(c[x])
                    for x in cols_c) + "\n")
        curve_figure(curve, out / "re2g_core_curve.svg", n_bio)
        key = ("subsampled_mito" if "subsampled_mito" in curve[0]
               else "observed_mito")
        print("\nshare of target genes with one element shared by "
              "k biosamples")
        marks = sorted({1, 2, 4, 8, 16, 24, 32, 40, 48, 56, n_bio}
                       & set(range(1, n_bio + 1)))
        for kk in marks:
            c = curve[kk - 1]
            line = f"  k = {kk:>3}   {c[key]:.3f}"
            if "observed_control" in c:
                ck = ("subsampled_control" if "subsampled_control" in c
                      else "observed_control")
                line += f"   control {c[ck]:.3f}"
            print(line)
        if n_bio > 8:
            lo, hi = curve[n_bio // 2 - 1][key], curve[n_bio - 1][key]
            slope = (lo - hi) / max(1, n_bio - n_bio // 2)
            print(f"  slope over the upper half: {slope:+.4f} per biosample "
                  f"({lo:.3f} at k = {n_bio//2} to {hi:.3f} at k = {n_bio})")
            print("  A slope near zero means the value at k = "
                  f"{n_bio} is a plateau and would survive more biosamples; a "
                  "slope that stays steep means it is a point on a decline.")

    txt = sentence(summarise(rows, "mito"), n_bio, n_types or "[X]",
                   a.threshold, a.window)
    (out / "re2g_sentence.txt").write_text(txt)
    print(txt)
    figure(rows, out / "re2g_elements_per_gene.svg", n_bio, a.window,
           labels=tuple(x.strip() for x in a.set_labels.split(",", 1)))
    print(f"Written to {out}/")


if __name__ == "__main__":
    main()
