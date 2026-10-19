#!/usr/bin/env python3
"""
promoter_features.py - sequence, architecture and DNA methylation of the
promoters of a gene set

Purpose
    Describe the promoters of a target gene set (nuclear-encoded mitochondrial
    genes) against control genes, from public resources only:

      1. sequence: GC content, CpG observed/expected ratio, overlap with a CpG
         island, promoter class (HCP, ICP, LCP) and occurrences of the
         motifs of the factors that drive broadly expressed and mitochondrial
         promoters (NRF1, ERR, GABPA/ETS, YY1, SP/KLF, CREB, NF-Y, CTCF,
         TATA);
      2. architecture: bidirectional promoters (another TSS on the opposite
         strand nearby) and, with --loops, the loops whose partner anchor
         holds another promoter (promoter-promoter contacts);
      3. DNA methylation, from whole-genome bisulfite (or long-read) tables
         of public tissues and cell types (ENCODE, Roadmap, the human
         methylome atlas): methylation of each promoter and, with
         --elements, of each replicable ENCODE-rE2G element, per sample and
         per group of samples; promoters whose methylation varies between
         groups; and, with --tpm, the relation between promoter methylation
         and expression across GTEx tissues.

    The question is whether mitochondrial promoters belong to the class of
    CpG-rich, constitutively unmethylated, motif-driven promoters that are
    intrinsically less responsive to distal enhancers, which would place
    their regulatory dependency on proximal and structural elements.

Inputs
    --gtf         GENCODE annotation (GTF, gzipped or not). The TSS of a gene
                  is that of its MANE Select transcript, else its gene record.
    --genes       target gene symbols, whitespace-separated. Symbols are
                  resolved as in re2g_panel.py (built-in HGNC alias table,
                  --aliases, --no-builtin-aliases).
    --controls    control gene symbols (make_controls.py). Optional; with
                  --controls-all, every gene of --gene-type not in --genes
                  is a control.
    --fasta       genome FASTA of the same assembly, uncompressed with a
                  .fai index (built if absent), or bgzipped when pysam is
                  installed. Needed for the sequence features.
    --cpg-islands CpG island BED (UCSC cpgIslandExt or cpgIslandExtUnmasked:
                  the first three columns are read; the UCSC table with a
                  leading bin column is recognised). Optional.
    --motifs      optional TSV, name and regular expression on the forward
                  strand (IUPAC written out as character classes), replacing
                  the built-in table MOTIFS.
    --loops       optional loops.tsv of v4c_loops.R (gene, chrom, vp_pos,
                  summit, called).
    --meth-table  optional TSV, one row per methylation file: sample, group
                  (tissue or cell type; samples of a group are pooled into
                  its mean), path, and optionally format and gtex_tissue
                  (the GTEx tissue matching the group, for --tpm).
                  Formats, guessed from the number of columns when absent:
                    bedmethyl  ENCODE WGBS bedMethyl or modkit pileup: chrom,
                               start, end, name, score, strand, ...,
                               coverage (10th), percent methylated (11th);
                               modkit rows other than 'm' are dropped
                    bismark    Bismark coverage file: chrom, start (1-based),
                               end, percent, methylated, unmethylated
                    bedgraph   chrom, start, end, methylation (fraction or
                               percent), optional coverage (5th); bigWig
                               files are converted first with
                               bigWigToBedGraph
    --elements    optional re2g_elements.tsv of re2g_panel.py; elements
                  replicable in >= --min-prevalence biosamples get their
                  methylation as well.
    --tpm         optional GTEx median TPM table (GCT, gzipped or not).
    --annotation  NAME=PATH[:MODE][:WINDOW], repeatable: experimental
                  annotation of promoters as a BED-like file (chrom, start,
                  end, any columns; gzipped or not), over the promoter, or
                  over TSS +/- WINDOW bp when given (0: the TSS base).
                  MODE:
                    overlap      (default) exp_NAME = 1 when an interval
                                 overlaps (ChIP-seq peaks of ReMap or
                                 ENCODE, EPDnew promoters)
                    width        expw_NAME = width of the overlapping
                                 interval nearest the TSS (FANTOM5 CAGE
                                 peaks: broad or sharp promoters)
                    label=REGEX  exp_NAME = 1 when an overlapping interval
                                 has a column beyond the third matching
                                 REGEX (ChromHMM state 'TssA', SCREEN cCRE
                                 class 'PLS')

Outputs (in --out)
    promoter_features.tsv        one row per gene: gene, set, chrom, tss,
                                 strand, promoter start and end, then the
                                 sequence features (gc, cpg_oe, cpg_island,
                                 promoter_class, motif_<name> counts),
                                 bidirectional, bidirectional_partner,
                                 bidirectional_distance, and with --loops
                                 n_loops, n_promoter_loops, loop_partners
    promoter_feature_tests.tsv   one row per feature: target and control
                                 values (share or median), odds ratio or
                                 probability of superiority, Fisher or
                                 Mann-Whitney p, and the p of a permutation
                                 stratified on the promoter class
    methylation_long.tsv         region, gene, set, kind (promoter or
                                 element), sample, group, n_cpg, coverage,
                                 beta
    methylation_by_group.tsv     promoters: gene, set, one column per group
                                 (mean beta of its samples)
    elements_methylation_by_group.tsv
                                 the same for the elements (region, gene,
                                 set, distance, prevalence, groups)
    methylation_summary.tsv      per kind, group and set: regions with data,
                                 median beta, share unmethylated, share
                                 methylated; target against control:
                                 probability of superiority and p
    methylation_variable.tsv     promoters whose group means span at least
                                 --variable-range: gene, set, min and max
                                 group, range, and with --focus the mean in
                                 the focus groups, in the others, and the
                                 difference
    methylation_expression.tsv   with --tpm: per gene whose promoter
                                 methylation spans at least --min-range
                                 across GTEx-matched groups, the Spearman
                                 correlation between methylation and
                                 log2(TPM + 1)
    genes_without_tss.txt        symbols left out for lack of a TSS

Method
    Promoter. --upstream bp upstream to --downstream bp downstream of the
    TSS, oriented on the gene strand.

    CpG features. GC is (C + G) / non-N bases; the CpG o/e ratio is
    CpG x length / (C x G) (Gardiner-Garden and Frommer). Promoter class as
    in Weber et al. 2007, on 500 bp windows sliding by 50 bp over the
    promoter: HCP when one window has CpG o/e > 0.75 and GC > 0.55, LCP when
    no window has CpG o/e > 0.48, ICP otherwise.

    Motifs. Each regular expression is searched on both strands, overlapping
    occurrences included; an occurrence found on both strands (palindrome)
    counts once. Consensus matches are a coarse proxy of binding; the class
    stratification below keeps CpG-rich motifs from reading as enrichment
    only because the target promoters are CpG-rich.

    Bidirectional promoters. Another gene of --partner-types with its TSS on
    the opposite strand within --bidirectional-distance (Trinklein et al.
    2004 used 1 kb).

    Promoter-promoter loops. A loop of a gene (loops.tsv) whose partner
    summit lies within --anchor-half-width of the TSS of another gene of
    --partner-types.

    Tests. Binary features: Fisher exact test and odds ratio; continuous
    ones: Mann-Whitney test and probability of superiority (P(target >
    control), 0.5 = interchangeable). For every feature, the difference of
    means is also tested by permuting the set labels within promoter
    classes (--perms), which asks whether a feature is enriched beyond what
    the CpG content of the promoters implies.

    Methylation. Per region and sample, beta = methylated calls / covered
    calls over the CpGs covered by at least --min-coverage reads (both
    strands pooled); a region needs --min-cpg such CpGs, else it is missing.
    Without a coverage column (bedGraph), beta is the mean of the CpG values.
    A group's value is the mean of its samples. A promoter is unmethylated
    below --unmethylated and methylated above --methylated. With --focus (a
    regular expression on group names, case ignored), the variable
    promoters report the focus groups (for instance skeletal muscle) against
    the others.

Usage
    python3 promoter_features.py --gtf data/gencode.gtf.gz \\
        --genes genes.txt --controls controls.txt \\
        --fasta data/GRCh38.fa --cpg-islands data/cpgIslandExt.txt.gz \\
        --loops results/loops/loops.tsv \\
        --meth-table config/methylation.tsv \\
        --elements results/re2g/re2g_elements.tsv --min-prevalence 51 \\
        --tpm data/gtex_median_tpm.gct.gz --focus "muscle" \\
        --out results/promoters

References
    Gardiner-Garden M, Frommer M. 1987. CpG islands in vertebrate genomes.
        J Mol Biol 196:261-282.
    Weber M et al. 2007. Distribution, silencing potential and evolutionary
        impact of promoter DNA methylation in the human genome. Nat Genet
        39:457-466.
    Trinklein ND et al. 2004. An abundance of bidirectional promoters in the
        human genome. Genome Res 14:62-66.
    Bergman DT et al. 2022. Compatibility rules of human enhancer and
        promoter sequences. Nature 607:176-184.
    Scarpulla RC. 2008. Transcriptional paradigms in mammalian mitochondrial
        biogenesis and function. Physiol Rev 88:611-638.
    Castro-Mondragon JA et al. 2022. JASPAR 2022. Nucleic Acids Res
        50:D165-D173 (motif consensus).
    ENCODE Project Consortium et al. 2020. Expanded encyclopaedias of DNA
        elements in the human and mouse genomes. Nature 583:699-710.
    Roadmap Epigenomics Consortium. 2015. Integrative analysis of 111
        reference human epigenomes. Nature 518:317-330.
    Loyfer N et al. 2023. A DNA methylation atlas of normal human cell
        types. Nature 613:355-364.
    GTEx Consortium. 2020. Science 369:1318-1330.
    Phipson B, Smyth GK. 2010. Permutation P-values should never be zero.
        Stat Appl Genet Mol Biol 9:39.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from re2g_panel import (NOT_GENE_PREFIXES, NOT_GENES, build_aliases,  # noqa: E402
                        opener, resolve_symbols)

# Consensus of the factors of broadly expressed and mitochondrial promoters,
# forward strand. The reverse strand is searched as well.
MOTIFS = {
    "NRF1": "GCGCA[CT]GCGC",                     # palindrome, NRF1 site
    "ERR": "T[ACGT]AAGGTCA",                     # ERRE, ESRRA/ESRRG
    "GABPA_ETS": "CCGGAAG[CT]",                  # NRF2/GABPA, ETS core
    "YY1": "GCCAT[ACGT]TT",
    "SP_KLF": "GGGGCGGGG",                       # GC box
    "CREB": "TGACGTCA",                          # full CRE
    "NFY": "CCAAT[CG][AG]",                      # CCAAT box
    "CTCF": "CCA[CG][CT]AG[AG][GT]GG[AC][AG]",   # core motif
    "TATA": "TATA[AT]A[AT][AG]",
}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


# ---------------------------------------------------------------------------
# annotation

def load_genes(gtf: Path) -> dict:
    """symbol -> (chrom, tss, strand, ensembl id, gene type).

    MANE Select transcript where available, otherwise the gene record. The
    TSS follows re2g_panel.load_tss: start - 1 on '+', end on '-'.
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
            if "_PAR_Y" in a:
                continue
            i = a.find('gene_name "')
            if i < 0:
                continue
            name = a[i + 11:a.find('"', i + 11)]
            rank = 2 if (f[2] == "transcript" and "MANE_Select" in a) else (
                1 if f[2] == "gene" else 0)
            if rank == 0:
                continue
            prev = best.get(name)
            if prev is not None and rank <= prev[5]:
                continue
            tss = int(f[3]) - 1 if f[6] == "+" else int(f[4])
            j = a.find('gene_id "')
            ensg = a[j + 9:a.find('"', j + 9)].split(".")[0] if j >= 0 else ""
            k = a.find('gene_type "')
            gtype = a[k + 11:a.find('"', k + 11)] if k >= 0 else ""
            best[name] = (f[0], tss, f[6], ensg, gtype, rank)
    return {k: v[:5] for k, v in best.items()}


def promoter_interval(tss, strand, up, down):
    return (tss - up, tss + down) if strand == "+" else (tss - down, tss + up)


def chrom_key(index: dict, chrom: str):
    for c in (chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom):
        if c in index:
            return c
    return None


# ---------------------------------------------------------------------------
# sequence

class Fasta:
    """Random access to a FASTA: pysam when installed, else a .fai reader."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.pysam = None
        try:
            import pysam
            self.pysam = pysam.FastaFile(str(self.path))
            self.names = set(self.pysam.references)
            return
        except ImportError:
            pass
        if str(self.path).endswith(".gz"):
            sys.exit("a gzipped FASTA needs pysam (bgzip + samtools faidx); "
                     "or give the uncompressed file")
        fai = Path(str(self.path) + ".fai")
        if not fai.exists():
            print(f"indexing {self.path}")
            build_fai(self.path, fai)
        self.index = {}
        with open(fai) as fh:
            for line in fh:
                n, ln, off, lb, lw = line.split("\t")[:5]
                self.index[n] = (int(ln), int(off), int(lb), int(lw))
        self.names = set(self.index)
        self.fh = open(self.path, "rb")

    def fetch(self, chrom, start, end) -> str:
        c = chrom_key({n: 1 for n in self.names}, chrom)
        if c is None:
            return ""
        if self.pysam is not None:
            return self.pysam.fetch(c, max(0, start), end).upper()
        ln, off, lb, lw = self.index[c]
        start, end = max(0, start), min(end, ln)
        if end <= start:
            return ""
        a = off + (start // lb) * lw + start % lb
        b = off + (end // lb) * lw + end % lb
        self.fh.seek(a)
        return self.fh.read(b - a).decode().replace("\n", "").replace(
            "\r", "").upper()


def build_fai(fasta: Path, fai: Path) -> None:
    rows = []
    with open(fasta, "rb") as fh:
        name, length, off, lb, lw = None, 0, 0, 0, 0
        pos = 0
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    rows.append((name, length, off, lb, lw))
                name = line[1:].split()[0].decode()
                length, lb, lw = 0, 0, 0
                off = pos + len(line)
            else:
                if lb == 0:
                    lb, lw = len(line.rstrip(b"\r\n")), len(line)
                length += len(line.rstrip(b"\r\n"))
            pos += len(line)
        if name is not None:
            rows.append((name, length, off, lb, lw))
    with open(fai, "w") as out:
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")


def gc_cpg(seq: str):
    s = seq.replace("N", "")
    n = len(s)
    if n == 0:
        return np.nan, np.nan
    c, g, cg = s.count("C"), s.count("G"), s.count("CG")
    gc = (c + g) / n
    oe = cg * n / (c * g) if c and g else 0.0
    return gc, oe


def promoter_class(seq: str, win=500, step=50) -> str:
    """Weber et al. 2007: HCP, ICP or LCP from 500 bp windows."""
    if len(seq.replace("N", "")) < win // 2:
        return "NA"
    starts = range(0, max(1, len(seq) - win + 1), step)
    hcp = icp_ok = False
    for s in starts:
        gc, oe = gc_cpg(seq[s:s + win])
        if np.isnan(gc):
            continue
        if oe > 0.75 and gc > 0.55:
            hcp = True
        if oe > 0.48:
            icp_ok = True
    return "HCP" if hcp else ("ICP" if icp_ok else "LCP")


def count_motif(seq: str, rx: re.Pattern, rx_rc: re.Pattern) -> int:
    """Occurrences on both strands; a palindrome counts once."""
    hits = {(m.start(), m.start() + len(m.group(1))) for m in rx.finditer(seq)}
    hits |= {(m.start(), m.start() + len(m.group(1)))
             for m in rx_rc.finditer(seq)}
    return len(hits)


def rc_regex(expr: str) -> str:
    """Reverse complement of a regular expression made of bases and
    character classes."""
    tokens = re.findall(r"\[[^\]]+\]|[ACGTN.]", expr.upper())
    if "".join(tokens) != expr.upper():
        sys.exit(f"motif {expr!r}: only bases and [..] classes are supported")
    out = []
    for t in reversed(tokens):
        if t.startswith("["):
            out.append("[" + "".join(sorted(revcomp(t[1:-1]))) + "]")
        elif t == ".":
            out.append(".")
        else:
            out.append(revcomp(t))
    return "".join(out)


def load_intervals(path: Path) -> dict:
    """BED-like intervals per chromosome, sorted; the UCSC table with a
    leading bin column is recognised."""
    iv = defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[0].isdigit() and not f[1].isdigit():
                f = f[1:]
            try:
                iv[f[0]].append((int(f[1]), int(f[2])))
            except (IndexError, ValueError):
                continue
    out = {}
    for c, v in iv.items():
        v.sort()
        out[c] = (np.array([x[0] for x in v]), np.array([x[1] for x in v]))
    return out


def overlaps(iv: dict, chrom, s, e) -> bool:
    c = chrom_key(iv, chrom)
    if c is None:
        return False
    st, en = iv[c]
    i = int(np.searchsorted(st, e, "left"))
    # intervals starting before e; CpG islands do not nest, so the ends of
    # the few before e are enough
    for k in range(i - 1, max(-1, i - 50), -1):
        if en[k] > s:
            return True
        if en[k] < s - 100_000:
            break
    return False


# ---------------------------------------------------------------------------
# experimental annotations (--annotation)

def load_labelled(path: Path) -> dict:
    """chrom -> (starts, ends, rest of the line), sorted by start"""
    iv = defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[0].isdigit() and not f[1].isdigit():
                f = f[1:]                               # UCSC bin column
            try:
                iv[f[0]].append((int(f[1]), int(f[2]), "\t".join(f[3:])))
            except (IndexError, ValueError):
                continue
    out = {}
    for c, v in iv.items():
        v.sort()
        st = np.array([x[0] for x in v])
        en = np.array([x[1] for x in v])
        out[c] = (st, en, [x[2] for x in v], int((en - st).max()))
    return out


def hits(iv: dict, chrom, s, e) -> list:
    """indices of the intervals overlapping [s, e) (nested ones included)"""
    c = chrom_key(iv, chrom)
    if c is None:
        return [], None
    st, en, _rest, maxlen = iv[c]
    lo = int(np.searchsorted(st, s - maxlen, "left"))
    hi = int(np.searchsorted(st, e, "left"))
    k = np.arange(lo, hi)
    return k[en[lo:hi] > s], c


def annotate(spec: str, rows: list):
    m = re.match(r"^([A-Za-z0-9_.-]+)=(.+?)(?::(overlap|width|label=[^:]+))?"
                 r"(?::(\d+))?$", spec)
    if not m:
        sys.exit(f"--annotation {spec!r}: expected NAME=PATH[:MODE][:WINDOW]")
    name, path, mode, win = m.groups()
    mode = mode or "overlap"
    path = Path(path)
    if not path.exists():
        print(f"--annotation {name}: {path} absent, skipped")
        return [], []
    iv = load_labelled(path)
    rx = re.compile(mode[6:], re.I) if mode.startswith("label=") else None
    col = f"expw_{name}" if mode == "width" else f"exp_{name}"
    n = 0
    for r in rows:
        if win is not None:
            s, e = r["tss"] - int(win), r["tss"] + int(win) + 1
        else:
            s, e = r["prom_start"], r["prom_end"]
        k, c = hits(iv, r["chrom"], s, e)
        if mode == "width":
            if len(k):
                st, en = iv[c][0][k], iv[c][1][k]
                mid = (st + en) / 2
                j = int(np.argmin(np.abs(mid - r["tss"])))
                r[col] = int(en[j] - st[j])
                n += 1
            else:
                r[col] = np.nan
        elif rx is not None:
            r[col] = int(any(rx.search(x) for j in k
                             for x in iv[c][2][j].split("\t"))) if len(k) else 0
            n += r[col]
        else:
            r[col] = int(len(k) > 0)
            n += r[col]
    print(f"annotation {name} ({mode}, "
          f"{'promoter' if win is None else f'TSS +/- {win} bp'}): "
          f"{n} of {len(rows)} promoters")
    return ([], [col]) if mode == "width" else ([col], [])


# ---------------------------------------------------------------------------
# methylation

class RegionIndex:
    """Regions per chromosome for vectorised point-in-interval assignment."""

    def __init__(self, regions):
        # regions: list of (chrom, start, end); the index of a region is its
        # position in the list
        per = defaultdict(list)
        for i, (c, s, e) in enumerate(regions):
            per[c].append((s, e, i))
        self.idx = {}
        for c, v in per.items():
            v.sort()
            st = np.array([x[0] for x in v], dtype=np.int64)
            en = np.array([x[1] for x in v], dtype=np.int64)
            ri = np.array([x[2] for x in v], dtype=np.int64)
            self.idx[c] = (st, en, ri, int((en - st).max()))
        self.n = len(regions)

    def assign(self, chrom, pos):
        """(positions index, region index) for every CpG inside a region."""
        c = chrom_key(self.idx, chrom)
        if c is None:
            return np.empty(0, np.int64), np.empty(0, np.int64)
        st, en, ri, maxlen = self.idx[c]
        lo = np.searchsorted(st, pos - maxlen, "right")
        hi = np.searchsorted(st, pos, "right")
        span = hi - lo
        if span.size == 0 or span.max() == 0:
            return np.empty(0, np.int64), np.empty(0, np.int64)
        pi, rr = [], []
        for k in range(int(span.max())):
            j = lo + k
            ok = np.nonzero(j < hi)[0]
            jj = j[ok]
            hit = en[jj] > pos[ok]
            pi.append(ok[hit])
            rr.append(ri[jj[hit]])
        return np.concatenate(pi), np.concatenate(rr)


def guess_format(path: Path) -> str:
    with opener(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            n = len(line.rstrip("\n").split("\t"))
            break
        else:
            sys.exit(f"{path}: empty")
    if n >= 11:
        return "bedmethyl"
    if n == 6:
        return "bismark"
    if n in (4, 5):
        return "bedgraph"
    sys.exit(f"{path}: {n} columns, format not recognised; give it in the "
             "format column of --meth-table")


def header_lines(path: Path) -> int:
    k = 0
    with opener(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                k += 1
            else:
                break
    return k


def read_methylation(path: Path, fmt: str, index: RegionIndex, min_cov: int,
                     chunksize: int = 2_000_000):
    """Methylated and covered calls and CpG counts per region."""
    import pandas as pd
    msum = np.zeros(index.n)
    csum = np.zeros(index.n)
    ncpg = np.zeros(index.n, dtype=np.int64)
    skip = header_lines(path)
    if fmt == "bedmethyl":
        cols, names = [0, 1, 3, 9, 10], ["chrom", "pos", "code", "cov", "pct"]
    elif fmt == "bismark":
        cols, names = [0, 1, 4, 5], ["chrom", "pos", "m", "u"]
    elif fmt == "bedgraph":
        ncol = len(pd.read_csv(path, sep="\t", header=None, skiprows=skip,
                               nrows=1).columns)
        cols = [0, 1, 3] + ([4] if ncol >= 5 else [])
        names = ["chrom", "pos", "value"] + (["cov"] if ncol >= 5 else [])
    else:
        sys.exit(f"unknown methylation format {fmt!r}")
    scale = None
    reader = pd.read_csv(path, sep="\t", header=None, skiprows=skip,
                         usecols=cols, names=None, chunksize=chunksize,
                         dtype={0: str, 3: str} if fmt == "bedmethyl"
                         else {0: str})
    for ch in reader:
        ch.columns = names
        if fmt == "bedmethyl":
            codes = ch["code"].astype(str)
            if codes.isin(["m", "h", "a"]).any():
                ch = ch[codes == "m"]
            cov = ch["cov"].to_numpy(float)
            meth = cov * ch["pct"].to_numpy(float) / 100.0
            pos = ch["pos"].to_numpy(np.int64)
        elif fmt == "bismark":
            m = ch["m"].to_numpy(float)
            cov = m + ch["u"].to_numpy(float)
            meth = m
            pos = ch["pos"].to_numpy(np.int64) - 1
        else:
            v = ch["value"].to_numpy(float)
            if scale is None:
                scale = 100.0 if np.nanmax(v) > 1.0 else 1.0
            v = v / scale
            if "cov" in ch:
                cov = ch["cov"].to_numpy(float)
                meth = v * cov
            else:
                cov = np.ones_like(v)
                meth = v
            pos = ch["pos"].to_numpy(np.int64)
        keep = cov >= (min_cov if fmt != "bedgraph" or "cov" in ch else 0)
        keep &= cov > 0
        chroms = ch["chrom"].to_numpy()[keep]
        pos, meth, cov = pos[keep], meth[keep], cov[keep]
        if pos.size == 0:
            continue
        # chunks are sorted by chromosome in practice, but grouping is safe
        order = np.argsort(chroms, kind="stable")
        chroms, pos, meth, cov = chroms[order], pos[order], meth[order], \
            cov[order]
        bounds = np.flatnonzero(chroms[1:] != chroms[:-1]) + 1
        for a, b in zip(np.r_[0, bounds], np.r_[bounds, len(chroms)]):
            pi, rr = index.assign(str(chroms[a]), pos[a:b])
            if pi.size == 0:
                continue
            np.add.at(msum, rr, meth[a:b][pi])
            np.add.at(csum, rr, cov[a:b][pi])
            np.add.at(ncpg, rr, 1)
    return msum, csum, ncpg


# ---------------------------------------------------------------------------
# statistics

def prob_superiority(a, b) -> float:
    """P(a > b) + 0.5 P(a = b), as in gtex_eqtl.py."""
    from scipy.stats import rankdata
    n, m = len(a), len(b)
    if not n or not m:
        return float("nan")
    r = rankdata(np.concatenate([a, b]))
    return float((r[:n].sum() - n * (n + 1) / 2) / (n * m))


def strat_perm_p(values, is_target, strata, perms, rng) -> float:
    """Two-sided p of the difference of means, labels permuted within
    strata."""
    v = np.asarray(values, float)
    t = np.asarray(is_target, bool)
    ok = ~np.isnan(v)
    v, t, s = v[ok], t[ok], np.asarray(strata)[ok]
    if t.sum() == 0 or (~t).sum() == 0 or perms <= 0:
        return float("nan")
    obs = v[t].mean() - v[~t].mean()
    groups = [np.flatnonzero(s == k) for k in np.unique(s)]
    nt = t.sum()
    tot = v.sum()
    hits = 0
    lab = t.copy()
    for _ in range(perms):
        for g in groups:
            lab[g] = rng.permutation(lab[g])
        st = v[lab].sum()
        d = st / nt - (tot - st) / (len(v) - nt)
        hits += abs(d) >= abs(obs) - 1e-12
    return (1 + hits) / (perms + 1)


def fmt(x):
    if isinstance(x, (float, np.floating)):
        return "NA" if np.isnan(x) else f"{x:.4g}"
    return str(x)


def write_tsv(path: Path, rows: list, cols: list | None = None) -> None:
    cols = cols or (list(rows[0]) if rows else [])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r.get(c, "")) for c in cols) + "\n")


def read_gene_list(path) -> list:
    return [g for g in Path(path).read_text().split() if g]


def read_table(path) -> list:
    with opener(path) as fh:
        head = fh.readline().rstrip("\n").lstrip("#").split("\t")
        rows = []
        for line in fh:
            if line.strip():
                f = line.rstrip("\n").split("\t")
                rows.append(dict(zip(head, f + [""] * (len(head) - len(f)))))
    return rows


def load_tpm(path: Path):
    """GTEx GCT: (Ensembl id -> row, symbol -> row, tissue columns)."""
    with opener(path) as fh:
        fh.readline()
        fh.readline()
        head = fh.readline().rstrip("\n").split("\t")
        by_id, by_sym = {}, {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            vals = np.array([float(x) if x else np.nan for x in f[2:]])
            by_id[f[0].split(".")[0]] = vals
            by_sym.setdefault(f[1], vals)
    return by_id, by_sym, head[2:]


# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtf", required=True, type=Path)
    p.add_argument("--genes", required=True, type=Path)
    p.add_argument("--controls", type=Path)
    p.add_argument("--controls-all", action="store_true",
                   help="every gene of --gene-type outside --genes is a "
                        "control")
    p.add_argument("--gene-type", default="protein_coding",
                   help="gene type of --controls-all (default %(default)s)")
    p.add_argument("--aliases", type=Path)
    p.add_argument("--no-builtin-aliases", action="store_true")
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--upstream", type=int, default=1000,
                   help="promoter length upstream of the TSS, bp "
                        "(default %(default)s)")
    p.add_argument("--downstream", type=int, default=500,
                   help="promoter length downstream of the TSS, bp "
                        "(default %(default)s)")
    p.add_argument("--fasta", type=Path)
    p.add_argument("--cpg-islands", type=Path)
    p.add_argument("--motifs", type=Path)
    p.add_argument("--partner-types", default="protein_coding,lncRNA",
                   help="gene types that count as a partner promoter "
                        "(bidirectional pairs, promoter-promoter loops; "
                        "default %(default)s)")
    p.add_argument("--bidirectional-distance", type=int, default=1000,
                   help="maximal distance between the two TSS of a "
                        "bidirectional pair, bp (default %(default)s)")
    p.add_argument("--loops", type=Path)
    p.add_argument("--anchor-half-width", type=int, default=5000,
                   help="half-width of a loop anchor around its summit, bp; "
                        "half a 10 kb bin (default %(default)s)")
    p.add_argument("--meth-table", type=Path)
    p.add_argument("--min-coverage", type=int, default=5,
                   help="reads covering a CpG for it to count "
                        "(default %(default)s)")
    p.add_argument("--min-cpg", type=int, default=3,
                   help="covered CpGs a region needs (default %(default)s)")
    p.add_argument("--unmethylated", type=float, default=0.2)
    p.add_argument("--methylated", type=float, default=0.6)
    p.add_argument("--variable-range", type=float, default=0.3,
                   help="span of group means for a promoter to be reported "
                        "as variable (default %(default)s)")
    p.add_argument("--focus", default="",
                   help="regular expression on group names, case ignored: "
                        "the groups compared with the others in "
                        "methylation_variable.tsv")
    p.add_argument("--elements", type=Path)
    p.add_argument("--min-prevalence", type=int, default=1,
                   help="biosamples an element must be linked in "
                        "(default %(default)s)")
    p.add_argument("--tpm", type=Path)
    p.add_argument("--min-tissues", type=int, default=5,
                   help="GTEx-matched groups needed for a correlation "
                        "(default %(default)s)")
    p.add_argument("--min-range", type=float, default=0.1,
                   help="span of methylation across the GTEx-matched groups "
                        "below which no correlation is computed "
                        "(default %(default)s)")
    p.add_argument("--annotation", action="append", default=[],
                   help="NAME=PATH[:MODE][:WINDOW], repeatable")
    p.add_argument("--perms", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()

    a.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    print("reading the annotation")
    allg = load_genes(a.gtf)
    names = set(allg)
    aliases = build_aliases(a.aliases, not a.no_builtin_aliases)

    def resolve(symbols):
        m, unresolved, _ = resolve_symbols(symbols, names, aliases,
                                           NOT_GENES, NOT_GENE_PREFIXES)
        return sorted(set(m.values())), unresolved

    target, miss = resolve(read_gene_list(a.genes))
    control, miss_c = [], []
    if a.controls:
        control, miss_c = resolve(read_gene_list(a.controls))
    elif a.controls_all:
        ts = set(target)
        control = sorted(g for g, v in allg.items()
                         if v[4] == a.gene_type and g not in ts
                         and not v[0].startswith("chrM"))
    # genes on the mitochondrial genome have no nuclear promoter
    chrm = [g for g in target + control if allg[g][0] in ("chrM", "MT", "M")]
    target = [g for g in target if g not in chrm]
    control = [g for g in control if g not in chrm and g not in set(target)]
    (a.out / "genes_without_tss.txt").write_text(
        "".join(f"{g}\n" for g in miss + miss_c + chrm))
    print(f"  {len(target)} target genes, {len(control)} control genes "
          f"({len(miss) + len(miss_c) + len(chrm)} left out)")
    genes = [(g, "mito") for g in target] + [(g, "control") for g in control]

    rows = []
    for g, s in genes:
        c, tss, strand, ensg, gtype = allg[g]
        ps, pe = promoter_interval(tss, strand, a.upstream, a.downstream)
        rows.append({"gene": g, "gene_id": ensg, "set": s, "chrom": c,
                     "tss": tss,
                     "strand": strand, "prom_start": ps, "prom_end": pe})

    # ------------------------------------------------------------ sequence
    feats_bin, feats_cont = [], []
    if a.fasta:
        print("sequence features")
        motifs = MOTIFS
        if a.motifs:
            motifs = {r[0]: r[1] for r in (
                line.rstrip("\n").split("\t")
                for line in a.motifs.read_text().splitlines()
                if line.strip() and not line.startswith("#"))
                if len(r) >= 2 and r[0] != "name"}
        cmp = {k: (re.compile(f"(?=({v.upper()}))"),
                   re.compile(f"(?=({rc_regex(v)}))"))
               for k, v in motifs.items()}
        fa = Fasta(a.fasta)
        for r in rows:
            seq = fa.fetch(r["chrom"], r["prom_start"], r["prom_end"])
            r["gc"], r["cpg_oe"] = gc_cpg(seq)
            r["promoter_class"] = promoter_class(seq)
            for k, (rx, rxr) in cmp.items():
                r[f"motif_{k}"] = count_motif(seq, rx, rxr) if seq else np.nan
        feats_cont += ["gc", "cpg_oe"]
        feats_bin += [f"motif_{k}" for k in motifs]
    if a.cpg_islands:
        cgi = load_intervals(a.cpg_islands)
        for r in rows:
            r["cpg_island"] = int(overlaps(cgi, r["chrom"], r["prom_start"],
                                           r["prom_end"]))
        feats_bin.append("cpg_island")

    # -------------------------------------------------------- architecture
    ptypes = set(a.partner_types.split(","))
    tss_by_chrom = defaultdict(list)
    for g, (c, t, st, _e, gt) in allg.items():
        if gt in ptypes:
            tss_by_chrom[c].append((t, st, g))
    for c in tss_by_chrom:
        tss_by_chrom[c].sort()
    tss_pos = {c: np.array([x[0] for x in v]) for c, v in tss_by_chrom.items()}
    for r in rows:
        v, pos = tss_by_chrom.get(r["chrom"], []), tss_pos.get(r["chrom"])
        best = None
        if v:
            lo = np.searchsorted(pos, r["tss"] - a.bidirectional_distance)
            hi = np.searchsorted(pos, r["tss"] + a.bidirectional_distance,
                                 "right")
            for t, st, g in v[lo:hi]:
                if g != r["gene"] and st != r["strand"]:
                    d = abs(t - r["tss"])
                    if best is None or d < best[0]:
                        best = (d, g)
        r["bidirectional"] = int(best is not None)
        r["bidirectional_partner"] = best[1] if best else ""
        r["bidirectional_distance"] = best[0] if best else ""
    feats_bin.append("bidirectional")

    if a.loops:
        loops = defaultdict(list)
        for lr in read_table(a.loops):
            try:
                loops[lr["gene"]].append((lr["chrom"], float(lr["summit"])))
            except (KeyError, ValueError):
                continue
        for r in rows:
            L = loops.get(r["gene"], [])
            partners, npp = set(), 0
            for c, sm in L:
                v, pos = tss_by_chrom.get(c, []), tss_pos.get(c)
                if not v:
                    continue
                lo = np.searchsorted(pos, sm - a.anchor_half_width)
                hi = np.searchsorted(pos, sm + a.anchor_half_width, "right")
                hit = {g for _t, _s, g in v[lo:hi] if g != r["gene"]}
                if hit:
                    npp += 1
                    partners |= hit
            r["n_loops"] = len(L)
            r["n_promoter_loops"] = npp
            r["loop_partners"] = ",".join(sorted(partners))
            r["has_promoter_loop"] = int(npp > 0) if L else np.nan
        feats_bin.append("has_promoter_loop")

    # ---------------------------------------- experimental annotations
    for spec in a.annotation:
        b, c = annotate(spec, rows)
        feats_bin += b
        feats_cont += c

    # -------------------------------------------------------------- tests
    tests = []
    if control:
        from scipy.stats import fisher_exact, mannwhitneyu
        is_t = np.array([r["set"] == "mito" for r in rows])
        strata = np.array([r.get("promoter_class", "NA") for r in rows])
        for f in feats_bin + feats_cont:
            v = np.array([float(r[f]) if r.get(f, "") != "" else np.nan
                          for r in rows], float)
            binary = f in feats_bin
            if binary:
                v = np.where(np.isnan(v), np.nan, (v > 0).astype(float))
            vt, vc = v[is_t & ~np.isnan(v)], v[~is_t & ~np.isnan(v)]
            if not len(vt) or not len(vc):
                continue
            t = {"feature": f, "type": "binary" if binary else "continuous",
                 "n_target": len(vt), "n_control": len(vc)}
            if binary:
                k1, k2 = int(vt.sum()), int(vc.sum())
                odds, pv = fisher_exact([[k1, len(vt) - k1],
                                         [k2, len(vc) - k2]])
                t.update(target=k1 / len(vt), control=k2 / len(vc),
                         effect=float(odds), effect_kind="odds_ratio", p=pv)
            else:
                t.update(target=float(np.median(vt)),
                         control=float(np.median(vc)),
                         effect=prob_superiority(vt, vc),
                         effect_kind="P(target>control)",
                         p=float(mannwhitneyu(vt, vc).pvalue))
            t["p_within_class"] = strat_perm_p(v, is_t, strata, a.perms, rng)
            tests.append(t)
        write_tsv(a.out / "promoter_feature_tests.tsv", tests,
                  ["feature", "type", "n_target", "n_control", "target",
                   "control", "effect", "effect_kind", "p", "p_within_class"])
        if "promoter_class" in rows[0]:
            print("\npromoter classes (share of genes)")
            for s in ("mito", "control"):
                sub = [r["promoter_class"] for r in rows if r["set"] == s]
                print(f"  {s:8s}" + "".join(
                    f"  {k} {100 * sub.count(k) / len(sub):5.1f}%"
                    for k in ("HCP", "ICP", "LCP")))
        print(f"\n{'feature':22s}{'target':>9s}{'control':>9s}"
              f"{'effect':>9s}{'p':>10s}{'p|class':>10s}")
        for t in tests:
            print(f"{t['feature']:22s}{t['target']:9.3f}{t['control']:9.3f}"
                  f"{t['effect']:9.2f}{t['p']:10.2g}"
                  f"{t['p_within_class']:10.2g}")

    cols = ["gene", "gene_id", "set", "chrom", "tss", "strand", "prom_start",
            "prom_end"]
    for r in rows:
        for c in r:
            if c not in cols:
                cols.append(c)
    write_tsv(a.out / "promoter_features.tsv", rows, cols)

    # --------------------------------------------------------- methylation
    if a.meth_table:
        methylation(a, rows, rng)
    print(f"\nWritten to {a.out}")


def methylation(a, rows, rng):
    from scipy.stats import mannwhitneyu, spearmanr
    samples = read_table(a.meth_table)
    for s in samples:
        for k in ("sample", "group", "path"):
            if not s.get(k):
                sys.exit(f"--meth-table needs the columns sample, group, path "
                         f"(row {s})")
        s["format"] = s.get("format") or guess_format(Path(s["path"]))
    regions = []                     # (chrom, start, end)
    meta = []                        # region id, gene, set, kind, extra
    for r in rows:
        regions.append((r["chrom"], r["prom_start"], r["prom_end"]))
        meta.append({"region": f"{r['gene']}:promoter", "gene": r["gene"],
                     "set": r["set"], "kind": "promoter"})
    sets = {r["gene"]: r["set"] for r in rows}
    if a.elements:
        for e in read_table(a.elements):
            try:
                if int(float(e["prevalence"])) < a.min_prevalence:
                    continue
                st, en = int(e["start"]), int(e["end"])
            except (KeyError, ValueError):
                continue
            if e["gene"] not in sets:
                continue
            regions.append((e["chrom"], st, en))
            meta.append({"region": f"{e['chrom']}:{st}-{en}",
                         "gene": e["gene"], "set": sets[e["gene"]],
                         "kind": "element", "distance": e.get("distance", ""),
                         "prevalence": e.get("prevalence", "")})
    index = RegionIndex(regions)

    groups = []
    for s in samples:
        if s["group"] not in groups:
            groups.append(s["group"])
    beta = np.full((len(regions), len(samples)), np.nan)
    long_rows = []
    for j, s in enumerate(samples):
        print(f"methylation: {s['sample']} ({s['format']})")
        m, c, n = read_methylation(Path(s["path"]), s["format"], index,
                                   a.min_coverage)
        ok = n >= a.min_cpg
        beta[ok, j] = m[ok] / c[ok]
        for i in np.flatnonzero(n > 0):
            long_rows.append({**meta[i], "sample": s["sample"],
                              "group": s["group"], "n_cpg": int(n[i]),
                              "coverage": float(c[i]),
                              "beta": beta[i, j]})
    write_tsv(a.out / "methylation_long.tsv", long_rows,
              ["region", "gene", "set", "kind", "sample", "group", "n_cpg",
               "coverage", "beta"])

    gmean = np.full((len(regions), len(groups)), np.nan)
    for k, g in enumerate(groups):
        cols = [j for j, s in enumerate(samples) if s["group"] == g]
        sub = beta[:, cols]
        with np.errstate(all="ignore"):
            cnt = (~np.isnan(sub)).sum(1)
            gmean[:, k] = np.where(cnt > 0, np.nansum(sub, 1) / np.maximum(
                cnt, 1), np.nan)

    kinds = np.array([m["kind"] for m in meta])
    isel = np.flatnonzero(kinds == "element")
    iprom = np.flatnonzero(kinds == "promoter")
    write_tsv(a.out / "methylation_by_group.tsv",
              [{"gene": meta[i]["gene"], "set": meta[i]["set"],
                **{g: gmean[i, k] for k, g in enumerate(groups)}}
               for i in iprom], ["gene", "set"] + groups)
    if len(isel):
        write_tsv(a.out / "elements_methylation_by_group.tsv",
                  [{**{c: meta[i].get(c, "") for c in
                       ("region", "gene", "set", "distance", "prevalence")},
                    **{g: gmean[i, k] for k, g in enumerate(groups)}}
                   for i in isel],
                  ["region", "gene", "set", "distance", "prevalence"] + groups)

    # summary per kind, group and set
    summ = []
    tset = np.array([m["set"] == "mito" for m in meta])
    print(f"\n{'kind':9s}{'group':28s}{'set':9s}{'n':>6s}{'median':>8s}"
          f"{'unmeth':>8s}{'meth':>8s}{'P(t>c)':>8s}{'p':>9s}")
    for kind in ("promoter", "element"):
        sel = kinds == kind
        if not sel.any():
            continue
        for k, g in enumerate(groups):
            vals = {}
            for s_ in ("mito", "control"):
                v = gmean[sel & (tset == (s_ == "mito")), k]
                vals[s_] = v[~np.isnan(v)]
            ps = (prob_superiority(vals["mito"], vals["control"])
                  if len(vals["control"]) else np.nan)
            pv = (float(mannwhitneyu(vals["mito"], vals["control"]).pvalue)
                  if len(vals["mito"]) and len(vals["control"]) else np.nan)
            for s_, v in vals.items():
                if not len(v):
                    continue
                r = {"kind": kind, "group": g, "set": s_, "n": len(v),
                     "median_beta": float(np.median(v)),
                     "share_unmethylated": float((v < a.unmethylated).mean()),
                     "share_methylated": float((v > a.methylated).mean()),
                     "P_target_gt_control": ps, "p": pv}
                summ.append(r)
                print(f"{kind:9s}{g[:27]:28s}{s_:9s}{r['n']:6d}"
                      f"{r['median_beta']:8.3f}"
                      f"{r['share_unmethylated']:8.2f}"
                      f"{r['share_methylated']:8.2f}{ps:8.3f}{pv:9.2g}")
    write_tsv(a.out / "methylation_summary.tsv", summ,
              ["kind", "group", "set", "n", "median_beta",
               "share_unmethylated", "share_methylated",
               "P_target_gt_control", "p"])

    # variable promoters
    focus = np.array([bool(a.focus) and bool(re.search(a.focus, g, re.I))
                      for g in groups])
    var = []
    for i in iprom:
        v = gmean[i]
        ok = ~np.isnan(v)
        if ok.sum() < 2:
            continue
        lo, hi = int(np.nanargmin(v)), int(np.nanargmax(v))
        rng_ = v[hi] - v[lo]
        if rng_ < a.variable_range:
            continue
        r = {"gene": meta[i]["gene"], "set": meta[i]["set"],
             "min_group": groups[lo], "min_beta": v[lo],
             "max_group": groups[hi], "max_beta": v[hi], "range": rng_}
        if focus.any():
            f_ = v[focus & ok]
            o_ = v[~focus & ok]
            r["focus_mean"] = float(f_.mean()) if len(f_) else np.nan
            r["other_mean"] = float(o_.mean()) if len(o_) else np.nan
            r["focus_minus_other"] = r["focus_mean"] - r["other_mean"]
        var.append(r)
    var.sort(key=lambda r: (r["set"] != "mito", -r["range"]))
    vcols = ["gene", "set", "min_group", "min_beta", "max_group", "max_beta",
             "range"] + (["focus_mean", "other_mean", "focus_minus_other"]
                         if focus.any() else [])
    write_tsv(a.out / "methylation_variable.tsv", var, vcols)
    nv = {s_: sum(r["set"] == s_ for r in var) for s_ in ("mito", "control")}
    nt = {s_: int(sum((tset[i] == (s_ == "mito")) for i in iprom))
          for s_ in ("mito", "control")}
    print(f"\npromoters with a methylation range >= {a.variable_range}: "
          + ", ".join(f"{s_} {nv[s_]}/{nt[s_]}" for s_ in nv if nt[s_]))

    # methylation against expression across GTEx tissues
    if a.tpm:
        gtex = {}
        for s in samples:
            if s.get("gtex_tissue") and s["group"] not in gtex:
                gtex[s["group"]] = s["gtex_tissue"]
        by_id, by_sym, tcols = load_tpm(a.tpm)
        tcol = {t: k for k, t in enumerate(tcols)}
        pairs = [(k, tcol[gtex[g]]) for k, g in enumerate(groups)
                 if g in gtex and gtex[g] in tcol]
        if len(pairs) < a.min_tissues:
            print(f"--tpm: {len(pairs)} groups matched to a GTEx tissue "
                  f"column, fewer than --min-tissues; skipped")
            return
        allg = {r["gene"]: r for r in rows}
        out = []
        for i in iprom:
            gname = meta[i]["gene"]
            e = by_id.get(allg[gname]["gene_id"], by_sym.get(gname))
            if e is None:
                continue
            x = np.array([gmean[i, k] for k, _ in pairs])
            y = np.log2(np.array([e[t] for _, t in pairs]) + 1)
            ok = ~np.isnan(x) & ~np.isnan(y)
            if ok.sum() < a.min_tissues or np.ptp(x[ok]) < a.min_range \
                    or np.ptp(y[ok]) == 0:
                continue
            rho, pv = spearmanr(x[ok], y[ok])
            out.append({"gene": gname, "set": allg[gname]["set"],
                        "n_tissues": int(ok.sum()),
                        "methylation_range": float(np.ptp(x[ok])),
                        "rho": float(rho), "p": float(pv)})
        write_tsv(a.out / "methylation_expression.tsv", out,
                  ["gene", "set", "n_tissues", "methylation_range", "rho",
                   "p"])
        for s_ in ("mito", "control"):
            v = np.array([r["rho"] for r in out if r["set"] == s_])
            if len(v):
                print(f"methylation-expression, {s_}: median rho "
                      f"{np.median(v):.2f} over {len(v)} promoters")


if __name__ == "__main__":
    main()
