#!/usr/bin/env python3
"""
l1_score.py - per-gene score of transcriptional control (L1): does the
protein follow its mRNA from one tissue to another? And its validation.

Purpose
    L1 asks, for each gene, how much of the variation of its protein across
    tissues is already set by the variation of its mRNA. A gene whose
    protein follows its mRNA is controlled at transcription, and a
    regulatory variant that lowers its mRNA is expected to lower its
    protein; a gene whose protein-to-mRNA ratio changes between tissues, or
    whose protein stays flat while its mRNA moves, is controlled after
    transcription (translation, degradation, assembly).

      score     L1 from matched RNA and protein measured in the same samples
                of several tissues (the GTEx tissue proteome of Jiang et al.
                2020 and the GTEx RNA-seq of the same samples), corrected for
                measurement noise, with an explicit status when a gene cannot
                be scored, and at the level of protein complexes.
      validate  agreement of L1 with independent data: L1 computed on
                another atlas, published per-gene scores (protein buffering
                of copy number, eQTL-to-pQTL transmission, half-lives, ...),
                expected high and low gene sets, and patient fibroblasts with
                RNA-seq and proteomics (Kopajtich et al.): when the mRNA of a
                gene collapses in one patient, does the protein follow?

Inputs
    score
      --rna        matrix, genes x samples: first column the gene (Ensembl id,
                   versioned or not, or symbol), then optionally Description,
                   then one column per sample; GCT accepted. Read counts or
                   TPM (--rna-units, guessed: integers are counts).
      --protein    matrix, genes x samples, same layout: log2 abundances or
                   linear intensities (--protein-scale, guessed).
      --samples    TSV, one row per sample measured for both: tissue,
                   rna_sample, protein_sample, optionally donor and the
                   columns named in --covariates (for instance the GTEx
                   ischemic time SMTSISCH).
      --gtf        GENCODE GTF: symbols, Ensembl ids, and the genes of the
                   mitochondrial genome (left out of the RNA normalisation).
      --genes      the gene set (symbols, or a table with a gene column such
                   as pathways.tsv); --controls optional.
      --complexes  optional TSV complex, gene; or --pathways pathways.tsv
                   with --complex-regex, which takes as complexes the
                   MitoPathways whose name matches (default subunits and
                   ribosome).
    validate
      --l1         l1.tsv of `score`.
      --replicate  optional l1.tsv computed on another atlas.
      --external   NAME=PATH[:COLUMN][:+|-], repeatable: per-gene table and
                   the expected sign of its correlation with L1.
      --expected-high, --expected-low   optional gene lists.
      --patients-rna, --patients-protein, --patients-samples
                   patient cohort: count matrix, protein matrix, and the
                   table pairing rna_sample and protein_sample (one row per
                   individual), with an optional batch column (TMT plex):
                   each protein is then centred on its median within batch. --rna-outliers optionally gives the RNA
                   outliers already called (columns rna_sample, gene), for
                   instance by OUTRIDER; otherwise they are called here.

Outputs
    score (in --out)
      l1.tsv            per gene: gene, gene_id, set, n_tissues, n_pairs,
                        sd_rna, sd_protein (across tissue means, log2),
                        rel_rna, rel_protein, rel_ratio (split-half
                        reliabilities), r and rho (raw Pearson and Spearman
                        across tissues), r_corrected, ci_low, ci_high,
                        share_mrna, status, score; with complexes: complex,
                        complex_score, protein_coupling, rna_coupling
      l1_complexes.tsv  the same statistics for each complex
      l1_summary.tsv    per set: genes per status, median score
      content_index.tsv per sample: mitochondrial content index, RNA and
                        protein
    validate (in --out)
      validation.tsv    one row per test: test, genes or events, statistic,
                        value, p, expected direction, agrees
      patient_events.tsv one row per RNA outlier: rna_sample, gene, rna_z,
                        protein_z, protein_missing, propagated, l1
      patient_genes.tsv per gene: fibroblast RNA-protein correlation across
                        individuals, events, propagated events, l1

Method
    Normalisation. RNA counts are scaled by median-of-ratios size factors
    computed on nuclear genes only, TPM are rescaled so that nuclear genes
    sum to one million: the mtDNA-encoded transcripts, up to 40% of the
    reads in heart, would otherwise shrink every other gene in the tissues
    rich in mitochondria. Values are log2(x + 1). Protein abundances are
    log2-transformed if linear and centred on their median per sample.

    Covariates. Each covariate is removed gene by gene with a slope
    estimated within tissues, so that a covariate that differs between
    tissues (ischemic time) does not take away tissue differences.

    Mitochondrial content. Tissues differ about 30-fold in mitochondrial
    mass, which moves every mitochondrial mRNA and protein together and
    would make them all look transcriptionally controlled. For the genes of
    the content set (--content-genes, default the gene set), the median of
    the set in each sample is subtracted, RNA and protein separately: L1
    then asks whether the share of the gene within mitochondria follows its
    mRNA share. --no-content-adjustment keeps raw values.

    Concordance. Per gene, the tissue means of RNA and protein (samples with
    both values) give two profiles; r is their Pearson correlation. The
    reliability of each profile is estimated by splitting the samples of
    every tissue into two halves, correlating the two half-profiles across
    tissues and applying the Spearman-Brown correction, averaged over
    --split-draws random splits. The corrected correlation is r /
    sqrt(rel_rna x rel_protein), bounded to [-1, 1] (Spearman's correction
    for attenuation; Franks et al. 2017). share_mrna is the noise-corrected
    share of the protein variation carried by the mRNA:
    rel_rna var(rna) / (rel_rna var(rna) + rel_ratio var(protein - rna)).
    A 95% interval comes from --boot resamplings of the tissues.

    Status. too_few_tissues: fewer than --min-tissues tissues. flat: neither
    the mRNA nor the protein varies across tissues by --min-sd (log2, noise
    corrected); such a gene cannot show concordance and is not called
    post-transcriptional. unreliable: a reliability under --min-reliability.
    estimable: score = r_corrected. When no tissue has two samples the
    reliabilities cannot be computed; the score is then the raw r and the
    status estimable_uncorrected.

    Complexes. Complex profiles are the mean of their subunits. A subunit
    controlled by assembly follows its complex at the protein level
    (protein_coupling: correlation of its protein with the mean of the
    other subunits) more than at the mRNA level (rna_coupling), and its own
    score is low even when the complex score is high.

    Validation. Spearman correlation of L1 with each external score, with
    the expected sign; probability that a gene of --expected-high scores
    above one of --expected-low. Patients: RNA and protein are normalised
    as above; a robust z-score (median, MAD) is computed per gene across
    individuals; an RNA outlier has rna_z <= -(--outlier-z); it propagates
    when protein_z <= -(--protein-z) or the protein is missing in this
    individual while measured in most others. L1 should be higher for
    events that propagate (Spearman of L1 with protein_z, expected
    negative; Fisher test of propagation between the top and bottom thirds
    of L1).

Usage
    python3 l1_score.py score --rna data/jiang2020_rna_samples.tsv \\
        --protein data/jiang2020_protein_samples.tsv \\
        --samples data/jiang2020_samples.tsv \\
        --covariates SMTSISCH --gtf data/gencode.gtf.gz \\
        --genes results/regulome/pathways.tsv \\
        --pathways results/regulome/pathways.tsv --out results/l1

    python3 l1_score.py validate --l1 results/l1/l1.tsv \\
        --replicate results/l1_hpa/l1.tsv \\
        --external buffering=data/schukken2022.tsv:buffering:- \\
        --patients-rna data/fibro_counts.tsv \\
        --patients-protein data/fibro_protein.tsv \\
        --patients-samples data/fibro_samples.tsv --gtf data/gencode.gtf.gz \\
        --out results/l1/validation

References
    Jiang L et al. 2020. A quantitative proteome map of the human body. Cell
        183:269-283.
    Wang D et al. 2019. A deep proteome and transcriptome abundance atlas of
        29 healthy human tissues. Mol Syst Biol 15:e8503.
    Franks A et al. 2017. Post-transcriptional regulation across human
        tissues. PLoS Comput Biol 13:e1005535.
    Upadhya SR, Ryan CJ. 2022. Experimental reproducibility limits the
        correlation between mRNA and protein abundances in tumor proteomic
        profiles. Cell Rep Methods 2:100288.
    Williams EG et al. 2018. Multiomic profiling of the liver across diets
        and age in a diverse mouse population. Mol Cell Proteomics 17:1890.
    Taggart JC, Li GW. 2018. Production of protein-complex components is
        stoichiometric and lacks general feedback regulation in eukaryotes.
        Cell Syst 7:580-589.
    Perez MF et al. 2021. Mito-nuclear co-expression across tissues and its
        artefacts. Genome Biol 22:327.
    Love MI et al. 2014. DESeq2 (median-of-ratios size factors). Genome
        Biol 15:550.
    Kopajtich R et al. Integration of proteomics with genomics and
        transcriptomics increases the diagnostic rate of Mendelian disorders.
    Spearman C. 1904. The proof and measurement of association between two
        things. Am J Psychol 15:72-101.
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
from regulome import (read_gene_list, read_table, to_float,  # noqa: E402
                      write_tsv)


# ---------------------------------------------------------------------------
# input

def load_annotation(gtf: Path):
    """Ensembl id -> symbol, symbol -> chromosome, for every gene."""
    ens, chrom = {}, {}
    with opener(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene" or "_PAR_Y" in f[8]:
                continue
            a = f[8]
            i = a.find('gene_id "')
            j = a.find('gene_name "')
            e = a[i + 9:a.find('"', i + 9)].split(".")[0]
            s = a[j + 11:a.find('"', j + 11)] if j >= 0 else e
            ens[e] = s
            chrom.setdefault(s, f[0])
    return ens, chrom


def is_mt(symbol, chrom):
    c = chrom.get(symbol, "")
    return c in ("chrM", "MT", "M") or (not c and symbol.startswith("MT-"))


def read_matrix(path: Path, ens2sym: dict):
    """(symbols, sample names, matrix); of several rows with one symbol,
    the first is kept."""
    with opener(path) as fh:
        first = fh.readline()
        if first.startswith("#1.2"):
            fh.readline()
            head = fh.readline().rstrip("\n").split("\t")
        else:
            head = first.rstrip("\n").split("\t")
        skip = 2 if len(head) > 1 and head[1].lower() in ("description",
                                                           "symbol",
                                                           "gene_name") \
            else 1
        samples = head[skip:]
        genes, rows, seen = [], [], set()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f[0]:
                continue
            g = f[0].split(".")[0] if f[0].startswith("ENS") else f[0]
            g = ens2sym.get(g, g)
            if g in seen:
                continue
            seen.add(g)
            genes.append(g)
            rows.append([to_float(x) for x in f[skip:]])
    return genes, samples, np.array(rows, float)


def normalise_rna(M, genes, chrom, units):
    """log2 normalised expression; mtDNA genes left out of the scaling."""
    nuc = np.array([not is_mt(g, chrom) for g in genes])
    if units == "auto":
        x = M[~np.isnan(M)]
        units = "counts" if x.size and np.allclose(x, np.round(x)) \
            and x.max() > 1000 else "tpm"
    if units == "counts":
        N = M[nuc]
        ok = np.all(N > 0, axis=1)
        if ok.sum() >= 50:
            logs = np.log(N[ok])
            ref = logs.mean(1, keepdims=True)
            sf = np.exp(np.median(logs - ref, axis=0))
        else:
            tot = np.nansum(N, 0)
            sf = tot / np.exp(np.mean(np.log(tot)))
        out = np.log2(M / sf[None, :] + 1)
    else:
        tot = np.nansum(M[nuc], 0)
        out = np.log2(M / tot[None, :] * 1e6 + 1)
    return out, units


def normalise_protein(M, scale):
    if scale == "auto":
        med = np.nanmedian(M)
        scale = "linear" if med > 100 else "log2"
    if scale == "linear":
        with np.errstate(divide="ignore", invalid="ignore"):
            M = np.where(M > 0, np.log2(M), np.nan)
    return M - np.nanmedian(M, axis=0, keepdims=True), scale


def resolve_set(path, names, aliases):
    if not path:
        return []
    m, _u, _n = resolve_symbols(read_gene_list(path), names, aliases,
                                NOT_GENES, NOT_GENE_PREFIXES)
    return list(dict.fromkeys(m.values()))


# ---------------------------------------------------------------------------
# statistics

def row_pearson(A, B):
    """Pearson correlation row by row, over the columns where both are
    finite; NaN with fewer than 3 such columns."""
    ok = ~np.isnan(A) & ~np.isnan(B)
    n = ok.sum(1)
    A0 = np.where(ok, A, 0.0)
    B0 = np.where(ok, B, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        ma = A0.sum(1) / n
        mb = B0.sum(1) / n
        da = np.where(ok, A - ma[:, None], 0.0)
        db = np.where(ok, B - mb[:, None], 0.0)
        r = (da * db).sum(1) / np.sqrt((da ** 2).sum(1) * (db ** 2).sum(1))
    r[n < 3] = np.nan
    return r


def tissue_means(X, groups):
    """genes x tissues: mean over the samples of each tissue."""
    out = np.full((X.shape[0], len(groups)), np.nan)
    for k, idx in enumerate(groups):
        if len(idx):
            with np.errstate(invalid="ignore"):
                sub = X[:, idx]
                cnt = (~np.isnan(sub)).sum(1)
                out[:, k] = np.where(cnt > 0, np.nansum(sub, 1) /
                                     np.maximum(cnt, 1), np.nan)
    return out


def split_half_reliability(X, groups, draws, rng, min_tissues):
    """Spearman-Brown corrected split-half reliability of the tissue
    profile, per row."""
    multi = [g for g in groups if len(g) >= 2]
    if len(multi) < min_tissues:
        return np.full(X.shape[0], np.nan)
    rs = []
    for _ in range(draws):
        ga, gb = [], []
        for g in multi:
            p = rng.permutation(g)
            h = len(p) // 2
            ga.append(p[:h])
            gb.append(p[h:2 * h] if len(p) % 2 == 0 else p[h:])
        rs.append(row_pearson(tissue_means(X, ga), tissue_means(X, gb)))
    with np.errstate(all="ignore"):
        r = np.nanmean(np.array(rs), 0)
    r = np.clip(r, 0, None)
    return 2 * r / (1 + r)


def remove_covariates(X, cov, tissue_of):
    """Remove each covariate with a slope estimated within tissues."""
    X = X.copy()
    for c in cov.T:
        cc = c.copy()
        for t in np.unique(tissue_of):
            sel = tissue_of == t
            cc[sel] = c[sel] - np.nanmean(c[sel])
        ok_c = ~np.isnan(cc)
        Xc = X.copy()
        for t in np.unique(tissue_of):
            sel = tissue_of == t
            with np.errstate(all="ignore"):
                Xc[:, sel] = X[:, sel] - np.nanmean(X[:, sel], 1,
                                                    keepdims=True)
        okx = ~np.isnan(Xc) & ok_c[None, :]
        num = np.where(okx, Xc * cc[None, :], 0).sum(1)
        den = np.where(okx, cc[None, :] ** 2, 0).sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            slope = np.where(den > 0, num / den, 0.0)
        centred = np.where(np.isnan(c), 0.0, c - np.nanmean(c))
        X = X - slope[:, None] * centred[None, :]
    return X


def concordance(R, P, tissue_of, a, rng, labels):
    """L1 statistics for every row of R and P (genes or complexes)."""
    tissues = sorted(set(tissue_of))
    both = ~np.isnan(R) & ~np.isnan(P)
    Rb = np.where(both, R, np.nan)
    Pb = np.where(both, P, np.nan)
    groups = [np.flatnonzero(tissue_of == t) for t in tissues]
    mR, mP = tissue_means(Rb, groups), tissue_means(Pb, groups)
    D = Pb - Rb
    n_t = (~np.isnan(mR) & ~np.isnan(mP)).sum(1)
    rel_r = split_half_reliability(Rb, groups, a.split_draws, rng,
                                   a.min_reliability_tissues)
    rel_p = split_half_reliability(Pb, groups, a.split_draws, rng,
                                   a.min_reliability_tissues)
    rel_d = split_half_reliability(D, groups, a.split_draws, rng,
                                   a.min_reliability_tissues)
    r = row_pearson(mR, mP)
    with np.errstate(all="ignore"):
        sd_r = np.nanstd(mR, 1, ddof=1)
        sd_p = np.nanstd(mP, 1, ddof=1)
        var_d = np.nanvar(tissue_means(D, groups), 1, ddof=1)
        factor = 1 / np.sqrt(rel_r * rel_p)
        rc = np.clip(r * factor, -1, 1)
        share = rel_r * sd_r ** 2 / (rel_r * sd_r ** 2 + rel_d * var_d)
        true_sd = np.fmax(np.sqrt(rel_r) * sd_r, np.sqrt(rel_p) * sd_p)
    no_rel = np.all(np.isnan(rel_r))
    if no_rel:
        true_sd = np.fmax(sd_r, sd_p)
    out = []
    for i, lab in enumerate(labels):
        o = {"label": lab, "n_tissues": int(n_t[i]),
             "n_pairs": int(both[i].sum()), "sd_rna": sd_r[i],
             "sd_protein": sd_p[i], "rel_rna": rel_r[i],
             "rel_protein": rel_p[i], "rel_ratio": rel_d[i], "r": r[i],
             "r_corrected": rc[i], "share_mrna": share[i],
             "ci_low": np.nan, "ci_high": np.nan}
        ok = ~np.isnan(mR[i]) & ~np.isnan(mP[i])
        if ok.sum() >= 3:
            from scipy.stats import spearmanr
            o["rho"] = float(spearmanr(mR[i, ok], mP[i, ok])[0]) \
                if np.ptp(mR[i, ok]) > 0 and np.ptp(mP[i, ok]) > 0 \
                else np.nan
        else:
            o["rho"] = np.nan
        if n_t[i] < a.min_tissues:
            o["status"] = "too_few_tissues"
        elif not true_sd[i] >= a.min_sd:
            o["status"] = "flat"
        elif no_rel:
            o["status"] = "estimable_uncorrected"
        elif not (rel_r[i] >= a.min_reliability and
                  rel_p[i] >= a.min_reliability):
            o["status"] = "unreliable"
        else:
            o["status"] = "estimable"
        if o["status"].startswith("estimable"):
            f = 1.0 if no_rel else factor[i]
            o["score"] = o["r"] * f if no_rel else o["r_corrected"]
            if a.boot > 0:
                x, y = mR[i, ok], mP[i, ok]
                idx = rng.integers(0, len(x), (a.boot, len(x)))
                bx, by = x[idx], y[idx]
                bx = bx - bx.mean(1, keepdims=True)
                by = by - by.mean(1, keepdims=True)
                with np.errstate(all="ignore"):
                    br = (bx * by).sum(1) / np.sqrt((bx ** 2).sum(1) *
                                                    (by ** 2).sum(1))
                br = np.clip(br[~np.isnan(br)] * f, -1, 1)
                if len(br):
                    o["ci_low"], o["ci_high"] = np.percentile(br, [2.5, 97.5])
        else:
            o["score"] = np.nan
        out.append(o)
    return out, mR, mP


# ---------------------------------------------------------------------------
# score

def cmd_score(a) -> None:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    ens2sym, chrom = load_annotation(a.gtf)
    names = set(chrom)
    aliases = build_aliases(a.aliases, not a.no_builtin_aliases)
    target = resolve_set(a.genes, names, aliases)
    control = resolve_set(a.controls, names, aliases)
    set_of = {g: "control" for g in control}
    set_of.update({g: "mito" for g in target})

    rg, rs, RM = read_matrix(a.rna, ens2sym)
    pg, ps, PM = read_matrix(a.protein, ens2sym)
    samples = read_table(a.samples)
    for k in ("tissue", "rna_sample", "protein_sample"):
        if samples and k not in samples[0]:
            sys.exit("--samples needs the columns tissue, rna_sample, "
                     "protein_sample")
    ri = {s: i for i, s in enumerate(rs)}
    pi = {s: i for i, s in enumerate(ps)}
    pairs = [s for s in samples if s["rna_sample"] in ri
             and s["protein_sample"] in pi]
    if len(pairs) < len(samples):
        print(f"{len(samples) - len(pairs)} rows of --samples without both "
              "matrices left out")
    RM, units = normalise_rna(RM, rg, chrom, a.rna_units)
    PM, scale = normalise_protein(PM, a.protein_scale)
    print(f"RNA: {len(rg)} genes ({units}); protein: {len(pg)} genes "
          f"({scale}); {len(pairs)} matched samples")

    genes = [g for g in pg if g in set(rg) and not is_mt(g, chrom)]
    rpos = {g: i for i, g in enumerate(rg)}
    ppos = {g: i for i, g in enumerate(pg)}
    R = RM[[rpos[g] for g in genes]][:, [ri[s["rna_sample"]] for s in pairs]]
    P = PM[[ppos[g] for g in genes]][:, [pi[s["protein_sample"]]
                                        for s in pairs]]
    tissue_of = np.array([s["tissue"] for s in pairs])

    if a.covariates:
        cols = a.covariates.split(",")
        cov = np.array([[to_float(s.get(c)) for c in cols] for s in pairs])
        R = remove_covariates(R, cov, tissue_of)
        P = remove_covariates(P, cov, tissue_of)
        print(f"covariates removed within tissues: {', '.join(cols)}")

    content = set(resolve_set(a.content_genes, names, aliases)) \
        if a.content_genes else set(target)
    cidx = np.array([g in content for g in genes])
    ci_rows = []
    if not a.no_content_adjustment and cidx.sum() >= 10:
        with np.errstate(all="ignore"):
            cR = np.nanmedian(R[cidx], 0)
            cP = np.nanmedian(P[cidx], 0)
        R[cidx] -= cR[None, :]
        P[cidx] -= cP[None, :]
        ci_rows = [{"rna_sample": s["rna_sample"],
                    "protein_sample": s["protein_sample"],
                    "tissue": s["tissue"], "content_rna": cR[k],
                    "content_protein": cP[k]} for k, s in enumerate(pairs)]
        write_tsv(out / "content_index.tsv", ci_rows,
                  ["rna_sample", "protein_sample", "tissue", "content_rna",
                   "content_protein"])
        print(f"mitochondrial content removed for {int(cidx.sum())} genes")

    stats, mR, mP = concordance(R, P, tissue_of, a, rng, genes)
    if all(np.isnan(s["rel_rna"]) for s in stats):
        print("no tissue with two samples: reliabilities not computed, "
              "scores are uncorrected correlations")

    # complexes
    members = defaultdict(list)
    if a.complexes:
        for r in read_table(a.complexes):
            v = list(r.values())
            if len(v) >= 2 and v[0] and v[1]:
                members[v[0]].append(v[1])
    elif a.pathways:
        rx = re.compile(a.complex_regex, re.I)
        for r in read_table(a.pathways):
            for pth in r.get("pathways", "").split("|"):
                leaf = pth.split(" > ")[-1].strip()
                if leaf and rx.search(leaf):
                    members[leaf].append(r["gene"])
    gpos = {g: i for i, g in enumerate(genes)}
    cstats, comp_of = [], {}
    if members:
        names_c, CR, CP, mem_idx = [], [], [], []
        for c, gs in members.items():
            idx = [gpos[g] for g in dict.fromkeys(gs) if g in gpos]
            if len(idx) < 3:
                continue
            with np.errstate(all="ignore"):
                CR.append(np.nanmean(R[idx], 0))
                CP.append(np.nanmean(P[idx], 0))
            names_c.append(c)
            mem_idx.append(idx)
        if names_c:
            cst, cmR, cmP = concordance(np.array(CR), np.array(CP),
                                        tissue_of, a, rng, names_c)
            for k, c in enumerate(names_c):
                cst[k]["n_subunits"] = len(mem_idx[k])
                cstats.append(cst[k])
                for i in mem_idx[k]:
                    others = [j for j in mem_idx[k] if j != i]
                    with np.errstate(all="ignore"):
                        oR = np.nanmean(mR[others], 0)
                        oP = np.nanmean(mP[others], 0)
                    pc = row_pearson(mP[i][None], oP[None])[0]
                    rc = row_pearson(mR[i][None], oR[None])[0]
                    prev = comp_of.get(i)
                    if prev is None or len(mem_idx[k]) > prev[2]:
                        comp_of[i] = (c, cst[k]["score"], len(mem_idx[k]),
                                      pc, rc)
            print(f"{len(names_c)} complexes scored")

    rows = []
    for i, (g, s) in enumerate(zip(genes, stats)):
        r = {"gene": g, "gene_id": "", "set": set_of.get(g, "other"), **s}
        if i in comp_of:
            c, cs, _n, pc, rc = comp_of[i]
            r.update(complex=c, complex_score=cs, protein_coupling=pc,
                     rna_coupling=rc)
        rows.append(r)
    sym2ens = {}
    for e, x in ens2sym.items():
        sym2ens.setdefault(x, e)
    for r in rows:
        r["gene_id"] = sym2ens.get(r["gene"], "")
    cols = ["gene", "gene_id", "set", "n_tissues", "n_pairs", "sd_rna",
            "sd_protein", "rel_rna", "rel_protein", "rel_ratio", "r", "rho",
            "r_corrected", "ci_low", "ci_high", "share_mrna", "status",
            "score"]
    if comp_of:
        cols += ["complex", "complex_score", "protein_coupling",
                 "rna_coupling"]
    write_tsv(out / "l1.tsv", rows, cols)
    if cstats:
        for c in cstats:
            c["complex"] = c.pop("label")
        write_tsv(out / "l1_complexes.tsv", cstats,
                  ["complex", "n_subunits"] + [c for c in cols[3:18]
                                               if c != "rho"] + ["rho"])

    summ = []
    statuses = ["estimable", "estimable_uncorrected", "unreliable", "flat",
                "too_few_tissues"]
    print(f"\n{'set':9s}{'genes':>7s}" + "".join(f"{s[:12]:>14s}"
                                                 for s in statuses)
          + f"{'median':>9s}")
    for lab in ("mito", "control", "other"):
        sub = [r for r in rows if r["set"] == lab]
        if not sub:
            continue
        sc = np.array([r["score"] for r in sub], float)
        sc = sc[~np.isnan(sc)]
        x = {"set": lab, "genes": len(sub),
             **{s: sum(r["status"] == s for r in sub) for s in statuses},
             "median_score": float(np.median(sc)) if len(sc) else np.nan}
        summ.append(x)
        print(f"{lab:9s}{len(sub):7d}" + "".join(f"{x[s]:14d}"
                                                 for s in statuses)
              + f"{x['median_score']:9.3f}")
    write_tsv(out / "l1_summary.tsv", summ,
              ["set", "genes"] + statuses + ["median_score"])
    print(f"Written to {out}")


# ---------------------------------------------------------------------------
# validate

def robust_z(X):
    med = np.nanmedian(X, 1, keepdims=True)
    mad = np.nanmedian(np.abs(X - med), 1, keepdims=True) * 1.4826
    with np.errstate(all="ignore"):
        return (X - med) / np.where(mad > 0, mad, np.nan)


def cmd_validate(a) -> None:
    from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

    from gtex_eqtl import prob_superiority
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    l1 = {r["gene"]: r for r in read_table(a.l1)}
    score = {g: to_float(r.get(a.score_column)) for g, r in l1.items()}
    score = {g: v for g, v in score.items() if not np.isnan(v)}
    mito = {g for g, r in l1.items() if r.get("set") == "mito"}
    tests = []

    def corr_test(name, other, expect, restrict=None):
        gs = [g for g in score if g in other and not np.isnan(other[g])
              and (restrict is None or g in restrict)]
        if len(gs) < 5:
            print(f"  {name}: {len(gs)} genes in common, skipped")
            return
        rho, p = spearmanr([score[g] for g in gs], [other[g] for g in gs])
        tests.append({"test": name, "n": len(gs), "statistic": "spearman",
                      "value": float(rho), "p": float(p), "expected": expect,
                      "agrees": int(np.sign(rho) == (1 if expect == "+"
                                                     else -1))})

    if a.replicate:
        rep = {r["gene"]: to_float(r.get(a.score_column))
               for r in read_table(a.replicate)}
        corr_test("replicate", rep, "+")
        corr_test("replicate_mito", rep, "+", mito)
    for spec in a.external or []:
        m = re.match(r"^([^=]+)=(.+?)(?::([^:]+))?(?::([+-]))?$", spec)
        if not m:
            sys.exit(f"--external {spec!r}: expected NAME=PATH[:COLUMN][:+|-]")
        name, path, col, sign = m.groups()
        if col in ("+", "-") and sign is None:
            col, sign = None, col
        rows = read_table(path)
        cols = list(rows[0])
        col = col or cols[-1]
        if col not in cols:
            sys.exit(f"--external {name}: no column {col!r} in {path}")
        other = {r[cols[0]]: to_float(r[col]) for r in rows}
        corr_test(name, other, sign or "+")
        corr_test(f"{name}_mito", other, sign or "+", mito)
    if a.expected_high and a.expected_low:
        hi = [score[g] for g in read_gene_list(a.expected_high) if g in score]
        lo = [score[g] for g in read_gene_list(a.expected_low) if g in score]
        if hi and lo:
            ps = prob_superiority(np.array(hi), np.array(lo))
            tests.append({"test": "expected_sets", "n": len(hi) + len(lo),
                          "statistic": "P(high>low)", "value": ps,
                          "p": float(mannwhitneyu(hi, lo).pvalue),
                          "expected": "+", "agrees": int(ps > 0.5)})

    if a.patients_rna and a.patients_protein and a.patients_samples:
        ens2sym, chrom = load_annotation(a.gtf) if a.gtf else ({}, {})
        rg, rs, RM = read_matrix(a.patients_rna, ens2sym)
        pg, ps_, PM = read_matrix(a.patients_protein, ens2sym)
        RM, _u = normalise_rna(RM, rg, chrom, "auto")
        PM, _s = normalise_protein(PM, "auto")
        ri = {s: i for i, s in enumerate(rs)}
        pi = {s: i for i, s in enumerate(ps_)}
        pairs = [r for r in read_table(a.patients_samples)
                 if r.get("rna_sample") in ri
                 and r.get("protein_sample") in pi]
        genes = [g for g in pg if g in set(rg)]
        rpos = {g: i for i, g in enumerate(rg)}
        ppos = {g: i for i, g in enumerate(pg)}
        R = RM[[rpos[g] for g in genes]][:, [ri[p["rna_sample"]]
                                            for p in pairs]]
        P = PM[[ppos[g] for g in genes]][:, [pi[p["protein_sample"]]
                                            for p in pairs]]
        batches = [p.get("batch", "") for p in pairs]
        if any(batches):
            # TMT plexes: each protein is centred within its batch
            for b in set(batches):
                sel = np.array([x == b for x in batches])
                with np.errstate(all="ignore"):
                    P[:, sel] -= np.nanmedian(P[:, sel], 1, keepdims=True)
            print(f"protein centred within {len(set(batches))} batches")
        print(f"patients: {len(pairs)} individuals with RNA and protein, "
              f"{len(genes)} genes")
        zR, zP = robust_z(R), robust_z(P)
        miss = np.isnan(P)
        frac_measured = 1 - miss.mean(1)
        called = None
        if a.rna_outliers:
            called = {(r["rna_sample"], r["gene"])
                      for r in read_table(a.rna_outliers)}
        events = []
        for i, g in enumerate(genes):
            for k, p in enumerate(pairs):
                if called is not None:
                    hit = (p["rna_sample"], g) in called
                else:
                    hit = zR[i, k] <= -a.outlier_z
                if not hit:
                    continue
                pm = bool(miss[i, k]) and frac_measured[i] >= 0.5
                prop = pm or (not np.isnan(zP[i, k])
                              and zP[i, k] <= -a.protein_z)
                events.append({"rna_sample": p["rna_sample"], "gene": g,
                               "rna_z": zR[i, k], "protein_z": zP[i, k],
                               "protein_missing": int(pm),
                               "propagated": int(prop),
                               "l1": score.get(g, np.nan)})
        write_tsv(out / "patient_events.tsv", events,
                  ["rna_sample", "gene", "rna_z", "protein_z",
                   "protein_missing", "propagated", "l1"])
        # per gene: correlation across individuals
        fib = {}
        rr = []
        for i, g in enumerate(genes):
            ok = ~np.isnan(R[i]) & ~np.isnan(P[i])
            c = np.nan
            if ok.sum() >= 10 and np.ptp(R[i, ok]) > 0 and \
                    np.ptp(P[i, ok]) > 0:
                c = float(np.corrcoef(R[i, ok], P[i, ok])[0, 1])
            fib[g] = c
            ev = [e for e in events if e["gene"] == g]
            rr.append({"gene": g, "fibroblast_r": c, "events": len(ev),
                       "propagated": sum(e["propagated"] for e in ev),
                       "l1": score.get(g, np.nan)})
        write_tsv(out / "patient_genes.tsv", rr,
                  ["gene", "fibroblast_r", "events", "propagated", "l1"])
        corr_test("fibroblast_rna_protein_r", fib, "+")
        corr_test("fibroblast_rna_protein_r_mito", fib, "+", mito)
        ev = [e for e in events if not np.isnan(e["l1"])]
        withz = [e for e in ev if not np.isnan(e["protein_z"])]
        if len(withz) >= 5:
            rho, p = spearmanr([e["l1"] for e in withz],
                               [e["protein_z"] for e in withz])
            tests.append({"test": "outlier_protein_z", "n": len(withz),
                          "statistic": "spearman", "value": float(rho),
                          "p": float(p), "expected": "-",
                          "agrees": int(rho < 0)})
        if len(ev) >= 6:
            vals = np.array([e["l1"] for e in ev])
            lo, hi = np.quantile(vals, [1 / 3, 2 / 3])
            top = [e["propagated"] for e in ev if e["l1"] >= hi]
            bot = [e["propagated"] for e in ev if e["l1"] <= lo]
            if top and bot:
                k1, k2 = sum(top), sum(bot)
                odds, p = fisher_exact([[k1, len(top) - k1],
                                        [k2, len(bot) - k2]])
                tests.append({"test": "outlier_propagation_top_vs_bottom",
                              "n": len(top) + len(bot),
                              "statistic": "odds_ratio", "value": float(odds),
                              "p": float(p), "expected": "+",
                              "agrees": int(odds > 1),
                              "detail": f"{k1}/{len(top)} vs "
                                        f"{k2}/{len(bot)}"})
        print(f"  {len(events)} RNA outliers, "
              f"{sum(e['propagated'] for e in events)} reaching the protein")

    write_tsv(out / "validation.tsv", tests,
              ["test", "n", "statistic", "value", "p", "expected", "agrees",
               "detail"])
    print(f"\n{'test':40s}{'n':>7s}{'value':>9s}{'p':>10s}  expected")
    for t in tests:
        print(f"{t['test']:40s}{t['n']:7d}{t['value']:9.3f}{t['p']:10.2g}  "
              f"{t['expected']}  {'ok' if t['agrees'] else 'NO'}")
    print(f"Written to {out}")


# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("score", help="L1 from matched RNA and protein")
    s.add_argument("--rna", required=True, type=Path)
    s.add_argument("--protein", required=True, type=Path)
    s.add_argument("--samples", required=True, type=Path)
    s.add_argument("--gtf", required=True, type=Path)
    s.add_argument("--genes", required=True, type=Path)
    s.add_argument("--controls", type=Path)
    s.add_argument("--aliases", type=Path)
    s.add_argument("--no-builtin-aliases", action="store_true")
    s.add_argument("--rna-units", choices=("auto", "counts", "tpm"),
                   default="auto")
    s.add_argument("--protein-scale", choices=("auto", "log2", "linear"),
                   default="auto")
    s.add_argument("--covariates", default="",
                   help="comma-separated columns of --samples removed "
                        "within tissues (for instance SMTSISCH)")
    s.add_argument("--content-genes", type=Path,
                   help="genes of the mitochondrial content index "
                        "(default: --genes)")
    s.add_argument("--no-content-adjustment", action="store_true")
    s.add_argument("--complexes", type=Path)
    s.add_argument("--pathways", type=Path)
    s.add_argument("--complex-regex", default="subunits|ribosome",
                   help="MitoPathways taken as complexes with --pathways "
                        "(default %(default)s)")
    s.add_argument("--min-tissues", type=int, default=10)
    s.add_argument("--min-sd", type=float, default=0.5,
                   help="noise-corrected sd across tissues (log2) below "
                        "which a gene is flat (default %(default)s)")
    s.add_argument("--min-reliability", type=float, default=0.3)
    s.add_argument("--min-reliability-tissues", type=int, default=5,
                   help="tissues with two samples needed for the "
                        "reliabilities (default %(default)s)")
    s.add_argument("--split-draws", type=int, default=20)
    s.add_argument("--boot", type=int, default=200)
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--out", required=True)

    v = sub.add_parser("validate", help="agreement with independent data")
    v.add_argument("--l1", required=True, type=Path)
    v.add_argument("--score-column", default="score")
    v.add_argument("--replicate", type=Path)
    v.add_argument("--external", action="append")
    v.add_argument("--expected-high", type=Path)
    v.add_argument("--expected-low", type=Path)
    v.add_argument("--patients-rna", type=Path)
    v.add_argument("--patients-protein", type=Path)
    v.add_argument("--patients-samples", type=Path)
    v.add_argument("--rna-outliers", type=Path)
    v.add_argument("--gtf", type=Path)
    v.add_argument("--outlier-z", type=float, default=3.0)
    v.add_argument("--protein-z", type=float, default=2.0)
    v.add_argument("--out", required=True)
    a = p.parse_args()
    {"score": cmd_score, "validate": cmd_validate}[a.cmd](a)


if __name__ == "__main__":
    main()
