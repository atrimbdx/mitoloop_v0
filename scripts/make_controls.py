#!/usr/bin/env python3
"""
make_controls.py - control genes matched on expression, cross-tissue
variability, gene length and gene density.

Purpose
-------
Draws, for each gene of a gene set (for example nuclear-encoded mitochondrial
genes), --controls control genes with comparable expression level,
comparable variability of expression across tissues, comparable gene length
and a comparable density of neighbouring genes. The control list is used as
the reference set in the downstream analyses of the pipeline.

Inputs
------
  --genes   plain text, one gene symbol per line. The set under study, for
            example the union of MitoCarta3.0 and the PanelApp
            "Mitochondrial disorders" panel. Symbols are matched against the
            GTF gene_name and the GCT Description column.
  --gtf     GENCODE GTF (gzip). Lines of feature type "gene" are read; the
            attributes gene_id, gene_name and gene_type are used.
  --tpm     GTEx median TPM per tissue, GCT format (gzip): two header lines,
            then a table with a symbol column (--tpm-symbol-column, default
            "Description") and one numeric column per tissue. Every numeric
            column is taken as a tissue.
  The GTEx and GENCODE files are the ones `gtex_eqtl.py fetch` downloads.

Outputs
-------
  --out      one control gene symbol per line, sorted (default: controls.txt)
  --pairs    CSV with mito_gene, control_gene, distance (distance in
             standardised feature space; control_gene is empty and distance
             NaN when the pool ran out for that gene) (default:
             control_pairs.csv)
  --report   TSV with set, n, median_log_tpm, median_tpm_sd,
             median_length_kb, median_density, one row for the gene set and
             one for the controls (default: controls_check.tsv)

Method
------
Genes of the selected gene type (--gene-type, protein_coding by default) are
taken from the GTF, with the Y pseudo-autosomal copies ("_PAR_Y") left out
and one entry kept per symbol. The TSS is the gene start on the + strand and
the gene end on the - strand. Four features are computed per gene: the
median across tissues of log10(TPM + 1), its standard deviation across
tissues, log10(length) and log10(1 + number of TSS within --density-window).
Genes whose maximum TPM over tissues does not exceed --min-tpm are dropped.

Cross-tissue variability is matched because a gene expressed at a similar
level in every tissue, as many mitochondrial genes are, would otherwise be
paired with genes expressed in a few tissues only. Their profiles would then
differ across tissues for reasons unrelated to the locus.

Candidate controls whose TSS lies within --exclusion of a gene-set TSS are
removed. Exclusion uses a short distance and not a whole viewpoint window:
gene sets of a thousand genes or more lie on average only a few megabases
apart, so requiring non-overlapping +/- --window windows would leave almost no
pool. Overlapping windows are accepted and quantified instead: the script
reports the share of control windows that contain at least one gene-set TSS,
which measures how much the controls are diluted by gene-set promoters.

Features are z-scored on the union of gene set and pool. Gene-set genes are
visited in random order (--seed), and each takes its --controls nearest
unused pool genes in Euclidean distance.

Usage
-----
    python3 make_controls.py --gtf data/gencode.gtf.gz \\
        --tpm data/gtex_median_tpm.gct.gz --genes genes.txt \\
        --out controls.txt --report controls_check.tsv
    python3 make_controls.py --gtf data/gencode.gtf.gz \\
        --tpm data/gtex_median_tpm.gct.gz --genes genes.txt \\
        --controls 3 --exclusion 250000 --seed 1

References
----------
  GTEx Consortium 2020. The GTEx Consortium atlas of genetic regulatory
    effects across human tissues. Science. doi:10.1126/science.aaz1776
  Frankish A et al. 2023. GENCODE: reference annotation for the human and
    mouse genomes in 2023. Nucleic Acids Research. doi:10.1093/nar/gkac1071
  Rath S et al. 2021. MitoCarta3.0: an updated mitochondrial proteome now
    with sub-organelle localization and pathway annotations. Nucleic Acids
    Research. doi:10.1093/nar/gkaa1011
  Martin AR et al. 2019. PanelApp crowdsources expert knowledge to establish
    consensus diagnostic gene panels. Nature Genetics.
    doi:10.1038/s41588-019-0528-2
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas required")


DEFAULT_CONTROLS_PER_GENE = 2
DEFAULT_EXCLUSION = 100_000
DEFAULT_DENSITY_WINDOW = 1_000_000
DEFAULT_VIEW_WINDOW = 1_000_000

# A pool smaller than this multiple of the number of controls requested is
# flagged, since late gene-set genes then get distant matches.
THIN_POOL_FACTOR = 3


def load_genes(gtf: Path, gene_type: str, density_window: int) -> pd.DataFrame:
    type_tag = f'gene_type "{gene_type}"'
    rows = []
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            a = f[8]
            if type_tag not in a or f[0].endswith("_PAR_Y"):
                continue
            i = a.find('gene_id "')
            gid = a[i + 9:a.find('"', i + 9)].split(".")[0]
            j = a.find('gene_name "')
            name = a[j + 11:a.find('"', j + 11)]
            start, end = int(f[3]) - 1, int(f[4])
            rows.append({"gene_id": gid, "gene_name": name, "chrom": f[0],
                         "start": start, "end": end,
                         "tss": start if f[6] == "+" else end,
                         "length": end - start})
    if not rows:
        sys.exit(f"no gene of type {gene_type!r} in {gtf}")
    g = pd.DataFrame(rows).drop_duplicates("gene_name").reset_index(drop=True)

    dens = []
    for _, sub in g.groupby("chrom"):
        t = np.sort(sub["tss"].to_numpy())
        v = sub["tss"].to_numpy()
        lo = np.searchsorted(t, v - density_window, "left")
        hi = np.searchsorted(t, v + density_window, "right")
        dens.append(pd.Series(hi - lo, index=sub.index))
    g["density"] = pd.concat(dens).reindex(g.index)
    return g


def load_tpm_matrix(path: Path, symbol_col: str):
    """Median TPM per gene per tissue, indexed by Ensembl id without version
    (the Name column) and by gene symbol: the ids are stable across GENCODE
    releases, the symbols are not (GTEx v8 is GENCODE v26)."""
    d = pd.read_csv(path, sep="\t", skiprows=2, low_memory=False)
    if symbol_col not in d.columns:
        sys.exit(f"expected a {symbol_col} column in the GCT")
    num = d.select_dtypes("number")
    by_sym = num.set_axis(d[symbol_col].astype(str))
    by_sym = by_sym[~by_sym.index.duplicated()]
    by_id = None
    if "Name" in d.columns:
        by_id = num.set_axis(d["Name"].astype(str).str.split(".").str[0])
        by_id = by_id[~by_id.index.duplicated()]
    return by_id, by_sym


def exclude_near(pool: pd.DataFrame, mito: pd.DataFrame,
                 span: int) -> pd.DataFrame:
    """Drop pool genes whose TSS is closer than span to a gene-set TSS."""
    keep = np.ones(len(pool), dtype=bool)
    for chrom, sub in mito.groupby("chrom"):
        idx = np.where(pool["chrom"].to_numpy() == chrom)[0]
        if not len(idx):
            continue
        mt = np.sort(sub["tss"].to_numpy())
        pt = pool["tss"].to_numpy()[idx]
        pos = np.searchsorted(mt, pt)
        near = np.zeros(len(idx), dtype=bool)
        # nearest gene-set TSS is either just before or just after
        for off in (-1, 0):
            j = np.clip(pos + off, 0, len(mt) - 1)
            near |= np.abs(mt[j] - pt) < span
        keep[idx[near]] = False
    return pool[keep]


def window_overlap(pool: pd.DataFrame, mito: pd.DataFrame,
                   window: int) -> pd.Series:
    """How many gene-set TSS fall inside each control window."""
    out = pd.Series(0, index=pool.index, dtype=int)
    for chrom, sub in mito.groupby("chrom"):
        idx = pool.index[pool["chrom"] == chrom]
        if not len(idx):
            continue
        mt = np.sort(sub["tss"].to_numpy())
        pt = pool.loc[idx, "tss"].to_numpy()
        lo = np.searchsorted(mt, pt - window, "left")
        hi = np.searchsorted(mt, pt + window, "right")
        out.loc[idx] = hi - lo
    return out


def match(mito: pd.DataFrame, pool: pd.DataFrame, feats: list[str],
          k: int, rng) -> pd.DataFrame:
    mu = pd.concat([mito[feats], pool[feats]]).mean()
    sd = pd.concat([mito[feats], pool[feats]]).std().replace(0, 1)
    zm = ((mito[feats] - mu) / sd).to_numpy()
    zc = ((pool[feats] - mu) / sd).to_numpy()
    mito_names = mito["gene_name"].to_numpy()
    pool_names = pool["gene_name"].to_numpy()

    used, pairs = set(), []
    for i in rng.permutation(len(zm)):
        dist = np.sqrt(((zc - zm[i]) ** 2).sum(1))
        taken = 0
        for j in np.argsort(dist):
            if pool_names[j] in used:
                continue
            used.add(pool_names[j])
            pairs.append({"mito_gene": mito_names[i],
                          "control_gene": pool_names[j],
                          "distance": float(dist[j])})
            taken += 1
            if taken >= k:
                break
        if taken < k:
            pairs.append({"mito_gene": mito_names[i], "control_gene": "",
                          "distance": float("nan")})
    return pd.DataFrame(pairs, columns=["mito_gene", "control_gene", "distance"])


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True, type=Path,
                    help="GENCODE GTF, gzip")
    ap.add_argument("--tpm", required=True, type=Path,
                    help="GTEx median TPM per tissue, GCT, gzip")
    ap.add_argument("--genes", required=True, type=Path,
                    help="gene-set symbols, one per line")
    ap.add_argument("--out", default="controls.txt",
                    help="control list, one symbol per line "
                         "(default: controls.txt)")
    ap.add_argument("--pairs", default="control_pairs.csv",
                    help="matched pairs (default: control_pairs.csv)")
    ap.add_argument("--report", default="controls_check.tsv",
                    help="matching check (default: controls_check.tsv)")
    ap.add_argument("--controls", type=int, default=DEFAULT_CONTROLS_PER_GENE,
                    help="controls drawn per gene-set gene "
                         f"(default: {DEFAULT_CONTROLS_PER_GENE})")
    ap.add_argument("--exclusion", type=int, default=DEFAULT_EXCLUSION,
                    help="minimum distance in bp between a control TSS and "
                         "the nearest gene-set TSS "
                         f"(default: {DEFAULT_EXCLUSION})")
    ap.add_argument("--window", type=int, default=DEFAULT_VIEW_WINDOW,
                    help="viewpoint half-window in bp, used only to report "
                         f"overlap (default: {DEFAULT_VIEW_WINDOW})")
    ap.add_argument("--density-window", type=int, default=DEFAULT_DENSITY_WINDOW,
                    help="half-window in bp for local gene density "
                         f"(default: {DEFAULT_DENSITY_WINDOW})")
    ap.add_argument("--gene-type", default="protein_coding",
                    help="GTF gene_type kept (default: protein_coding)")
    ap.add_argument("--tpm-symbol-column", default="Description",
                    help="gene symbol column of the GCT "
                         "(default: Description)")
    ap.add_argument("--min-tpm", type=float, default=0.1,
                    help="genes need a maximum TPM over tissues above this "
                         "(default: 0.1)")
    ap.add_argument("--seed", type=int, default=0,
                    help="seed of the matching order (default: 0)")
    args = ap.parse_args()

    wanted = set(Path(args.genes).read_text().split())
    genes = load_genes(args.gtf, args.gene_type, args.density_window)
    tpm_id, tpm_sym = load_tpm_matrix(args.tpm, args.tpm_symbol_column)
    print(f"{len(genes):,} {args.gene_type} genes, "
          f"{tpm_sym.shape[1]} tissues in the TPM matrix")

    # per gene: by Ensembl id, else by symbol
    def per_gene(f):
        v = genes["gene_name"].map(f(tpm_sym))
        if tpm_id is not None:
            v = genes["gene_id"].map(f(tpm_id)).fillna(v)
        return v
    genes["log_tpm"] = per_gene(lambda t: np.log10(t + 1).median(axis=1))
    genes["tpm_sd"] = per_gene(lambda t: np.log10(t + 1).std(axis=1))
    genes["tpm_max"] = per_gene(lambda t: t.max(axis=1))
    genes["log_len"] = np.log10(genes["length"].clip(lower=1))
    genes["log_dens"] = np.log10(genes["density"] + 1)
    feats = ["log_tpm", "tpm_sd", "log_len", "log_dens"]

    genes = genes.dropna(subset=feats)
    genes = genes[genes["tpm_max"] > args.min_tpm]

    is_mito = genes["gene_name"].isin(wanted)
    mito = genes[is_mito].copy()
    print(f"{len(mito):,} gene-set genes matched "
          f"(of {len(wanted):,} symbols given)")
    missing = wanted - set(mito["gene_name"])
    if missing:
        print(f"  {len(missing)} symbols absent from the annotation or below "
              f"the expression floor, e.g. {sorted(missing)[:6]}")
    if mito.empty:
        sys.exit("no gene of the set passed the filters; nothing to match")

    pool = genes[~is_mito].copy()
    before = len(pool)
    pool = exclude_near(pool, mito, args.exclusion)
    print(f"control pool {before:,} -> {len(pool):,} after removing genes "
          f"within {args.exclusion/1000:g} kb of a gene-set TSS")
    if len(pool) < THIN_POOL_FACTOR * args.controls * len(mito):
        print("  !! pool is thin relative to the number of controls requested")

    rng = np.random.default_rng(args.seed)
    pairs = match(mito, pool, feats, args.controls, rng)
    for target in (args.out, args.pairs, args.report):
        Path(target).parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.pairs, index=False)
    ctrl = sorted({c for c in pairs["control_gene"] if c})
    Path(args.out).write_text("\n".join(ctrl) + "\n")

    unmatched = int((pairs["control_gene"] == "").sum())
    print(f"\n{len(ctrl):,} control genes -> {args.out}")
    if unmatched:
        print(f"  {unmatched} gene-set genes short of {args.controls} "
              f"controls; the pool ran out")

    ov = window_overlap(pool[pool["gene_name"].isin(ctrl)], mito, args.window)
    print(f"\ncontrol windows of +/-{args.window/1000:g} kb containing at least "
          f"one gene-set TSS: {(ov > 0).mean():.1%}, median count "
          f"{ov.median():.0f}")
    print("This measures the dilution of the control set by gene-set "
          "promoters. A larger --exclusion reduces it.")

    rep = []
    for lab, sub in (("mito", mito),
                     ("control", pool[pool["gene_name"].isin(ctrl)])):
        rep.append({"set": lab, "n": len(sub),
                    "median_log_tpm": sub["log_tpm"].median(),
                    "median_tpm_sd": sub["tpm_sd"].median(),
                    "median_length_kb": sub["length"].median() / 1e3,
                    "median_density": sub["density"].median()})
    rep = pd.DataFrame(rep)
    rep.to_csv(args.report, sep="\t", index=False)
    print("\n" + rep.round(3).to_string(index=False))
    print("\nThe two rows should be close on every column; otherwise the "
          "gene set and the controls are not comparable.")


if __name__ == "__main__":
    main()
