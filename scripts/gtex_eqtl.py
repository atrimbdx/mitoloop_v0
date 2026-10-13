#!/usr/bin/env python3
"""
gtex_eqtl.py - effect size and TSS distance of independent cis-eQTL in a gene
set against matched control genes.

Purpose
-------
Asks whether a gene set (for example nuclear-encoded mitochondrial genes)
carries distal regulation of ordinary strength and reach in one tissue.
Conditionally independent GTEx cis-eQTL signals of the gene set are compared
with those of control genes matched on expression, gene length and local gene
density, on two quantities:

  effect size   |slope| of each independent signal
  distance      |distance to the TSS| of each independent signal

Inputs
------
  --genes   plain text, one gene symbol per line (whitespace-separated
            symbols are also accepted). The set under study, for example the
            union of MitoCarta3.0 and the PanelApp "Mitochondrial disorders"
            panel. Symbols are matched against the GTF gene_name.
  --gtf     GENCODE GTF (gzip). Lines of feature type "gene" are read; the
            attributes gene_id, gene_name and gene_type are used. Gene ids are
            compared without their version suffix.
  --tpm     GTEx median TPM per tissue, GCT format (gzip): two header lines,
            then a table with the id column (--tpm-id-column, default "Name",
            versioned Ensembl ids), "Description" and one column per tissue.
            The tissue column is chosen with --tpm-column.
  --eqtl    GTEx conditionally independent cis-eQTL. Either the release tar
            holding one gzipped table per tissue (the member whose name
            contains --tissue and ends in ".gz" is read), or one gzipped
            per-tissue table. Columns used, matched case-insensitively:
            gene id (--eqtl-gene-column), signed distance to the TSS
            (--eqtl-distance-column) and slope (--eqtl-slope-column).
  --datadir optional; supplies any of --gtf, --tpm, --eqtl left unset, under
            the file names written by `fetch`.

Outputs (in --outdir)
---------------------
  gtex_matched_pairs.csv   mito_gene, control_gene, distance (Ensembl ids;
                           distance in standardised feature space)
  gtex_matching_check.csv  set, n, median_tpm, median_length_kb,
                           median_density
  gtex_eqtl_signals.csv    gene_id, abs_slope, abs_dist, tss_distance, slope,
                           set (one row per independent signal)
  gtex_summary.csv         quantity, n_mito, n_control, median_mito,
                           median_control, cohens_d, prob_superiority, overlap
  gtex_per_gene.csv        set, gene_id, n_signals, max_dist, frac_distal
                           (share of signals beyond --distal-cutoff)

Method
------
Genes of the selected gene type (--gene-type, protein_coding by default)
are taken from the GTF, with genes on the Y pseudo-autosomal copies ("_PAR_Y"
chromosomes) left out. The TSS is the gene start on the + strand and the gene
end on the - strand. Local density is the number of TSS of the same gene type
within --density-window of each TSS, a proxy for how crowded the regulatory
neighbourhood is. Genes with median TPM above --min-tpm in the tissue are
kept.

eQTL discovery power increases with expression level, and mitochondrial genes
are more highly and more uniformly expressed than the genomic average. A
comparison with all other genes would therefore reflect a difference in
power. Controls are drawn by greedy nearest-neighbour matching on z-scored
log10(TPM + 1), log10(length) and log10(density + 1): gene-set genes are
visited in random order (--seed) and each takes its --controls nearest
unused genes in Euclidean distance.

The question is whether the two sets differ, so the comparison is reported
as effect sizes and distribution overlap. With thousands of signals per set,
p values would reach significance for shifts too small to matter. Reported
per quantity: the probability of superiority, P(mito > control) +
0.5 P(equal), equal to the Mann-Whitney U statistic divided by n*m; the
overlap coefficient of the two histograms on a common grid (--overlap-bins);
Cohen's d. Distances are log10(d + 1) transformed for d and the overlap.

`fetch` downloads the three inputs and is the only step that uses the
network. URLs are options; each list is tried in order.

Usage
-----
    python3 gtex_eqtl.py fetch --datadir data
    python3 gtex_eqtl.py run --gtf data/gencode.gtf.gz \\
        --tpm data/gtex_median_tpm.gct.gz \\
        --eqtl data/GTEx_eQTL_independent.tar \\
        --genes genes.txt --outdir results/eqtl
    python3 gtex_eqtl.py run --datadir data --genes genes.txt \\
        --tissue Heart_Left_Ventricle --tpm-column "Heart - Left Ventricle"

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
  Mann HB, Whitney DR 1947. On a test of whether one of two random variables
    is stochastically larger than the other. Annals of Mathematical
    Statistics.
  McGraw KO, Wong SP 1992. A common language effect size statistic.
    Psychological Bulletin.
  Inman HF, Bradley EL 1989. The overlapping coefficient as a measure of
    agreement between probability distributions and point estimation of the
    overlap of two normal densities. Communications in Statistics - Theory
    and Methods.
  Cohen J 1988. Statistical Power Analysis for the Behavioral Sciences.
    Lawrence Erlbaum.
"""

from __future__ import annotations

import argparse
import gzip
import sys
import tarfile
from pathlib import Path

import numpy as np

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas required")


DEFAULT_TISSUE = "Muscle_Skeletal"
DEFAULT_TISSUE_TPM_COLUMN = "Muscle - Skeletal"
DEFAULT_CONTROLS_PER_GENE = 3
DEFAULT_DISTAL_CUTOFF = 20_000
DEFAULT_DENSITY_WINDOW = 1_000_000

# GTEx v8 release files; the second URL of each pair is the older bucket.
DEFAULT_EQTL_URLS = [
    "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/"
    "GTEx_Analysis_v8_eQTL_independent.tar",
    "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/"
    "GTEx_Analysis_v8_eQTL_independent.tar",
]
DEFAULT_TPM_URLS = [
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
    "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
]
DEFAULT_GENCODE_RELEASES = [49, 48, 47, 46, 45]
GENCODE_URL_TEMPLATE = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{v}/"
    "gencode.v{v}.basic.annotation.gtf.gz")

# File names written by `fetch` and looked up under --datadir by `run`.
EQTL_FILE = "GTEx_eQTL_independent.tar"
TPM_FILE = "gtex_median_tpm.gct.gz"
GTF_FILE = "gencode.gtf.gz"


# ----------------------------------------------------------------------------

def download(urls, dest: Path, label: str, timeout: int) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  [skip] {dest.name} ({dest.stat().st_size/1e6:.0f} MB)")
        return dest
    import requests
    dest.parent.mkdir(parents=True, exist_ok=True)
    for url in urls:
        print(f"  [try ] {url}")
        try:
            r = requests.get(url, stream=True, timeout=timeout)
            if r.status_code != 200:
                print(f"         HTTP {r.status_code}")
                continue
            tmp = dest.with_suffix(dest.suffix + ".part")
            with open(tmp, "wb") as fh:
                for c in r.iter_content(1 << 20):
                    fh.write(c)
            tmp.rename(dest)
            print(f"  [ok  ] {dest.name}")
            return dest
        except Exception as exc:
            print(f"         {type(exc).__name__}")
    sys.exit(f"Could not download {label}; fetch it manually into {dest.parent}.")


def load_gene_table(gtf: Path, gene_type: str,
                    density_window: int) -> pd.DataFrame:
    """Genes of one type with length and local density, keyed by Ensembl id."""
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
    g = pd.DataFrame(rows).drop_duplicates("gene_id")

    # local density: TSS within the window, a proxy for how crowded the
    # regulatory neighbourhood is
    dens = []
    for chrom, sub in g.groupby("chrom"):
        t = np.sort(sub["tss"].to_numpy())
        lo = np.searchsorted(t, sub["tss"].to_numpy() - density_window, "left")
        hi = np.searchsorted(t, sub["tss"].to_numpy() + density_window, "right")
        dens.append(pd.Series(hi - lo, index=sub.index))
    g["density"] = pd.concat(dens).reindex(g.index)
    return g


def read_eqtl_table(path: Path, tissue: str) -> pd.DataFrame:
    """Read the per-tissue table, from the release tar or a single file."""
    if tarfile.is_tarfile(path):
        with tarfile.open(path) as tf:
            member = next((m for m in tf.getmembers()
                           if tissue in m.name and m.name.endswith(".gz")), None)
            if member is None:
                sys.exit(f"{tissue} not found in {path.name}. Members:\n  " +
                         "\n  ".join(m.name for m in tf.getmembers()[:10]))
            print(f"  reading {member.name}")
            return pd.read_csv(tf.extractfile(member), sep="\t",
                               compression="gzip", low_memory=False)
    print(f"  reading {path.name}")
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def load_eqtl(path: Path, tissue: str, gene_col: str, dist_col: str,
              slope_col: str) -> pd.DataFrame:
    """Conditionally independent signals for one tissue."""
    d = read_eqtl_table(path, tissue)
    cols = {c.lower(): c for c in d.columns}
    wanted = {"gene_id": gene_col, "tss_distance": dist_col, "slope": slope_col}
    need = {k: cols[v.lower()] for k, v in wanted.items() if v.lower() in cols}
    if len(need) < 3:
        sys.exit(f"expected columns {list(wanted.values())}; "
                 f"found {list(d.columns)[:12]}")
    d = d[list(need.values())].rename(columns={v: k for k, v in need.items()})
    d["gene_id"] = d["gene_id"].astype(str).str.split(".").str[0]
    d["abs_slope"] = d["slope"].abs()
    d["abs_dist"] = d["tss_distance"].abs()
    return d[["gene_id", "abs_slope", "abs_dist", "tss_distance", "slope"]]


def load_tpm(path: Path, col: str, id_col: str) -> pd.Series:
    d = pd.read_csv(path, sep="\t", skiprows=2, low_memory=False)
    name = next((c for c in d.columns if c.strip() == col), None)
    if name is None:
        sys.exit(f"column {col!r} absent. Available: {list(d.columns)[:8]}")
    if id_col not in d.columns:
        sys.exit(f"id column {id_col!r} absent from {path}")
    s = pd.Series(d[name].to_numpy(),
                  index=d[id_col].astype(str).str.split(".").str[0])
    return s[~s.index.duplicated()]


def match_controls(genes: pd.DataFrame, mito: set, k: int,
                   rng) -> pd.DataFrame:
    """Nearest-neighbour controls on standardised expression, length, density.

    Matching makes an absence of difference interpretable: the comparison is
    otherwise dominated by the higher expression of mitochondrial genes and
    therefore by eQTL discovery power.
    """
    feats = ["log_tpm", "log_len", "log_dens"]
    d = genes.dropna(subset=feats).copy()
    # a constant feature carries no information; a unit scale keeps it finite
    z = (d[feats] - d[feats].mean()) / d[feats].std().replace(0, 1)
    is_mito = d["gene_id"].isin(mito).to_numpy()
    zm, zc = z[is_mito].to_numpy(), z[~is_mito].to_numpy()
    cid = d.loc[~is_mito, "gene_id"].to_numpy()
    if not len(zm) or not len(zc):
        sys.exit("empty mito or control pool after filtering")

    mito_ids = d.loc[is_mito, "gene_id"].to_numpy()
    used, pairs = set(), []
    for i in rng.permutation(len(zm)):
        dist = np.sqrt(((zc - zm[i]) ** 2).sum(1))
        taken = 0
        for j in np.argsort(dist):
            if cid[j] in used:
                continue
            used.add(cid[j])
            pairs.append({"mito_gene": mito_ids[i], "control_gene": cid[j],
                          "distance": float(dist[j])})
            taken += 1
            if taken >= k:
                break
    return pd.DataFrame(pairs, columns=["mito_gene", "control_gene", "distance"])


def overlap_coefficient(a: np.ndarray, b: np.ndarray, bins: int = 80) -> float:
    """Area shared by two densities on a common grid. 1 = identical."""
    lo, hi = min(a.min(), b.min()), max(a.max(), b.max())
    e = np.linspace(lo, hi, bins + 1)
    ha = np.histogram(a, bins=e, density=True)[0]
    hb = np.histogram(b, bins=e, density=True)[0]
    return float(np.minimum(ha, hb).sum() * (e[1] - e[0]))


def prob_superiority(a: np.ndarray, b: np.ndarray) -> float:
    """P(a > b) + 0.5 P(a = b): the chance that a signal drawn from the first
    set exceeds one drawn from the second. 0.5 means the two are
    interchangeable. It is the Mann-Whitney statistic divided by n*m, and it
    needs no reference table to be read."""
    from scipy.stats import rankdata
    n, m = len(a), len(b)
    if not n or not m:
        return float("nan")
    r = rankdata(np.concatenate([a, b]))
    return float((r[:n].sum() - n * (n + 1) / 2) / (n * m))


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    s = np.sqrt(((len(a) - 1) * a.var(ddof=1) +
                 (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / s) if s else 0.0


# ----------------------------------------------------------------------------

def cmd_fetch(args) -> None:
    d = args.datadir
    d.mkdir(parents=True, exist_ok=True)
    eqtl_urls = args.eqtl_url or DEFAULT_EQTL_URLS
    tpm_urls = args.tpm_url or DEFAULT_TPM_URLS
    gencode_urls = args.gencode_url or [
        GENCODE_URL_TEMPLATE.format(v=v)
        for v in (args.gencode_release or DEFAULT_GENCODE_RELEASES)]
    print("GTEx independent eQTL")
    download(eqtl_urls, d / EQTL_FILE, "the eQTL archive", args.timeout)
    print("GTEx median TPM")
    download(tpm_urls, d / TPM_FILE, "the TPM matrix", args.timeout)
    print("GENCODE")
    download(gencode_urls, d / GTF_FILE, "GENCODE", args.timeout)


def resolve_input(value: Path | None, datadir: Path | None, name: str,
                  option: str) -> Path:
    path = value if value is not None else (datadir / name if datadir else None)
    if path is None:
        sys.exit(f"{option} is required (or --datadir holding {name}).")
    if not path.exists():
        sys.exit(f"Missing {path}. Run `fetch` first or pass {option}.")
    return path


def cmd_run(args) -> None:
    gtf = resolve_input(args.gtf, args.datadir, GTF_FILE, "--gtf")
    tpm_path = resolve_input(args.tpm, args.datadir, TPM_FILE, "--tpm")
    eqtl_path = resolve_input(args.eqtl, args.datadir, EQTL_FILE, "--eqtl")

    wanted = set(Path(args.genes).read_text().split())
    if not wanted:
        sys.exit(f"{args.genes} contains no gene symbol.")

    print("Loading")
    genes = load_gene_table(gtf, args.gene_type, args.density_window)
    tpm = load_tpm(tpm_path, args.tpm_column, args.tpm_id_column)
    genes["tpm"] = genes["gene_id"].map(tpm)
    genes = genes[genes["tpm"].notna() & (genes["tpm"] > args.min_tpm)].copy()
    genes["log_tpm"] = np.log10(genes["tpm"] + 1)
    genes["log_len"] = np.log10(genes["length"])
    genes["log_dens"] = np.log10(genes["density"] + 1)

    mito_ids = set(genes.loc[genes["gene_name"].isin(wanted), "gene_id"])
    print(f"  {len(genes):,} expressed {args.gene_type} genes")
    print(f"  {len(mito_ids):,} gene-set genes matched "
          f"(of {len(wanted):,} symbols given)")

    eqtl = load_eqtl(eqtl_path, args.tissue, args.eqtl_gene_column,
                     args.eqtl_distance_column, args.eqtl_slope_column)
    print(f"  {len(eqtl):,} independent signals in {args.tissue}")

    rng = np.random.default_rng(args.seed)
    pairs = match_controls(genes, mito_ids, args.controls, rng)
    ctrl_ids = set(pairs["control_gene"])
    print(f"  {len(ctrl_ids):,} matched control genes "
          f"({args.controls} per gene-set gene requested)")

    args.outdir.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.outdir / "gtex_matched_pairs.csv", index=False)

    # matching quality: if these rows differ, the comparison below is not
    # interpretable
    chk = []
    for lab, ids in (("mito", mito_ids), ("control", ctrl_ids)):
        sub = genes[genes["gene_id"].isin(ids)]
        chk.append({"set": lab, "n": len(sub),
                    "median_tpm": sub["tpm"].median(),
                    "median_length_kb": sub["length"].median() / 1e3,
                    "median_density": sub["density"].median()})
    chk = pd.DataFrame(chk)
    chk.to_csv(args.outdir / "gtex_matching_check.csv", index=False)
    print("\n" + chk.round(2).to_string(index=False))

    eqtl["set"] = np.where(eqtl["gene_id"].isin(mito_ids), "mito",
                  np.where(eqtl["gene_id"].isin(ctrl_ids), "control", None))
    sig = eqtl[eqtl["set"].notna()].copy()
    sig.to_csv(args.outdir / "gtex_eqtl_signals.csv", index=False)

    m = sig[sig.set == "mito"]
    c = sig[sig.set == "control"]
    if not len(m) or not len(c):
        sys.exit("one of the sets has no eQTL signal")

    summ = []
    for lab, col in (("abs_slope", "abs_slope"), ("abs_dist", "abs_dist")):
        a, b = m[col].to_numpy(), c[col].to_numpy()
        # distances span orders of magnitude: compare them on a log scale
        ta = np.log10(a + 1) if lab == "abs_dist" else a
        tb = np.log10(b + 1) if lab == "abs_dist" else b
        summ.append({
            "quantity": lab, "n_mito": len(a), "n_control": len(b),
            "median_mito": float(np.median(a)), "median_control": float(np.median(b)),
            "cohens_d": cohens_d(ta, tb),
            "prob_superiority": prob_superiority(a, b),
            "overlap": overlap_coefficient(ta, tb, args.overlap_bins),
        })
    summ = pd.DataFrame(summ)
    summ.to_csv(args.outdir / "gtex_summary.csv", index=False)

    cutoff = args.distal_cutoff
    per_gene = (sig.groupby(["set", "gene_id"])
                .agg(n_signals=("abs_slope", "size"),
                     max_dist=("abs_dist", "max"),
                     frac_distal=("abs_dist", lambda x: (x > cutoff).mean()))
                .reset_index())
    per_gene.to_csv(args.outdir / "gtex_per_gene.csv", index=False)

    print("\n" + "=" * 70)
    print(f"{'quantity':<12}{'median mito':>14}{'median ctrl':>14}"
          f"{'P(mito>ctrl)':>14}{'overlap':>10}")
    for _, r in summ.iterrows():
        print(f"{r['quantity']:<12}{r.median_mito:>14.4g}{r.median_control:>14.4g}"
              f"{r.prob_superiority:>14.3f}{r.overlap:>10.2f}")
    print("=" * 70)
    print("A probability of superiority near 0.5 and an overlap near 1 indicate")
    print("similar distributions. Cohen's d is also given in gtex_summary.csv.")
    print("P values are not reported: with many signals per set they reach")
    print("significance for negligible shifts.")
    for lab, s in (("mito", m), ("control", c)):
        print(f"  {lab:<8} beyond {cutoff/1000:g} kb: "
              f"{(s.abs_dist > cutoff).mean():.1%}   "
              f"max {s.abs_dist.max()/1e3:.0f} kb")
    print(f"\nWritten to {args.outdir}/")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="download the inputs (needs network)")
    f.add_argument("--datadir", type=Path, default=Path("data"),
                   help="output directory (default: data)")
    f.add_argument("--eqtl-url", action="append",
                   help="URL of the independent cis-eQTL tar; repeat to give "
                        "fallbacks (default: GTEx v8, two mirrors)")
    f.add_argument("--tpm-url", action="append",
                   help="URL of the median TPM GCT; repeat to give fallbacks "
                        "(default: GTEx v8, two mirrors)")
    f.add_argument("--gencode-url", action="append",
                   help="URL of the GENCODE GTF; repeat to give fallbacks "
                        "(overrides --gencode-release)")
    f.add_argument("--gencode-release", action="append", type=int,
                   help="GENCODE human release to try, basic annotation; "
                        "repeat for fallbacks (default: "
                        f"{' '.join(map(str, DEFAULT_GENCODE_RELEASES))})")
    f.add_argument("--timeout", type=int, default=600,
                   help="HTTP timeout in seconds (default: 600)")

    p = sub.add_parser("run", help="compare gene-set genes to matched controls "
                                   "(offline)")
    p.add_argument("--genes", required=True, type=Path,
                   help="file of gene symbols, one per line")
    p.add_argument("--gtf", type=Path, default=None,
                   help="GENCODE GTF, gzip")
    p.add_argument("--tpm", type=Path, default=None,
                   help="GTEx median TPM GCT, gzip")
    p.add_argument("--eqtl", type=Path, default=None,
                   help="GTEx independent cis-eQTL: release tar or one "
                        "per-tissue table")
    p.add_argument("--datadir", type=Path, default=None,
                   help=f"directory holding {GTF_FILE}, {TPM_FILE} and "
                        f"{EQTL_FILE}, used for inputs not given explicitly")
    p.add_argument("--outdir", type=Path, default=Path("results"),
                   help="output directory (default: results)")
    p.add_argument("--tissue", default=DEFAULT_TISSUE,
                   help="tissue name as it appears in the eQTL file names "
                        f"(default: {DEFAULT_TISSUE})")
    p.add_argument("--tpm-column", default=DEFAULT_TISSUE_TPM_COLUMN,
                   help="tissue column of the TPM GCT "
                        f"(default: {DEFAULT_TISSUE_TPM_COLUMN!r})")
    p.add_argument("--tpm-id-column", default="Name",
                   help="Ensembl id column of the TPM GCT (default: Name)")
    p.add_argument("--eqtl-gene-column", default="gene_id",
                   help="gene id column of the eQTL table (default: gene_id)")
    p.add_argument("--eqtl-distance-column", default="tss_distance",
                   help="signed distance-to-TSS column of the eQTL table "
                        "(default: tss_distance)")
    p.add_argument("--eqtl-slope-column", default="slope",
                   help="effect size column of the eQTL table (default: slope)")
    p.add_argument("--gene-type", default="protein_coding",
                   help="GTF gene_type kept (default: protein_coding)")
    p.add_argument("--density-window", type=int, default=DEFAULT_DENSITY_WINDOW,
                   help="half-window in bp for local gene density "
                        f"(default: {DEFAULT_DENSITY_WINDOW})")
    p.add_argument("--controls", type=int, default=DEFAULT_CONTROLS_PER_GENE,
                   help="controls drawn per gene-set gene "
                        f"(default: {DEFAULT_CONTROLS_PER_GENE})")
    p.add_argument("--min-tpm", type=float, default=0.1,
                   help="genes need median TPM above this in the tissue "
                        "(default: 0.1)")
    p.add_argument("--distal-cutoff", type=int, default=DEFAULT_DISTAL_CUTOFF,
                   help="distance to the TSS in bp beyond which a signal is "
                        f"called distal (default: {DEFAULT_DISTAL_CUTOFF})")
    p.add_argument("--overlap-bins", type=int, default=80,
                   help="histogram bins for the overlap coefficient "
                        "(default: 80)")
    p.add_argument("--seed", type=int, default=0,
                   help="seed of the matching order (default: 0)")
    args = ap.parse_args()
    {"fetch": cmd_fetch, "run": cmd_run}[args.cmd](args)


if __name__ == "__main__":
    main()
