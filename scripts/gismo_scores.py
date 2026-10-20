#!/usr/bin/env python3
"""
gismo_scores.py - gene table of GISMO and GISMO-mis (Liao et al. 2026) from the
Zenodo record 20547722, which holds the code and its inputs but no score table.

Purpose      One row per gene with GISMO (loss of the ortholog across 462
             mammals) and GISMO-mis (missense over synonymous substitutions
             across orthologs), and their deciles, for mitoloop_paquet.py.
Inputs       the record's archive GISMO-GISMO.zip (or its extracted folder):
             GenerateScore/inputs/unmerged-species_combined_matrix_*.tsv
             (one2one, one2zero... per gene and species),
             GenerateScore/inputs/gen.dist.withassembly.txt (species,
             assembly, lin.num, div.time, contig_N50_bp),
             partitioned-heritability/inputs/gismo-mis_consensus_genenames_*.txt
             (per transcript: mean.mis, mean.syn, mis.cal, dec, genename);
             a GENCODE GTF for gene symbols and the MANE Select transcript.
Outputs      gismo_genes.tsv: gene_id, gene, gismo_obs, gismo, gismo_decile,
             gismo_mis, gismo_mis_decile, gismo_mis_transcript.
Method       GISMO as in GenerateScore/01_GISMOscores.R: one2zero = 1, any
             other state = 0; one genome per species (largest contig N50);
             loss rate p of each gene per lineage (lin.num); weights
             div.time / max(div.time); observed statistic = mean over species
             of weight x loss; null = the same statistic on Bernoulli(p)
             draws (--iterations, seed 1337), quantiles 2.5, 50 and 97.5 %.
             gismo = the 97.5 % quantile (upper.ci in the script), decile =
             dplyr::ntile of it (1 = fewest losses, most conserved).
             GISMO-mis = mis.cal and the decile 'dec' of the consensus table,
             on the MANE Select transcript of the gene when it is listed,
             else on its first transcript.
Usage        python3 scripts/gismo_scores.py --zip data/gismo/cliao5/GISMO-GISMO.zip \
                 --gtf data/gencode.gtf.gz --out data/gismo/gismo_genes.tsv
"""

from __future__ import annotations

import argparse
import gzip
import io
import re
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


class Source:
    """members of the archive, or files of the extracted folder, by pattern"""

    def __init__(self, path: Path):
        self.zip = zipfile.ZipFile(path) if path.is_file() else None
        self.names = (self.zip.namelist() if self.zip else
                      [str(p.relative_to(path)) for p in path.rglob("*") if p.is_file()])
        self.root = path

    def find(self, pattern: str) -> str:
        hits = sorted(n for n in self.names if re.search(pattern, n))
        if not hits:
            sys.exit(f"{pattern} absent de {self.root}")
        return hits[-1]

    def open(self, name: str):
        raw = self.zip.open(name) if self.zip else open(self.root / name, "rb")
        return gzip.open(raw) if name.endswith(".gz") else raw


def ntile(x, n=10):
    """dplyr::ntile: ranks by value, ties in order of appearance"""
    order = np.argsort(x, kind="stable")
    rank = np.empty(len(x), int)
    rank[order] = np.arange(len(x))
    return (rank * n // len(x)) + 1


def gismo(src: Source, iterations: int, seed: int) -> pd.DataFrame:
    stat = pd.read_csv(src.open(src.find(r"GenerateScore/inputs/gen\.dist\.withassembly\.txt$")),
                       sep="\t")
    for c in ("species", "assembly_name", "contig_N50_bp", "lin.num", "div.time"):
        if c not in stat.columns:
            sys.exit(f"gen.dist.withassembly.txt : colonne {c} absente ({list(stat.columns)[:12]}...)")
    stat = stat.assign(_n50=stat["contig_N50_bp"].abs()).sort_values(
        ["species", "_n50"], ascending=[True, False]).drop_duplicates("species")
    mfile = src.find(r"GenerateScore/inputs/unmerged-species_combined_matrix.*\.tsv$")
    with src.open(mfile) as fh:
        head = io.TextIOWrapper(fh).readline().rstrip("\n").split("\t")
    o2o = pd.read_csv(src.open(mfile), sep="\t", skiprows=1, header=None, index_col=0,
                      dtype=str, na_filter=False)
    if o2o.shape[1] != len(head):          # header with a first empty column
        head = head[1:] if head and head[0] == "" else head
    o2o.columns = [re.sub(r".*_", "", h) for h in head[-o2o.shape[1]:]]
    o2o = o2o.loc[:, ~o2o.columns.duplicated()]
    use = set(stat["assembly_name"]) & set(o2o.columns)
    stat = stat[stat["assembly_name"].isin(use)].sort_values("lin.num", kind="stable")
    x = (o2o[stat["assembly_name"].tolist()].to_numpy() == "one2zero").astype(np.float32)
    lin = stat["lin.num"].to_numpy()
    w = (stat["div.time"] / stat["div.time"].max()).to_numpy(np.float32)
    ns = x.shape[1]
    print(f"GISMO : {x.shape[0]} gènes, {ns} espèces, {len(set(lin))} lignées, "
          f"{iterations} tirages")
    # loss rate of each gene per lineage, then per species column
    p = np.empty_like(x)
    for L in np.unique(lin):
        cols = lin == L
        p[:, cols] = x[:, cols].mean(axis=1, keepdims=True)
    obs = (x * w).mean(axis=1)
    rng = np.random.default_rng(seed)
    sims = np.empty((x.shape[0], iterations), np.float32)
    for i in range(iterations):
        sims[:, i] = (rng.random(p.shape, dtype=np.float32) < p) @ w / ns
    lo, med, up = np.quantile(sims, [0.025, 0.5, 0.975], axis=1)
    return pd.DataFrame({"gene_id": [g.split(".")[0] for g in o2o.index], "gismo_obs": obs,
                         "gismo_lower": lo, "gismo_median": med, "gismo": up,
                         "gismo_decile": ntile(up)})


def gismo_mis(src: Source, mane: set) -> pd.DataFrame:
    t = pd.read_csv(src.open(src.find(r"gismo-mis_consensus_genenames.*\.txt$")), sep="\t")
    for c in ("genes", "mis.cal", "dec", "genename"):
        if c not in t.columns:
            sys.exit(f"table GISMO-mis : colonne {c} absente ({list(t.columns)})")
    t["tx"] = t["genes"].astype(str).str.split(".").str[0]
    t["_mane"] = t["tx"].isin(mane)
    t = t.sort_values(["genename", "_mane"], ascending=[True, False], kind="stable")
    t = t.drop_duplicates("genename")
    print(f"GISMO-mis : {len(t)} gènes ({int(t['_mane'].sum())} sur leur MANE Select)")
    return pd.DataFrame({"gene": t["genename"].astype(str), "gismo_mis": t["mis.cal"],
                         "gismo_mis_decile": t["dec"].astype("Int64"),
                         "gismo_mis_transcript": t["tx"]})


def gtf_genes(gtf: Path):
    """ENSG -> symbol, and the MANE Select transcripts"""
    sym, mane = {}, set()
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if f[2] == "gene":
                g = re.search(r'gene_id "(ENSG\d+)', f[8])
                n = re.search(r'gene_name "([^"]+)"', f[8])
                if g and n:
                    sym[g.group(1)] = n.group(1)
            elif f[2] == "transcript" and "MANE_Select" in f[8]:
                m = re.search(r'transcript_id "(ENST\d+)', f[8])
                if m:
                    mane.add(m.group(1))
    return sym, mane


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--zip", type=Path, required=True,
                   help="GISMO-GISMO.zip of the Zenodo record, or its extracted folder")
    p.add_argument("--gtf", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--iterations", type=int, default=1000,
                   help="null draws (10000 in the original script; 1000 gives the same deciles "
                        "up to genes at a boundary)")
    p.add_argument("--seed", type=int, default=1337)
    a = p.parse_args()
    src = Source(a.zip)
    sym, mane = gtf_genes(a.gtf)
    g = gismo(src, a.iterations, a.seed)
    g["gene"] = g["gene_id"].map(sym)
    m = gismo_mis(src, mane)
    out = g.merge(m, on="gene", how="outer")
    out = out[["gene_id", "gene", "gismo_obs", "gismo_lower", "gismo_median", "gismo",
               "gismo_decile", "gismo_mis", "gismo_mis_decile", "gismo_mis_transcript"]]
    out["gismo_decile"] = out["gismo_decile"].astype("Int64")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False, float_format="%.5g")
    print(f"{a.out} : {len(out)} lignes ; GISMO {out['gismo'].notna().sum()}, "
          f"GISMO-mis {out['gismo_mis'].notna().sum()}, sans symbole "
          f"{out['gene'].isna().sum()}")


if __name__ == "__main__":
    main()
