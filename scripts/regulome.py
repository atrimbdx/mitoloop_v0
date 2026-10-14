#!/usr/bin/env python3
"""
regulome.py - functional layer and mode of regulation of the genes of the
set, and the regulatory architecture of the most transcriptionally
controlled ones

Purpose
    Three steps towards a mitochondrial regulome:

      pathways  the layer on which each gene acts, taken as published: the
                MitoPathways hierarchy of MitoCarta3.0 (sections, then
                pathways). Genes of the set absent from MitoCarta (PanelApp
                genes, for instance) or without a MitoPathway get one
                inferred from Gene Ontology: their nearest MitoCarta genes in
                GO annotation, with a leave-one-out estimate of how often
                this inference is right. The section "Mitochondrial central
                dogma" (mtDNA maintenance, mtRNA metabolism, translation) is
                kept as its own layer, mtDNA.
      modes     supporting per-gene indicators, ranked among all
                protein-coding genes within bins of expression level. The
                built-in layer, variability, says how much the mRNA varies
                between tissues and between individuals (tissue sd, eQTL,
                Vg, protein-mRNA correlation of tissue medians); further
                modes (splicing, translation, turnover, ...) are added as
                indicator tables. The L1 score of transcriptional control
                itself is computed by l1_score.py from matched samples.
      compare   the promoters, rE2G elements, loops and methylation of the
                genes of the set along a score (l1.tsv, or a layer of
                modes.tsv), and the top against the bottom of it.

Inputs
    pathways
      --mitocarta     Human.MitoCarta3.0.xls (needs xlrd), or the sheet "A
                      Human MitoCarta3.0" exported as TSV. Columns used:
                      Symbol, Synonyms ('|'-separated), the column whose name
                      holds MitoPathways and the one holding
                      SubMitoLocalization when present.
      --mitopathways  optional sheet "C MitoPathways" as TSV (MitoPathways
                      Hierarchy, MitoPathway, Genes), read from the .xls
                      otherwise. Pathways given in sheet A as leaf names are
                      placed in the hierarchy through it.
      --genes         the gene set (MitoCarta3.0 and PanelApp symbols).
      --gaf, --obo    GO annotations (goa_human.gaf, gzipped or not) and the
                      ontology (go-basic.obo).
      --gtf           optional; symbols are then resolved as in
                      re2g_panel.py (alias table, --aliases).
    modes
      --gtf           GENCODE GTF: protein-coding genes, symbol and Ensembl id.
      --genes         the gene set (set = mito in the output), as symbols or
                      as a table with a gene column (pathways.tsv, whose
                      symbols are already matched to MitoCarta); --controls
                      optional (set = control); other genes are 'other'.
      --tpm           GTEx median TPM per tissue (GCT).
      --protein       optional protein abundance per tissue: a gene column
                      (symbol or Ensembl id) then one column per tissue, log
                      scale or not (rank correlation), for instance the
                      GTEx-matched proteome of Jiang et al. 2020.
                      --tissue-map (protein column, GTEx tissue) when the
                      column names are not GTEx tissue names.
      --eqtl          optional GTEx independent cis-eQTL (release tar or one
                      tissue table, as in gtex_eqtl.py) for --tissue.
      --vg            optional per-gene variance of genetic effects on
                      expression (ANEVA Vg, Mohammadi et al. 2019): gene
                      column then value column (--vg-column).
      --indicator     LAYER:NAME=PATH[:COLUMN][:-], repeatable: any other
                      per-gene table; '-' reverses the direction (high value
                      = less controlled at that layer).
    compare
      --scores        l1.tsv of l1_score.py (--score-column, default score);
                      or --modes modes.tsv with --layer. --expression
                      (modes.tsv) supplies the expression level of the
                      partial correlations when the score table lacks it.
      --pathways      optional pathways.tsv: the permutation is stratified on
                      the first section of each gene.
      --promoters     promoter_features.tsv (promoter_features.py)
      --methylation   methylation_by_group.tsv (promoter_features.py)
      --re2g          re2g_per_gene.tsv (re2g_panel.py)
      --loops         loops.tsv (v4c_loops.R); genes of the set missing from
                      it count as viewpoints without a loop.
      --assoc         assoc_loops.tsv (v4c_loops.R assoc)

Outputs (in --out)
    pathways   pathways.tsv: gene, source (mitocarta, go_inferred or
               unassigned), layers, sections, level2, pathways (full
               hierarchies, '|'-separated), sub-mitochondrial localisation,
               and for inferred genes support, expected_precision and the
               nearest MitoCarta genes; layer_matrix.tsv (gene x layer, 0/1);
               loo.tsv (leave-one-out on MitoCarta genes: gene, true
               sections, predicted section, support, correct);
               calibration.tsv (support bin, genes, share correct);
               genes_unmatched.txt
    modes      modes.tsv: gene, gene_id, set, mean_log_expression, then per
               indicator its value and its expression-matched percentile
               (pct_<name>), then per layer score_<layer> (mean percentile)
               and n_<layer> (indicators available); modes_summary.tsv
    compare    compare.tsv: one row per feature: source, feature, kind,
               genes, Spearman rho with the score and its p, partial rho
               given expression level and its p, p of the permutation
               within sections, high and low group values, P(high > low)
               or odds ratio and p, q (Benjamini-Hochberg on the partial p);
               groups.tsv: gene, score, group (high, middle, low), section

Method
    GO inference. Annotations are propagated to the ancestors (is_a,
    part_of); NOT qualifiers and --exclude-evidence codes are dropped. The
    information content of a term is -ln of the share of annotated genes
    carrying it. Two genes are compared by simGIC, the summed information
    content of the terms they share over that of the terms either carries
    (Pesquita et al. 2008). A gene takes the sections and pathways of its
    --k most similar MitoCarta genes, each weighted by its similarity; the
    support of a section is its share of that weight. The same inference
    run on every MitoCarta gene with its own annotation hidden gives the
    share of correct sections per support bin, reported as the expected
    precision of the inferred genes. MitoCarta itself used GO among its
    sources, so this estimate is an upper bound.

    Variability indicators. sd of log2(TPM + 1) across GTEx
    tissues (tissue_sd), number of independent eQTL signals and largest
    |slope| in --tissue (eqtl_signals, eqtl_effect; genes expressed above
    --min-tpm without a signal count 0), log10 Vg (vg), and the Spearman
    correlation between protein and mRNA across tissues (mrna_protein_rho;
    a protein that follows its mRNA is set before translation). Each
    indicator is turned into a percentile among the genes of the same
    expression bin (--bins quantile bins of mean log2(TPM + 1)), since
    eQTL power and protein detection grow with expression. The score of a
    layer is the mean of its percentiles over at least --min-indicators of
    them. Neither these indicators nor L1 use regulatory-element data, so
    that `compare` is not circular.

    Comparison. Spearman correlation of each feature with the score over
    the genes of --set; partial correlation on the residuals of the ranks
    of both after a linear fit of the rank of expression level; a
    permutation of the score within MitoPathways sections (--perms), which
    asks whether the relation holds beyond differences between sections;
    and the top against the bottom --fraction of genes (Mann-Whitney or
    Fisher test).

Usage
    python3 regulome.py pathways --mitocarta data/Human.MitoCarta3.0.xls \\
        --genes genes.txt --gaf data/goa_human.gaf.gz \\
        --obo data/go-basic.obo --gtf data/gencode.gtf.gz \\
        --out results/regulome

    python3 regulome.py modes --gtf data/gencode.gtf.gz --genes genes.txt \\
        --controls controls.txt --tpm data/gtex_median_tpm.gct.gz \\
        --protein data/jiang2020_protein.tsv \\
        --tissue-map config/protein_tissues.tsv \\
        --eqtl data/GTEx_Analysis_v8_eQTL_independent.tar \\
        --tissue Muscle_Skeletal --vg data/aneva_vg.tsv \\
        --out results/regulome

    python3 regulome.py compare --scores results/l1/l1.tsv \\
        --expression results/regulome/modes.tsv \\
        --pathways results/regulome/pathways.tsv \\
        --promoters results/promoters/promoter_features.tsv \\
        --methylation results/promoters/methylation_by_group.tsv \\
        --re2g results/re2g/re2g_per_gene.tsv \\
        --loops results/loops/loops.tsv \\
        --assoc results/loops/assoc_loops.tsv --out results/regulome

References
    Rath S et al. 2021. MitoCarta3.0: an updated mitochondrial proteome now
        with sub-organelle localization and pathway annotations. Nucleic
        Acids Res 49:D1541-D1547.
    Gene Ontology Consortium. 2023. The Gene Ontology knowledgebase in 2023.
        Genetics 224:iyad031.
    Pesquita C et al. 2008. Metrics for GO based protein semantic
        similarity: a systematic evaluation. BMC Bioinformatics 9:S4.
    GTEx Consortium. 2020. Science 369:1318-1330.
    Jiang L et al. 2020. A quantitative proteome map of the human body. Cell
        183:269-283.
    Mohammadi P et al. 2019. Genetic regulatory variation in populations
        informs transcriptome analysis in rare disease. Science 366:351-356.
    Martin AR et al. 2019. PanelApp. Nat Genet 51:1560-1565.
    Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate. J R
        Stat Soc B 57:289-300.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from re2g_panel import (ALIASES, NOT_GENE_PREFIXES, NOT_GENES,  # noqa: E402
                        build_aliases, gtf_gene_names, opener,
                        resolve_symbols)

DNA_SECTION = "Mitochondrial central dogma"
DNA_LAYER = "mtDNA"


def fmt(x):
    if isinstance(x, (float, np.floating)):
        return "NA" if np.isnan(x) else f"{x:.4g}"
    return str(x)


def write_tsv(path: Path, rows: list, cols: list) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r.get(c, "")) for c in cols) + "\n")


def read_table(path) -> list:
    with opener(path) as fh:
        head = fh.readline().rstrip("\n").lstrip("#").split("\t")
        rows = []
        for line in fh:
            if line.strip():
                f = line.rstrip("\n").split("\t")
                rows.append(dict(zip(head, f + [""] * (len(head) - len(f)))))
    return rows


def read_gene_list(path) -> list:
    """symbols separated by white space, or the gene column of a table
    with a header (pathways.tsv, for instance)"""
    text = Path(path).read_text()
    first = text.split("\n", 1)[0].split("\t")
    if "gene" in first and len(first) > 1:
        return [r["gene"] for r in read_table(path) if r.get("gene")]
    return [g for g in text.split() if g]


def to_float(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return np.nan
    return v


# ===========================================================================
# pathways

def read_mitocarta(path: Path, c_path: Path | None):
    """(gene -> record, leaf name -> hierarchy list) from MitoCarta3.0."""
    import pandas as pd
    a_rows, c_rows = None, None
    if str(path).endswith((".xls", ".xlsx")):
        try:
            sheets = pd.read_excel(path, sheet_name=None)
        except ImportError:
            sys.exit("reading the .xls needs xlrd; or export the sheets 'A "
                     "Human MitoCarta3.0' and 'C MitoPathways' as TSV")
        for name, d in sheets.items():
            if name.startswith("A") and a_rows is None:
                a_rows = d.fillna("").astype(str).to_dict("records")
            elif name.startswith("C") and c_rows is None:
                c_rows = d.fillna("").astype(str).to_dict("records")
    else:
        a_rows = read_table(path)
    if c_path:
        c_rows = read_table(c_path)
    if not a_rows:
        sys.exit(f"{path}: no gene table")

    hier = {}                  # leaf -> list of levels
    c_members = defaultdict(set)
    for r in c_rows or []:
        keys = {k.lower().strip(): k for k in r}
        hk = next((keys[k] for k in keys if "hierarchy" in k), None)
        lk = next((keys[k] for k in keys if k in ("mitopathway",
                                                   "mitopathways")), None)
        gk = next((keys[k] for k in keys if k == "genes"), None)
        if hk is None:
            continue
        levels = [x.strip() for x in r[hk].split(">") if x.strip()]
        leaf = r[lk].strip() if lk else levels[-1]
        if not levels:
            continue
        if levels[-1] != leaf:
            levels.append(leaf)
        hier[leaf] = levels
        if gk:
            for g in re.split(r"[,;|]\s*", r[gk]):
                if g.strip():
                    c_members[g.strip()].add(leaf)

    cols = list(a_rows[0])
    low = {c.lower(): c for c in cols}
    sym = low.get("symbol") or cols[0]
    syn = next((c for c in cols if c.lower().startswith("synonym")), None)
    pw = next((c for c in cols if "mitopathways" in c.lower()), None)
    loc = next((c for c in cols if "submitolocalization" in c.lower()), None)
    genes = {}
    for r in a_rows:
        g = str(r[sym]).strip()
        if not g:
            continue
        paths = []
        raw = str(r.get(pw, "")) if pw else ""
        for item in raw.split("|"):
            item = item.strip()
            if item in ("", "0", "nan"):
                continue
            if ">" in item:
                paths.append([x.strip() for x in item.split(">")
                              if x.strip()])
            else:
                paths.append(hier.get(item, [item]))
        for leaf in c_members.get(g, ()):
            p = hier[leaf]
            if p not in paths:
                paths.append(p)
        synonyms = [s.strip() for s in str(r.get(syn, "")).split("|")
                    if s.strip() and s.strip() not in ("-", "0")] if syn \
            else []
        genes[g] = {"paths": paths, "synonyms": synonyms,
                    "loc": str(r.get(loc, "")) if loc else ""}
    return genes, hier


def read_obo(path: Path):
    """term -> (namespace, parents); obsolete terms left out."""
    terms = {}
    cur = None
    with opener(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = {"id": None, "ns": "", "parents": [], "obs": False}
                continue
            if line.startswith("[") and line.endswith("]"):
                if cur and cur["id"] and not cur["obs"]:
                    terms[cur["id"]] = (cur["ns"], cur["parents"])
                cur = None
                continue
            if cur is None:
                continue
            if not line:
                if cur["id"] and not cur["obs"]:
                    terms[cur["id"]] = (cur["ns"], cur["parents"])
                cur = None
            elif line.startswith("id: "):
                cur["id"] = line[4:].strip()
            elif line.startswith("namespace: "):
                cur["ns"] = line[11:].strip()
            elif line.startswith("is_a: "):
                cur["parents"].append(line[6:].split("!")[0].strip())
            elif line.startswith("relationship: part_of "):
                cur["parents"].append(line[22:].split("!")[0].strip())
            elif line.startswith("is_obsolete: true"):
                cur["obs"] = True
    if cur and cur["id"] and not cur["obs"]:
        terms[cur["id"]] = (cur["ns"], cur["parents"])
    return terms


def ancestors_fn(terms):
    memo = {}

    def anc(t):
        if t in memo:
            return memo[t]
        memo[t] = {t}                       # guards against cycles
        out = {t}
        for p in terms.get(t, ("", []))[1]:
            if p in terms:
                out |= anc(p)
        memo[t] = out
        return out
    return anc


def read_gaf(path: Path, terms, exclude: set, namespaces: set):
    """symbol -> propagated term set; synonym -> symbol."""
    anc = ancestors_fn(terms)
    direct = defaultdict(set)
    synonyms = {}
    ns_code = {"P": "biological_process", "F": "molecular_function",
               "C": "cellular_component"}
    with opener(path) as fh:
        for line in fh:
            if line[0] == "!":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            if "NOT" in f[3] or f[6] in exclude:
                continue
            if namespaces and ns_code.get(f[8], f[8]) not in namespaces:
                continue
            if f[4] not in terms:
                continue
            direct[f[2]].add(f[4])
            for s in f[10].split("|"):
                if s and s != f[2]:
                    synonyms.setdefault(s, f[2])
    ann = {g: set().union(*(anc(t) for t in ts)) for g, ts in direct.items()}
    return ann, synonyms


def information_content(ann: dict):
    """term -> column index, and -ln(share of annotated genes) per term"""
    n = len(ann)
    freq = defaultdict(int)
    for ts in ann.values():
        for t in ts:
            freq[t] += 1
    order = sorted(freq)
    return ({t: i for i, t in enumerate(order)},
            np.array([-math.log(freq[t] / n) for t in order]))


def simgic(Xa, Xb, ic):
    """simGIC between every row of Xa and every row of Xb."""
    from scipy import sparse
    W = sparse.diags(ic)
    shared = (Xa @ W @ Xb.T).toarray()
    sa = np.asarray(Xa @ ic).ravel()
    sb = np.asarray(Xb @ ic).ravel()
    union = sa[:, None] + sb[None, :] - shared
    with np.errstate(invalid="ignore", divide="ignore"):
        s = np.where(union > 0, shared / union, 0.0)
    return s


def knn_vote(sim_row, ref_paths, k, exclude=None):
    """section, level-2 and pathway weights from the k nearest genes."""
    order = np.argsort(-sim_row)
    nb = []
    for j in order:
        if j == exclude:
            continue
        if sim_row[j] <= 0:
            break
        nb.append(j)
        if len(nb) == k:
            break
    tot = sum(sim_row[j] for j in nb)
    votes = [defaultdict(float) for _ in range(3)]
    if tot <= 0:
        return nb, votes
    for j in nb:
        w = sim_row[j] / tot
        seen = [set(), set(), set()]
        for p in ref_paths[j]:
            keys = (p[0], " > ".join(p[:2]), " > ".join(p))
            for lvl, key in enumerate(keys):
                if key not in seen[lvl]:
                    votes[lvl][key] += w
                    seen[lvl].add(key)
    return nb, votes


def layer_of(section: str) -> str:
    return DNA_LAYER if section == DNA_SECTION else section


def cmd_pathways(a) -> None:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    mc, _hier = read_mitocarta(a.mitocarta, a.mitopathways)
    print(f"MitoCarta: {len(mc)} genes, "
          f"{sum(bool(v['paths']) for v in mc.values())} with a MitoPathway")

    symbols = read_gene_list(a.genes)
    if a.gtf:
        names = gtf_gene_names(a.gtf)
        aliases = build_aliases(a.aliases, not a.no_builtin_aliases)
        mapping, unresolved, _ = resolve_symbols(symbols, names, aliases,
                                                 NOT_GENES, NOT_GENE_PREFIXES)
        mc_map, _u, _n = resolve_symbols(list(mc), names | set(mc), aliases,
                                         NOT_GENES, NOT_GENE_PREFIXES)
        mc = {mc_map.get(g, g): v for g, v in mc.items()}
    else:
        mapping, unresolved = {}, []
        for g in symbols:
            if any(t in g for t in NOT_GENES) or g.startswith(
                    NOT_GENE_PREFIXES):
                continue
            mapping[g] = g
    # symbols of the set that MitoCarta knows under another name
    syn_of = {}
    for g, v in mc.items():
        for s in v["synonyms"]:
            syn_of.setdefault(s, g)
    for old, new in ALIASES.items():
        for c in new:
            if c in mc:
                syn_of.setdefault(old, c)
    genes = []
    for g in symbols:
        if g not in mapping:
            continue
        x = mapping[g]
        if x not in mc and g in mc:
            x = g
        elif x not in mc and g in syn_of:
            x = syn_of[g]
        if x not in genes:
            genes.append(x)

    terms = read_obo(a.obo)
    namespaces = set(a.namespaces.split(",")) if a.namespaces else set()
    ann, gaf_syn = read_gaf(a.gaf, terms, set(a.exclude_evidence.split(",")),
                            namespaces)
    print(f"GO: {len(terms)} terms, {len(ann)} annotated genes")

    def go_key(g):
        if g in ann:
            return g
        for s in [g] + mc.get(g, {}).get("synonyms", []):
            if s in gaf_syn:
                return gaf_syn[s]
        return None

    ref = [g for g, v in mc.items() if v["paths"] and go_key(g)]
    ref_paths = [mc[g]["paths"] for g in ref]
    queries = [g for g in genes if not (g in mc and mc[g]["paths"])]
    qkeys = [go_key(g) for g in queries]
    # information content over every annotated gene, not only these
    from scipy import sparse
    tid, ic = information_content(ann)

    def mat(keys):
        r, c = [], []
        for i, kk in enumerate(keys):
            for t in ann.get(kk, ()):
                r.append(i)
                c.append(tid[t])
        return sparse.csr_matrix((np.ones(len(r)), (r, c)),
                                 shape=(len(keys), len(tid)))
    Xr = mat([go_key(g) for g in ref])
    print(f"{len(ref)} MitoCarta genes with a MitoPathway and GO terms as "
          f"reference; {len(queries)} genes of the set to infer")

    # leave-one-out on the reference
    S = simgic(Xr, Xr, ic)
    loo = []
    for i, g in enumerate(ref):
        _nb, votes = knn_vote(S[i], ref_paths, a.k, exclude=i)
        if not votes[0]:
            continue
        sec, sup = max(votes[0].items(), key=lambda x: x[1])
        true = {p[0] for p in ref_paths[i]}
        loo.append({"gene": g, "true_sections": "|".join(sorted(true)),
                    "predicted_section": sec, "support": sup,
                    "correct": int(sec in true)})
    write_tsv(out / "loo.tsv", loo, ["gene", "true_sections",
                                     "predicted_section", "support",
                                     "correct"])
    edges = np.array([0, 0.25, 0.5, 0.75, 1.0001])
    calib = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sub = [r["correct"] for r in loo if lo <= r["support"] < hi]
        calib.append({"support_min": lo, "support_max": min(hi, 1.0),
                      "genes": len(sub),
                      "share_correct": float(np.mean(sub)) if sub
                      else np.nan})
    write_tsv(out / "calibration.tsv", calib,
              ["support_min", "support_max", "genes", "share_correct"])
    acc = np.mean([r["correct"] for r in loo]) if loo else np.nan
    print(f"leave-one-out: section right for {100 * acc:.1f}% of "
          f"{len(loo)} MitoCarta genes")

    def precision(sup):
        for c in calib:
            if c["support_min"] <= sup < c["support_max"] + 1e-9:
                return c["share_correct"]
        return np.nan

    # inference
    inferred = {}
    if queries:
        ok = [i for i, kk in enumerate(qkeys) if kk]
        if ok:
            Xq = mat([qkeys[i] for i in ok])
            Sq = simgic(Xq, Xr, ic)
            for row, i in enumerate(ok):
                nb, votes = knn_vote(Sq[row], ref_paths, a.k)
                if not votes[0]:
                    continue
                inferred[queries[i]] = (nb, votes, Sq[row])

    rows = []
    for g in genes:
        r = {"gene": g, "in_mitocarta": int(g in mc),
             "sub_localization": mc.get(g, {}).get("loc", "")}
        if g in mc and mc[g]["paths"]:
            paths = mc[g]["paths"]
            r["source"] = "mitocarta"
            secs = list(dict.fromkeys(p[0] for p in paths))
            lv2 = list(dict.fromkeys(" > ".join(p[:2]) for p in paths))
            r["pathways"] = "|".join(" > ".join(p) for p in paths)
        elif g in inferred:
            nb, votes, srow = inferred[g]
            r["source"] = "go_inferred"
            best = max(votes[0].values())
            secs = [s for s, v in sorted(votes[0].items(),
                                         key=lambda x: -x[1])
                    if v >= a.min_support or v == best]
            lv2 = [s for s, v in sorted(votes[1].items(),
                                        key=lambda x: -x[1])
                   if v >= a.min_support and s.split(" > ")[0] in secs]
            pws = [s for s, v in sorted(votes[2].items(),
                                        key=lambda x: -x[1])
                   if v >= a.min_support and s.split(" > ")[0] in secs]
            r["pathways"] = "|".join(pws)
            r["support"] = best
            r["expected_precision"] = precision(best)
            r["neighbours"] = ",".join(f"{ref[j]}:{srow[j]:.2f}"
                                       for j in nb[:5])
        else:
            r["source"] = "unassigned"
            secs, lv2 = [], []
            r["pathways"] = ""
        r["sections"] = "|".join(secs)
        r["layers"] = "|".join(dict.fromkeys(layer_of(s) for s in secs))
        r["level2"] = "|".join(lv2)
        rows.append(r)
    write_tsv(out / "pathways.tsv", rows,
              ["gene", "source", "layers", "sections", "level2", "pathways",
               "in_mitocarta", "sub_localization", "support",
               "expected_precision", "neighbours"])
    layers = sorted({x for r in rows for x in r["layers"].split("|") if x})
    write_tsv(out / "layer_matrix.tsv",
              [{"gene": r["gene"],
                **{L: int(L in r["layers"].split("|")) for L in layers}}
               for r in rows], ["gene"] + layers)
    (out / "genes_unmatched.txt").write_text(
        "".join(f"{g}\n" for g in unresolved))
    tally = defaultdict(int)
    for r in rows:
        tally[r["source"]] += 1
    print("genes of the set: " + ", ".join(f"{k} {v}" for k, v in
                                           sorted(tally.items())))
    print("genes per layer: " + ", ".join(
        f"{L} {sum(L in r['layers'].split('|') for r in rows)}"
        for L in layers))
    print(f"Written to {out}")


# ===========================================================================
# modes

def load_gtf_genes(gtf: Path, gene_type: str) -> dict:
    """Ensembl id -> symbol for one gene type."""
    out = {}
    tag = f'gene_type "{gene_type}"'
    with opener(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene" or tag not in f[8] \
                    or "_PAR_Y" in f[8]:
                continue
            a = f[8]
            i = a.find('gene_id "')
            j = a.find('gene_name "')
            out[a[i + 9:a.find('"', i + 9)].split(".")[0]] = \
                a[j + 11:a.find('"', j + 11)]
    return out


def read_gct(path: Path):
    """(ids without version, symbols, tissues, matrix)"""
    with opener(path) as fh:
        fh.readline()
        fh.readline()
        head = fh.readline().rstrip("\n").split("\t")
        ids, syms, vals = [], [], []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            ids.append(f[0].split(".")[0])
            syms.append(f[1])
            vals.append([to_float(x) for x in f[2:]])
    return ids, syms, head[2:], np.array(vals, float)


def gene_value_table(path, column=None):
    """first column -> value of `column` (header name, 1-based index, or
    the last column)."""
    rows = read_table(path)
    if not rows:
        return {}
    cols = list(rows[0])
    if column is None:
        col = cols[-1]
    elif column.isdigit():
        col = cols[int(column) - 1]
    elif column in cols:
        col = column
    else:
        sys.exit(f"{path}: no column {column!r}")
    return {r[cols[0]].split(".")[0]: to_float(r[col]) for r in rows}


def binned_percentile(values, covariate, bins):
    """percentile of each value among the genes of its covariate bin."""
    from scipy.stats import rankdata
    v = np.asarray(values, float)
    out = np.full(len(v), np.nan)
    ok = ~np.isnan(v)
    if covariate is None or bins <= 1:
        idx = np.flatnonzero(ok)
        if len(idx):
            out[idx] = (rankdata(v[idx]) - 0.5) / len(idx)
        return out
    c = np.asarray(covariate, float)
    ok &= ~np.isnan(c)
    idx = np.flatnonzero(ok)
    if not len(idx):
        return out
    edges = np.unique(np.quantile(c[idx], np.linspace(0, 1, bins + 1)))
    b = np.clip(np.searchsorted(edges, c[idx], "right") - 1, 0,
                len(edges) - 2)
    for k in np.unique(b):
        sel = idx[b == k]
        out[sel] = (rankdata(v[sel]) - 0.5) / len(sel)
    return out


def cmd_modes(a) -> None:
    from scipy.stats import spearmanr
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ens2sym = load_gtf_genes(a.gtf, a.gene_type)
    sym2ens = {}
    for e, s in ens2sym.items():
        sym2ens.setdefault(s, e)
    names = set(sym2ens)
    aliases = build_aliases(a.aliases, not a.no_builtin_aliases)
    gene_set = {}
    for path, label in ((a.genes, "mito"), (a.controls, "control")):
        if not path:
            continue
        m, unresolved, _n = resolve_symbols(read_gene_list(path), names,
                                            aliases, NOT_GENES,
                                            NOT_GENE_PREFIXES)
        for g in m.values():
            gene_set.setdefault(g, label)
        if unresolved:
            print(f"{label}: {len(unresolved)} symbols not in the GTF: "
                  + " ".join(unresolved[:20]))
    ids = sorted(ens2sym)
    n = len(ids)
    pos = {e: i for i, e in enumerate(ids)}

    def key_of(x):
        x = x.split(".")[0]
        if x in pos:
            return pos[x]
        e = sym2ens.get(x)
        return pos.get(e) if e else None

    ind = {}                      # name -> (layer, values, reversed)
    mean_expr = None
    if a.tpm:
        gid, _gs, tissues, M = read_gct(a.tpm)
        keep_t = [i for i, t in enumerate(tissues)
                  if not (a.exclude_tissues and
                          re.search(a.exclude_tissues, t))]
        L = np.log2(M[:, keep_t] + 1)
        mean_expr = np.full(n, np.nan)
        sd = np.full(n, np.nan)
        rows_of = {}
        for r, e in enumerate(gid):
            if e in pos and e not in rows_of:
                rows_of[e] = r
        for e, r in rows_of.items():
            mean_expr[pos[e]] = L[r].mean()
            sd[pos[e]] = L[r].std(ddof=1)
        expressed = np.zeros(n, bool)
        for e, r in rows_of.items():
            expressed[pos[e]] = np.nanmax(M[r, keep_t]) >= a.min_tpm
        sd[~expressed] = np.nan
        ind["tissue_sd"] = ("variability", sd, False)
        print(f"GTEx: {len(rows_of)} genes, {len(keep_t)} tissues")
        tcol = {t: i for i, t in enumerate(tissues)}

        if a.protein:
            prot = read_table(a.protein)
            pcols = list(prot[0])[1:]
            tmap = {}
            if a.tissue_map:
                for r in read_table(a.tissue_map):
                    v = list(r.values())
                    tmap[v[0]] = v[1]
            pairs = [(c, tcol[tmap.get(c, c)]) for c in pcols
                     if tmap.get(c, c) in tcol]
            if len(pairs) < a.min_tissues:
                sys.exit(f"--protein: {len(pairs)} columns matched to GTEx "
                         "tissues; give --tissue-map")
            rho = np.full(n, np.nan)
            first = list(prot[0])[0]
            for r in prot:
                k = key_of(r[first])
                if k is None or ids[k] not in rows_of:
                    continue
                x = np.array([to_float(r[c]) for c, _ in pairs])
                y = M[rows_of[ids[k]], [t for _, t in pairs]]
                ok = ~np.isnan(x) & ~np.isnan(y)
                if ok.sum() >= a.min_tissues and np.ptp(x[ok]) > 0 \
                        and np.ptp(y[ok]) > 0:
                    rho[k] = spearmanr(x[ok], y[ok])[0]
            ind["mrna_protein_rho"] = ("variability", rho, False)
            print(f"protein: {np.sum(~np.isnan(rho))} genes over "
                  f"{len(pairs)} tissues")

    if a.eqtl:
        from gtex_eqtl import load_eqtl
        d = load_eqtl(Path(a.eqtl), a.tissue, "gene_id", "tss_distance",
                      "slope")
        nsig = np.zeros(n)
        eff = np.zeros(n)
        seen = np.zeros(n, bool)
        for e, s in zip(d["gene_id"], d["abs_slope"]):
            if e in pos:
                k = pos[e]
                nsig[k] += 1
                eff[k] = max(eff[k], s)
                seen[k] = True
        tested = seen.copy()
        if a.tpm and a.tpm_column:
            if a.tpm_column not in tcol:
                sys.exit(f"--tpm-column {a.tpm_column!r} not in the GCT")
            for e, r in rows_of.items():
                if M[r, tcol[a.tpm_column]] >= a.min_tpm:
                    tested[pos[e]] = True
        nsig[~tested] = np.nan
        eff[~tested] = np.nan
        ind["eqtl_signals"] = ("variability", nsig, False)
        ind["eqtl_effect"] = ("variability", eff, False)
        print(f"eQTL: {int(seen.sum())} genes with a signal, "
              f"{int(tested.sum())} tested")

    if a.vg:
        tab = gene_value_table(a.vg, a.vg_column)
        v = np.full(n, np.nan)
        for g, x in tab.items():
            k = key_of(g)
            if k is not None and x > 0:
                v[k] = math.log10(x)
        ind["vg"] = ("variability", v, False)

    for spec in a.indicator or []:
        m = re.match(r"^([^:=]+):([^=]+)=(.+)$", spec)
        if not m:
            sys.exit(f"--indicator {spec!r}: expected "
                     "LAYER:NAME=PATH[:COLUMN][:-]")
        layer, name, rest = m.groups()
        rev = rest.endswith(":-")
        rest = rest[:-2] if rev else rest
        col = None
        mm = re.match(r"^(.*\.(?:tsv|txt|gz))(?::(.+))?$", rest)
        if mm:
            rest, col = mm.group(1), mm.group(2)
        tab = gene_value_table(rest, col)
        v = np.full(n, np.nan)
        for g, x in tab.items():
            k = key_of(g)
            if k is not None:
                v[k] = x
        ind[name] = (layer, v, rev)
    if not ind:
        sys.exit("no indicator: give --tpm, --protein, --eqtl, --vg or "
                 "--indicator")

    pct = {}
    for name, (layer, v, rev) in ind.items():
        p = binned_percentile(-v if rev else v, mean_expr, a.bins)
        pct[name] = p
    layers = list(dict.fromkeys(x[0] for x in ind.values()))
    score, count = {}, {}
    for L in layers:
        P = np.column_stack([pct[k] for k, x in ind.items() if x[0] == L])
        c = (~np.isnan(P)).sum(1)
        with np.errstate(all="ignore"):
            s = np.nanmean(np.where(np.isnan(P), np.nan, P), 1)
        s[c < min(a.min_indicators, P.shape[1])] = np.nan
        score[L], count[L] = s, c

    rows = []
    for i, e in enumerate(ids):
        g = ens2sym[e]
        r = {"gene": g, "gene_id": e, "set": gene_set.get(g, "other"),
             "mean_log_expression": mean_expr[i] if mean_expr is not None
             else np.nan}
        for k, (_l, v, _r) in ind.items():
            r[k] = v[i]
            r[f"pct_{k}"] = pct[k][i]
        for L in layers:
            r[f"score_{L}"] = score[L][i]
            r[f"n_{L}"] = int(count[L][i])
        rows.append(r)
    cols = (["gene", "gene_id", "set", "mean_log_expression"]
            + [x for k in ind for x in (k, f"pct_{k}")]
            + [x for L in layers for x in (f"score_{L}", f"n_{L}")])
    write_tsv(out / "modes.tsv", rows, cols)

    from gtex_eqtl import prob_superiority
    summ = []
    sets = np.array([r["set"] for r in rows])
    print(f"\n{'layer':16s}{'set':9s}{'genes':>7s}{'median':>8s}"
          f"{'P(set>other)':>14s}")
    for L in layers:
        s = score[L]
        other = s[(sets == "other") & ~np.isnan(s)]
        for lab in ("mito", "control"):
            x = s[(sets == lab) & ~np.isnan(s)]
            if not len(x):
                continue
            ps = prob_superiority(x, other) if len(other) else np.nan
            summ.append({"layer": L, "set": lab, "genes": len(x),
                         "median_score": float(np.median(x)),
                         "P_set_gt_other": ps})
            print(f"{L:16s}{lab:9s}{len(x):7d}{np.median(x):8.3f}{ps:14.3f}")
    write_tsv(out / "modes_summary.tsv", summ,
              ["layer", "set", "genes", "median_score", "P_set_gt_other"])
    print(f"Written to {out}")


# ===========================================================================
# compare

def bh(p):
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = np.flatnonzero(~np.isnan(p))
    if not len(ok):
        return q
    o = ok[np.argsort(p[ok])]
    m = len(o)
    vals = p[o] * m / np.arange(1, m + 1)
    vals = np.minimum.accumulate(vals[::-1])[::-1]
    q[o] = np.minimum(vals, 1.0)
    return q


def rank_resid(y, x):
    """residuals of rank(y) after a linear fit on rank(x)"""
    from scipy.stats import rankdata
    ry, rx = rankdata(y), rankdata(x)
    A = np.column_stack([np.ones(len(rx)), rx])
    beta = np.linalg.lstsq(A, ry, rcond=None)[0]
    return ry - A @ beta


def cmd_compare(a) -> None:
    from scipy.stats import fisher_exact, mannwhitneyu, rankdata, spearmanr

    from gtex_eqtl import prob_superiority
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.scores:
        modes = read_table(a.scores)
        col = a.score_column
        if col not in modes[0]:
            sys.exit(f"{a.scores} has no column {col}")
    elif a.modes:
        modes = read_table(a.modes)
        col = f"score_{a.layer}"
        if col not in modes[0]:
            sys.exit(f"{a.modes} has no {col}; layers: " + ", ".join(
                c[6:] for c in modes[0] if c.startswith("score_")))
    else:
        sys.exit("give --scores (l1.tsv of l1_score.py) or --modes")
    expr_of = {}
    if a.expression:
        expr_of = {r["gene"]: r.get("mean_log_expression")
                   for r in read_table(a.expression)}
    genes, score, expr = [], [], []
    for r in modes:
        if a.set and r.get("set", a.set) != a.set:
            continue
        s = to_float(r[col])
        if np.isnan(s):
            continue
        genes.append(r["gene"])
        score.append(s)
        expr.append(to_float(r.get("mean_log_expression",
                                   expr_of.get(r["gene"]))))
    score, expr = np.array(score), np.array(expr)
    gi = {g: i for i, g in enumerate(genes)}
    n = len(genes)
    if n < 10:
        sys.exit(f"{n} genes of set {a.set!r} with a {a.layer} score")
    order = np.argsort(score)
    k = max(1, int(round(a.fraction * n)))
    group = np.array(["middle"] * n, dtype=object)
    group[order[:k]] = "low"
    group[order[-k:]] = "high"

    section = np.array(["none"] * n, dtype=object)
    if a.pathways:
        for r in read_table(a.pathways):
            if r["gene"] in gi and r.get("sections"):
                section[gi[r["gene"]]] = r["sections"].split("|")[0]

    feats = {}          # name -> (source, kind, values)

    def add(source, name, kind, values):
        feats[f"{source}.{name}"] = (source, kind, values)

    def per_gene(path, cols, source, binary=()):
        rows = read_table(path)
        if not rows:
            return
        avail = list(rows[0])
        use = [c for c in avail if (cols is None or c in cols
                                    or any(c.startswith(p) for p in cols
                                           if p.endswith("_")))]
        for c in use:
            if c in ("gene", "set", "chrom", "tss", "strand", "gene_id",
                     "prom_start", "prom_end"):
                continue
            v = np.full(n, np.nan)
            for r in rows:
                i = gi.get(r.get("gene", ""))
                if i is not None:
                    v[i] = to_float(r[c])
            if np.all(np.isnan(v)):
                continue
            is_bin = c in binary or c.startswith(binary) if binary else False
            if is_bin:
                v = np.where(np.isnan(v), np.nan, (v > 0).astype(float))
            add(source, c, "binary" if is_bin else "continuous", v)

    if a.promoters:
        per_gene(a.promoters, ["gc", "cpg_oe", "cpg_island", "bidirectional",
                               "n_promoter_loops", "motif_"], "promoter",
                 binary=("motif_", "cpg_island", "bidirectional"))
        hcp = np.full(n, np.nan)
        for r in read_table(a.promoters):
            i = gi.get(r["gene"])
            if i is not None and r.get("promoter_class") not in (None, "",
                                                                 "NA"):
                hcp[i] = float(r["promoter_class"] == "HCP")
        if not np.all(np.isnan(hcp)):
            add("promoter", "HCP", "binary", hcp)
    if a.methylation:
        rows = read_table(a.methylation)
        groups = [c for c in rows[0] if c not in ("gene", "set")]
        mean = np.full(n, np.nan)
        rng_ = np.full(n, np.nan)
        for r in rows:
            i = gi.get(r["gene"])
            if i is None:
                continue
            v = np.array([to_float(r[c]) for c in groups])
            v = v[~np.isnan(v)]
            if len(v):
                mean[i], rng_[i] = v.mean(), np.ptp(v)
        add("methylation", "mean_beta", "continuous", mean)
        add("methylation", "range_beta", "continuous", rng_)
    if a.re2g:
        per_gene(a.re2g, ["n_ccre_window", "links_per_biosample_median",
                          "n_union", "n_intersection", "max_prevalence",
                          "frac_never_linked"], "re2g")
    if a.loops:
        nl = np.zeros(n)
        terms = defaultdict(list)
        for r in read_table(a.loops):
            i = gi.get(r.get("gene", ""))
            if i is None:
                continue
            nl[i] += 1
            terms[i].append(to_float(r.get("n_terms")))
        add("loops", "n_loops", "continuous", nl)
        add("loops", "has_loop", "binary", (nl > 0).astype(float))
        mt = np.full(n, np.nan)
        for i, v in terms.items():
            v = [x for x in v if not np.isnan(x)]
            if v:
                mt[i] = float(np.mean(v))
        if not np.all(np.isnan(mt)):
            add("loops", "mean_terms_per_loop", "continuous", mt)
    if a.assoc:
        acc = defaultdict(list)
        for r in read_table(a.assoc):
            i = gi.get(r.get("gene", ""))
            if i is not None:
                acc[i].append(to_float(r.get("log2fc_per_doubling")))
        v = np.full(n, np.nan)
        for i, x in acc.items():
            x = [y for y in x if not np.isnan(y)]
            if x:
                v[i] = float(np.mean(x))
        add("loops", "mtcn_log2fc_per_doubling", "continuous", v)
    if not feats:
        sys.exit("no feature table given")

    rng = np.random.default_rng(a.seed)
    strata = [np.flatnonzero(section == s) for s in np.unique(section)]
    have_expr = np.sum(~np.isnan(expr)) >= 10
    res = []
    for name, (source, kind, v) in feats.items():
        ok = ~np.isnan(v)
        if ok.sum() < 10 or np.ptp(v[ok]) == 0:
            continue
        r = {"source": source, "feature": name.split(".", 1)[1],
             "kind": kind, "genes": int(ok.sum())}
        rho, p = spearmanr(score[ok], v[ok])
        r["rho"], r["p_rho"] = float(rho), float(p)
        if have_expr:
            ok2 = ok & ~np.isnan(expr)
            ry = rank_resid(v[ok2], expr[ok2])
            rx = rank_resid(score[ok2], expr[ok2])
            pr, pp = spearmanr(rx, ry)
            r["partial_rho"], r["p_partial"] = float(pr), float(pp)
        # permutation of the score within sections
        if a.perms > 0:
            rv = rankdata(v[ok])
            idx_ok = np.flatnonzero(ok)
            loc = {j: t for t, j in enumerate(idx_ok)}
            st = [np.array([loc[j] for j in s if j in loc]) for s in strata]
            st = [s for s in st if len(s)]
            rs = rankdata(score[ok])
            obs = np.corrcoef(rs, rv)[0, 1]
            hits = 0
            perm = rs.copy()
            for _ in range(a.perms):
                for s in st:
                    perm[s] = rng.permutation(rs[s])
                hits += abs(np.corrcoef(perm, rv)[0, 1]) >= abs(obs) - 1e-12
            r["p_within_sections"] = (1 + hits) / (a.perms + 1)
        hi = v[(group == "high") & ok]
        lo = v[(group == "low") & ok]
        if len(hi) and len(lo):
            if kind == "binary":
                k1, k2 = int(hi.sum()), int(lo.sum())
                odds, pg = fisher_exact([[k1, len(hi) - k1],
                                         [k2, len(lo) - k2]])
                r.update(high=k1 / len(hi), low=k2 / len(lo),
                         effect=float(odds), effect_kind="odds_ratio",
                         p_groups=float(pg))
            else:
                r.update(high=float(np.median(hi)), low=float(np.median(lo)),
                         effect=prob_superiority(hi, lo),
                         effect_kind="P(high>low)",
                         p_groups=float(mannwhitneyu(hi, lo).pvalue))
        res.append(r)
    key = "p_partial" if have_expr else "p_rho"
    q = bh([r.get(key, np.nan) for r in res])
    for r, x in zip(res, q):
        r["q"] = x
    res.sort(key=lambda r: (np.nan_to_num(r.get(key, 1), nan=1.0)))
    write_tsv(out / "compare.tsv", res,
              ["source", "feature", "kind", "genes", "rho", "p_rho",
               "partial_rho", "p_partial", "p_within_sections", "high",
               "low", "effect", "effect_kind", "p_groups", "q"])
    write_tsv(out / "groups.tsv",
              [{"gene": g, "score": score[i], "group": group[i],
                "section": section[i],
                "mean_log_expression": expr[i]}
               for i, g in enumerate(genes)],
              ["gene", "score", "group", "section", "mean_log_expression"])
    print(f"{n} genes of set {a.set!r} with a {a.layer} score; high and low "
          f"groups of {k} genes")
    print(f"\n{'feature':38s}{'rho':>7s}{'partial':>9s}{'p':>10s}"
          f"{'p|sect':>10s}{'high':>10s}{'low':>10s}{'q':>10s}")
    for r in res:
        print(f"{r['source'] + '.' + r['feature']:38s}{r['rho']:7.2f}"
              f"{r.get('partial_rho', np.nan):9.2f}"
              f"{r.get(key, np.nan):10.2g}"
              f"{r.get('p_within_sections', np.nan):10.2g}"
              f"{r.get('high', np.nan):10.4g}{r.get('low', np.nan):10.4g}"
              f"{r['q']:10.2g}")
    print(f"Written to {out}")


# ===========================================================================

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    w = sub.add_parser("pathways", help="MitoPathways, inferred from GO "
                                        "where missing")
    w.add_argument("--mitocarta", required=True, type=Path)
    w.add_argument("--mitopathways", type=Path)
    w.add_argument("--genes", required=True, type=Path)
    w.add_argument("--gaf", required=True, type=Path)
    w.add_argument("--obo", required=True, type=Path)
    w.add_argument("--gtf", type=Path)
    w.add_argument("--aliases", type=Path)
    w.add_argument("--no-builtin-aliases", action="store_true")
    w.add_argument("--k", type=int, default=10,
                   help="nearest MitoCarta genes voting (default %(default)s)")
    w.add_argument("--min-support", type=float, default=0.3,
                   help="share of the vote for a section or pathway to be "
                        "kept beside the best one (default %(default)s)")
    w.add_argument("--exclude-evidence", default="ND",
                   help="comma-separated GO evidence codes left out "
                        "(default %(default)s)")
    w.add_argument("--namespaces", default="",
                   help="comma-separated GO namespaces kept "
                        "(biological_process, molecular_function, "
                        "cellular_component); default all")
    w.add_argument("--out", required=True)

    m = sub.add_parser("modes", help="score per gene and mode of regulation")
    m.add_argument("--gtf", required=True, type=Path)
    m.add_argument("--gene-type", default="protein_coding")
    m.add_argument("--genes", required=True, type=Path)
    m.add_argument("--controls", type=Path)
    m.add_argument("--aliases", type=Path)
    m.add_argument("--no-builtin-aliases", action="store_true")
    m.add_argument("--tpm", type=Path)
    m.add_argument("--exclude-tissues", default="^Cells",
                   help="regular expression on GTEx tissues left out "
                        "(default %(default)s: cultured cells)")
    m.add_argument("--min-tpm", type=float, default=1.0,
                   help="median TPM a gene needs in at least one tissue "
                        "(default %(default)s)")
    m.add_argument("--protein", type=Path)
    m.add_argument("--tissue-map", type=Path)
    m.add_argument("--min-tissues", type=int, default=8,
                   help="tissues needed for a protein-mRNA correlation "
                        "(default %(default)s)")
    m.add_argument("--eqtl", type=Path)
    m.add_argument("--tissue", default="Muscle_Skeletal",
                   help="eQTL tissue (file name) (default %(default)s)")
    m.add_argument("--tpm-column", default="Muscle - Skeletal",
                   help="TPM column of the eQTL tissue, for the genes "
                        "tested without a signal (default %(default)s)")
    m.add_argument("--vg", type=Path)
    m.add_argument("--vg-column", default=None)
    m.add_argument("--indicator", action="append")
    m.add_argument("--bins", type=int, default=20,
                   help="expression bins of the percentiles "
                        "(default %(default)s)")
    m.add_argument("--min-indicators", type=int, default=2,
                   help="indicators a gene needs for a layer score "
                        "(default %(default)s)")
    m.add_argument("--out", required=True)

    c = sub.add_parser("compare", help="architecture of high- against "
                                       "low-score genes")
    c.add_argument("--scores", type=Path,
                   help="per-gene score table, such as l1.tsv")
    c.add_argument("--score-column", default="score")
    c.add_argument("--modes", type=Path,
                   help="modes.tsv, with --layer, instead of --scores")
    c.add_argument("--layer", default="variability")
    c.add_argument("--expression", type=Path,
                   help="table with gene and mean_log_expression (modes.tsv) "
                        "when --scores lacks it")
    c.add_argument("--set", default="mito",
                   help="genes compared (set column of modes.tsv); '' for "
                        "all (default %(default)s)")
    c.add_argument("--pathways", type=Path)
    c.add_argument("--promoters", type=Path)
    c.add_argument("--methylation", type=Path)
    c.add_argument("--re2g", type=Path)
    c.add_argument("--loops", type=Path)
    c.add_argument("--assoc", type=Path)
    c.add_argument("--fraction", type=float, default=1 / 3,
                   help="share of genes in the high and in the low group "
                        "(default 1/3)")
    c.add_argument("--perms", type=int, default=1000)
    c.add_argument("--seed", type=int, default=0)
    c.add_argument("--out", required=True)

    a = p.parse_args()
    {"pathways": cmd_pathways, "modes": cmd_modes,
     "compare": cmd_compare}[a.cmd](a)


if __name__ == "__main__":
    main()
