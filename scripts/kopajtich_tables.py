#!/usr/bin/env python3
"""
kopajtich_tables.py - the patient fibroblasts of Kopajtich et al. (RNA-seq
and TMT proteomics of the same individuals) as the inputs of
`l1_score.py validate`

Purpose
    Kopajtich et al. measured RNA-seq and proteomics in fibroblasts of about
    145 individuals with suspected mitochondrial disease (Zenodo 4501904,
    raw_data.zip: raw_counts.tsv, proteomics_not_normalized.tsv,
    proteomics_annotation.tsv). This script finds the gene column and the
    sample columns of the two matrices and the columns of the annotation
    that pair them (RNA id, proteome id, TMT batch), and writes the three
    tables that `l1_score.py validate --patients-*` reads. It also writes a
    description of every table (first lines and columns), so that a layout
    it does not recognise can be read from the run archive.

Inputs
    --dir      directory holding the extracted tables (default
               data/kopajtich/raw_data); --counts, --protein, --annotation
               to give them otherwise
    --gtf      GENCODE GTF: recognises the gene columns

Outputs (in --out, default data/kopajtich)
    patients_rna.tsv      gene, then one column per RNA sample: read counts
    patients_protein.tsv  gene, then one column per proteome: intensities
    patients_samples.tsv  rna_sample, protein_sample, batch
    and in --report (default results/kopajtich): layout.txt, conversion.txt

Usage
    python3 scripts/kopajtich_tables.py --dir data/kopajtich/raw_data \\
        --gtf data/gencode.gtf.gz

Method
    The gene column of a matrix is, of its non-numeric columns, the one
    with most values that are Ensembl gene ids or GTF symbols (protein
    groups 'A;B' are read as their first gene); the sample columns are
    those whose filled cells ('NA' counts as empty) are numeric. In the
    annotation, the proteome column is the one holding most proteome
    names, the RNA column the one holding most RNA names (both compared
    with case, '.', '-' and '_' ignored), and the batch the column named
    like batch, plex or run with between 2 and n/2 values. Without an RNA
    column, the samples are paired by name.

References
    Kopajtich R et al. 2021. Integration of proteomics with genomics and
    transcriptomics increases the diagnostic rate of Mendelian disorders.
    medRxiv 2021.03.09.21253187; Zenodo 10.5281/zenodo.4501904.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from jiang2020_tables import ENSG, gtf_names, is_num, log, missing  # noqa: E402
from re2g_panel import opener  # noqa: E402


def rows_of(path: Path):
    with opener(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.strip():
                f = line.rstrip("\n").split("\t")
                yield head, f + [""] * (len(head) - len(f))


def describe(path: Path, fh, lines=5, cols=12):
    with opener(path) as f:
        first = [next(f, "") for _ in range(lines)]
        n = 4 + sum(1 for _ in f)
    width = len(first[0].split("\t"))
    fh.write(f"### {path.name}: {width} columns, about {n} lines\n")
    for k, ln in enumerate(first):
        cells = [c[:30] for c in ln.rstrip("\n").split("\t")[:cols]]
        fh.write(f"  line {k + 1}: " + " | ".join(cells) + "\n")


def norm(s: str) -> str:
    return re.sub(r"[-._\s]+", "", str(s).upper())


def gene_of(v: str) -> str:
    v = str(v).strip().split(";")[0].strip()
    return v.split(".")[0] if ENSG.match(v) else v


def read_matrix(path: Path, names: set, rep, label: str):
    data = list(rows_of(path))
    if not data:
        sys.exit(f"{path}: empty")
    head = data[0][0]
    body = [f for _h, f in data]
    probe = body[:300]
    numeric, best, best_n = [], None, 0
    for j in range(len(head)):
        filled = [r[j] for r in probe if not missing(r[j])]
        num = len(filled) >= 0.1 * len(probe) and \
            sum(is_num(v) for v in filled) >= 0.9 * len(filled)
        numeric.append(num)
        if not num:
            n = sum(bool(ENSG.match(gene_of(v))) or gene_of(v) in names
                    for v in filled)
            if n > best_n:
                best, best_n = j, n
    if best is None:
        sys.exit(f"{label}: no column of Ensembl ids or gene symbols; see layout.txt")
    samples = [j for j in range(len(head)) if numeric[j]]
    log(rep, f"{label}: {path.name}, gene column {head[best]!r} "
             f"({best_n}/{len(probe)} recognised), {len(samples)} sample "
             f"columns, {len(body)} rows; first samples "
             f"{[head[j] for j in samples[:4]]}")
    return head, best, samples, body


def write_matrix(path: Path, head, gcol, samples, body):
    seen = set()
    with open(path, "w") as fh:
        fh.write("gene\t" + "\t".join(head[j] for j in samples) + "\n")
        for r in body:
            g = gene_of(r[gcol])
            if not g or g in seen:
                continue
            seen.add(g)
            fh.write(g + "\t" + "\t".join(
                "" if missing(r[j]) else r[j] for j in samples) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[1].strip(),
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", type=Path, default=Path("data/kopajtich/raw_data"))
    p.add_argument("--counts", type=Path)
    p.add_argument("--protein", type=Path)
    p.add_argument("--annotation", type=Path)
    p.add_argument("--gtf", type=Path, required=True)
    p.add_argument("--out", type=Path, default=Path("data/kopajtich"))
    p.add_argument("--report", type=Path, default=Path("results/kopajtich"))
    a = p.parse_args()
    a.report.mkdir(parents=True, exist_ok=True)
    a.out.mkdir(parents=True, exist_ok=True)
    counts = a.counts or a.dir / "raw_counts.tsv"
    prot = a.protein or a.dir / "proteomics_not_normalized.tsv"
    ann = a.annotation or a.dir / "proteomics_annotation.tsv"
    for f in (counts, prot):
        if not f.exists():
            sys.exit(f"{f} not found")
    with open(a.report / "layout.txt", "w") as fh:
        for f in sorted(a.dir.glob("*.tsv")) or [counts, prot, ann]:
            if f.exists():
                describe(f, fh)
    rep = open(a.report / "conversion.txt", "w")
    names, _coding = gtf_names(a.gtf)
    rh, rg, rs, rb = read_matrix(counts, names, rep, "RNA")
    ph, pgc, pss, pb = read_matrix(prot, names, rep, "protein")
    rnames = {norm(rh[j]): rh[j] for j in rs}
    pnames = {norm(ph[j]): ph[j] for j in pss}

    pairs, batch_of = [], {}
    if ann.exists():
        arows = list(rows_of(ann))
        ahead = arows[0][0] if arows else []
        abody = [f for _h, f in arows]
        hit = lambda j, ref: sum(norm(r[j]) in ref for r in abody)
        pc = max(range(len(ahead)), key=lambda j: hit(j, pnames), default=None)
        rc = max(range(len(ahead)), key=lambda j: hit(j, rnames), default=None)
        if pc is not None and hit(pc, pnames) == 0:
            pc = None
        if rc is not None and hit(rc, rnames) == 0:
            rc = None
        bc = None
        for j, h in enumerate(ahead):
            vals = {r[j] for r in abody if not missing(r[j])}
            if re.search(r"batch|plex|run", h, re.I) and 2 <= len(vals) <= len(abody) / 2:
                bc = j
                break
        log(rep, f"annotation {ann.name}: {len(abody)} rows; proteome column "
                 f"{ahead[pc] if pc is not None else None!r}, RNA column "
                 f"{ahead[rc] if rc is not None else None!r}, batch column "
                 f"{ahead[bc] if bc is not None else None!r}")
        if pc is not None and rc is not None:
            for r in abody:
                q, x = pnames.get(norm(r[pc])), rnames.get(norm(r[rc]))
                if q and x:
                    pairs.append((x, q))
                    batch_of[q] = r[bc] if bc is not None else ""
        elif pc is not None and bc is not None:
            for r in abody:
                q = pnames.get(norm(r[pc]))
                if q:
                    batch_of[q] = r[bc]
    if not pairs:
        pairs = [(rnames[k], q) for k, q in pnames.items() if k in rnames]
        log(rep, f"paired by sample name: {len(pairs)}")
    log(rep, f"{len(pairs)} individuals with RNA and protein "
             f"({len(set(batch_of.values()) - {''})} batches)")
    if len(pairs) < 10:
        log(rep, f"too few pairs: see layout.txt; RNA {[rh[j] for j in rs[:5]]}, "
                 f"protein {[ph[j] for j in pss[:5]]}")
        sys.exit(1)
    write_matrix(a.out / "patients_rna.tsv", rh, rg, rs, rb)
    write_matrix(a.out / "patients_protein.tsv", ph, pgc, pss, pb)
    with open(a.out / "patients_samples.tsv", "w") as fh:
        fh.write("rna_sample\tprotein_sample\tbatch\n")
        for x, q in pairs:
            fh.write(f"{x}\t{q}\t{batch_of.get(q, '')}\n")
    log(rep, f"written: {a.out}/patients_rna.tsv, patients_protein.tsv, "
             "patients_samples.tsv")
    rep.close()


if __name__ == "__main__":
    main()
