#!/usr/bin/env python3
"""
jiang2020_tables.py - the matched RNA and protein of the GTEx tissue
proteome (Jiang et al. 2020, Cell 183:269) as the inputs of l1_score.py

Purpose
    Jiang et al. measured the proteome of 201 GTEx samples (32 tissues, TMT)
    and give, in their supplementary tables, the per-sample protein
    abundances (Table S2, sheet D: cleaned relative abundances, log2) and
    the RNA of the same samples (Table S3, sheet A: cleaned and normalised
    TPM, log scale). This script finds these sheets, the gene column, the
    sample columns and the tissue of each sample, and writes the three
    tables that `l1_score.py score` reads. It also writes a description of
    every sheet (first rows and columns), so that a layout it does not
    recognise can be read from the run archive without the 150 MB of
    spreadsheets.

Inputs
    --dir      directory holding the supplementary files (names ending in
               mmc3.xlsx = Table S2, mmc4.xlsx = Table S3, mmc5.xlsx = Table
               S4; --s2, --s3 to give them otherwise)
    --gtf      GENCODE GTF: recognises the gene column (Ensembl ids or
               symbols)
    --tpm      optional GTEx median TPM (GCT): its tissue names, for the
               per-tissue protein table of `regulome.py modes`
    --attributes  optional GTEx v8 SampleAttributesDS.txt: tissue (SMTSD)
               and ischemic time (SMTSISCH) of samples named by GTEx id;
               downloaded with --download-attributes
    --protein-sheet, --rna-sheet   sheet names, when the letters D and A of
               the tables are not found in the sheet names

Outputs (in --out, default data/)
    jiang2020_protein_samples.tsv  gene, then one column per sample: log2
                           relative abundance as given
    jiang2020_rna_samples.tsv  gene, then one column per sample: TPM (the
                           log scale undone, base found from the column sums)
    jiang2020_samples.tsv  tissue, rna_sample, protein_sample, donor, and
                           SMTSISCH when known
    jiang2020_protein.tsv  with --tpm: gene, then the median protein of each
                           tissue, under its GTEx name (`regulome.py modes
                           --protein`); written when at least 10 tissues
                           match a column of --tpm
    and in --report (default results/jiang2020):
    layout.txt             every sheet of every file: first rows, first
                           columns
    conversion.txt         what was found: sheets, gene column, samples,
                           tissues, log base

Usage
    python3 scripts/jiang2020_tables.py --dir data --gtf data/gencode.gtf.gz \\
        --attributes data/GTEx_v8_SampleAttributesDS.txt

Method
    Sheets are read with openpyxl in read-only mode (rows streamed). The
    header is the row just above the first row holding a gene (Ensembl id
    or GTF symbol in one of its first five cells); the sample columns are
    those whose filled cells ('NA' counts as empty) are numeric for at least
    80%; the gene column is, of the other columns, the one with most values
    that are Ensembl gene ids or GTF symbols. A row above the header whose
    values repeat over the sample columns is read as a tissue row (as
    published: GTEx tissue names above the GTEx sample ids). RNA and protein
    samples are paired by name, by GTEx tissue sample (GTEX-donor-site: the
    RNA-seq and proteomics aliquots differ in their SM- suffix) or by donor
    and tissue, whichever pairs most. The tissue of a sample is, in this
    order: the tissue row, the GTEx attributes (exact sample id, or donor
    and tissue-site code), the abbreviations of Table S2 sheet A found as a
    token of the sample name. The RNA values are log_b(TPM + c), c = 0 when
    some are negative (as published: log2 TPM down to about -10) and 1
    otherwise, and b (2, e or 10) brings the median sample sum of
    protein-coding TPM closest to 10^6.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from re2g_panel import opener  # noqa: E402

ATTR_URL = ("https://storage.googleapis.com/adult-gtex/annotations/v8/"
            "metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
ENSG = re.compile(r"^ENSG\d{11}")
DONOR = re.compile(r"GTEX-[A-Z0-9]+", re.I)
SITE = re.compile(r"GTEX-[A-Z0-9]+-\d{4}", re.I)


def log(fh, *x):
    s = " ".join(str(v) for v in x)
    print(s, flush=True)
    fh.write(s + "\n")


def sheet_letter(name: str) -> str | None:
    """'D', from 'D', 'D. Cleaned ...', 'Table S2D', 'S2-D', 'S2_D ...'"""
    m = re.match(r"^\s*(?:table\s*)?(?:s\d+\s*[-_.]?\s*)?([A-Ha-h])(?:\b|[-_.\s])",
                 name.strip())
    return m.group(1).upper() if m else None


def open_book(path: Path):
    import openpyxl
    return openpyxl.load_workbook(path, read_only=True, data_only=True)


def pick_sheet(wb, letter: str, explicit: str | None):
    if explicit:
        return wb[explicit]
    for ws in wb.worksheets:
        if sheet_letter(ws.title) == letter:
            return ws
    i = ord(letter) - ord("A")
    names = [ws.title for ws in wb.worksheets]
    # sheet order A, B, C...; a leading "readme" sheet shifts it by one
    if names and not sheet_letter(names[0]) and re.search(
            r"readme|legend|content|descr", names[0], re.I):
        i += 1
    return wb.worksheets[i] if i < len(wb.worksheets) else None


def describe(path: Path, fh, rows=6, cols=12):
    wb = open_book(path)
    fh.write(f"### {path.name}\n")
    for ws in wb.worksheets:
        try:
            dim = ws.calculate_dimension()
        except Exception:  # unsized read-only sheet
            dim = "?"
        fh.write(f"## sheet {ws.title!r} (dimensions {dim})\n")
        for k, r in enumerate(ws.iter_rows(max_row=rows, values_only=True)):
            cells = ["" if v is None else str(v)[:30] for v in r[:cols]]
            fh.write(f"  row {k + 1} ({len(r)} cells): " + " | ".join(cells) + "\n")
    wb.close()


MISSING = ("", "na", "nan", "n/a", "none", "null", "-")


def missing(v) -> bool:
    return v is None or (isinstance(v, str) and v.strip().lower() in MISSING)


def is_num(v) -> bool:
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return True
    if isinstance(v, str):
        try:
            float(v)
            return v.strip().lower() not in ("nan", "")
        except ValueError:
            return False
    return False


def to_num(v) -> float:
    try:
        x = float(v)
        return x if math.isfinite(x) else math.nan
    except (TypeError, ValueError):
        return math.nan


def read_sheet(ws, names: set, rep, label: str):
    """(gene column values, sample names, tissue row or None, data rows)"""
    head_rows, data = [], []
    it = ws.iter_rows(values_only=True)
    for r in it:
        head_rows.append(list(r))
        if len(head_rows) >= 30:
            break
    width = max(len(r) for r in head_rows)
    # the header is the row just above the first row that starts with a
    # gene (Ensembl id or GTF symbol in one of its first 5 cells); else the
    # first row with half of its cells filled
    first = next((i for i, r in enumerate(head_rows) if i > 0 and any(
        isinstance(v, str) and (ENSG.match(v.strip()) or v.strip() in names)
        for v in r[:5])), None)
    hi = first - 1 if first is not None else next(
        i for i, r in enumerate(head_rows)
        if sum(v not in (None, "") for v in r) >= width / 2)
    header = ["" if v is None else str(v).strip() for v in head_rows[hi]]
    header += [""] * (width - len(header))
    above = head_rows[hi - 1] if hi > 0 else None
    data = [r for r in head_rows[hi + 1:]]
    for r in it:
        data.append(r)
    data = [list(r) + [None] * (width - len(r)) for r in data
            if any(v not in (None, "") for v in r)]
    probe = data[:300]
    # numeric: 80% of the filled cells ('NA' is empty), and some filled
    numeric = []
    for j in range(width):
        filled = [r[j] for r in probe if not missing(r[j])]
        numeric.append(len(filled) >= 0.1 * len(probe) and
                       sum(is_num(v) for v in filled) >= 0.8 * len(filled))
    # the gene column: Ensembl ids or GTF symbols
    best, best_n = None, 0
    for j in range(width):
        if numeric[j]:
            continue
        vals = [str(r[j]).strip() for r in probe if not missing(r[j])]
        n = sum(bool(ENSG.match(v)) or v in names for v in vals)
        if n > best_n:
            best, best_n = j, n
    if best is None:
        sys.exit(f"{label}: no column of Ensembl ids or gene symbols; see "
                 "the layout file")
    samples = [j for j in range(width) if numeric[j] and header[j]]
    log(rep, f"{label}: sheet {ws.title!r}, header row {hi + 1}, "
             f"gene column {header[best]!r} ({best_n}/{len(probe)} recognised), "
             f"{len(samples)} sample columns, {len(data)} rows; "
             f"first samples {[header[j] for j in samples[:4]]}")
    tissue_row = None
    if above is not None:
        t = [above[j] if j < len(above) else None for j in samples]
        filled = [str(v).strip() for v in t if v not in (None, "")]
        if len(filled) >= 0.9 * len(samples) and len(set(filled)) < len(filled):
            tissue_row = {header[j]: str(above[j]).strip() for j in samples}
            log(rep, f"{label}: tissue row above the header "
                     f"({len(set(filled))} values)")
    genes = []
    for r in data:
        g = str(r[best]).strip() if r[best] is not None else ""
        genes.append(g.split(".")[0] if ENSG.match(g) else g)
    M = [[to_num(r[j]) for j in samples] for r in data]
    return genes, [header[j] for j in samples], tissue_row, M


def key(s: str) -> str:
    return re.sub(r"[._\s]+", "-", s.strip().upper())


def gtf_names(gtf: Path):
    names, coding = set(), set()
    with opener(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "gene":
                continue
            a = f[8]
            i = a.find('gene_name "')
            if i < 0:
                continue
            n = a[i + 11:a.find('"', i + 11)]
            names.add(n)
            j = a.find('gene_id "')
            e = a[j + 9:a.find('"', j + 9)].split(".")[0]
            if 'gene_type "protein_coding"' in a:
                coding.update((n, e))
    return names, coding


def abbreviations(wb, rep) -> dict:
    """abbreviation -> tissue, from Table S2 sheet A"""
    ws = pick_sheet(wb, "A", None)
    out = {}
    if ws is None:
        return out
    for r in ws.iter_rows(max_row=200, max_col=10, values_only=True):
        cells = [str(v).strip() for v in r if v not in (None, "")]
        if len(cells) < 2:
            continue
        ab, full = sorted(cells[:2], key=len)       # either order
        if re.fullmatch(r"[A-Za-z0-9]{2,12}", ab) and len(full) > len(ab) \
                and not re.search(r"abbreviation", ab + full, re.I):
            out[ab.upper()] = full
    log(rep, f"Table S2 sheet {ws.title!r}: {len(out)} tissue abbreviations "
             f"{list(out.items())[:3]}")
    return out


def read_attributes(path: Path | None) -> dict:
    out = {}
    if not path or not path.exists():
        return out
    with open(path, encoding="utf-8", errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(head)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            s = f[c["SAMPID"]]
            out[s.upper()] = (f[c["SMTSD"]], f[c.get("SMTSISCH", 0)]
                              if "SMTSISCH" in c else "")
    return out


def tissue_of(name, tissue_row, attrs, by_site, abbr):
    if tissue_row and tissue_row.get(name):
        return tissue_row[name], "tissue row"
    k = key(name)
    if k in attrs:
        return attrs[k][0], "GTEx id"
    m = SITE.search(k)
    if m and m.group(0) in by_site:
        return by_site[m.group(0)], "GTEx donor and site"
    toks = re.split(r"[^A-Z0-9]+", k)
    hits = [abbr[t] for t in toks if t in abbr]
    if hits:
        return hits[0], "abbreviation"
    return None, None


def log_base(M, genes, coding, rep):
    """(b, c) of values log_b(TPM + c): c = 1, or 0 when some values are
    negative; b the base bringing the median column sum of protein-coding
    TPM closest to 1e6"""
    rows = [i for i, g in enumerate(genes) if g in coding]
    if len(rows) < 1000:
        rows = list(range(len(genes)))
    vmax = max((x for i in rows for x in M[i] if not math.isnan(x)), default=0)
    if vmax > 50:
        log(rep, f"RNA: maximum {vmax:.0f}, values taken as TPM (not log)")
        return None
    vmin = min((x for i in rows for x in M[i] if not math.isnan(x)), default=0)
    c = 0.0 if vmin < 0 else 1.0
    best, err = 2.0, math.inf
    for b in (2.0, math.e, 10.0):
        sums = sorted(sum(b ** M[i][j] - c for i in rows
                          if not math.isnan(M[i][j]))
                      for j in range(len(M[0])))
        med = sums[len(sums) // 2]
        e = abs(math.log10(max(med, 1e-9)) - 6)
        log(rep, f"RNA log base {b:.3g} (+{c:g}): median sample sum of TPM {med:,.0f}")
        if e < err:
            best, err = b, e
    log(rep, f"RNA: log base {best:.3g} taken, values log(TPM + {c:g}) "
             f"(minimum {vmin:.2f})")
    return best, c


def write_matrix(path: Path, genes, samples, M, f=None):
    seen = set()
    with open(path, "w") as fh:
        fh.write("gene\t" + "\t".join(samples) + "\n")
        for g, row in zip(genes, M):
            if not g or g in seen:
                continue
            seen.add(g)
            vals = [(f(x) if f else x) for x in row]
            fh.write(g + "\t" + "\t".join("" if math.isnan(x) else f"{x:.5g}"
                                          for x in vals) + "\n")


def tissue_table(path: Path, genes, samples, M, tissue, tpm: Path, rep):
    """median protein per tissue, columns named as in the GTEx TPM table"""
    with opener(tpm) as fh:
        fh.readline()
        fh.readline()
        cols = fh.readline().rstrip("\n").split("\t")[2:]
    norm = lambda x: re.sub(r"[^a-z0-9]", "", x.lower())
    gtex = {norm(c): c for c in cols}
    idx = defaultdict(list)
    for j, s in enumerate(samples):
        if tissue.get(s) and norm(tissue[s]) in gtex:
            idx[gtex[norm(tissue[s])]].append(j)
    miss = sorted({t for t in tissue.values() if t and norm(t) not in gtex})
    log(rep, f"protein per tissue: {len(idx)} tissues named as in {tpm.name}"
             + (f"; not matched: {miss[:8]}" if miss else ""))
    if len(idx) < 10:
        log(rep, "  fewer than 10: jiang2020_protein.tsv not written (modes "
                 "takes --tissue-map otherwise)")
        return
    ts = sorted(idx)
    seen = set()
    with open(path, "w") as fh:
        fh.write("gene\t" + "\t".join(ts) + "\n")
        for g, row in zip(genes, M):
            if not g or g in seen:
                continue
            seen.add(g)
            out = []
            for t in ts:
                v = sorted(row[j] for j in idx[t] if not math.isnan(row[j]))
                out.append(f"{(v[(len(v) - 1) // 2] + v[len(v) // 2]) / 2:.5g}"
                           if v else "")
            fh.write(g + "\t" + "\t".join(out) + "\n")
    log(rep, f"written: {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[1].strip(),
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", type=Path, default=Path("data"))
    p.add_argument("--s2", type=Path)
    p.add_argument("--s3", type=Path)
    p.add_argument("--gtf", type=Path, required=True)
    p.add_argument("--tpm", type=Path)
    p.add_argument("--attributes", type=Path)
    p.add_argument("--download-attributes", action="store_true")
    p.add_argument("--protein-sheet")
    p.add_argument("--rna-sheet")
    p.add_argument("--out", type=Path, default=Path("data"))
    p.add_argument("--report", type=Path, default=Path("results/jiang2020"))
    a = p.parse_args()
    a.report.mkdir(parents=True, exist_ok=True)
    a.out.mkdir(parents=True, exist_ok=True)
    s2 = a.s2 or next(iter(sorted(a.dir.glob("*mmc3.xlsx"))), None)
    s3 = a.s3 or next(iter(sorted(a.dir.glob("*mmc4.xlsx"))), None)
    if not s2 or not s3:
        sys.exit(f"Table S2 (*mmc3.xlsx) and S3 (*mmc4.xlsx) not found in {a.dir}")
    with open(a.report / "layout.txt", "w") as fh:
        for f in [s2, s3] + sorted(a.dir.glob("*mmc5.xlsx")):
            print(f"describing {f.name}", flush=True)
            describe(f, fh)
    rep = open(a.report / "conversion.txt", "w")
    if a.download_attributes and a.attributes and not a.attributes.exists():
        try:
            req = urllib.request.Request(ATTR_URL, headers={"User-Agent": "mitoloop_v0"})
            with urllib.request.urlopen(req, timeout=300) as r:
                a.attributes.write_bytes(r.read())
            log(rep, f"GTEx sample attributes -> {a.attributes}")
        except Exception as e:  # network errors of every kind
            log(rep, f"GTEx sample attributes not downloaded: {e}")
    names, coding = gtf_names(a.gtf)
    wb2, wb3 = open_book(s2), open_book(s3)
    abbr = abbreviations(wb2, rep)
    ws = pick_sheet(wb2, "D", a.protein_sheet)
    pg, ps, ptr, PM = read_sheet(ws, names, rep, "protein (S2 D)")
    ws = pick_sheet(wb3, "A", a.rna_sheet)
    rg, rs, rtr, RM = read_sheet(ws, names, rep, "RNA (S3 A)")
    attrs = read_attributes(a.attributes)
    by_site = {}
    for s, (t, _i) in attrs.items():
        m = SITE.search(s)
        if m:
            by_site.setdefault(m.group(0), t)
    if attrs:
        log(rep, f"GTEx attributes: {len(attrs)} samples")

    def annotate(samples, tr):
        out = {}
        how = Counter()
        for s in samples:
            t, h = tissue_of(s, tr, attrs, by_site, abbr)
            d = DONOR.search(s)
            isch = attrs.get(key(s), ("", ""))[1]
            out[s] = (t, d.group(0).upper() if d else "", isch)
            how[h] += 1
        return out, how

    pa, ph = annotate(ps, ptr)
    ra, rh = annotate(rs, rtr or ptr)
    log(rep, f"tissue of protein samples: {dict(ph)}; of RNA samples: {dict(rh)}")
    # pairing: same name; same GTEx tissue sample (GTEX-donor-site, the
    # RNA-seq and proteomics aliquots differ in their SM- suffix); same
    # donor and tissue. The rule pairing most samples is kept.
    def by(keyf):
        idx = {}
        for s in rs:
            k = keyf(s, ra[s])
            if k:
                idx.setdefault(k, s)
        out = []
        for s in ps:
            k = keyf(s, pa[s])
            if k and k in idx:
                out.append((idx[k], s))
        return out

    def site(s, _a):
        m = SITE.search(key(s))
        return m.group(0) if m else None

    rules = [("same name", by(lambda s, _a: key(s))),
             ("same GTEx tissue sample", by(site)),
             ("donor and tissue", by(lambda s, a: (a[1], a[0]) if a[0] and a[1] else None))]
    for name, pr in rules:
        log(rep, f"  pairing by {name}: {len(pr)} pairs")
    how, pairs = max(rules, key=lambda x: len(x[1]))
    log(rep, f"{len(pairs)} RNA-protein pairs ({how}) of {len(ps)} protein and "
             f"{len(rs)} RNA samples")
    pairs = [(r, q) for r, q in pairs if pa[q][0] or ra[r][0]]
    if len(pairs) < 20:
        log(rep, "too few pairs with a tissue: see layout.txt; first names: "
                 f"protein {ps[:6]}, RNA {rs[:6]}")
        sys.exit(1)
    ntis = len({pa[q][0] or ra[r][0] for r, q in pairs})
    log(rep, f"{len(pairs)} pairs in {ntis} tissues")
    with open(a.out / "jiang2020_samples.tsv", "w") as fh:
        fh.write("tissue\trna_sample\tprotein_sample\tdonor\tSMTSISCH\n")
        for r, q in pairs:
            t = pa[q][0] or ra[r][0]
            fh.write(f"{t}\t{r}\t{q}\t{pa[q][1] or ra[r][1]}\t{pa[q][2] or ra[r][2]}\n")
    b = log_base(RM, rg, coding, rep)
    write_matrix(a.out / "jiang2020_rna_samples.tsv", rg, rs, RM,
                 (lambda x: max(b[0] ** x - b[1], 0.0)) if b else None)
    write_matrix(a.out / "jiang2020_protein_samples.tsv", pg, ps, PM)
    log(rep, f"written: {a.out}/jiang2020_rna_samples.tsv ({len(set(rg))} genes), "
             f"jiang2020_protein_samples.tsv ({len(set(pg))} genes), "
             "jiang2020_samples.tsv")
    if a.tpm and a.tpm.exists():
        paired = {q for _r, q in pairs}
        tissue = {q: pa[q][0] or ra[r][0] for r, q in pairs}
        tissue.update({s: pa[s][0] for s in ps if s not in paired and pa[s][0]})
        tissue_table(a.out / "jiang2020_protein.tsv", pg, ps, PM, tissue,
                     a.tpm, rep)
    rep.close()


if __name__ == "__main__":
    main()
