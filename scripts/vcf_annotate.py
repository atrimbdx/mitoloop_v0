#!/usr/bin/env python3
"""
vcf_annotate.py - rare non-coding variants of WGS data in the regulatory
elements of a gene set, with frequency, scores, inheritance and second hits

Purpose
    Prepare the interpretation of non-coding variants in unsolved genomes
    (MitoLoop WP3). Two subcommands:

      regions   build one annotated BED of the regions attributed to each
                gene of the set: promoter, 5' and 3' UTR, coding exons
                (for second hits), gene body (for structural variants),
                replicable ENCODE-rE2G elements, loop anchors of
                v4c_loops.R, and any other BED (capture Hi-C or Micro-C
                interactions, elements validated by CRISPRi, ...).
      annotate  read one or more VCFs (single genomes, trios or a joint
                call), keep the rare variants of the probands that fall in
                these regions, add population frequency, variant scores
                (CADD, ReMM, AlphaGenome or any table) and the parental
                origin or phase, then rank each regulatory variant by the
                genetic model: a regulatory variant in trans of a rare
                coding variant of the same gene, homozygous or hemizygous
                regulatory variants, two regulatory variants in trans, de
                novo variants in dominant genes. Structural variants are
                read from the same or another VCF and reported when they
                delete, duplicate or break an element, or separate it from
                its promoter.

Inputs
    regions
      --gtf        GENCODE annotation (GTF). The transcript used per gene is
                   MANE Select, else Ensembl_canonical, else the longest.
      --genes      gene symbols, resolved as in re2g_panel.py (--aliases,
                   --no-builtin-aliases).
      --elements   optional re2g_elements.tsv (re2g_panel.py); elements
                   linked in >= --min-prevalence biosamples are kept.
      --loops      optional loops.tsv (v4c_loops.R); the partner anchor,
                   summit +/- --anchor-half-width, restricted with
                   --called-in.
      --bed        optional NAME=PATH, repeatable: chrom, start, end, and
                   in the 4th column the gene; without a gene column the
                   region goes to every gene of the set whose TSS lies
                   within --window. Extra columns are kept as detail.
    annotate
      --vcf        VCF(.gz), repeatable: small variants and/or structural
                   variants (SVTYPE, END or SVLEN in INFO, or symbolic
                   ALT). Multi-allelic records are handled allele by allele;
                   normalise first (bcftools norm -m -any -f ref.fa) so that
                   positions and alleles match the score tables.
      --regions    the BED written by `regions`.
      --ped        optional PED (family, individual, father, mother, sex,
                   phenotype). Probands are the affected individuals present
                   in a VCF; without a PED, every sample is a proband
                   without parents. --samples restricts the probands.
      --inheritance optional gene table: two columns (gene, mode with mode
                   among AR, AD, XL, AD/AR, MT) or a PanelApp download (Gene
                   Symbol or Entity Name, Model_Of_Inheritance).
      --score      NAME=PATH[:COLUMN], repeatable: a variant score table,
                   bgzipped and tabix-indexed (looked up by position) or
                   plain (read once). Columns: chrom, pos (1-based), then
                   ref and alt when the table is allele-specific (CADD,
                   AlphaGenome exports), or not (ReMM: chrom, pos, score).
                   COLUMN is a header name or a 1-based index; default
                   PHRED when present, else the last column.
      --info-score NAME=KEY, repeatable: a score already in the VCF INFO or
                   in the VEP CSQ fields (for instance CADD=CADD_PHRED).
      --gnomad     optional gnomAD sites VCF, tabix-indexed, local or URL;
                   {chrom} in the path is replaced by the chromosome
                   without prefix (write chr{chrom} for gnomAD v4 files).
                   --gnomad-keys are the INFO keys read (maximum taken).
      --af-keys    INFO or CSQ keys holding a population frequency already
                   annotated in the VCF (maximum taken).

Outputs
    regions     --out, a BED with header: chrom, start, end, gene, class
                (promoter, utr5, utr3, cds, exon_nc, gene_body, element,
                anchor, or the NAME of --bed), source, tss, strand, detail.
    annotate    in --out:
      candidates.tsv  one row per proband, variant and gene: tier, category,
                      gene, inheritance mode, role (regulatory or coding),
                      region classes and coordinates, variant, genotype,
                      zygosity, parental origin, phase block, frequency and
                      its source, carriers among the probands, VEP or
                      snpEff consequence, scores, score_flags (the
                      --min-score thresholds reached), partners (the other
                      rare variants of the gene in the proband with their
                      relation: trans, cis or unknown). Coding variants are
                      listed only in genes holding a regulatory candidate.
      genes.tsv       proband, gene, inheritance, regulatory and coding
                      variant counts, best tier and category
      cohort_genes.tsv  gene, probands with a tier 1 candidate, probands
                      with a tier 1 or 2 candidate, regulatory variants
      log.txt         counts at every filter

Method
    Region filter. A small variant is kept when its reference span overlaps
    a region; a structural variant when it overlaps a region or has a
    breakpoint between an element and the TSS of its gene.

    Genotype filter. FILTER is PASS or '.' with --pass-only; the proband's
    GQ and DP (when present) reach --min-gq and --min-dp.

    Frequency. The maximum over the --af-keys found in INFO or CSQ and the
    --gnomad lookup. Regulatory variants are rare when this maximum is at
    most --max-af, coding ones at most --max-af-coding; a variant without
    any frequency is kept and its af_source is empty. Carriers among the
    probands are counted across the VCFs, a recurrent variant in unrelated
    unsolved patients pointing to an artefact or a frequent allele missing
    from gnomAD.

    Roles. A variant is coding for a gene when VEP (CSQ) or snpEff (ANN)
    gives it an impact of --coding-impacts on that gene, or, without such an
    annotation, when it overlaps a 'cds' region of the gene (coding exon +/-
    --splice bp); it is regulatory when it overlaps any other region class
    except gene_body. A structural variant is regulatory for every gene
    whose element or promoter it hits, coding when it deletes or breaks a
    coding exon.

    Origin and phase. With both parents: paternal or maternal when one
    parent carries the allele, de novo when both are genotyped without it
    (each with DP >= --min-dp when DP is present), unknown otherwise.
    Without parents, phased genotypes ('|') of the same PS block (or of the
    chromosome when PS is absent) give the haplotype. Two variants are in
    trans when they come from different parents or sit on different
    haplotypes of one block, in cis when the same, unknown otherwise.

    Tiers, per regulatory variant of a gene:
      1  homozygous; hemizygous (male, chrX outside the GRCh38 PARs);
         in trans of a rare coding variant (trans_with_coding); in trans of
         another regulatory variant (compound_regulatory); de novo in a
         gene with a dominant mode (de_novo_dominant)
      2  with a coding or regulatory partner of unknown phase
         (possible_trans_with_coding, possible_compound_regulatory); de
         novo in another gene (de_novo); heterozygous in a dominant gene
         (dominant_gene)
      3  in cis only (cis_with_coding), single heterozygous variant in a
         recessive or unclassified gene (single_heterozygous), heterozygous
         in an X-linked gene in a female (x_linked_heterozygous)
    Within a tier, rows are sorted by the number of score thresholds
    reached, then by the first score.

    Structural variants. Effects per region: deleted or partially_deleted
    (DEL), duplicated or partially_duplicated (DUP), inverted (INV
    containing the region), breakpoint (a breakpoint inside the region),
    and, for an element or anchor and the TSS of its gene: separated (one
    breakpoint between them, INV or BND) or spans_between (a DEL, DUP or
    INV lying between them, which changes their distance).

Usage
    python3 vcf_annotate.py regions --gtf data/gencode.gtf.gz \\
        --genes genes.txt --elements results/re2g/re2g_elements.tsv \\
        --min-prevalence 51 --loops results/loops/loops.tsv \\
        --called-in "psoas muscle|gastrocnemius" \\
        --bed chic=results/chic/other_ends.bed --out results/wgs/regions.bed

    python3 vcf_annotate.py annotate --vcf cohort.snv.vcf.gz \\
        --vcf cohort.sv.vcf.gz --ped cohort.ped \\
        --regions results/wgs/regions.bed --inheritance panelapp_mito.tsv \\
        --gnomad "gnomad.genomes.v4.1.sites.{chrom}.vcf.bgz" \\
        --score CADD=whole_genome_SNVs.tsv.gz:PHRED \\
        --score ReMM=ReMM.v0.4.hg38.tsv.gz \\
        --score AlphaGenome=alphagenome_scores.tsv.gz:score \\
        --min-score CADD=15 --min-score ReMM=0.9 \\
        --out results/wgs/annot

References
    Rentzsch P et al. 2019. CADD: predicting the deleteriousness of variants
        throughout the human genome. Nucleic Acids Res 47:D886-D894.
    Smedley D et al. 2016. A whole-genome analysis framework for effective
        identification of pathogenic regulatory variants in Mendelian
        disease (ReMM). Am J Hum Genet 99:595-606.
    Avsec Z et al. 2026. AlphaGenome. Nature 649:1206-1218.
    Chen S et al. 2024. A genomic mutational constraint map using variation
        in 76,156 human genomes (gnomAD). Nature 625:92-100.
    McLaren W et al. 2016. The Ensembl Variant Effect Predictor. Genome Biol
        17:122.
    Cingolani P et al. 2012. A program for annotating and predicting the
        effects of single nucleotide polymorphisms, SnpEff. Fly 6:80-92.
    Martin AR et al. 2019. PanelApp. Nat Genet 51:1560-1565.
    Morales J et al. 2022. MANE. Nature 604:310-315.
    Martin M et al. 2016. WhatsHap: fast and accurate read-based phasing.
        bioRxiv 085050.
    Danecek P et al. 2021. Twelve years of SAMtools and BCFtools.
        GigaScience 10:giab008.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from re2g_panel import (NOT_GENE_PREFIXES, NOT_GENES, build_aliases,  # noqa: E402
                        opener, resolve_symbols)

# GRCh38 pseudo-autosomal regions of chrX, 0-based half-open
PAR_X = ((10_000, 2_781_479), (155_701_382, 156_030_895))
IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
# when one variant holds several region classes of a gene, the first listed
# is reported first
CLASS_ORDER = ["promoter", "utr5", "cds", "element", "anchor", "utr3",
               "exon_nc", "gene_body"]
REGION_COLS = ["chrom", "start", "end", "gene", "class", "source", "tss",
               "strand", "detail"]


def norm_chrom(c: str) -> str:
    c = c[3:] if c.startswith("chr") else c
    return "M" if c == "MT" else c


def read_table(path) -> list:
    with opener(path) as fh:
        head = fh.readline().rstrip("\n").lstrip("#").split("\t")
        rows = []
        for line in fh:
            if line.strip():
                f = line.rstrip("\n").split("\t")
                rows.append(dict(zip(head, f + [""] * (len(head) - len(f)))))
    return rows


# ===========================================================================
# regions

def load_models(gtf: Path, wanted: set) -> dict:
    """gene -> chosen transcript: chrom, strand, tss, exons, cds, gene span,
    ensembl id. MANE Select, else Ensembl_canonical, else the longest."""
    genes, tx = {}, {}
    with opener(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] not in ("gene", "transcript", "exon",
                                          "CDS"):
                continue
            a = f[8]
            if "_PAR_Y" in a:
                continue
            i = a.find('gene_name "')
            if i < 0:
                continue
            name = a[i + 11:a.find('"', i + 11)]
            if name not in wanted:
                continue
            s, e = int(f[3]) - 1, int(f[4])
            if f[2] == "gene":
                j = a.find('gene_id "')
                genes[name] = (f[0], s, e, f[6],
                               a[j + 9:a.find('"', j + 9)].split(".")[0])
                continue
            j = a.find('transcript_id "')
            if j < 0:
                continue
            tid = a[j + 15:a.find('"', j + 15)]
            t = tx.setdefault(tid, {"gene": name, "chrom": f[0],
                                    "strand": f[6], "exons": [], "cds": [],
                                    "rank": 0, "span": (s, e)})
            if f[2] == "transcript":
                t["span"] = (s, e)
                t["rank"] = 2 if "MANE_Select" in a else (
                    1 if "Ensembl_canonical" in a else 0)
            elif f[2] == "exon":
                t["exons"].append((s, e))
            else:
                t["cds"].append((s, e))
    best = {}
    for tid, t in tx.items():
        if not t["exons"]:
            continue
        length = sum(e - s for s, e in t["exons"])
        key = (t["rank"], bool(t["cds"]), length, tid)
        if t["gene"] not in best or key > best[t["gene"]][0]:
            best[t["gene"]] = (key, t)
    out = {}
    for g, (_k, t) in best.items():
        s, e = t["span"]
        gs = genes.get(g)
        out[g] = {"chrom": t["chrom"], "strand": t["strand"],
                  "tss": s if t["strand"] == "+" else e,
                  "exons": sorted(t["exons"]), "cds": sorted(t["cds"]),
                  "body": (gs[1], gs[2]) if gs else (s, e),
                  "ensg": gs[4] if gs else ""}
    return out


def cmd_regions(a) -> None:
    from re2g_panel import gtf_gene_names
    names = gtf_gene_names(a.gtf)
    aliases = build_aliases(a.aliases, not a.no_builtin_aliases)
    symbols = [g for g in Path(a.genes).read_text().split() if g]
    mapping, unresolved, _ = resolve_symbols(symbols, names, aliases,
                                             NOT_GENES, NOT_GENE_PREFIXES)
    wanted = set(mapping.values())
    models = load_models(a.gtf, wanted)
    for g in sorted(wanted - set(models)):
        unresolved.append(g)
    models = {g: m for g, m in models.items()
              if norm_chrom(m["chrom"]) != "M"}
    print(f"{len(models)} genes with a transcript model"
          + (f", {len(unresolved)} left out: {' '.join(sorted(unresolved))}"
             if unresolved else ""))

    rows = []

    def add(chrom, s, e, gene, cls, source, detail=""):
        m = models[gene]
        rows.append({"chrom": chrom, "start": max(0, int(s)), "end": int(e),
                     "gene": gene, "class": cls, "source": source,
                     "tss": m["tss"], "strand": m["strand"],
                     "detail": detail})

    for g, m in models.items():
        c, t, st = m["chrom"], m["tss"], m["strand"]
        if st == "+":
            add(c, t - a.upstream, t + a.downstream, g, "promoter", "gtf")
        else:
            add(c, t - a.downstream, t + a.upstream, g, "promoter", "gtf")
        add(c, m["body"][0], m["body"][1], g, "gene_body", "gtf")
        if m["cds"]:
            lo = min(s for s, _ in m["cds"])
            hi = max(e for _, e in m["cds"])
            for s, e in m["cds"]:
                add(c, s - a.splice, e + a.splice, g, "cds", "gtf")
            for s, e in m["exons"]:
                for ps, pe, side in ((s, min(e, lo), "low"),
                                     (max(s, hi), e, "high")):
                    if pe > ps:
                        five = (side == "low") == (st == "+")
                        add(c, ps, pe, g, "utr5" if five else "utr3", "gtf")
        else:
            for s, e in m["exons"]:
                add(c, s, e, g, "exon_nc", "gtf")

    if a.elements:
        n = 0
        for r in read_table(a.elements):
            g = mapping.get(r.get("gene", ""), r.get("gene", ""))
            if g not in models:
                continue
            try:
                if int(float(r["prevalence"])) < a.min_prevalence:
                    continue
                add(r["chrom"], int(r["start"]), int(r["end"]), g, "element",
                    "re2g", f"prevalence={r['prevalence']}")
                n += 1
            except (KeyError, ValueError):
                continue
        print(f"  {n} elements linked in >= {a.min_prevalence} biosamples")

    if a.loops:
        rx = re.compile(a.called_in, re.I) if a.called_in else None
        n = 0
        for r in read_table(a.loops):
            g = mapping.get(r.get("gene", ""), r.get("gene", ""))
            if g not in models:
                continue
            called = r.get("called", "")
            if rx and not any(rx.search(t) for t in called.split("|") if t):
                continue
            try:
                sm = int(float(r["summit"]))
            except (KeyError, ValueError):
                continue
            add(r["chrom"], sm - a.anchor_half_width,
                sm + a.anchor_half_width, g, "anchor", "loops",
                f"called={called}" if called else "")
            n += 1
        print(f"  {n} loop anchors")

    for spec in a.bed or []:
        if "=" not in spec:
            sys.exit(f"--bed {spec!r}: expected NAME=PATH")
        name, path = spec.split("=", 1)
        n = 0
        by_chrom = defaultdict(list)
        for g, m in models.items():
            by_chrom[norm_chrom(m["chrom"])].append((m["tss"], g))
        with opener(path) as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")) or \
                        not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                try:
                    s, e = int(f[1]), int(f[2])
                except (IndexError, ValueError):
                    continue
                g = f[3] if len(f) > 3 else ""
                g = mapping.get(g, g)
                detail = ";".join(f[4:])
                if g in models:
                    targets = [g]
                else:
                    mid = (s + e) // 2
                    targets = [x for t, x in by_chrom.get(norm_chrom(f[0]),
                                                          [])
                               if abs(t - mid) <= a.window]
                for t in targets:
                    add(f[0], s, e, t, name, path.split("/")[-1], detail)
                    n += 1
        print(f"  {n} regions from {name}")

    rows.sort(key=lambda r: (r["chrom"], r["start"], r["end"], r["gene"]))
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("#" + "\t".join(REGION_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in REGION_COLS) + "\n")
    tally = defaultdict(int)
    for r in rows:
        tally[r["class"]] += 1
    print(f"{len(rows)} regions written to {out}: "
          + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))


# ===========================================================================
# annotate: regions, scores, frequencies

class Regions:
    def __init__(self, path: Path):
        self.rows = read_table(path)
        if not self.rows or "class" not in self.rows[0]:
            sys.exit(f"{path}: not a regions file of `vcf_annotate.py "
                     "regions`")
        per = defaultdict(list)
        for i, r in enumerate(self.rows):
            r["start"], r["end"] = int(r["start"]), int(r["end"])
            r["tss"] = int(r["tss"]) if r["tss"] else None
            per[norm_chrom(r["chrom"])].append(i)
        self.idx = {}
        self.span = {}
        for c, ii in per.items():
            ii.sort(key=lambda i: self.rows[i]["start"])
            st = np.array([self.rows[i]["start"] for i in ii])
            en = np.array([self.rows[i]["end"] for i in ii])
            self.idx[c] = (st, en, np.array(ii), int((en - st).max()))
            # element-to-TSS spans, for structural variants
            jj = [i for i in ii if self.rows[i]["class"] not in
                  ("promoter", "cds", "utr5", "utr3", "exon_nc", "gene_body")
                  and self.rows[i]["tss"] is not None]
            if jj:
                mids = np.array([(self.rows[i]["start"] + self.rows[i]["end"])
                                 // 2 for i in jj])
                tss = np.array([self.rows[i]["tss"] for i in jj])
                self.span[c] = (np.minimum(mids, tss), np.maximum(mids, tss),
                                np.array(jj))

    def query(self, chrom, s, e) -> list:
        """indices of the regions overlapping [s, e)"""
        x = self.idx.get(norm_chrom(chrom))
        if x is None:
            return []
        st, en, ii, maxlen = x
        lo = int(np.searchsorted(st, s - maxlen, "right"))
        hi = int(np.searchsorted(st, e, "left"))
        return [int(ii[k]) for k in range(lo, hi) if en[k] > s]

    def between(self, chrom, bp) -> list:
        """element regions whose element-to-TSS span holds position bp"""
        x = self.span.get(norm_chrom(chrom))
        if x is None:
            return []
        lo, hi, ii = x
        return [int(i) for i in ii[(lo < bp) & (bp < hi)]]


def tabix_rows(path: str, chrom: str, pos: int):
    """Text rows of a tabix-indexed file at one 1-based position, trying the
    chromosome with and without 'chr'."""
    names = [chrom, "chr" + norm_chrom(chrom), norm_chrom(chrom)]
    if norm_chrom(chrom) == "M":
        names += ["chrM", "MT"]
    try:
        import pysam
        tb = _TABIX.get(path)
        if tb is None:
            tb = _TABIX[path] = pysam.TabixFile(path)
        contigs = set(tb.contigs)
        for c in dict.fromkeys(names):
            if c in contigs:
                return list(tb.fetch(c, pos - 1, pos))
        return []
    except ImportError:
        pass
    if not shutil.which("tabix"):
        sys.exit(f"{path}: an indexed lookup needs pysam or tabix")
    for c in dict.fromkeys(names):
        r = subprocess.run(["tabix", path, f"{c}:{pos}-{pos}"],
                           capture_output=True, text=True)
        if r.returncode == 0 and r.stdout:
            return r.stdout.splitlines()
    return []


_TABIX = {}


def is_indexed(path: str) -> bool:
    if "://" in path:
        return True
    return Path(path + ".tbi").exists() or Path(path + ".csi").exists()


class ScoreTable:
    """A variant score table, looked up by position (tabix) or read once."""

    def __init__(self, name, spec):
        self.name = name
        path, col = spec, None
        m = re.match(r"^(.*\.(?:gz|bgz|tsv|txt|bed))[:](.+)$", spec)
        if m:
            path, col = m.group(1), m.group(2)
        self.path = path
        head = None
        with opener(path) as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                f = line.rstrip("\n").lstrip("#").split("\t")
                if line.startswith("#") or not f[1].isdigit():
                    head = f
                    continue
                first = f
                break
            else:
                sys.exit(f"{path}: no data line")
        n = len(first)
        low = [h.lower() for h in head] if head else []
        self.ref = self.alt = None
        if head and "ref" in low and "alt" in low:
            self.ref, self.alt = low.index("ref"), low.index("alt")
        elif n >= 5 and re.fullmatch(r"[ACGTN]+", first[2] or "") and \
                re.fullmatch(r"[ACGTN]+", first[3] or ""):
            self.ref, self.alt = 2, 3
        if col is None:
            self.col = low.index("phred") if "phred" in low else n - 1
        elif col.isdigit():
            self.col = int(col) - 1
        elif head and col in head:
            self.col = head.index(col)
        else:
            sys.exit(f"--score {name}: column {col!r} not in {path}")
        self.indexed = is_indexed(path)
        self.mem = None

    def preload(self, keys: set) -> None:
        if self.indexed:
            return
        want = {(norm_chrom(c), p) for c, p, _r, _a in keys}
        self.mem = defaultdict(list)
        with opener(self.path) as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                try:
                    k = (norm_chrom(f[0]), int(f[1]))
                except (IndexError, ValueError):
                    continue
                if k in want:
                    self.mem[k].append(f)

    def lookup(self, chrom, pos, ref, alt):
        if self.mem is not None:
            rows = self.mem.get((norm_chrom(chrom), pos), [])
        else:
            rows = [r.split("\t") for r in tabix_rows(self.path, chrom, pos)]
        best = None
        for f in rows:
            try:
                if int(f[1]) != pos:
                    continue
                if self.ref is not None and (f[self.ref] != ref or
                                             f[self.alt] != alt):
                    continue
                v = float(f[self.col])
            except (IndexError, ValueError):
                continue
            best = v if best is None else max(best, v)
        return best


class GnomadVcf:
    """Allele frequency from a gnomAD sites VCF, by position."""

    def __init__(self, spec, keys):
        self.spec, self.keys = spec, keys
        self.mem = None

    def path(self, chrom):
        return self.spec.replace("{chrom}", norm_chrom(chrom))

    def preload(self, keys: set) -> None:
        if "{chrom}" in self.spec or is_indexed(self.spec):
            return
        want = {(norm_chrom(c), p) for c, p, _r, _a in keys}
        self.mem = defaultdict(list)
        with opener(self.spec) as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.split("\t", 8)
                try:
                    k = (norm_chrom(f[0]), int(f[1]))
                except (IndexError, ValueError):
                    continue
                if k in want:
                    self.mem[k].append(f)

    def lookup(self, chrom, pos, ref, alt):
        if self.mem is not None:
            rows = self.mem.get((norm_chrom(chrom), pos), [])
        else:
            p = self.path(chrom)
            if "://" not in p and not Path(p).exists():
                return None
            rows = [r.split("\t", 8) for r in tabix_rows(p, chrom, pos)]
        best = None
        for f in rows:
            if len(f) < 8 or int(f[1]) != pos or f[3] != ref:
                continue
            alts = f[4].split(",")
            if alt not in alts:
                continue
            k = alts.index(alt)
            info = parse_info(f[7])
            for key in self.keys:
                v = info.get(key)
                if v is None:
                    continue
                try:
                    x = float(v.split(",")[k] if "," in v else v)
                except (ValueError, IndexError):
                    continue
                best = x if best is None else max(best, x)
        return best


def parse_info(s: str) -> dict:
    out = {}
    if s == ".":
        return out
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        else:
            out[kv] = True
    return out


def load_inheritance(path) -> dict:
    """gene -> AR, AD, XL, AD/AR, MT or '' from a two-column table or a
    PanelApp download."""
    out = {}
    with opener(path) as fh:
        first = fh.readline().rstrip("\n").split("\t")
        low = [h.lower() for h in first]
        gcol = mcol = None
        for i, h in enumerate(low):
            if h in ("gene symbol", "gene_symbol", "entity name",
                     "entity_name", "gene", "symbol"):
                gcol = i if gcol is None else gcol
            if "inheritance" in h or h == "mode":
                mcol = i
        rows = []
        if gcol is None or mcol is None:
            gcol, mcol = 0, 1
            rows.append(first)
        rows += [line.rstrip("\n").split("\t") for line in fh]
    for f in rows:
        if len(f) <= max(gcol, mcol) or not f[gcol]:
            continue
        m = f[mcol].upper()
        if m in ("AR", "AD", "XL", "AD/AR", "MT", "XLR", "XLD"):
            mode = "XL" if m.startswith("XL") else m
        elif "X-LINKED" in m or "X LINKED" in m:
            mode = "XL"
        elif m.startswith("BOTH") or ("MONOALLELIC" in m and "BIALLELIC" in m):
            mode = "AD/AR"
        elif "MONOALLELIC" in m:
            mode = "AD"
        elif "BIALLELIC" in m:
            mode = "AR"
        elif "MITOCHONDRIAL" in m:
            mode = "MT"
        else:
            mode = ""
        out[f[gcol].strip()] = mode
    return out


def read_ped(path) -> dict:
    """individual -> family, father, mother, sex (1 male, 2 female), affected"""
    out = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 6:
                continue
            out[f[1]] = {"family": f[0],
                         "father": f[2] if f[2] not in ("0", ".") else "",
                         "mother": f[3] if f[3] not in ("0", ".") else "",
                         "sex": f[4], "affected": f[5] == "2"}
    return out


# ===========================================================================
# annotate: VCF records

def gt_alleles(gt: str):
    if gt in (".", "./.", ".|."):
        return None, False
    phased = "|" in gt
    al = []
    for x in re.split(r"[/|]", gt):
        al.append(None if x == "." else int(x))
    return al, phased


def csq_fields(header_lines) -> tuple:
    csq = ann = None
    for h in header_lines:
        if h.startswith("##INFO=<ID=CSQ"):
            m = re.search(r"Format: ([^\"]+)", h)
            if m:
                csq = m.group(1).strip().split("|")
        elif h.startswith("##INFO=<ID=ANN"):
            ann = True
    return csq, ann


def vep_allele(ref, alts, alt):
    """VEP's allele string: the shared first base is removed from indels."""
    if all(len(x) == 1 for x in [ref] + alts):
        return alt
    if all(x[:1] == ref[:1] for x in alts if not x.startswith("<")):
        return alt[1:] or "-"
    return alt


def allele_consequences(info, csq_names, has_ann, ref, alts, k):
    """For ALT allele k (1-based): gene -> (impact, consequences), and the
    CSQ field values of this allele (for frequencies and scores)."""
    genes = defaultdict(lambda: ["", set()])
    fields = defaultdict(list)
    alt = alts[k - 1]
    if csq_names and "CSQ" in info:
        col = {n: i for i, n in enumerate(csq_names)}
        va = vep_allele(ref, alts, alt)
        entries = [e.split("|") for e in info["CSQ"].split(",")]
        if "ALLELE_NUM" in col:
            mine = [e for e in entries if e[col["ALLELE_NUM"]] == str(k)]
        else:
            mine = [e for e in entries if e[col.get("Allele", 0)] in (va, alt)]
            if not mine and len(alts) == 1:
                mine = entries
        for e in mine:
            g = e[col["SYMBOL"]] if "SYMBOL" in col else ""
            imp = e[col["IMPACT"]] if "IMPACT" in col else ""
            cons = e[col["Consequence"]] if "Consequence" in col else ""
            if g:
                cur = genes[g]
                if IMPACT_RANK.get(imp, -1) > IMPACT_RANK.get(cur[0], -1):
                    cur[0] = imp
                cur[1].update(c for c in cons.split("&") if c)
            for n, i in col.items():
                if i < len(e) and e[i] not in ("", "."):
                    fields[n].append(e[i])
    if has_ann and "ANN" in info:
        for e in info["ANN"].split(","):
            f = e.split("|")
            if len(f) < 4 or f[0] != alt:
                continue
            cur = genes[f[3]]
            if IMPACT_RANK.get(f[2], -1) > IMPACT_RANK.get(cur[0], -1):
                cur[0] = f[2]
            cur[1].update(c for c in f[1].split("&") if c)
    return genes, fields


def max_float(values):
    best = None
    for v in values:
        for x in str(v).replace("&", ",").split(","):
            try:
                y = float(x)
            except ValueError:
                continue
            if not np.isnan(y):
                best = y if best is None else max(best, y)
    return best


def sv_interval(pos, ref, alt, info):
    """(svtype, start0, end) of a structural variant, or None."""
    svtype = info.get("SVTYPE")
    if svtype is None and alt.startswith("<"):
        svtype = alt.strip("<>").split(":")[0]
    if svtype is None and ("[" in alt or "]" in alt):
        svtype = "BND"
    if svtype is None:
        return None
    s = pos - 1
    if "END" in info and svtype != "BND":
        e = int(info["END"])
    elif "SVLEN" in info and svtype not in ("BND", "INS"):
        e = pos + abs(int(str(info["SVLEN"]).split(",")[0]))
    else:
        e = pos
    return svtype, s, max(e, pos)


def in_par(chrom, pos0):
    return norm_chrom(chrom) == "X" and any(a <= pos0 < b for a, b in PAR_X)


def relation(v1, v2) -> str:
    """trans, cis or unknown, from parental origin, else from phase."""
    o1, o2 = v1["origin"], v2["origin"]
    if o1 in ("paternal", "maternal") and o2 in ("paternal", "maternal"):
        return "trans" if o1 != o2 else "cis"
    p1, p2 = v1["phase"], v2["phase"]
    if p1 and p2 and p1[0] == p2[0]:
        return "trans" if p1[1] != p2[1] else "cis"
    return "unknown"


def cmd_annotate(a) -> None:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    log = []

    def say(msg):
        print(msg)
        log.append(msg)

    regions = Regions(a.regions)
    say(f"{len(regions.rows)} regions from {a.regions}")
    modes = load_inheritance(a.inheritance) if a.inheritance else {}
    ped = read_ped(a.ped) if a.ped else {}
    only = set(a.samples.split(",")) if a.samples else None
    coding_imp = set(a.coding_impacts.split(","))
    af_keys = [k for k in a.af_keys.split(",") if k]
    info_scores = []
    for spec in a.info_score or []:
        n, k = spec.split("=", 1)
        info_scores.append((n, k))
    tables = []
    for spec in a.score or []:
        if "=" not in spec:
            sys.exit(f"--score {spec!r}: expected NAME=PATH[:COLUMN]")
        n, p_ = spec.split("=", 1)
        tables.append(ScoreTable(n, p_))
    thresholds = {}
    for spec in a.min_score or []:
        n, v = spec.split("=", 1)
        thresholds[n] = float(v)
    score_names = [n for n, _ in info_scores] + [t.name for t in tables]
    for n in thresholds:
        if n not in score_names:
            sys.exit(f"--min-score {n}: no --score or --info-score of that "
                     "name")
    gnomad = (GnomadVcf(a.gnomad, a.gnomad_keys.split(","))
              if a.gnomad else None)

    calls = []          # one dict per proband x allele kept
    carriers = defaultdict(set)
    counts = defaultdict(int)
    probands_all = set()
    for vcf in a.vcf:
        head_lines, samples = [], []
        with opener(vcf) as fh:
            for line in fh:
                if line.startswith("##"):
                    head_lines.append(line)
                    continue
                if line.startswith("#CHROM"):
                    samples = line.rstrip("\n").split("\t")[9:]
                    break
            csq_names, has_ann = csq_fields(head_lines)
            sidx = {s: i for i, s in enumerate(samples)}
            fam = []
            for s in samples:
                if only is not None and s not in only:
                    continue
                p = ped.get(s)
                if ped and (p is None or not p["affected"]):
                    if only is None:
                        continue
                fam.append({"id": s, "i": sidx[s],
                            "family": p["family"] if p else s,
                            "father": sidx.get(p["father"]) if p and
                            p["father"] else None,
                            "mother": sidx.get(p["mother"]) if p and
                            p["mother"] else None,
                            "sex": p["sex"] if p else "0"})
            if not fam:
                say(f"{vcf}: no proband among its {len(samples)} samples")
                continue
            probands_all.update(f["id"] for f in fam)
            say(f"{vcf}: {len(fam)} probands ("
                f"{sum(f['father'] is not None and f['mother'] is not None for f in fam)}"
                " with both parents)")
            for line in fh:
                counts["records"] += 1
                # the sample columns are split only for records kept
                f = line.split("\t", 8)
                chrom, pos, ref = f[0], int(f[1]), f[3]
                if norm_chrom(chrom) == "M":
                    continue
                alts = f[4].split(",")
                info = None
                sv = None
                if "SVTYPE" in f[7] or f[4].startswith("<") or \
                        "[" in f[4] or "]" in f[4]:
                    info = parse_info(f[7])
                    sv = sv_interval(pos, ref, alts[0], info)
                if sv:
                    svtype, s0, e0 = sv
                    bps = [s0] if svtype in ("BND", "INS") else [s0, e0]
                    hits = regions.query(chrom, s0, max(e0, s0 + 1))
                    btw = set()
                    for bp in bps:
                        btw.update(regions.between(chrom, bp))
                    if not hits and not btw:
                        continue
                else:
                    s0, e0 = pos - 1, pos - 1 + len(ref)
                    hits = regions.query(chrom, s0, e0)
                    btw = set()
                    if not hits:
                        continue
                counts["records_in_regions"] += 1
                f = line.rstrip("\n").split("\t")
                if a.pass_only and f[6] not in ("PASS", "."):
                    counts["filtered_FILTER"] += 1
                    continue
                if info is None:
                    info = parse_info(f[7])
                fmt = f[8].split(":") if len(f) > 8 else []
                fi = {k: i for i, k in enumerate(fmt)}

                def sample_field(si, key):
                    if si is None or key not in fi:
                        return None
                    v = f[9 + si].split(":")
                    return v[fi[key]] if fi[key] < len(v) else None

                for k in range(1, len(alts) + 1):
                    alt = alts[k - 1]
                    if alt in ("*", "."):
                        continue
                    genes_csq, csqf = allele_consequences(
                        info, csq_names, has_ann, ref, alts, k)
                    for pr in fam:
                        gt = sample_field(pr["i"], "GT")
                        al, phased = gt_alleles(gt or ".")
                        if not al or k not in al:
                            continue
                        gq = sample_field(pr["i"], "GQ")
                        dp = sample_field(pr["i"], "DP")
                        try:
                            if gq not in (None, ".") and float(gq) < a.min_gq:
                                counts["filtered_GQ"] += 1
                                continue
                            if dp not in (None, ".") and float(dp) < a.min_dp:
                                counts["filtered_DP"] += 1
                                continue
                        except ValueError:
                            pass
                        n_alt = sum(x == k for x in al)
                        male = pr["sex"] == "1"
                        hemi = (norm_chrom(chrom) == "X" and male and
                                not in_par(chrom, s0)) or \
                            (len(al) == 1 and norm_chrom(chrom) in ("X", "Y"))
                        zyg = "hemizygous" if hemi else (
                            "homozygous" if n_alt >= 2 else "heterozygous")
                        par = {}
                        for role in ("father", "mother"):
                            si = pr[role]
                            if si is None:
                                par[role] = None
                                continue
                            pal, _ = gt_alleles(sample_field(si, "GT") or ".")
                            pdp = sample_field(si, "DP")
                            low = False
                            try:
                                low = pdp not in (None, ".") and \
                                    float(pdp) < a.min_dp
                            except ValueError:
                                pass
                            if pal is None or None in pal:
                                par[role] = None
                            elif k in pal:
                                par[role] = True
                            else:
                                par[role] = None if low else False
                        if par["father"] and par["mother"]:
                            origin = "both_parents"
                        elif par["father"]:
                            origin = "paternal"
                        elif par["mother"]:
                            origin = "maternal"
                        elif par["father"] is False and par["mother"] is False:
                            origin = "de_novo"
                        else:
                            origin = "unknown"
                        phase = None
                        if phased and len(al) == 2 and n_alt == 1:
                            ps = sample_field(pr["i"], "PS")
                            block = f"{chrom}:{ps}" if ps not in (None, ".") \
                                else f"{chrom}:phased"
                            phase = (block, al.index(k))
                        key = (chrom, pos, ref, alt)
                        carriers[key].add(pr["id"])
                        calls.append({
                            "proband": pr["id"], "family": pr["family"],
                            "male": male, "chrom": chrom, "pos": pos,
                            "ref": ref, "alt": alt, "id": f[2],
                            "filter": f[6], "gt": gt, "zygosity": zyg,
                            "origin": origin, "phase": phase,
                            "gq": gq or "", "dp": dp or "", "sv": sv,
                            "hits": hits, "between": sorted(btw),
                            "csq": genes_csq, "csqf": csqf, "info": info,
                            "key": key, "k": k, "n_alts": len(alts)})
    say(f"{counts['records']} records read, {counts['records_in_regions']} "
        f"in the regions; {len(calls)} proband alleles kept "
        f"(FILTER {counts['filtered_FILTER']}, GQ {counts['filtered_GQ']}, "
        f"DP {counts['filtered_DP']} removed)")

    # ------------------------------------------------ frequency and scores
    keys = {c["key"] for c in calls if not c["sv"]}
    for t in tables:
        t.preload(keys)
    if gnomad:
        gnomad.preload(keys)
    cache_score, cache_af = {}, {}
    for c in calls:
        afs = []
        for k_ in af_keys:
            v = None
            if k_ in c["info"]:
                v = c["info"][k_]
                parts = str(v).split(",")
                # Number=A values: the one of this allele
                if len(parts) == c["n_alts"] > 1:
                    v = parts[c["k"] - 1]
                v = max_float([v])
            if v is None and k_ in c["csqf"]:
                v = max_float(c["csqf"][k_])
            if v is not None:
                afs.append((v, k_))
        if gnomad and not c["sv"]:
            if c["key"] not in cache_af:
                cache_af[c["key"]] = gnomad.lookup(*c["key"])
            v = cache_af[c["key"]]
            if v is not None:
                afs.append((v, "gnomad"))
        if afs:
            c["af"], c["af_source"] = max(afs)
        else:
            c["af"], c["af_source"] = None, ""
        sc = {}
        for n, k_ in info_scores:
            v = max_float([c["info"][k_]]) if k_ in c["info"] else None
            if v is None and k_ in c["csqf"]:
                v = max_float(c["csqf"][k_])
            sc[n] = v
        if not c["sv"]:
            if c["key"] not in cache_score:
                cache_score[c["key"]] = {t.name: t.lookup(*c["key"])
                                         for t in tables}
            sc.update(cache_score[c["key"]])
        c["scores"] = sc
        c["flags"] = [n for n, thr in thresholds.items()
                      if sc.get(n) is not None and sc[n] >= thr]

    # ------------------------------------------------------- gene hits
    rows_by = defaultdict(list)          # (proband, gene) -> variant hits
    n_common = n_noscore = 0
    for c in calls:
        per_gene = defaultdict(lambda: {"classes": [], "regions": [],
                                        "sources": [], "details": [],
                                        "effects": [], "tss": None})
        if c["sv"]:
            svtype, s0, e0 = c["sv"]
            bps = [s0] if svtype in ("BND", "INS") else [s0, e0]
            for i in c["hits"]:
                r = regions.rows[i]
                inside = s0 <= r["start"] and r["end"] <= e0
                if svtype in ("DEL", "DUP", "INV", "CNV") and inside:
                    eff = {"DEL": "deleted", "DUP": "duplicated",
                           "INV": "inverted", "CNV": "copy_number"}[svtype]
                elif any(r["start"] <= b < r["end"] for b in bps):
                    eff = "breakpoint"
                    if svtype in ("DEL", "DUP"):
                        eff = "partially_" + ("deleted" if svtype == "DEL"
                                              else "duplicated")
                elif svtype in ("DEL", "DUP", "INV", "CNV"):
                    # the region lies inside the variant without a
                    # breakpoint: inside covers it; nothing else to say
                    continue
                else:
                    continue
                if r["class"] == "gene_body" and eff not in ("deleted",
                                                             "breakpoint"):
                    continue
                _add_hit(per_gene, r, eff)
            for i in c["between"]:
                r = regions.rows[i]
                lo = min((r["start"] + r["end"]) // 2, r["tss"])
                hi = max((r["start"] + r["end"]) // 2, r["tss"])
                inside = [lo < b < hi for b in bps]
                if svtype in ("BND", "INS"):
                    eff = "separated" if svtype == "BND" else "spans_between"
                elif all(inside):
                    eff = "spans_between"
                elif svtype == "INV":
                    eff = "separated"
                else:
                    # a DEL or DUP with one end between the element and the
                    # TSS overlaps one of them, reported above
                    continue
                _add_hit(per_gene, r, eff)
        else:
            for i in c["hits"]:
                r = regions.rows[i]
                if r["class"] == "gene_body":
                    continue
                _add_hit(per_gene, r, "in_region")
        if not per_gene:
            continue
        for g, h in per_gene.items():
            imp = c["csq"].get(g, ["", set()])
            if c["csq"]:
                coding = imp[0] in coding_imp
            else:
                coding = "cds" in h["classes"]
            if c["sv"]:
                coding = "cds" in h["classes"] and any(
                    e in ("deleted", "partially_deleted", "breakpoint",
                          "inverted") for e in h["effects"])
            reg = any(x not in ("cds", "gene_body") for x in h["classes"])
            if not coding and not reg:
                continue
            role = "coding" if coding else "regulatory"
            limit = a.max_af_coding if coding else a.max_af
            if c["af"] is not None and c["af"] > limit:
                n_common += 1
                continue
            if role == "regulatory" and a.require_score and not c["sv"] \
                    and not c["flags"]:
                n_noscore += 1
                continue
            rows_by[(c["proband"], g)].append({
                **c, "gene": g, "role": role, "impact": imp[0],
                "consequence": "&".join(sorted(imp[1])), **h})
    say(f"{n_common} gene hits above the frequency limits removed"
        + (f", {n_noscore} regulatory hits without a score threshold "
           f"reached removed" if a.require_score else ""))

    # ------------------------------------------------------------ tiers
    out_rows, gene_rows = [], []
    for (pr, g), vs in rows_by.items():
        mode = modes.get(g, "")
        dominant = mode in ("AD", "AD/AR")
        regs = [v for v in vs if v["role"] == "regulatory"]
        cods = [v for v in vs if v["role"] == "coding"]
        if not regs:
            continue
        for v in regs:
            parts = []
            for w in vs:
                if w is v:
                    continue
                parts.append((w, relation(v, w)))
            cat, tier = None, 3
            rel_c = [r for w, r in parts if w["role"] == "coding"]
            rel_r = [r for w, r in parts if w["role"] == "regulatory"]
            if v["zygosity"] == "homozygous":
                cat, tier = "homozygous", 1
            elif v["zygosity"] == "hemizygous":
                cat, tier = "hemizygous", 1
            elif "trans" in rel_c:
                cat, tier = "trans_with_coding", 1
            elif "trans" in rel_r:
                cat, tier = "compound_regulatory", 1
            elif v["origin"] == "de_novo" and dominant:
                cat, tier = "de_novo_dominant", 1
            elif "unknown" in rel_c:
                cat, tier = "possible_trans_with_coding", 2
            elif "unknown" in rel_r:
                cat, tier = "possible_compound_regulatory", 2
            elif v["origin"] == "de_novo":
                cat, tier = "de_novo", 2
            elif dominant:
                cat, tier = "dominant_gene", 2
            elif "cis" in rel_c:
                cat, tier = "cis_with_coding", 3
            elif mode == "XL":
                cat, tier = "x_linked_heterozygous", 3
            else:
                cat, tier = "single_heterozygous", 3
            v["tier"], v["category"] = tier, cat
            v["partners"] = ";".join(
                f"{label(w)}({w['role']},{r})" for w, r in parts)
            out_rows.append(v)
        for w in cods:
            w["tier"], w["category"] = "", "coding_partner"
            w["partners"] = ";".join(
                f"{label(v)}(regulatory,{relation(w, v)})" for v in regs)
            out_rows.append(w)
        best = min(regs, key=lambda v: v["tier"])
        gene_rows.append({"proband": pr, "gene": g, "inheritance": mode,
                          "n_regulatory": len(regs), "n_coding": len(cods),
                          "best_tier": best["tier"],
                          "best_category": best["category"]})

    def sort_key(v):
        first = v["scores"].get(score_names[0]) if score_names else None
        return (v["tier"] if v["tier"] != "" else 9, v["proband"], v["gene"],
                v["role"] != "regulatory", -len(v["flags"]),
                -(first if first is not None else -1e9))
    out_rows.sort(key=sort_key)

    cols = (["proband", "family", "tier", "category", "gene", "inheritance",
             "role", "region_class", "region", "region_source",
             "region_detail", "distance_to_tss", "chrom", "pos", "ref", "alt",
             "id", "svtype", "sv_end", "effect", "filter", "genotype",
             "zygosity", "origin", "phase_block", "gq", "dp", "af",
             "af_source", "cohort_carriers", "impact", "consequence"]
            + score_names + ["score_flags", "partners"])
    with open(out / "candidates.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for v in out_rows:
            order = sorted(range(len(v["classes"])),
                           key=lambda i: CLASS_ORDER.index(v["classes"][i])
                           if v["classes"][i] in CLASS_ORDER else 5.5)
            pick = lambda key: "|".join(dict.fromkeys(  # noqa: E731
                str(v[key][i]) for i in order))
            mid = v["pos"]
            row = {
                "proband": v["proband"], "family": v["family"],
                "tier": v["tier"], "category": v["category"],
                "gene": v["gene"], "inheritance": modes.get(v["gene"], ""),
                "role": v["role"], "region_class": pick("classes"),
                "region": pick("regions"), "region_source": pick("sources"),
                "region_detail": pick("details"),
                "distance_to_tss": "" if v["tss"] is None else
                mid - 1 - v["tss"],
                "chrom": v["chrom"], "pos": v["pos"], "ref": v["ref"],
                "alt": v["alt"], "id": v["id"],
                "svtype": v["sv"][0] if v["sv"] else "",
                "sv_end": v["sv"][2] if v["sv"] else "",
                "effect": pick("effects"), "filter": v["filter"],
                "genotype": v["gt"], "zygosity": v["zygosity"],
                "origin": v["origin"],
                "phase_block": f"{v['phase'][0]}#{v['phase'][1]}"
                if v["phase"] else "",
                "gq": v["gq"], "dp": v["dp"],
                "af": "" if v["af"] is None else f"{v['af']:.3g}",
                "af_source": v["af_source"],
                "cohort_carriers": len(carriers[v["key"]]),
                "impact": v["impact"], "consequence": v["consequence"],
                "score_flags": ",".join(v["flags"]),
                "partners": v["partners"]}
            for n in score_names:
                x = v["scores"].get(n)
                row[n] = "" if x is None else f"{x:.4g}"
            fh.write("\t".join(str(row[c_]) for c_ in cols) + "\n")

    gcols = ["proband", "gene", "inheritance", "n_regulatory", "n_coding",
             "best_tier", "best_category"]
    gene_rows.sort(key=lambda r: (r["best_tier"], r["proband"], r["gene"]))
    with open(out / "genes.tsv", "w") as fh:
        fh.write("\t".join(gcols) + "\n")
        for r in gene_rows:
            fh.write("\t".join(str(r[c_]) for c_ in gcols) + "\n")

    coh = defaultdict(lambda: [set(), set(), 0])
    for r in gene_rows:
        x = coh[r["gene"]]
        x[2] += r["n_regulatory"]
        if r["best_tier"] == 1:
            x[0].add(r["proband"])
        if r["best_tier"] in (1, 2):
            x[1].add(r["proband"])
    with open(out / "cohort_genes.tsv", "w") as fh:
        fh.write("gene\tinheritance\tprobands_tier1\tprobands_tier1_2\t"
                 "regulatory_variants\n")
        for g, (t1, t12, n) in sorted(coh.items(),
                                      key=lambda x: (-len(x[1][0]),
                                                     -len(x[1][1]), x[0])):
            fh.write(f"{g}\t{modes.get(g, '')}\t{len(t1)}\t{len(t12)}\t{n}\n")

    tally = defaultdict(int)
    for v in out_rows:
        if v["role"] == "regulatory":
            tally[(v["tier"], v["category"])] += 1
    say(f"{len(probands_all)} probands, {sum(tally.values())} regulatory "
        f"candidates:")
    for (t, cat), n in sorted(tally.items()):
        say(f"  tier {t}  {cat:30s} {n}")
    (out / "log.txt").write_text("\n".join(log) + "\n")
    print(f"Written to {out}")


def label(v) -> str:
    if v["sv"]:
        return f"{v['chrom']}:{v['pos']}-{v['sv'][2]}:{v['sv'][0]}"
    return f"{v['chrom']}:{v['pos']}{v['ref']}>{v['alt']}"


def _add_hit(per_gene, r, effect):
    h = per_gene[r["gene"]]
    h["classes"].append(r["class"])
    h["regions"].append(f"{r['chrom']}:{r['start']}-{r['end']}")
    h["sources"].append(r["source"])
    h["details"].append(r["detail"])
    h["effects"].append(effect)
    if h["tss"] is None:
        h["tss"] = r["tss"]


# ===========================================================================

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("regions", help="build the annotated region BED")
    r.add_argument("--gtf", required=True, type=Path)
    r.add_argument("--genes", required=True, type=Path)
    r.add_argument("--aliases", type=Path)
    r.add_argument("--no-builtin-aliases", action="store_true")
    r.add_argument("--upstream", type=int, default=1000,
                   help="promoter length upstream of the TSS, bp "
                        "(default %(default)s)")
    r.add_argument("--downstream", type=int, default=500,
                   help="promoter length downstream of the TSS, bp "
                        "(default %(default)s)")
    r.add_argument("--splice", type=int, default=10,
                   help="bp added on each side of a coding exon "
                        "(default %(default)s)")
    r.add_argument("--elements", type=Path)
    r.add_argument("--min-prevalence", type=int, default=1,
                   help="biosamples an rE2G element must be linked in "
                        "(default %(default)s)")
    r.add_argument("--loops", type=Path)
    r.add_argument("--called-in", default="",
                   help="regular expression on the terms that called a loop, "
                        "case ignored")
    r.add_argument("--anchor-half-width", type=int, default=5000,
                   help="half-width of a loop anchor around its summit, bp "
                        "(default %(default)s)")
    r.add_argument("--bed", action="append",
                   help="NAME=PATH, repeatable: other regions")
    r.add_argument("--window", type=int, default=1_000_000,
                   help="for --bed rows without a gene: genes whose TSS lies "
                        "within this distance (default %(default)s)")
    r.add_argument("--out", required=True)

    n = sub.add_parser("annotate", help="annotate and rank the variants")
    n.add_argument("--vcf", required=True, action="append")
    n.add_argument("--regions", required=True, type=Path)
    n.add_argument("--ped")
    n.add_argument("--samples", default="",
                   help="comma-separated probands; default: the affected "
                        "individuals of --ped, or every sample")
    n.add_argument("--inheritance")
    n.add_argument("--af-keys",
                   default="gnomADg_AF,gnomAD_AF,gnomADe_AF,MAX_AF,AF_joint,"
                           "gnomad_AF",
                   help="INFO or CSQ keys of a population frequency "
                        "(default %(default)s)")
    n.add_argument("--gnomad", default="")
    n.add_argument("--gnomad-keys", default="AF_joint,AF",
                   help="INFO keys read from --gnomad (default %(default)s)")
    n.add_argument("--score", action="append",
                   help="NAME=PATH[:COLUMN], repeatable")
    n.add_argument("--info-score", action="append",
                   help="NAME=KEY, repeatable")
    n.add_argument("--min-score", action="append",
                   help="NAME=VALUE, repeatable: threshold of a score")
    n.add_argument("--require-score", action="store_true",
                   help="drop the regulatory small variants that reach no "
                        "--min-score threshold")
    n.add_argument("--max-af", type=float, default=0.01,
                   help="maximal frequency of a regulatory variant "
                        "(default %(default)s)")
    n.add_argument("--max-af-coding", type=float, default=0.01,
                   help="maximal frequency of a coding second hit "
                        "(default %(default)s)")
    n.add_argument("--coding-impacts", default="HIGH,MODERATE",
                   help="VEP/snpEff impacts that make a variant coding "
                        "(default %(default)s)")
    n.add_argument("--pass-only", action="store_true")
    n.add_argument("--min-gq", type=float, default=20)
    n.add_argument("--min-dp", type=float, default=8)
    n.add_argument("--out", required=True)
    a = p.parse_args()
    if a.cmd == "regions":
        cmd_regions(a)
    else:
        cmd_annotate(a)


if __name__ == "__main__":
    main()
