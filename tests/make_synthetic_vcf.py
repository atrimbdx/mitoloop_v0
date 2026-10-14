#!/usr/bin/env python3
"""
make_synthetic_vcf.py - synthetic inputs for vcf_annotate.py

Purpose
    Build a GTF, element and loop tables, a capture BED, a PED of two trios
    and one singleton, small-variant and structural-variant VCFs, score
    tables, a gnomAD sites VCF and a PanelApp-like inheritance table, with
    one planted case per genetic model that vcf_annotate.py ranks.

Inputs
    none

Outputs (in --out)
    gencode.gtf.gz, genes.txt, re2g_elements.tsv, loops.tsv, chic.bed,
    cohort.ped, snv.vcf.gz, sv.vcf, cadd.tsv.gz (plain), remm.tsv.gz
    (bgzipped and tabix-indexed when pysam is installed, plain otherwise),
    gnomad.vcf, panelapp.tsv, truth.tsv (proband, gene, pos, tier, category)

Method
    Genes: GENEA (+, AR) at chr1:1 Mb, GENEB (-, AR) at chr1:3 Mb, GENEC (+,
    AD) at chr1:5 Mb, TAFAZZIN (+, X-linked) at chrX:50 Mb, listed in
    genes.txt under its former symbol TAZ. Each has three exons, coding from
    TSS + 200 to TSS + 10,300. Planted cases:
      P1 (male, trio)  GENEA element variant from the mother, missense from
                       the father: trans_with_coding; TAFAZZIN element
                       variant hemizygous; GENEC element variant with GQ 5
                       (removed)
      P2 (female, trio) GENEA loop-anchor variant homozygous; GENEB promoter
                       variant and stop-gain both from the father:
                       cis_with_coding; GENEC element variant de novo in a
                       dominant gene; TAFAZZIN element variant from the
                       mother: x_linked_heterozygous; de novo inversion with
                       one breakpoint between the GENEA element and its TSS
      P3 (female, alone) GENEB capture-region variant and coding variant on
                       opposite haplotypes of one PS block:
                       trans_with_coding; second allele of a multi-allelic
                       GENEA promoter record (the first is common);
                       recurrent GENEA element variant (also in P1); GENEC
                       element deleted by a DEL: dominant_gene
    Removed on purpose: a common promoter variant (INFO frequency), one
    common in the gnomAD VCF, one with FILTER LowQual (--pass-only), one at
    a loop anchor called only in liver, one in an element of low
    prevalence.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

GENES = [  # name, chrom, tss, strand, panelapp mode
    ("GENEA", "chr1", 1_000_000, "+",
     "BIALLELIC, autosomal or pseudoautosomal"),
    ("GENEB", "chr1", 3_000_000, "-",
     "BIALLELIC, autosomal or pseudoautosomal"),
    ("GENEC", "chr1", 5_000_000, "+",
     "MONOALLELIC, autosomal or pseudoautosomal, imprinted status unknown"),
    ("TAFAZZIN", "chrX", 50_000_000, "+",
     "X-LINKED: hemizygous mutation in males, biallelic mutations in "
     "females"),
]
EXONS = [(0, 500), (5000, 5300), (10000, 11000)]
CDS = [(200, 500), (5000, 5300), (10000, 10300)]
SAMPLES = ["P1", "F1", "M1", "P2", "F2", "M2", "P3"]


def iv(tss, strand, a, b):
    return (tss + a, tss + b) if strand == "+" else (tss - b, tss - a)


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    with gzip.open(out / "gencode.gtf.gz", "wt") as fh:
        for n, (g, c, tss, st, _m) in enumerate(GENES):
            attr = (f'gene_id "ENSG{n:011d}.1"; gene_type "protein_coding"; '
                    f'gene_name "{g}";')
            s, e = iv(tss, st, 0, 11000)
            fh.write(f"{c}\tSYN\tgene\t{s + 1}\t{e}\t.\t{st}\t.\t{attr}\n")
            for tid, tag in ((f"ENST{n:011d}", ' tag "MANE_Select";'),
                             (f"ENST{n:011d}B", "")):
                # a shorter, non-MANE transcript that must be ignored
                exons = EXONS if tag else EXONS[1:]
                ts = iv(tss, st, exons[0][0], exons[-1][1])
                ta = f'{attr} transcript_id "{tid}";{tag}'
                fh.write(f"{c}\tSYN\ttranscript\t{ts[0] + 1}\t{ts[1]}\t.\t{st}"
                         f"\t.\t{ta}\n")
                for x in exons:
                    s, e = iv(tss, st, *x)
                    fh.write(f"{c}\tSYN\texon\t{s + 1}\t{e}\t.\t{st}\t.\t{ta}\n")
                for x in CDS[3 - len(exons):]:
                    s, e = iv(tss, st, *x)
                    fh.write(f"{c}\tSYN\tCDS\t{s + 1}\t{e}\t.\t{st}\t0\t{ta}\n")
    (out / "genes.txt").write_text("GENEA\nGENEB\nGENEC\nTAZ\n")
    with open(out / "panelapp.tsv", "w") as fh:
        fh.write("Entity Name\tModel_Of_Inheritance\tPhenotypes\n")
        for g, _c, _t, _s, m in GENES:
            fh.write(f"{g}\t{m}\tsynthetic\n")

    with open(out / "re2g_elements.tsv", "w") as fh:
        fh.write("gene\tset\tchrom\tstart\tend\tdistance\tprevalence\n")
        for g, c, tss, off, prev in (("GENEA", "chr1", 1_000_000, 50_000, 6),
                                     ("GENEA", "chr1", 1_000_000, 80_000, 1),
                                     ("GENEC", "chr1", 5_000_000, 40_000, 6),
                                     ("TAZ", "chrX", 50_000_000, 30_000, 5)):
            m = tss + off
            fh.write(f"{g}\tmito\t{c}\t{m - 200}\t{m + 200}\t{off}\t{prev}\n")
    with open(out / "loops.tsv", "w") as fh:
        fh.write("gene\tchrom\tvp_pos\tsummit\tstart\tend\tn_terms\tcalled\t"
                 "loop_id\n")
        fh.write("GENEA\tchr1\t1000000\t800000\t795000\t805000\t2\t"
                 "psoas muscle|liver\tGENEA:800000\n")
        fh.write("GENEA\tchr1\t1000000\t700000\t695000\t705000\t1\t"
                 "liver\tGENEA:700000\n")
    with open(out / "chic.bed", "w") as fh:
        fh.write("chr1\t3119500\t3120500\tGENEB\tscore=7.1\n")
        fh.write("chr1\t5200000\t5201000\n")

    with open(out / "cohort.ped", "w") as fh:
        fh.write("FAM1\tP1\tF1\tM1\t1\t2\nFAM1\tF1\t0\t0\t1\t1\n"
                 "FAM1\tM1\t0\t0\t2\t1\nFAM2\tP2\tF2\tM2\t2\t2\n"
                 "FAM2\tF2\t0\t0\t1\t1\nFAM2\tM2\t0\t0\t2\t1\n"
                 "FAM3\tP3\t0\t0\t2\t2\n")

    ref = "0/0:30:40"
    # chrom, pos, ref, alt, filter, info, {sample: GT:GQ:DP[:PS]}
    recs = [
        ("chr1", 700001, "A", "G", "PASS", ".", {"P2": "0/1:50:30"}),
        ("chr1", 800001, "C", "T", "PASS", "gnomADg_AF=0.002",
         {"P2": "1/1:50:30", "F2": "0/1:50:30", "M2": "0/1:50:30"}),
        ("chr1", 999701, "G", "A", "PASS", "gnomADg_AF=0.2",
         {"P1": "0/1:50:30", "M1": "0/1:50:30"}),
        ("chr1", 999801, "C", "T,G", "PASS", "gnomADg_AF=0.5,0.0001",
         {"P3": "0/2:50:30"}),
        ("chr1", 1005101, "C", "A", "PASS",
         "CSQ=A|missense_variant|MODERATE|GENEA|0.00002",
         {"P1": "0/1:50:30", "F1": "0/1:50:30"}),
        ("chr1", 1050001, "C", "T", "PASS",
         "gnomADg_AF=0.0001;CSQ=T|upstream_gene_variant|MODIFIER|GENEA|",
         {"P1": "0/1:50:30", "M1": "0/1:50:30", "P3": "0/1:50:30"}),
        ("chr1", 1080001, "G", "T", "PASS", ".", {"P1": "0/1:50:30"}),
        ("chr1", 2989851, "G", "C", "PASS", ".", {"P3": "1|0:50:30:100"}),
        ("chr1", 2994851, "A", "G", "PASS",
         "CSQ=G|stop_gained|HIGH|GENEB|",
         {"P2": "0/1:50:30", "F2": "0/1:50:30"}),
        ("chr1", 3000301, "T", "C", "PASS", ".",
         {"P2": "0/1:50:30", "F2": "0/1:50:30"}),
        ("chr1", 3120001, "G", "A", "PASS", ".", {"P3": "0|1:50:30:100"}),
        ("chr1", 4999901, "C", "G", "PASS", ".", {"P1": "0/1:50:30"}),
        ("chr1", 4999951, "A", "C", "LowQual", ".", {"P3": "0/1:50:30"}),
        ("chr1", 5040001, "T", "A", "PASS", "AG=0.8", {"P2": "0/1:50:30"}),
        ("chr1", 5040101, "G", "A", "PASS", ".", {"P1": "0/1:5:30"}),
        ("chrX", 50030001, "A", "T", "PASS", ".",
         {"P1": "1/1:50:30", "M1": "0/1:50:30", "P2": "0/1:50:30",
          "M2": "0/1:50:30"}),
    ]
    with gzip.open(out / "snv.vcf.gz", "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="'
                 'Consequence annotations from Ensembl VEP. Format: '
                 'Allele|Consequence|IMPACT|SYMBOL|gnomADg_AF">\n')
        fh.write('##INFO=<ID=gnomADg_AF,Number=A,Type=Float,Description="">\n')
        fh.write('##INFO=<ID=AG,Number=A,Type=Float,Description="">\n')
        fh.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 + "\t".join(SAMPLES) + "\n")
        for c, pos, r, alt, flt, info, gts in recs:
            cols = []
            for s in SAMPLES:
                g = gts.get(s, ref)
                cols.append(g if g.count(":") == 3 else g + ":.")
            fh.write(f"{c}\t{pos}\tv{pos}\t{r}\t{alt}\t50\t{flt}\t{info}\t"
                     "GT:GQ:DP:PS\t" + "\t".join(cols) + "\n")

    with open(out / "sv.vcf", "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 + "\t".join(SAMPLES) + "\n")
        for c, pos, sid, alt, info, gts in (
                ("chr1", 1020001, "inv1", "<INV>", "SVTYPE=INV;END=1070000",
                 {"P2": "0/1"}),
                ("chr1", 5039001, "del1", "<DEL>",
                 "SVTYPE=DEL;END=5041000;SVLEN=-2000", {"P3": "0/1"}),
                ("chr1", 9000001, "del2", "<DEL>", "SVTYPE=DEL;END=9001000",
                 {"P3": "0/1"})):
            cols = [gts.get(s, "0/0") for s in SAMPLES]
            fh.write(f"{c}\t{pos}\t{sid}\tN\t{alt}\t50\tPASS\t{info}\tGT\t"
                     + "\t".join(cols) + "\n")

    with gzip.open(out / "cadd.tsv.gz", "wt") as fh:
        fh.write("## CADD GRCh38-v1.7 (c) University of Washington\n")
        fh.write("#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n")
        for pos, r, alt, ph in ((800001, "C", "T", 22.0),
                                (1050001, "C", "T", 18.0),
                                (1050001, "C", "G", 30.0),
                                (3120001, "G", "A", 12.0),
                                (5040001, "T", "A", 25.0)):
            fh.write(f"1\t{pos}\t{r}\t{alt}\t1.5\t{ph}\n")
    remm = out / "remm.tsv"
    with open(remm, "w") as fh:
        fh.write("#chr\tpos\tscore\n")
        for pos, v in ((800001, 0.5), (1050001, 0.95), (3120001, 0.92)):
            fh.write(f"chr1\t{pos}\t{v}\n")
    try:
        import pysam
        pysam.tabix_index(str(remm), seq_col=0, start_col=1, end_col=1,
                          force=True, meta_char="#")
    except ImportError:
        with open(remm) as src, gzip.open(str(remm) + ".gz", "wt") as dst:
            dst.write(src.read())
        remm.unlink()

    with open(out / "gnomad.vcf", "w") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                 "FILTER\tINFO\n")
        fh.write("chr1\t3120001\t.\tG\tA\t.\tPASS\tAF=0.0002;AF_joint=0.0001\n")
        fh.write("chr1\t4999901\t.\tC\tG\t.\tPASS\tAF=0.05\n")

    with open(out / "truth.tsv", "w") as fh:
        fh.write("proband\tgene\tpos\ttier\tcategory\n")
        for row in (("P1", "GENEA", 1050001, 1, "trans_with_coding"),
                    ("P1", "TAFAZZIN", 50030001, 1, "hemizygous"),
                    ("P2", "GENEA", 800001, 1, "homozygous"),
                    ("P2", "GENEB", 3000301, 3, "cis_with_coding"),
                    ("P2", "GENEC", 5040001, 1, "de_novo_dominant"),
                    ("P2", "GENEA", 1020001, 2,
                     "possible_compound_regulatory"),
                    ("P2", "TAFAZZIN", 50030001, 3, "x_linked_heterozygous"),
                    ("P3", "GENEB", 3120001, 1, "trans_with_coding"),
                    ("P3", "GENEA", 999801, 2,
                     "possible_compound_regulatory"),
                    ("P3", "GENEA", 1050001, 2,
                     "possible_compound_regulatory"),
                    ("P3", "GENEC", 5039001, 2, "dominant_gene")):
            fh.write("\t".join(str(x) for x in row) + "\n")


if __name__ == "__main__":
    main()
