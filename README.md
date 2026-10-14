# mitoloop_v0

A pipeline that combines public resources to describe the regulatory
landscape of nuclear-encoded mitochondrial genes, with a focus on skeletal
muscle:

- **GTEx** independent cis-eQTLs, to compare the effect sizes and the
  distances to the transcription start site (TSS) of the eQTLs of a gene set
  with those of matched control genes;
- **ENCODE-rE2G** element-gene links, to count the regulatory elements
  attributed to each gene and measure how reproducibly they are attributed
  across biosamples;
- **ENCODE intact Hi-C** maps of human tissues, converted into virtual 4C
  profiles at the promoters of the gene set, to call promoter loops with
  peakC, keep those confirmed from both anchors, and relate their strength to
  the mitochondrial DNA copy number (mtCN) of each tissue.

Two modules extend it:

- **promoter features**: sequence (CpG content, promoter class, motifs of
  NRF1, ERR, GABPA and other factors), bidirectional promoters,
  promoter-promoter loops, and DNA methylation of the promoters and of the
  replicable elements from public whole-genome bisulfite data;
- **WGS annotation**: a region set built from the outputs above and the
  rare variants of unsolved genomes that fall in it, ranked by frequency,
  scores and genetic model (second coding hit in trans, homozygous,
  de novo, structural variants), to prepare the cohort analysis.

The gene set is an input (for instance the union of MitoCarta3.0 and of the
PanelApp "Mitochondrial disorders" panel), as are the control genes, which
`make_controls.py` can draw.

## Contents

| Script | Language | Role |
|---|---|---|
| `scripts/gtex_eqtl.py` | Python | download GTEx and GENCODE files; compare eQTL effect sizes and TSS distances between the gene set and matched controls |
| `scripts/make_controls.py` | Python | draw control genes matched on expression, cross-tissue expression variability, gene length and gene density |
| `scripts/re2g_panel.py` | Python | element-gene links per gene and their replicability across ENCODE-rE2G biosamples |
| `scripts/v4c_pipeline.py` | Python | virtual 4C extraction from `.hic`, `.cool` or `.mcool` libraries, library selection, export for R, pairwise comparison of terms |
| `scripts/v4c_loops.R` | R | peakC window choice, loop calling per ENCODE term, reciprocity filter, association of loop strength with mtCN (limma), one loop-strength value per GTEx tissue |
| `scripts/peaks_vs_elements.py` | Python | overlap between loop anchors and replicable rE2G elements, with permutation nulls |
| `scripts/locus_figure.py` | Python | one locus: rotated Hi-C map, loop arcs, gene track, virtual 4C profile with the called loops |
| `scripts/folding_figure.py` | Python | schematic chromatin folding of a locus drawn from the called loops |
| `scripts/svg_to_pptx.py` | Python | convert a matplotlib SVG into editable PowerPoint shapes |
| `scripts/promoter_features.py` | Python | promoter sequence, motifs, bidirectional pairs, promoter-promoter loops and DNA methylation (promoters and elements) against controls; methylation against GTEx expression |
| `scripts/vcf_annotate.py` | Python | `regions`: annotated BED of promoters, UTRs, coding exons, elements, loop anchors and other BEDs per gene; `annotate`: rare variants of WGS VCFs in these regions with frequency, scores, parental origin or phase, second hits and tiers |
| `scripts/v4c_pairs.py`, `scripts/v4c_pairs.R` | Python, R | optional module: loops at mitochondrial promoters in two metabolic states of one cell type |
| `config/` | | example library tables |
| `slurm/` | | job templates for SLURM clusters |
| `tests/` | | synthetic end-to-end tests |

## Installation

Two conda environments are used, one for Python and one for R.

```bash
git clone https://github.com/atrimbdx/mitoloop_v0.git
cd mitoloop_v0
conda env create -f environment.yml        # environment mitoloop_v0
conda env create -f environment-r.yml      # environment mitoloop_v0-r
conda activate mitoloop_v0-r
Rscript -e 'remotes::install_github("deWitLab/peakC")'
```

peakC is distributed on GitHub only. `hic-straw` is installed from PyPI by
`environment.yml`; it compiles against libcurl, which conda provides.

Python scripts run in `mitoloop_v0`, R scripts in `mitoloop_v0-r`.

## Input data

| Resource | Used by | Source |
|---|---|---|
| GENCODE annotation, GTF (the basic set downloaded by `gtex_eqtl.py fetch` is sufficient) | all | <https://www.gencodegenes.org/human/> |
| GTEx v8 independent cis-eQTLs and median TPM per tissue | `gtex_eqtl.py`, `make_controls.py` | GTEx portal, <https://gtexportal.org> (`gtex_eqtl.py fetch` downloads them) |
| Gene set | all | user supplied, one symbol per line; for instance MitoCarta3.0 (<https://www.broadinstitute.org/mitocarta>) and the PanelApp "Mitochondrial disorders" panel (<https://panelapp.genomicsengland.co.uk>) |
| ENCODE-rE2G element-gene predictions, one file per biosample | `re2g_panel.py` | ENCODE portal, <https://www.encodeproject.org> |
| ENCODE registry of candidate cis-regulatory elements (cCREs), BED | `re2g_panel.py` | SCREEN, <https://screen.wenglab.org> |
| ENCODE intact Hi-C maps of human tissues, `.hic` | `v4c_pipeline.py`, `locus_figure.py` | ENCODE portal |
| mtDNA copy number per GTEx tissue | `v4c_loops.R` (through the library table) | Rath et al. 2024, Datasets S1 to S4 |
| Genome sequence (GRCh38 FASTA) and CpG islands (`cpgIslandExt`) | `promoter_features.py` | UCSC, <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/> |
| Whole-genome bisulfite methylation (bedMethyl, Bismark coverage or bedGraph) of tissues and cell types | `promoter_features.py` (through `--meth-table`) | ENCODE portal (WGBS); Roadmap Epigenomics; Loyfer et al. 2023 atlas, GEO GSE186458 |
| WGS VCFs (small variants, structural variants) and PED | `vcf_annotate.py` | diagnostic cohort |
| Variant scores: CADD, ReMM, AlphaGenome or any per-position table | `vcf_annotate.py` | <https://cadd.gs.washington.edu/download>; ReMM (Smedley et al. 2016); AlphaGenome variant scoring exported to a table |
| gnomAD sites VCFs | `vcf_annotate.py` | <https://gnomad.broadinstitute.org/downloads> |
| Mode of inheritance per gene | `vcf_annotate.py` | PanelApp panel download, or a two-column table |

### Library table

`v4c_pipeline.py` and `locus_figure.py` read a tab-separated table with one
row per Hi-C library (`config/tissues.example.tsv` shows the format;
its accessions, locations and values are placeholders):

| Column | Content |
|---|---|
| `accession` | file accession, used to name the outputs |
| `term_name` | ENCODE biosample term; libraries sharing a term are replicates |
| `gtex_tissue` | GTEx tissue matching the term; the unit of the mtCN analysis |
| `mtcn` | mtDNA copy number of that GTEx tissue |
| `s3_url`, `download_url`, `path` | location of the file: an http(s) URL or a local path; tried in this order (`--location-columns`) |
| `format` | optional: `hic`, `cool` or `mcool`, otherwise guessed from the extension |

Column names can be changed with `--accession-column`, `--term-column`,
`--gtex-column` and `--mtcn-column`.

## Workflow

The commands below write into `results/` and `figures/`, read inputs from
`data/`, and assume a gene set in `genes.txt`. Every option is described by
`--help` (`Rscript scripts/v4c_loops.R --help` for the R script).

### 1. Annotation, GTEx files and control genes

```bash
mkdir -p results
python3 scripts/gtex_eqtl.py fetch --datadir data
python3 scripts/make_controls.py --gtf data/gencode.gtf.gz \
    --tpm data/gtex_median_tpm.gct.gz --genes genes.txt \
    --out controls.txt --report results/controls_check.tsv
```

### 2. GTEx eQTLs

```bash
python3 scripts/gtex_eqtl.py run --datadir data --genes genes.txt \
    --outdir results/eqtl
```

The tissue is set with `--tissue` (eQTL file name) and `--tpm-column` (TPM
table column).

### 3. ENCODE-rE2G links

```bash
# list the biosamples available
python3 scripts/re2g_panel.py --links-dir data/re2g --list-biosamples \
    --gtf data/gencode.gtf.gz --genes genes.txt --out results/re2g

python3 scripts/re2g_panel.py --links-dir data/re2g \
    --ccre data/cCREs.bed.gz --gtf data/gencode.gtf.gz \
    --genes genes.txt --controls controls.txt \
    --keep "skeletal muscle|psoas muscle|gastrocnemius" \
    --out results/re2g --cores 8
```

Links are binarised at the published operating point of the ENCODE-rE2G model
that produced the files, detected from their feature columns. Several models
share the same columns; `--model` or `--threshold` sets the operating point
explicitly. `re2g_elements.tsv` lists every element with its prevalence, the
number of biosamples in which it is linked to the gene.

### 4. Virtual 4C profiles

```bash
python3 scripts/v4c_pipeline.py extract --tissues tissues.tsv \
    --gtf data/gencode.gtf.gz --genes genes.txt --controls controls.txt \
    --outdir results/raw
python3 scripts/v4c_pipeline.py select --tissues tissues.tsv \
    --indir results/raw --out results/libraries.tsv
python3 scripts/v4c_pipeline.py export --indir results/raw \
    --libraries results/libraries.tsv --set mito,control \
    --outdir results/export/both
```

Each profile counts the contacts between the TSS anchor (TSS ± 10 kb by
default) and every 10 kb bin within ± 1 Mb, as raw observed counts. Remote
`.hic` files are read over their URL; only the profiles are stored.

`select` drops the libraries whose depth is below the median minus
`--k-sd` standard deviations of log10 depth, but always keeps the
`--min-keep` deepest; adjust `--min-keep` to the size of the library table.

### 5. Loops

```bash
conda activate mitoloop_v0-r
# choose the peakC window size
Rscript scripts/v4c_loops.R scan --libraries results/libraries.tsv \
    --export results/export/both --out results/loops
# call peaks per ENCODE term, with the window proposed by scan
Rscript scripts/v4c_loops.R call --libraries results/libraries.tsv \
    --export results/export/both --wsize 3 --out results/loops

conda activate mitoloop_v0
# profiles at the partner anchors
python3 scripts/v4c_pipeline.py extract --tissues tissues.tsv \
    --anchors results/loops/partners.tsv --outdir results/raw_partners
python3 scripts/v4c_pipeline.py export --indir results/raw_partners \
    --libraries results/libraries.tsv --set partner \
    --outdir results/export/partners

conda activate mitoloop_v0-r
# keep the peaks found again from the partner anchor
Rscript scripts/v4c_loops.R reciprocal --libraries results/libraries.tsv \
    --export results/export/partners --wsize 3 --out results/loops
```

`loops.tsv` holds one row per loop, with the ENCODE terms in which it was
called (`called`, separated by `|`). The same `--wsize` must be used in
`call`, `reciprocal`, `assoc` and `tissue`.

### 6. Loop strength and mtDNA copy number

```bash
Rscript scripts/v4c_loops.R assoc --libraries results/libraries.tsv \
    --export results/export/both --wsize 3 --out results/loops \
    --reference controls.txt --min-terms 3
Rscript scripts/v4c_loops.R tissue --libraries results/libraries.tsv \
    --export results/export/both --wsize 3 --out results/loops \
    --reference controls.txt --min-terms 3
```

`assoc` models the log contact of each loop against log2 mtCN with limma
(vooma weights, libraries of one GTEx tissue treated as correlated through
duplicateCorrelation), tests the mean effect by permuting mtCN between GTEx
tissues, and compares the loops of the gene set with those of the reference
genes. `tissue` writes one loop-strength value per GTEx tissue, for plotting
against mtCN.

### 7. Loop anchors and rE2G elements

```bash
python3 scripts/peaks_vs_elements.py \
    --elements results/re2g/re2g_elements.tsv \
    --loops results/loops/loops.tsv --list-terms
python3 scripts/peaks_vs_elements.py \
    --elements results/re2g/re2g_elements.tsv \
    --loops results/loops/loops.tsv --min-prevalence 51 \
    --min-distance 50000 --called-in "psoas muscle|gastrocnemius" \
    --out results/overlap
```

`--called-in` keeps the loops called in the tissues the elements come from.
`--min-prevalence` is a number of biosamples: set it relative to the number
kept by `re2g_panel.py --keep`.
Two nulls are reported: elements placed uniformly in the window, and element
distances to the TSS redrawn from the pooled distribution of the same gene
set. The second keeps the concentration of both elements and anchors near
the promoter and is the one to read.

### 8. Figures of one locus

```bash
python3 scripts/locus_figure.py --npz results/raw/ACCESSION.npz \
    --gene GENE1 --tissues tissues.tsv --term "gastrocnemius medialis" \
    --loops results/loops/loops.tsv --gtf data/gencode.gtf.gz \
    --arcs --genes-track --window 1000000 --max-distance 1000000 \
    --norm SCALE --out figures/locus_GENE1
python3 scripts/folding_figure.py --loops results/loops/loops.tsv \
    --gene GENE1 --window 1000000 \
    --elements results/re2g/re2g_elements.tsv --gtf data/gencode.gtf.gz \
    --called-in "psoas muscle|gastrocnemius" --out figures/folding_GENE1
python3 scripts/svg_to_pptx.py figures/locus_GENE1.svg \
    figures/locus_GENE1.pptx
```

The virtual 4C reads one row of the contact matrix per viewpoint. A loop
between two positions of which neither is a viewpoint is therefore absent
from `loops.tsv`, even when it is visible on the map. To include the other
genes of a locus as viewpoints, run steps 4 and 5 with a gene list holding
them.

### 9. Promoter features and DNA methylation

```bash
python3 scripts/promoter_features.py --gtf data/gencode.gtf.gz \
    --genes genes.txt --controls controls.txt \
    --fasta data/GRCh38.fa --cpg-islands data/cpgIslandExt.txt.gz \
    --loops results/loops/loops.tsv \
    --meth-table config/methylation.tsv \
    --elements results/re2g/re2g_elements.tsv --min-prevalence 51 \
    --tpm data/gtex_median_tpm.gct.gz --focus "muscle" \
    --out results/promoters
```

`config/methylation.example.tsv` shows the methylation table (one row per
file: sample, group, path, optional format and GTEx tissue; its paths are
placeholders). Every feature is compared between the gene set and the
controls, by a Fisher or Mann-Whitney test and by a permutation of the set
labels within promoter classes (HCP, ICP, LCP), which separates a motif
enrichment from the CpG richness of the promoters. `--controls-all` takes
every other protein-coding gene as control. Only the CpGs inside the
promoters and elements are kept from each methylation file.

### 10. Annotation of WGS variants

```bash
# regions attributed to each gene of the set
python3 scripts/vcf_annotate.py regions --gtf data/gencode.gtf.gz \
    --genes genes.txt --elements results/re2g/re2g_elements.tsv \
    --min-prevalence 51 --loops results/loops/loops.tsv \
    --called-in "psoas muscle|gastrocnemius" \
    --bed chic=results/chic/other_ends.bed --out results/wgs/regions.bed

# normalised, multi-allelic records split
bcftools norm -m -any -f data/GRCh38.fa cohort.vcf.gz -Oz -o cohort.norm.vcf.gz

python3 scripts/vcf_annotate.py annotate --vcf cohort.norm.vcf.gz \
    --vcf cohort.sv.vcf.gz --ped cohort.ped \
    --regions results/wgs/regions.bed --inheritance panelapp_mito.tsv \
    --gnomad "data/gnomad/gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz" \
    --score CADD=data/cadd/whole_genome_SNVs.tsv.gz:PHRED \
    --score ReMM=data/ReMM.v0.4.hg38.tsv.gz \
    --min-score CADD=15 --min-score ReMM=0.9 --pass-only \
    --out results/wgs/annot
```

`--bed` adds any other region set with the gene in its fourth column
(capture Hi-C or Micro-C interactions, elements validated by perturbation);
rows without a gene go to the genes of the set within `--window`. Score
tables indexed with tabix are looked up by position; plain ones are read
once. Frequencies and scores already annotated by VEP (`CSQ`) are read with
`--af-keys` and `--info-score`. Each regulatory variant receives a tier
(1: homozygous, hemizygous, in trans of a rare coding variant or of
another regulatory variant, de novo in a dominant gene; 2: partner of
unknown phase, de novo, dominant gene; 3: cis, single heterozygous); the
method is described in the header of the script.

### 11. Optional: two metabolic states of one cell type

```bash
cp config/pairs.example.tsv config/pairs.tsv     # then edit the paths
python3 scripts/v4c_pairs.py extract --pairs config/pairs.tsv \
    --gtf hg38=data/gencode.v44.gtf.gz --gtf hg19=data/gencode.v19.gtf.gz \
    --gtf mm10=data/gencode.vM25.gtf.gz \
    --genes genes.txt --controls controls.txt --outdir results/raw_pairs
python3 scripts/v4c_pairs.py export --pairs config/pairs.tsv \
    --indir results/raw_pairs --outdir results/export_pairs
Rscript scripts/v4c_pairs.R compare \
    --libraries results/export_pairs/libraries_pairs.tsv \
    --export results/export_pairs --wsize 3 --out results/pairs_out
```

`config/pairs.example.tsv` documents the table and lists candidate public
datasets. One `--gtf` is needed per genome build named in the table.

## Running on a SLURM cluster

`slurm/extract.sbatch` runs one library per array task;
`slurm/run_r.sbatch` runs one `v4c_loops.R` subcommand. Both read three
environment variables: `CONDA_PROFILE` (path to `etc/profile.d/conda.sh`),
`MITOLOOP_DIR` (working directory), and optionally `CONDA_ENV` /
`CONDA_ENV_R`. Partition and account directives are site specific and are
to be added to the templates.

```bash
export CONDA_PROFILE=path/to/conda/etc/profile.d/conda.sh
export MITOLOOP_DIR=$PWD
mkdir -p logs
N=$(( $(wc -l < tissues.tsv) - 1 ))
sbatch --array=1-${N}%10 slurm/extract.sbatch --gtf data/gencode.gtf.gz \
    --genes genes.txt --controls controls.txt --outdir results/raw
sbatch slurm/run_r.sbatch call --libraries results/libraries.tsv \
    --export results/export/both --wsize 3 --out results/loops
```

## Tests

```bash
bash tests/run_all.sh
```

Each test builds small synthetic inputs with planted signal, runs the
scripts end to end and checks their outputs. The `.hic` route is tested when
`hicstraw`, java and Juicer Tools are available (`JUICER_TOOLS` pointing to
the jar), and skipped otherwise. The `deseq` subcommand of `v4c_loops.R`
needs DESeq2 and ashr, which the tests do not require.

## Citation

See `CITATION.cff`. Please also cite the data resources and methods listed
below.

## References

Data resources and annotation

- ENCODE Project Consortium, Moore JE, Purcaro MJ, et al. 2020. Expanded encyclopaedias of DNA elements in the human and mouse genomes. *Nature* 583:699-710. doi:10.1038/s41586-020-2493-4
- Frankish A, Carbonell-Sala S, Diekhans M, et al. 2023. GENCODE: reference annotation for the human and mouse genomes in 2023. *Nucleic Acids Research* 51:D942-D949. doi:10.1093/nar/gkac1071
- Gschwind AR, Mualim KS, Karbalayghareh A, et al. 2023. An encyclopedia of enhancer-gene regulatory interactions in the human genome. *bioRxiv*. doi:10.1101/2023.11.09.563812
- GTEx Consortium. 2020. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science* 369:1318-1330. doi:10.1126/science.aaz1776
- Martin AR, Williams E, Foulger RE, et al. 2019. PanelApp crowdsources expert knowledge to establish consensus diagnostic gene panels. *Nature Genetics* 51:1560-1565. doi:10.1038/s41588-019-0528-2
- Morales J, Pujar S, Loveland JE, et al. 2022. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research. *Nature* 604:310-315. doi:10.1038/s41586-022-04558-8
- Rath S, Sharma R, Gupta R, et al. 2021. MitoCarta3.0: an updated mitochondrial proteome now with sub-organelle localization and pathway annotations. *Nucleic Acids Research* 49:D1541-D1547. doi:10.1093/nar/gkaa1011
- Rath SP, Gupta R, Todres E, et al. 2024. Mitochondrial genome copy number variation across tissues in mice and humans. *Proceedings of the National Academy of Sciences USA* 121:e2402291121. doi:10.1073/pnas.2402291121
- Seal RL, Braschi B, Gray K, et al. 2023. Genenames.org: the HGNC resources in 2023. *Nucleic Acids Research* 51:D1003-D1009. doi:10.1093/nar/gkac888

- Chen S, Francioli LC, Goodrich JK, et al. 2024. A genomic mutational constraint map using variation in 76,156 human genomes (gnomAD). *Nature* 625:92-100. doi:10.1038/s41586-023-06045-0
- Loyfer N, Magenheim J, Peretz A, et al. 2023. A DNA methylation atlas of normal human cell types. *Nature* 613:355-364. doi:10.1038/s41586-022-05580-6
- Roadmap Epigenomics Consortium. 2015. Integrative analysis of 111 reference human epigenomes. *Nature* 518:317-330. doi:10.1038/nature14248

Methods

- Gardiner-Garden M, Frommer M. 1987. CpG islands in vertebrate genomes. *Journal of Molecular Biology* 196:261-282. doi:10.1016/0022-2836(87)90689-9
- Rentzsch P, Witten D, Cooper GM, et al. 2019. CADD: predicting the deleteriousness of variants throughout the human genome. *Nucleic Acids Research* 47:D886-D894. doi:10.1093/nar/gky1016
- Smedley D, Schubach M, Jacobsen JOB, et al. 2016. A whole-genome analysis framework for effective identification of pathogenic regulatory variants in Mendelian disease. *American Journal of Human Genetics* 99:595-606. doi:10.1016/j.ajhg.2016.07.005
- Trinklein ND, Aldred SF, Hartman SJ, et al. 2004. An abundance of bidirectional promoters in the human genome. *Genome Research* 14:62-66. doi:10.1101/gr.1982804
- Weber M, Hellmann I, Stadler MB, et al. 2007. Distribution, silencing potential and evolutionary impact of promoter DNA methylation in the human genome. *Nature Genetics* 39:457-466. doi:10.1038/ng1990
- Bianchi V, Geeven G, Tucker N, et al. 2019. Detailed regulatory interaction map of the human heart facilitates gene discovery for cardiovascular disease (peakHiC). *bioRxiv*. doi:10.1101/705715
- Breitling R, Armengaud P, Amtmann A, Herzyk P. 2004. Rank products: a simple, yet powerful, new method to detect differentially regulated genes in replicated microarray experiments. *FEBS Letters* 573:83-92. doi:10.1016/j.febslet.2004.07.055
- Cohen J. 1988. *Statistical Power Analysis for the Behavioral Sciences*, 2nd edition. Lawrence Erlbaum Associates.
- de Leeuw J, Hornik K, Mair P. 2009. Isotone optimization in R: pool-adjacent-violators algorithm (PAVA) and active set methods. *Journal of Statistical Software* 32(5). doi:10.18637/jss.v032.i05
- Fulco CP, Nasser J, Jones TR, et al. 2019. Activity-by-contact model of enhancer-promoter regulation from thousands of CRISPR perturbations. *Nature Genetics* 51:1664-1669. doi:10.1038/s41588-019-0538-0
- Geeven G, Teunissen H, de Laat W, de Wit E. 2018. peakC: a flexible, non-parametric peak calling package for 4C and Capture-C data. *Nucleic Acids Research* 46:e91. doi:10.1093/nar/gky443
- Inman HF, Bradley EL. 1989. The overlapping coefficient as a measure of agreement between probability distributions and point estimation of the overlap of two normal densities. *Communications in Statistics - Theory and Methods* 18:3851-3874. doi:10.1080/03610928908830127
- Jaccard P. 1912. The distribution of the flora in the alpine zone. *New Phytologist* 11:37-50. doi:10.1111/j.1469-8137.1912.tb05611.x
- Law CW, Chen Y, Shi W, Smyth GK. 2014. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology* 15:R29. doi:10.1186/gb-2014-15-2-r29
- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550. doi:10.1186/s13059-014-0550-8
- Mann HB, Whitney DR. 1947. On a test of whether one of two random variables is stochastically larger than the other. *Annals of Mathematical Statistics* 18:50-60. doi:10.1214/aoms/1177730491
- McGraw KO, Wong SP. 1992. A common language effect size statistic. *Psychological Bulletin* 111:361-365. doi:10.1037/0033-2909.111.2.361
- Phipson B, Smyth GK. 2010. Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. *Statistical Applications in Genetics and Molecular Biology* 9:39. doi:10.2202/1544-6115.1585
- Ritchie ME, Phipson B, Wu D, et al. 2015. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research* 43:e47. doi:10.1093/nar/gkv007
- Sen PK. 1968. Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association* 63:1379-1389. doi:10.1080/01621459.1968.10480934
- Smyth GK, Michaud J, Scott HS. 2005. Use of within-array replicate spots for assessing differential expression in microarray experiments. *Bioinformatics* 21:2067-2075. doi:10.1093/bioinformatics/bti270
- Stephens M. 2017. False discovery rates: a new deal. *Biostatistics* 18:275-294. doi:10.1093/biostatistics/kxw041

Software

- Danecek P, Bonfield JK, Liddle J, et al. 2021. Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008. doi:10.1093/gigascience/giab008
- McLaren W, Gil L, Hunt SE, et al. 2016. The Ensembl Variant Effect Predictor. *Genome Biology* 17:122. doi:10.1186/s13059-016-0974-4
- Abdennur N, Mirny LA. 2020. Cooler: scalable storage for Hi-C data and other genomically labeled arrays. *Bioinformatics* 36:311-316. doi:10.1093/bioinformatics/btz540
- Durand NC, Shamim MS, Machol I, et al. 2016. Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments. *Cell Systems* 3:95-98. doi:10.1016/j.cels.2016.07.002
- Hunter JD. 2007. Matplotlib: a 2D graphics environment. *Computing in Science and Engineering* 9:90-95. doi:10.1109/MCSE.2007.55
- Virtanen P, Gommers R, Oliphant TE, et al. 2020. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods* 17:261-272. doi:10.1038/s41592-019-0686-2
- ENCODE-rE2G models and thresholds: <https://github.com/EngreitzLab/ENCODE_rE2G>
- caTools, R package used by peakC: <https://CRAN.R-project.org/package=caTools>
- peakC: <https://github.com/deWitLab/peakC>
- python-pptx: <https://github.com/scanny/python-pptx>

Public datasets of the optional paired-state module are cited in
`config/pairs.example.tsv`.

## License

MIT, see `LICENSE`.
