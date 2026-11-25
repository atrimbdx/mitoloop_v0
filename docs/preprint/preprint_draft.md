# Preprint draft (bioRxiv) — first full draft

*Working title (French, from the authors): « prédiction de l'implication de régions non codantes dans les maladies mitochondriales ».*
*Draft of 25 November 2026. Every number comes from the AG2 v2 campaign handover note (reannot `docs/ag2/ETAT_CAMPAGNE_AG2_V2.md`), the code of the two repositories, or the facts list supplied with the brief. Anything not yet computed or checked is marked [TODO]; any missing reference is marked [REF NEEDED] or [DOI TO CHECK].*

## Title options

1. **Mapping the non-coding regulatory regions of nuclear mitochondrial disease genes with tissue Hi-C loops, ENCODE-rE2G enhancers and AlphaGenome in silico mutagenesis**
2. **In silico deletion scanning and saturation mutagenesis of 1,220 nuclear mitochondrial genes in four tissues point to proximal promoters as the main predictable non-coding target**
3. **Which non-coding regions matter for mitochondrial disease genes? A sequence-model campaign across promoters, enhancers and chromatin loop anchors in muscle, heart, brain and kidney**

**Chosen title: option 1.** It names what was done and the three data layers, without claiming a diagnostic yield that has not been measured. Option 2 states the main result but leans on a model that is known to underestimate distal effects (see Discussion), so the claim fits the abstract better than the title.

**Authors:** [authors]

**Affiliations:** [authors]

**Correspondence:** [authors]

---

## Abstract

Many patients with suspected mitochondrial disease remain undiagnosed after genome sequencing, and non-coding variants of nuclear mitochondrial genes are rarely interpreted. We predicted which non-coding bases and regions of these genes matter for their expression in four affected tissues. We assembled 1,220 genes (MitoCarta3.0 and the PanelApp "Mitochondrial disorders" panel), called promoter-anchored chromatin loops per tissue from pooled ENCODE intact Hi-C libraries (virtual 4C, peakC, reciprocal confirmation; 653 loops in skeletal muscle, 2,122 in heart, 1,494 in brain, 612 in kidney), and kept ENCODE-rE2G enhancers linked to the gene in at least 80 % of a tissue's biosamples. We then ran an in silico campaign with AlphaGenome (API version 2), on the reference sequence only: 58,089 10-bp deletion windows over 1,183 MANE promoters, 173 tissue-dominant promoters and 707 enhancer elements, multi-scale deletions, and base-by-base mutagenesis of the 8,272 most sensitive windows (about 248,000 alleles). Predicted effects concentrated in proximal promoters: about 28 % of MANE promoter alleles lowered the host RNA by at least 0.2 log2, critical positions were enriched from −125 to +25 bp and in YY1, ETS/GABPA, SP1, NRF1 and NF-Y sites. Single-nucleotide effects in enhancers were weak, and beyond 100 kb the model predicted contact loss without RNA change. Bases with predicted effects showed no clear depletion in gnomAD, consistent with recessive inheritance. These predictions, not yet validated experimentally, define a prioritisation filter of 14,249 promoter and element positions for unsolved patients, to be interpreted with a second allele and functional evidence.

*(Abstract: about 240 words.)*

---

## Introduction

Mitochondrial diseases are a clinically and genetically heterogeneous group of disorders. Most are caused by variants in nuclear genes that encode mitochondrial proteins, and diagnostic gene panels such as the PanelApp "Mitochondrial disorders" panel (Martin et al. 2019) now hold several hundred genes. Despite exome and genome sequencing, a substantial share of patients remains without a molecular diagnosis [REF NEEDED: diagnostic yield of WGS in mitochondrial disease]. RNA sequencing of patient fibroblasts has shown that some of these patients carry variants that reduce the expression of a known disease gene without altering its coding sequence (Kopajtich et al. 2021; Yépez et al. [year and DOI TO CHECK]). Such expression outliers point to the non-coding sequence that controls these genes, but they are found in one accessible cell type, while the disease often affects muscle, heart, brain or kidney.

Interpreting a non-coding variant requires two things: knowing which regions act on the gene in the relevant tissue, and knowing which bases inside those regions matter. For the first, the ENCODE-rE2G model links candidate enhancers to their target genes across hundreds of biosamples (Gschwind et al. 2026), building on the activity-by-contact model validated by CRISPR perturbations (Fulco et al. 2019), and tissue Hi-C maps identify the chromatin loops that bring distal elements close to promoters (Rao et al. 2014 [DOI TO CHECK]). For the second, sequence-to-function models now predict the effect of any sequence change on RNA, transcription initiation, accessibility and chromatin contacts from the DNA sequence alone. AlphaGenome (Avsec et al. 2026) predicts these outputs at base resolution over a megabase of context. Earlier models were shown to underestimate the contribution of distal enhancers to gene expression (Karollus et al. 2023), which is a known limit for this kind of analysis.

Population data offer a third, independent line of evidence: regions under purifying selection carry fewer variants in gnomAD than expected from their mutability (Samocha et al. 2014; Lek et al. 2016; Karczewski et al. 2020 [DOIs TO CHECK]; Chen et al. 2024), and base-level mutation-rate models such as Roulette (Seplyarskiy et al. 2023 [DOI TO CHECK]) or the non-coding constraint map of di Iulio et al. (2018 [DOI TO CHECK]) refine this expectation. Whether such constraint can be detected for the regulatory regions of mostly recessive disease genes is not known.

Here we combine these resources for the nuclear-encoded mitochondrial genes. We call tissue-specific promoter loops from ENCODE intact Hi-C, select reproducible ENCODE-rE2G enhancers per tissue, and run a systematic in silico campaign with AlphaGenome (API version 2, hereafter AG2) in four cell types: deletion scanning of promoters and elements, multi-scale deletions of elements and loop anchors, and saturation mutagenesis of the most sensitive windows. We then ask whether the bases predicted to matter are constrained in gnomAD, and derive a list of positions to prioritise in unsolved patients. Throughout, AG2 outputs are used to describe and rank regions; they are never used to train a model or to tune a threshold. We report predictions only: no experimental validation has been done yet, and we discuss where the model is expected to fail.

---

## Results

### 1. Gene set and tissues

The gene set is the union of MitoCarta3.0 (Rath et al. 2021) and of all confidence levels of the PanelApp "Mitochondrial disorders" panel (Martin et al. 2019), with symbols updated to GENCODE: 1,220 genes. Of these, 1,183 had a MANE Select promoter usable in the campaign [TODO: list the reasons for the 37 genes left out, e.g. no MANE transcript, mtDNA-encoded genes, non-gene panel entries].

We studied four cell types or tissues chosen for their relevance to mitochondrial disease and for the availability of matched data: skeletal muscle, heart (cardiomyocyte), brain (neurone) and kidney. Each was matched across three data sources: the ENCODE Hi-C libraries used to call loops, the ENCODE-rE2G biosamples used to select enhancers, and the AG2 output condition (gastrocnemius medialis, UBERON:0011907; heart left ventricle, UBERON:0002084; dorsolateral prefrontal cortex, UBERON:0009834; kidney, UBERON:0001225) (Supplementary Table S1).

### 2. Tissue-specific promoter loops from Hi-C

For each gene promoter we extracted a virtual 4C profile from every ENCODE intact Hi-C library (contacts of the TSS anchor, TSS ± 10 kb, with every 10-kb bin within ± 1 Mb), pooled the libraries of each tissue as replicates and called peaks with peakC (Geeven et al. 2018) in combined analysis with a window of three bins. A peak was kept as a loop only if a peak at the mitochondrial TSS was found again in the virtual 4C profile of the partner anchor (reciprocity). Loops were called separately per tissue, so that coordinates never merge across tissues. This gave 653 loops in skeletal muscle, 2,122 in heart, 1,494 in brain and 612 in kidney (Fig. 2a) [TODO: number of genes with at least one loop per tissue; distribution of loop lengths; number of libraries per tissue; share of loops also found in single libraries (`n_presents`)].

As a control, each library was also called on its own (`v4c_loops.R call_lib`), and the number of a tissue's libraries in which each pooled loop is found alone is reported without being used as a filter [TODO: summary of `n_presents`].

### 3. Selection of reproducible enhancers per tissue

ENCODE-rE2G predictions were binarised at the published operating point of the model that produced each file and mapped to the ENCODE cCRE registry so that biosamples share one coordinate system. For each tissue we grouped the biosamples matching that tissue (groups of at least 10 biosamples: 64 for muscle, 48 for heart, 60 for brain, 33 for kidney) and kept an element when it was linked to the gene in at least 80 % of the group's biosamples. This stringent rule keeps elements that are reproducibly attributed to the gene in the tissue, at the cost of missing elements active in a subset of cell states.

Elements were then classed by their position relative to the loops of the same gene in the same tissue. Under the initial anchor rule (loop summit ± 5 kb), 349 elements lay on a loop anchor. For disease genes, we added one off-anchor element per gene and per cell type (no anchor overlap, no GENCODE TSS within 1 kb, the best-linked then the closest element): 358 elements. The campaign therefore covered 707 elements (Fig. 1). The anchor rule was later revised (Section 6).

### 4. Design and calibration of the AlphaGenome campaign

All AG2 calls used the reference GRCh38 sequence only; no patient sequence was ever sent to the API. Before the campaign, two calls fixed the conventions of the API: the alternative sequence after a deletion is not realigned, and the RNA score is returned in natural-log units, which the code converts to log2 (`controle_convention`, `verifier_convention`).

**Targets.** MANE promoters were defined from −250 to +100 bp around the TSS, the downstream end stopping 10 bp before the end of exon 1. Because the MANE transcript is not always the one used in a tissue, we also identified, from GTEx v8 transcript-level medians, 173 tissue-dominant TSS that lie outside the MANE promoter deletion and are confirmed by a FANTOM5 CAGE peak (42 in heart, 37 in muscle, 49 in brain, 45 in kidney; for example VDAC1, whose dominant TSS lies 233 bp downstream of the MANE TSS).

**Calibration.** As a positive control of each condition, the whole promoter (−2,500/+100 bp) of every gene, and the dominant TSS where relevant, was deleted and the effect read on transcription initiation (CAGE) and RNA (`calibrer`; about 1,220 deletions per cell type). A condition that does not respond to a promoter deletion cannot reveal smaller effects, so a gene whose whole-promoter deletion had no effect was flagged rather than read as insensitive [TODO: number of genes per type with a coherent, incoherent or non-evaluable calibration status].

**Deletion scan.** Every 10-bp window of the MANE promoters, tissue promoters and on-anchor elements was deleted in each cell type: 58,089 windows and 178,605 calls (`balayage`). Multi-scale deletions (promoters at three scales, whole elements and elements ± 500 bp, 10-kb anchor blocks and their 1-kb slices) added 9,156 deletions (19,803 calls), and a signed reading of RNA, CAGE, accessibility and contact maps added 1,390 calls (`multi`, `lecture_multi`).

**Thresholds fixed beforehand.** An effect of at most −0.2 log2 on the RNA of the host gene was called a decrease; ≤ −1 log2 a strong decrease. The top 5 % of windows by absolute effect, per cell type and per class (promoter or element), were sent to mutagenesis. These values were set before the results were read and were never adjusted on AG2 outputs.

### 5. Promoter deletion scan and saturation mutagenesis

Before ranking windows for mutagenesis, we removed those that overlap a coding exon (13,599 windows) or a splice junction (587 windows) of the target gene, since deleting them damages the transcript itself rather than its regulation. The top 5 % of the remaining windows per type and class gave 8,272 windows, each mutated only in the cell type where it was sensitive. Each base was mutated to its three alternative alleles, for about 248,000 alleles; mutagenesis was completed for the four cell types.

**Promoters carry most predicted effects.** About 28 % of MANE promoter alleles lowered the RNA of the host gene by at least 0.2 log2, with the strongest alleles reaching −5 to −7 log2 (Fig. 3). Critical positions were enriched between −125 and +25 bp relative to the TSS (Fig. 3b). In total, 14,249 positions, corresponding to 51,539 allele × cell type combinations, had a mean predicted effect of at most −0.2 log2.

**Motif grammar.** Bases with predicted effects were enriched in predicted binding sites (JASPAR) of YY1 (4.3-fold), ETS/GABPA (3.2-fold), SP1, NRF1 and NF-Y (Fig. 4) [TODO: fold enrichments and intervals for SP1, NRF1 and NF-Y; state the background used for the enrichment]. These factors are the known regulators of nuclear mitochondrial promoters [REF NEEDED: NRF1/GABPA/YY1 regulation of OXPHOS gene promoters], so the model recovers the expected promoter grammar without being told about it.

**Tissue sharing.** Effects were highly correlated between skeletal muscle and heart (Spearman ρ = 0.92 on shared alleles). Differences between these two types mostly reflected which regions were retained for mutagenesis in each type (UCP3 and FBP2 on the muscle side; HMGCS2, COX6A2 and CPT1B on the heart side) rather than opposite effects of the same allele [TODO: correlations involving brain and kidney now that their mutagenesis is complete].

### 6. Enhancer elements and loop anchors

**Point mutations in elements have weak predicted effects.** No single-nucleotide change in an element lowered the host RNA by 1 log2 or more. Whole-element deletion was more informative: 36 of the 707 elements lowered the host RNA by at least 0.2 log2 ("active" elements). All 36 were off-anchor elements, 35 of them 10–100 kb from the TSS and one beyond 100 kb. Their 1,241 10-bp windows were added to the scan.

**Distance limits what the model reads as an RNA effect.** Almost all on-anchor elements lie more than 100 kb from the promoter. At these distances, AG2 predicted a loss of contact when the anchor or the element was deleted (189 of 334 deleted 10-kb anchors lost contact; some elements also did, e.g. at PRSS35, AHCYL1 and LYRM9) but almost never a decrease of RNA. This is consistent with the known underestimation of distal enhancer effects by sequence models (Karollus et al. 2023; Avsec et al. 2026), and we read it as a limit of the model rather than as evidence that these anchors are dispensable.

**Revision of the anchor rule.** Inspection of loci where an element sat at the end of a loop yet was classed as off-anchor (COQ9, Section 8) led us to revisit the anchor definition, from the loop-calling method rather than from AG2 outputs. peakC sums three 10-kb bins around each bin (window of three), so the resolution of a loop summit is the peakC window, summit ± 15 kb. The start–end region reported for each loop (median 30–40 kb, up to 230 kb) is broadened by this same smoothing and is not the anchor itself. Under the revised rule, an element is on-anchor when it overlaps summit ± 15 kb. This raised the number of on-anchor elements from 349 to 755, including 406 elements in a bin neighbouring the summit ("voisin" class). Anchors are now also deleted as intervals: the whole 30-kb anchor (never touching the target promoter, TSS ± 1 kb) and 10-kb windows tiling it. The complementary computation (new on-anchor elements, whole-anchor deletions and tiling) is running. [TODO: results of the anchor complement: number of active elements among the 406 new on-anchor elements; whole-anchor RNA and contact effects; check from `comparer_cibles` that no element of the campaign was lost.]

**Orthogonal check of contact effects.** To separate architectural from activity effects, a polymer simulation driven only by the position and orientation of CTCF sites (loop extrusion; it never sees the Hi-C map nor the rE2G elements) is being set up on selected loci, with deletion of the CTCF sites inside on-anchor elements or near distal summits, and uniform and randomised CTCF controls [TODO: results; first loci].

### 7. Population constraint in gnomAD

We asked whether alleles predicted to lower expression are depleted in gnomAD v4.1 genomes. For each mutated allele, the expected probability of being observed was derived from its mutability (trinucleotide context, 96 classes), calibrated on the unmutated bases of the campaign targets, and compared with its observed presence. Alleles were grouped into effect classes fixed beforehand (strong decrease ≤ −1; decrease ]−1, −0.2]; neutral ]−0.2, 0.2[; increase; strong increase) and into equal-size groups of predicted effect.

Observed/expected ratios were close to 1 across effect classes, with 95 % confidence intervals crossing 1 (Fig. 5a). A logistic regression of observation on the amplitude of the predicted decrease and increase, with the logit of the expected probability as offset, gave no clear association [TODO: odds ratios per log2 and their intervals, per cell type]. Homozygous genotypes compared with their Hardy–Weinberg expectation were dominated by common alleles and were not conclusive (Fig. 5b). We interpret the absence of clear constraint as expected for genes whose diseases are mostly recessive: a heterozygous allele that halves expression is not under strong selection. The absence of signal also means that gnomAD will not provide a constraint argument for an individual regulatory variant of these genes [TODO: result of the intergenic calibration, `ism_positions --neutre-gtf`].

### 8. Integrated locus examples

For each gene of the PanelApp panel we drew one plate per tissue that combines the Hi-C map, called loops, rE2G elements, AG2 tracks (whole deletion of each element and 10-kb anchor, predicted contact loss), genes and the virtual 4C profile on the left, with the AG2 zooms of the promoter and scanned elements (10-bp deletions and mutagenesis) on the right (Fig. 6) [TODO: number of plates per tissue, and number of genes drawn without AG2 zoom]. Three loci illustrate the main patterns.

**COX6A2 — a strong promoter effect in a highly expressed gene.** COX6A2 is highly expressed in muscle (median 2,284 TPM in GTEx). Its strongest promoter allele lowered the predicted RNA by 1.81 log2, and 33 positions of its promoter had a predicted effect of at most −0.2 log2 (Fig. 6a) [TODO: confirm the cell type of the −1.81 allele (COX6A2 is listed among the heart-side regions) and the exact meaning of "33 positions"].

**BOLA3 — contact loss without RNA change at a distal anchor.** BOLA3 has a loop of about 330 kb. Deleting its distal anchor reduced the predicted contact (−0.59, change in log observed/expected) without changing the predicted RNA (Fig. 6b), the pattern seen for most anchors beyond 100 kb.

**COQ9 — an element at the end of a loop.** At COQ9, an rE2G element lies at the end of a called loop but fell outside the initial anchor window (summit ± 5 kb). This case motivated the revised anchor rule (Section 6) (Fig. 2b, Fig. 6c) [TODO: AG2 results for this element under the revised rule].

### 9. Estimate of a clinical filter

The campaign yields a list of positions where a single-nucleotide change is predicted to lower the expression of a mitochondrial gene: 14,249 positions (51,539 allele × cell type combinations) with a mean predicted effect of at most −0.2 log2. The main contributors are proximal promoter positions, often in binding sites of the factors listed above. From these results we propose three rules for patients without a diagnosis after genome sequencing:

1. **SNVs** in proximal promoters, at positions of the list, especially in predicted binding sites;
2. **deletions and indels** for promoters that are predicted to be sensitive to deletion but show few sensitive single bases;
3. **structural variants** for enhancer elements and loop anchors, whose single-base effects are predicted to be weak.

Because these diseases are mostly recessive and gnomAD shows no constraint on these positions, a regulatory variant should be considered only with a second allele in the same gene (coding or regulatory) and functional evidence (RNA-seq of an accessible tissue, allelic imbalance, or biochemical tests of mitochondrial function). [TODO: expected number of rare variants per genome that fall on the 14,249 positions, from gnomAD allele frequencies; number of candidates in the unsolved cohort, computed locally.] The list is integrated as an annotation layer in the local variant re-annotation pipeline (MitoLoop module of reannot) and does not change any score threshold of other predictors.

---

## Discussion

We combined tissue Hi-C loops, reproducible ENCODE-rE2G enhancers and a systematic AlphaGenome campaign to predict which non-coding bases and regions of 1,220 nuclear mitochondrial genes matter for their expression in four tissues. The main result is that predicted effects concentrate in proximal promoters, in a window of about 150 bp around the TSS and in binding sites of YY1, ETS/GABPA, SP1, NRF1 and NF-Y, the known grammar of these promoters. Effects are largely shared between muscle and heart. Single-base changes in distal elements are predicted to have weak effects, and beyond 100 kb the model predicts loss of chromatin contact without a change in RNA. Bases with predicted effects are not depleted in gnomAD.

**Limits.** First, the weakness of predicted distal effects must not be read as an absence of distal regulation. Sequence models have been shown to underestimate the contribution of distal enhancers (Karollus et al. 2023), and AlphaGenome shares this limit at long range (Avsec et al. 2026). In our data, the dissociation between predicted contact loss and unchanged RNA beyond 100 kb is most likely a property of the model. A planned test will check, on a few CRISPRi-validated enhancer–gene pairs beyond 100 kb from the ENCODE-rE2G benchmark (Gschwind et al. 2026), whether the API reproduces the published behaviour; the BCL11A +58 enhancer (Canver et al. 2015) is a possible positive control, although it lies within 100 kb. Second, the kidney contact maps of AG2 are imputed and should be read as exploratory; the muscle, heart and brain maps are observed (gastrocnemius medialis, heart left ventricle, dorsolateral prefrontal cortex). Third, no prediction has been validated experimentally; the effect sizes are model outputs, and the −0.2 log2 threshold is a prioritisation threshold, not a pathogenicity threshold. Fourth, the 80 % reproducibility rule for enhancers is stringent and excludes elements active in a subset of cell states. Fifth, the loops are called from virtual 4C at promoter viewpoints at 10-kb resolution, which cannot detect loops between two non-viewpoint positions and places anchors to within ± 15 kb. Sixth, the population analysis is limited by recessivity: heterozygous loss of half the expression of most of these genes is not expected to be selected against, so the absence of constraint is uninformative about function. Finally, the tissue-dominant promoters were identified from GTEx isoform quantification, which is an estimate.

**Clinical implications.** The results suggest a practical division of labour for non-coding variants in these genes: SNVs in proximal promoter binding sites, deletions or indels for promoters sensitive only to deletion, and structural variants for enhancers and anchors. In all cases a second allele and functional evidence are required. Non-coding variants in rare disease are increasingly reported [REF NEEDED; possibly Varshney et al. 2026, doi:10.1038/s41588-026-02786-6 — TODO: confirm that its content supports this statement], and our list provides a prior for where to look in mitochondrial disease genes.

**Next steps.** (i) Extend the campaign to blood and fibroblasts, the tissues sampled in patients, and compare predictions with RNA-seq expression outliers in patient fibroblasts (Kopajtich et al. 2021; Yépez et al. [DOI TO CHECK]). (ii) Study exonic and intronic enhancers: exonic sequences can act as enhancers (Mouren et al. 2026), and we have counted the rE2G elements lying in introns of their own gene (`introniques` subcommand) [TODO: counts]. (iii) Cross-check the elements with independent enhancer resources such as VISTA and FANTOM5 enhancers [REF NEEDED for both]. (iv) Precompute, on the reference sequence, the contact effect of CTCF sites in anchors that lose contact, since patient variants are never sent to the API and such effects can only be looked up. (v) Complete the polymer simulation check. In parallel, a separate workstream examines non-coding RNA effector genes (snRNAs, miRNAs, lncRNAs) whose variants could produce a mitochondrial phenotype; it is not reported here.

---

## Methods

### Gene set
The gene list (`genes.txt` in mitoloop_v0) is the union of MitoCarta3.0 and of the PanelApp "Mitochondrial disorders" panel (all confidence levels), with PanelApp symbols mapped to the GENCODE GTF through an alias table of HGNC renames (`re2g_panel.py` alias table; renames listed in `data/genes_renamed.tsv`). Entries that are not genes (repeat expansions, copy-number regions) and mtDNA-encoded genes are dropped where a step requires a nuclear TSS. The gene set contains 1,220 genes; 1,183 had a MANE Select promoter in the campaign.

### Hi-C data and virtual 4C
ENCODE intact Hi-C libraries of human tissues were read remotely (`.hic`, hicstraw) with `scripts/v4c_pipeline.py extract`. For each promoter viewpoint (TSS ± 10 kb), the raw observed contacts with every 10-kb bin within ± 1 Mb were stored. Libraries were selected by depth (`v4c_pipeline.py select`: libraries below the median minus `--k-sd` standard deviations of log10 depth are dropped, the `--min-keep` deepest always kept) and exported for R (`export`) [TODO: number of libraries per tissue, accessions in Supplementary Table S2].

### Loop calling per tissue
Libraries were grouped by tissue with regular expressions on the ENCODE biosample term (`scripts/loops_par_tissu.py librairies`: muscle, heart, brain, kidney). Each group was called on its own with `scripts/v4c_loops.R call` (peakC `combined.analysis`, rank product over replicates; `single.analysis` for a group with one library), window size 3 bins (`--wsize 3`, chosen with `v4c_loops.R scan` on the reproducibility between libraries of one term, Jaccard index of peaks). Profiles at partner anchors were extracted (`loops_par_tissu.py partenaires`, `v4c_pipeline.py extract --anchors`, `fusionner`), and `v4c_loops.R reciprocal` kept a peak when the partner's profile had a peak at the mitochondrial TSS. `loops_par_tissu.py separer` wrote one file per tissue (`loops_<tissue>.tsv`: gene, chrom, viewpoint, summit, start, end). As a control only, `v4c_loops.R call_lib` called each library alone and `separer --peaks-lib` counted, for each loop, the libraries in which it is also found (`n_presents`).

### Loop anchor definition
The anchor of a loop is the peakC window centred on its summit: three 10-kb bins, summit ± 15 kb (`ancre_region` in `mcia/ag2_balayage.py`; option `--ancre-case`, default 10 kb). The start–end interval of a loop is widened by peakC smoothing and is not used as the anchor. The initial campaign used summit ± 5 kb. `comparer_cibles` compares the targets of the two rules without API call (elements reclassed, new, lost; `ancre_precision`: "sommet" inside the central 10-kb anchor deletion, "voisin" on the anchor outside it).

### ENCODE-rE2G elements
Per-biosample ENCODE-rE2G link files were processed with `scripts/re2g_panel.py`: links binarised at the published operating point of the detected model (score giving 70 % recall on the CRISPR benchmark), distal links only, elements mapped onto the ENCODE cCRE registry (SCREEN) to share one coordinate system, and each element listed with the biosamples linking it (`re2g_elements.tsv`). For the campaign, biosamples were grouped by tissue with regular expressions on their cell type (`CELLTYPES` in `ag2_balayage.py`), groups of at least 10 biosamples were kept (muscle 64, heart 48, brain 60, kidney 33), and an element was kept for a gene and a tissue when linked in at least 80 % of the group's biosamples.

### Tissue-dominant promoters
`scripts/gtex_promoteurs.py` grouped the transcripts of each gene by TSS (groups spanning at most 1 kb) using the GENCODE v26 GTF of GTEx v8, read the median TPM of every transcript per tissue from the GTEx portal API, and computed the share of expression initiated inside the campaign's MANE promoter deletion. A TSS was retained as a tissue target (`--tss-dominants`) when the gene was expressed in the tissue, less than half of its expression came from TSSs within the MANE calibration deletion, and the dominant TSS was confirmed by a FANTOM5 CAGE peak and not tagged `mRNA_start_NF`. EPDnew is reported but not used as confirmation, being derived from CAGE. The rule was fixed beforehand and never judged on AG2 outputs.

### AlphaGenome campaign
All calls were made with the AlphaGenome Python client (version ≥ 2.0) to the API (non-commercial research use) from a connected node of the MCIA Curta cluster, one cell type at a time, through SLURM (`mcia/ag2_balayage.sbatch`). Code: `mcia/ag2_balayage.py` in the reannot repository; results in `results/ag2_v2`. Subcommands, in order:

- `conditions`: lists the API conditions and chooses one per cell type by ontology identifier, preferring observed over imputed RNA-seq tracks (`choix_conditions.tsv`).
- `controle_convention`, `verifier_convention`: two calls fixing the coordinate convention of the alternative allele (no realignment after a deletion) and the RNA unit (natural log, converted to log2) (`convention_ag2.tsv`).
- `cibles`: targets and 10-bp windows without API call: MANE promoters (−250/+100 bp, strand-aware, downstream end capped 10 bp before the end of exon 1), tissue promoters (`--tss-dominants`), on-anchor elements (`--boucles-tissus`, `--elements`) and one off-anchor element per disease gene and cell type (`--panel-hors-ancre`).
- `annoter`: splice junctions (`jonctions.tsv`) and coding exons (`cds.tsv`) of the target genes; TSSs of other genes and CTCF peaks near each element.
- `calibrer`: whole-promoter deletion (−2,500/+100 bp) of every gene in every type as a positive control; initiation decomposed by annotated TSS (CAGE), regional initiation, mechanical share of the RNA effect, and a status (coherent, incoherent, non-evaluable).
- `balayage`: deletion of every 10-bp window in every type; RNA (mean of the tracks of the expression output of the condition) and accessibility.
- `multi`, `lecture_multi`: multi-scale deletions (promoter central, 1 kb upstream, ± 1.5 kb; whole elements and ± 500 bp; 10-kb anchor blocks and 1-kb slices; under the revised rule, the whole anchor and 10-kb windows tiling it, never touching TSS ± 1 kb), with signed readouts of RNA, CAGE, accessibility and contact maps (change in log observed/expected; for an element, the 3 × 3 pixel between promoter and summit; for a deleted anchor, contact towards the deletion flanks).
- `elements_actifs`: elements whose whole deletion lowers the host RNA in their tissue by at least 0.2 log2 (threshold fixed beforehand); their 10-bp windows are added to the scan.
- `sites --fenetres`: per type and class, the top 5 % of windows by absolute RNA effect, after removing windows overlapping a coding exon or a splice junction of the target gene (and, optionally, windows of genes not expressed in the tissue); each window is mutated only in the type where it is sensitive (`sites_ism.tsv`).
- `ism`: each base of the selected windows mutated to its three alternative alleles, in the selected type.
- `bilan`: one row per gene and type (strongest promoter window and its position, windows ≤ −0.2 log2, whole-promoter deletion and status, elements, strongest mutagenesis allele, GTEx TPM and expression flag, predicted contact loss of elements and anchors).
- `comparer_cibles`, `introniques`: call-free comparison of anchor rules and count of intronic rE2G elements (with distance to TSS: < 30 kb, 30–100 kb, > 100 kb).

A deletion or allele effect is reported as log2(ALT/REF) of the host gene RNA, averaged across the tracks of the condition's expression output. Thresholds: −0.2 log2 (decrease), −1 log2 (strong decrease), top 5 % (mutagenesis). Measured call time was about 0.75 s per deletion call (about 50,000 calls in about 10 h per type).

### Motif analysis
Predicted binding sites of the factors of mitochondrial promoters (NRF1, GABPA and ETS factors, YY1, SP1/KLF, NF-Y, CREB/ATF, ERRα/γ, THAP11, ZNF143, MAX/USF, ZFX) were found with JASPAR 2024 CORE vertebrate non-redundant matrices on both strands (`mcia/contrainte_motifs.py sites`; log-odds with a 0.25 pseudo-count, relative score ≥ `--min-rel`), with observed binding from ReMap 2022 where available. The factor list and the site threshold were fixed before reading AG2 results. Enrichment of positions with predicted effect in the sites of each factor [TODO: exact procedure and background; script and option used].

### gnomAD analysis
`ag2_balayage.py ism_positions` prepared the targets and mutated bases; `contrainte_motifs.py gnomad` read PASS alleles of gnomAD v4.1 genomes on each base (tabix) and computed the expected probability of observing each allele from its trinucleotide context (96 classes), calibrated on the unmutated bases of the campaign targets (Roulette rates may replace the context model; `--couverture` can exclude poorly covered bases). `ag2_balayage.py ism_gnomad` compared, per type and per effect class fixed beforehand, observed and expected counts (O/E with 95 % interval, ratio to the neutral class, share of common alleles), fitted a logistic regression of observation on decrease and increase amplitude with logit(expected) as offset (odds ratio per log2), and compared homozygote counts with Hardy–Weinberg expectation. An intergenic calibration far from genes is available (`ism_positions --neutre-gtf`). Plots: `mcia/ag2_profils.py --gnomad`.

### Figures
Deletion profiles: `mcia/ag2_profils.py` (`--top N`, `--tous-types`); mutagenesis maps: `ag2_profils.py --ism` (`--seulement-elements`, `--sans-elements`); single-gene zooms: `mcia/ag2_figure_locus.py --gene`; locus figures: `scripts/locus_figure.py` (rotated observed/expected Hi-C map, loop arcs, rE2G elements with on-anchor elements in gold, AG2 tracks with `--ag2-dir`, gene track, virtual 4C); choice of illustrative loci: `scripts/choisir_locus.py` (a stated ranking, not tuned on any output); plates: `scripts/planches_panelapp.py` (one plate per PanelApp gene and tissue).

### Polymer simulation control
`polymere/` uses path A of the `simpol` package: loop extrusion stopped by oriented CTCF sites, followed by a coarse-grained polymer simulation (OpenMM). Its only input is the position and orientation of ENCODE CTCF sites, so its predicted contacts are independent of the Hi-C maps and of the rE2G elements. For each locus, `preparer_locus.py` defines the window, anchors, deletion conditions (CTCF sites inside on-anchor elements; sites within 10 kb of each distal summit) and controls (uniform and random CTCF sites); `evaluer_4c.py` requires the wild-type simulation to beat the random control before a deletion is read; `contact_deletion.py` reports the change in promoter–anchor contact with bootstrap intervals.

### Rules for data use and ethics
- AG2 and SpliceAI outputs are never used to train a model or to tune a threshold. Thresholds were fixed beforehand (−0.2 log2, top 5 %) or derived from independent data (GTEx, CAGE, ENCODE-rE2G, Hi-C, gnomAD).
- Only reference sequences were sent to the AG2 API. No patient data leaves the local environment; patient variants are annotated offline from precomputed values.
- OMIM is licensed; it is used locally and not redistributed.
- Supplementary data from publishers are used internally only and are not redistributed.
- [TODO: ethics statement and consent for any patient cohort data, if patient results are added.]

### Software and versions
Python 3 [TODO: version], AlphaGenome client ≥ 2.0, pysam, hicstraw, matplotlib, numpy; R [TODO: version] with peakC (GitHub deWitLab/peakC), limma; GENCODE [TODO: release of the project GTF], GTEx v8, gnomAD v4.1, JASPAR 2024, ReMap 2022, FANTOM5 CAGE phases 1 and 2, EPDnew 006, GRCh38. Tests on synthetic data: `bash tests/run_all.sh` in both repositories.

---

## Data and code availability

- Code: mitoloop_v0 [URL/DOI TODO; Zenodo snapshot TODO], reannot [URL/DOI TODO; Zenodo snapshot TODO; note that files derived from licensed or publisher data must be removed before release].
- Loops per tissue, selected rE2G elements, campaign targets and the per-gene summary (`bilan_genes.tsv`): [TODO: repository/DOI].
- Precomputed AG2 predictions: [TODO: check AlphaGenome terms of use before redistributing model outputs].
- Public inputs: ENCODE (Hi-C, rE2G, cCRE), GTEx v8, gnomAD v4.1, MitoCarta3.0, PanelApp, JASPAR, ReMap, FANTOM5, EPDnew.

---

## Figure legends

**Figure 1. Overview of the campaign.** Flowchart from the gene set (1,220 genes; 1,183 MANE promoters) to tissue loops (muscle 653, heart 2,122, brain 1,494, kidney 612), reproducible rE2G elements (≥ 80 % of a tissue's biosamples), AG2 targets (MANE promoters, 173 tissue promoters, 349 on-anchor and 358 off-anchor elements), calibration, 10-bp deletion scan (58,089 windows, 178,605 calls), multi-scale deletions (19,803 calls) and signed readouts (1,390 calls), active elements (36/707), selection of the top 5 % windows (8,272) and mutagenesis (about 248,000 alleles), then gnomAD comparison and the clinical filter (14,249 positions). Thresholds fixed beforehand are shown in boxes.

**Figure 2. Tissue loops and anchors.** (a) Number of loops per tissue and distribution of loop lengths [TODO]. (b) Example locus (COQ9): rotated observed/expected Hi-C map, loop arcs, rE2G elements of the tissue (gold: on anchor), genes and virtual 4C profile with called loops. (c) Anchor rule: summit ± 5 kb (initial), summit ± 15 kb (peakC window, revised) and the start–end interval broadened by smoothing; counts of on-anchor elements under each rule (349, 755).

**Figure 3. Promoter deletion profiles and mutagenesis maps.** (a) Effect on host RNA (log2) of each 10-bp deletion along the MANE promoter for the top genes of the summary, small multiples on a common scale; dark bars ≤ −0.2 log2, hatched windows excluded from mutagenesis (coding exon or splice junction), dots for mutated windows. (b) Distribution of critical positions relative to the TSS (enrichment from −125 to +25 bp). (c) Mutagenesis maps (bases × alternative allele, diverging scale centred on 0, blue = decrease) for selected genes.

**Figure 4. Motif grammar of predicted critical bases.** Enrichment of positions with predicted effect ≤ −0.2 log2 in predicted binding sites of YY1 (4.3-fold), ETS/GABPA (3.2-fold), SP1, NRF1 and NF-Y [TODO: remaining values and intervals], and example sites with per-base effects.

**Figure 5. Predicted effect and gnomAD.** (a) Observed/expected ratio of mutated alleles in gnomAD v4.1 genomes per class of predicted effect, per cell type, with 95 % intervals. (b) Homozygotes observed against Hardy–Weinberg expectation per class.

**Figure 6. Integrated locus plates.** (a) COX6A2; (b) BOLA3 (loop of about 330 kb; anchor deletion: contact −0.59, no RNA change); (c) COQ9. Left: locus figure with AG2 tracks; right: AG2 zooms of the promoter and scanned elements.

**Supplementary Figure S1.** Muscle versus heart allele effects (Spearman 0.92) [TODO: pairs including brain and kidney].
**Supplementary Figure S2.** Active elements (36) by distance to TSS; contact loss versus RNA change by distance (189/334 anchors losing contact).
**Supplementary Figure S3.** Calibration of the four conditions (whole-promoter deletions) [TODO].

---

## Supplementary tables

- **Table S1.** Gene set (1,220 genes) with source (MitoCarta3.0, PanelApp level), MANE transcript, renamed symbols.
- **Table S2.** Hi-C libraries per tissue (ENCODE accessions, depth, kept) and rE2G biosamples per tissue group.
- **Table S3.** Loops per tissue (`loops_<tissue>.tsv`), with `n_presents`.
- **Table S4.** Campaign targets (`cibles.tsv`): promoters, tissue promoters, elements on/off anchor, with anchor precision under the revised rule.
- **Table S5.** Calibration per gene and type (`calibration.tsv`).
- **Table S6.** Per-gene summary (`bilan_genes.tsv`).
- **Table S7.** Active elements (`elements_actifs.tsv`).
- **Table S8.** Positions and alleles with predicted effect ≤ −0.2 log2 (14,249 positions, 51,539 allele × type combinations) [TODO: check redistribution terms of model outputs].
- **Table S9.** gnomAD comparison (`ism_gnomad_resume.tsv`, homozygotes).
- **Table S10.** Intronic rE2G elements (`elements_introniques.tsv`) [TODO].

---

## References

References carry a DOI only when it was supplied as verified or appears in the repository files; others are marked [DOI TO CHECK].

- Avsec Ž et al. 2026. AlphaGenome. *Nature*. doi:10.1038/s41586-025-10014-0
- Canver MC et al. 2015. *Nature*. doi:10.1038/nature15521
- Chen S et al. 2024. A genomic mutational constraint map using variation in 76,156 human genomes (Gnocchi). *Nature* 625:92-100. doi:10.1038/s41586-023-06045-0 (from the mitoloop_v0 README)
- di Iulio J et al. 2018. *Nat Genet*. [DOI TO CHECK]
- Fulco CP et al. 2019. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. *Nat Genet*. doi:10.1038/s41588-019-0538-0
- Geeven G et al. 2018. peakC. *Nucleic Acids Res* 46:e91. doi:10.1093/nar/gky443 (from the mitoloop_v0 README)
- Gschwind AR et al. 2026. ENCODE-rE2G. *Nature*. doi:10.1038/s41586-026-10781-4
- Karczewski KJ et al. 2020. gnomAD. *Nature*. [DOI TO CHECK]
- Karollus A et al. 2023. *Genome Biol*. doi:10.1186/s13059-023-02899-9
- Kopajtich R et al. 2021. RNA-seq and proteomics of patient fibroblasts in mitochondrial disease. *medRxiv*. doi:10.1101/2021.03.09.21253187 (from the mitoloop_v0 README; check for a journal version [DOI TO CHECK])
- Lek M et al. 2016. *Nature*. [DOI TO CHECK]
- Martin AR et al. 2019. PanelApp. *Nat Genet*. doi:10.1038/s41588-019-0528-2
- Mouren et al. 2026. Exonic enhancers. *Nat Commun*. doi:10.1038/s41467-026-71220-6
- Rao SSP et al. 2014. *Cell*. [DOI TO CHECK]
- Rath S et al. 2021. MitoCarta3.0. *Nucleic Acids Res* 49:D1541-D1547. doi:10.1093/nar/gkaa1011 (from the mitoloop_v0 README)
- Samocha KE et al. 2014. *Nat Genet*. [DOI TO CHECK]
- Seplyarskiy V et al. 2023. Roulette. [DOI TO CHECK]
- Varshney et al. 2026. *Nat Genet*. doi:10.1038/s41588-026-02786-6 [TODO: confirm content supports the sentence where it is cited]
- Yépez VA et al. Fibroblast RNA-seq diagnostics. [year and DOI TO CHECK]
- [REF NEEDED] Diagnostic yield of genome sequencing in mitochondrial disease.
- [REF NEEDED] NRF1/GABPA/YY1 regulation of nuclear OXPHOS promoters.
- [REF NEEDED] VISTA enhancer browser; FANTOM5 enhancer atlas.
- Data resources also to cite [REF NEEDED or DOI from mitoloop_v0 README]: ENCODE 2020, GENCODE 2023, GTEx 2020, MANE (Morales 2022), JASPAR 2024, ReMap 2022, FANTOM5 CAGE, EPDnew.
