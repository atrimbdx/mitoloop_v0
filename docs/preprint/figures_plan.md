# Figures plan for the preprint draft

Each figure is drawn from existing campaign outputs, without any new AlphaGenome call.
Paths: reannot results in `~/reannot_git/results/ag2_v2` (Curta), mitoloop_v0 loops in `results/loops_tissus`.

| Figure | Content | Script / command | Inputs | Status |
|---|---|---|---|---|
| 1 | Campaign flowchart with counts at each step | to draw by hand (SVG), or `scripts/svg_to_pptx.py` for an editable version | numbers of the draft | TODO |
| 2a | Loops per tissue, loop-length distribution | small new plot from `loops_<tissu>.tsv` (no script yet) | `results/loops_tissus/loops_*.tsv` | TODO |
| 2b | COQ9 locus: Hi-C map, arcs, rE2G elements, virtual 4C | `scripts/locus_figure.py` (call written by `scripts/choisir_locus.py`) | `.hic` of the tissue, `loops_<tissu>.tsv`, `re2g_elements.tsv` | TODO |
| 2c | Anchor rule (±5 kb, ±15 kb, start-end) | schematic | `ag2_balayage.py comparer_cibles` output | TODO |
| 3a | 10-bp deletion profiles along MANE promoters | `mcia/ag2_profils.py --type <t> --top 20` | `bilan_genes.tsv`, `balayage_<t>.tsv`, `sites_ism.tsv`, `sites_exclus.tsv` | done for muscle/heart? [TODO check] |
| 3b | Critical positions vs TSS | from `ism_<t>.tsv` (histogram; no script yet) | `ism_*.tsv` | TODO |
| 3c | Mutagenesis maps | `mcia/ag2_profils.py --ism --sans-elements` | `ism_*.tsv`, `bilan_genes.tsv` | TODO for 4 types |
| 4 | Motif enrichment | `mcia/contrainte_motifs.py` sites + enrichment table [TODO: script] | `sites.tsv`, `ism_*.tsv` | TODO |
| 5a | Predicted effect vs gnomAD O/E | `mcia/ag2_profils.py --gnomad` | `ism_gnomad_resume.tsv` | rerun on 4 types |
| 5b | Homozygotes vs Hardy-Weinberg | from `ism_gnomad_homozygotes.tsv` | | TODO |
| 6 | Integrated plates (COX6A2, BOLA3, COQ9) | `scripts/planches_panelapp.py --tissu <t>`; zooms `mcia/ag2_figure_locus.py --gene` | `.hic`, loops, elements, `ag2_v2` | TODO |
| S1 | Muscle vs heart allele effects | new scatter from `ism_*.tsv` | | TODO |
| S2 | Active elements by distance; contact loss vs RNA | from `elements_actifs.tsv`, `lecture_multi_<t>.tsv`, `bilan_genes.tsv` | | TODO |
| S3 | Calibration | from `calibration.tsv` | | TODO |

Style: blue = decrease of RNA in all AG2 panels (as in `ag2_profils.py`); coding-exon and splice windows hatched; kidney contact maps flagged "imputed".
