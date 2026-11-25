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
| 7a | Non-coding effector classes (47) by level of mitochondrial link and status | small new bar plot (no script yet) | `docs/ncrna_effecteur/inventaire/inventaire_ncrna_effecteurs.tsv` (column `niveau_lien_mito`); reannot `donnees/nc_types.tsv` (column `statut`) | TODO |
| 7b | Three paths and the U12 cascade, with counts of the 24/11 run | schematic by hand (SVG) | counts from `ETAT_CHANTIER_2026-11-24.md`; `elements_nc.rapport.tsv`, `chemin1.rapport.tsv`, `chemin1_cibles.tsv`, `regions_effet.bed`, path 3 `chemin3` report (Curta, `results/nc_effecteur`) | TODO |
| 7c | Base pairs per pre-filter level | small new bar plot (no script yet) | `prefiltre_mito.resume.tsv` from `mcia/elements_nc.py prefiltre` | TODO: rerun (pre-filter not yet run on real data) |
| S1 | Muscle vs heart allele effects | new scatter from `ism_*.tsv` | | TODO |
| S2 | Active elements by distance; contact loss vs RNA | from `elements_actifs.tsv`, `lecture_multi_<t>.tsv`, `bilan_genes.tsv` | | TODO |
| S3 | Calibration | from `calibration.tsv` | | TODO |

Style: blue = decrease of RNA in all AG2 panels (as in `ag2_profils.py`); coding-exon and splice windows hatched; kidney contact maps flagged "imputed".

Figure 7 (non-coding effector axis) uses no AlphaGenome output: it is drawn from the reannot non-coding effector run (`mcia/run_paquet.sh nc_effecteur`, results on Curta under `REANNOT_BASE`) and the mitoloop_v0 inventory.

- **Supplementary Fig. S4** (enhancer threshold sensitivity): reannot `mcia/seuil_re2g.py --elements results/paquet/re2g_tous/re2g_elements.tsv --genes ~/mitoloop_v0/genes.txt --out results/seuil_re2g/seuil_re2g` → `seuil_re2g.png/.svg`, `seuil_re2g.tsv`, `seuil_re2g_fractions.tsv`. Status: done (25 Nov 2026), figure on Curta.
