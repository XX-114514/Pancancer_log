# Immune atlas: 45 figure legends

## Fig01_atlas_and_cohorts — An immune-treatment-associated atlas

a, All 2,670,057 archived cells at unchanged atlas UMAP coordinates, shuffled drawing order. b, 199,998 study-balanced display cells, seed 20261014; counts in c and eligibility in d use every cell. Study colours match cohort bars. Endpoint eligibility is source-specific and does not count independent patients. Major labels include proposals, unresolved and boundary categories. UMAP fit used 100,000 landmarks and transformed the remaining 2,570,057 cells.

## Fig02_composition_and_identity_evidence — Composition and multigene identity evidence

All archived cells contribute to composition and marker statistics. Dot area is the positive fraction among cells in sources that measured that gene; colour is log1p(mean native chosen-count CP10K), scaled to each gene maximum across displayed rows. Unmeasured genes remain missing. Capture-biased composition is descriptive. Proposals require multigene review; naming markers alone are not definitive identities.

## Fig03_T_descriptive_subclusters — T: descriptive subclusters

1,571,352 archived T cells and 16 descriptive nodes; 188,453 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## Fig04_T_continuous_states — Continuous T-cell state axes

Five saved native panel-score axes, 145,326 study-balanced T display cells. Each score colour range runs from zero to its all-T-cell 99th percentile; saturation is visual only. The heatmap uses all 1,571,352 T cells and descriptive node means standardized separately within each axis, clipped visually at ±2. Scores are state axes, not mutually exclusive cell types or treatment-response tests.

## Fig05_CNV_evidence — Malignancy candidates and evaluation scope

All-cell atlas CNV states and source-specific cell counts; original-library methods, not cell counts, define panel c. Panel d shows each evaluable library median candidate burden divided by its own reference median (log2), with library median across points as black line; no new threshold or hypothesis test. 52,657 single-method supported malignant candidates, 57,958 uncertain cells. Skips and negative candidates do not establish normal identity.

## Fig06_whole_TME_covariance — CM scope: whole_TME

Actual saved whole_TME nsNMF, 80 nodes, rank 2, 2 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED01_CM_scope_coverage — Complete CM discovery scope

All 15 completed composition scopes: 10 have supported network-member modules, five do not; another six prepared scopes were insufficient. Selected NMF rank and supported-member module count are distinct. Scopes overlap and their ranks must not be summed as independent modules. No unsupported result is hidden.

## ED02_PCA_Harmony — Batch and biological structure diagnostics

Identical 199,998 study-balanced display cells before and after Harmony in the first two component coordinates, not a before/after UMAP. All RNA counts remain unchanged. The original 10,000-cell neighbor diagnostic reported entropy 0.4127→0.8077 and marker-proposal purity 0.8223→0.8136; these are descriptive indicators, not a ground-truth annotation or proof of absent overcorrection.

## ED03_depth_and_provenance — Available depth and correction provenance

Left: per-original-library medians of the complete native chosen-count denominator; points are libraries and black segments are study medians. Right: full-cell correction-adopted and primary-fit exclusion fractions. This is available provenance, not a replacement for source QC; original mitochondrial/nGene fields are not universally available and have not been fabricated. Prior-QC sources without new ambient assessment retain their documented status.

## ED04_all_cohort_UMAP — Every cohort in the shared atlas coordinates

All nine sources, each with unchanged common UMAP limits. Up to 60,000 major-balanced display points per source, seed 20261014, preserving rare observed classes; this sampling is not used for composition proportions. Titles map to GEO IDs and the CAR-T CXG source in the cohort registry.

## ED05_capture_scope_UMAP — Capture-specific atlas structure

Every recorded capture scope in the unchanged atlas UMAP, using up to 60,000 major-balanced display cells per scope, fixed seed20261014. Full scope counts use all 2,670,057 cells. Enriched immune, stimulated lymphocyte, selected nuclei, serial biopsy and infusion-product captures are distinct observation contexts; these facets do not estimate comparable population prevalences.

## ED06_assay_and_uncertainty — Technology and annotation uncertainty

Recorded assays are retained exactly, including any source-specific or unknown labels. UMAP uses 199,998 study-balanced display points, fixed seed20261014. Assay proportions and unresolved/boundary fractions are full-cell counts. The 171,341 unresolved and 12,409 T/NK-boundary cells remain archived; their labels have not been forced or equated to malignant or normal truth.

## ED_CM_01_whole_TME — CM scope: whole_TME

Actual saved whole_TME nsNMF, 80 nodes, rank 2, 2 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_02_immune_enriched — CM scope: immune_enriched

Actual saved immune_enriched nsNMF, 40 nodes, rank 9, 9 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_03_CAR_product — CM scope: CAR_product

Actual saved CAR_product nsNMF, 23 nodes, rank 12, 12 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_04_T_lineage_biopsy — CM scope: T_lineage_biopsy

Actual saved T_lineage_biopsy nsNMF, 16 nodes, rank 12, 12 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_05_whole_TME_primary_tumor — CM scope: whole_TME_primary_tumor

Actual saved whole_TME_primary_tumor nsNMF, 102 nodes, rank 12, 12 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_06_T_whole_TME_primary_tumor — CM scope: T_whole_TME_primary_tumor

Actual saved T_whole_TME_primary_tumor nsNMF, 15 nodes, rank 4, 4 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_07_whole_TME_metastatic_tumor — CM scope: whole_TME_metastatic_tumor

Actual saved whole_TME_metastatic_tumor nsNMF, 21 nodes, rank 2, 0 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_08_immune_enriched_primary_tumor — CM scope: immune_enriched_primary_tumor

Actual saved immune_enriched_primary_tumor nsNMF, 23 nodes, rank 3, 3 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_09_T_immune_enriched_primary_tumor — CM scope: T_immune_enriched_primary_tumor

Actual saved T_immune_enriched_primary_tumor nsNMF, 14 nodes, rank 3, 3 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_10_immune_enriched_blood — CM scope: immune_enriched_blood

Actual saved immune_enriched_blood nsNMF, 27 nodes, rank 5, 0 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_11_T_immune_enriched_blood — CM scope: T_immune_enriched_blood

Actual saved T_immune_enriched_blood nsNMF, 11 nodes, rank 5, 0 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_12_immune_enriched_normal_adjacent — CM scope: immune_enriched_normal_adjacent

Actual saved immune_enriched_normal_adjacent nsNMF, 22 nodes, rank 3, 3 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_13_T_immune_enriched_normal_adjacent — CM scope: T_immune_enriched_normal_adjacent

Actual saved T_immune_enriched_normal_adjacent nsNMF, 12 nodes, rank 3, 2 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_14_immune_enriched_unknown_material — CM scope: immune_enriched_unknown_material

Actual saved immune_enriched_unknown_material nsNMF, 38 nodes, rank 2, 0 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_CM_15_T_immune_enriched_unknown_material — CM scope: T_immune_enriched_unknown_material

Actual saved T_immune_enriched_unknown_material nsNMF, 14 nodes, rank 2, 0 factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.

## ED_lineage_01_T — T: descriptive subclusters

1,571,352 archived T cells and 16 descriptive nodes; 188,453 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_02_NK — NK: descriptive subclusters

149,285 archived NK cells and 10 descriptive nodes; 40,000 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_03_B — B: descriptive subclusters

247,517 archived B cells and 16 descriptive nodes; 39,596 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_04_Plasma — Plasma: descriptive subclusters

59,898 archived Plasma cells and 31 descriptive nodes; 31,995 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_05_Myeloid — Myeloid: descriptive subclusters

173,846 archived Myeloid cells and 13 descriptive nodes; 36,979 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_06_Neutrophil — Neutrophil: descriptive subclusters

41,960 archived Neutrophil cells and 15 descriptive nodes; 27,927 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_07_Mast — Mast: descriptive subclusters

8,060 archived Mast cells and 10 descriptive nodes; 8,060 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_08_Epithelial — Epithelial: descriptive subclusters

130,111 archived Epithelial cells and 33 descriptive nodes; 36,096 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_09_Fibroblast — Fibroblast: descriptive subclusters

38,476 archived Fibroblast cells and 14 descriptive nodes; 30,150 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_10_Endothelial — Endothelial: descriptive subclusters

16,907 archived Endothelial cells and 15 descriptive nodes; 16,773 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_11_Pericyte_SMC — Pericyte SMC: descriptive subclusters

7,019 archived Pericyte_SMC cells and 10 descriptive nodes; 7,019 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_12_Erythroid — Erythroid: descriptive subclusters

2,586 archived Erythroid cells and 13 descriptive nodes; 2,586 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_13_Glial — Glial: descriptive subclusters

18,675 archived Glial cells and 14 descriptive nodes; 15,690 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_14_Neuron — Neuron: descriptive subclusters

12,458 archived Neuron cells and 24 descriptive nodes; 12,257 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_15_Oligodendrocyte — Oligodendrocyte: descriptive subclusters

8,062 archived Oligodendrocyte cells and 11 descriptive nodes; 8,062 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_16_Melanocytic — Melanocytic: descriptive subclusters

95 archived Melanocytic cells and 1 descriptive nodes; 95 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_17_T_NK_boundary — T NK boundary: descriptive subclusters

12,409 archived T_NK_boundary cells and 9 descriptive nodes; 12,409 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.

## ED_lineage_18_Unresolved — Unresolved: descriptive subclusters

171,341 archived Unresolved cells and 15 descriptive nodes; 33,282 node-balanced display points, seed 20261014, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.
