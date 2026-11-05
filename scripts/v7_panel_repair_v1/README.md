# Literature-guided V7 evidence workflow

This version produces auditable evidence and review queues. It does not assign final cell labels. The authoritative panel is [registry.json](../../resources/annotation/v7_panel_repair_v1/registry.json); biological scope and adaptation decisions are in the [method](../../methods/V7-panel-repair-literature-workflow-20261105.md).

## Aggregate evidence

Requirements: Python3, numpy, pandas. Expression TSV has `cluster,gene,n_cells,n_measured,pct_expressing_measured`; fractions are 0–1 among measured cells. Optional `mean_cp10k_measured` must be linear CP10K, not mean log expression. Each cluster–gene pair is unique. Missing genes remain unknown.

```bash
python scripts/v7_panel_repair_v1/evidence_engine.py \
  --registry resources/annotation/v7_panel_repair_v1/registry.json \
  --expression measured_marker_expression.tsv \
  --de RNA_DEG_common_measured.tsv \
  --out NEW_EVIDENCE_DIRECTORY
```

Optional DE requires Seurat RNA columns `cluster,gene,avg_log2FC,pct.1,pct.2,p_val_adj`. Do not substitute legacy mean-log differences. Source-specific contrasts must be supplied separately. `--de-complete` declares an untruncated supplied contrast, not that every possible gene was tested; never use it with top10. Optional contexts JSON maps cluster IDs to `reviewed_parent_ids` and `organ_contexts`. It must encode reviewed evidence, not a copy of old labels.

The report retains all panel axes and never computes a winner/margin. `each_core_component_has_measured_gene` means at least one gene per component is measured. `complete_core_measurement` requires all core genes measured across the full cluster. Neither means expression or same-cell coexpression. DEG/naming support and coverage are independent fields. Expected-lower genes require a matched comparator and are not automatically interpreted as negative markers.

## Full raw-count lineage rebuilding

Requirements tested: Seurat5.3.1, Matrix1.6.5, jsonlite2.0.0. Harmony optional (installed1.2.4, optional branch not exercised in this smoke); official UCell optional (not installed in tested environment, scoring branch not exercised).

Input RDS is a list with:

```r
list(counts = sparse_gene_by_cell_raw_counts,
     metadata = cell_metadata_in_exact_counts_column_order,
     availability = logical_gene_by_source_measurement_mask,
     full_library_counts = TRUE,
     source_column = "source")
```

Gene names must be unique safe feature names, with an upstream stable-ID mapping for any necessary conversion. Source-mask absence must agree with raw counts. Normalisation uses the full library, not selected panel genes. `full_library_counts` is an upstream contract, not something the script can infer from a matrix alone.

Example config (values are illustrative project choices, not universal literature thresholds):

```json
{
  "registry": "resources/annotation/v7_panel_repair_v1/registry.json",
  "seed": 191,
  "n_hvg": 2000,
  "n_pcs": 50,
  "graph_pcs": 30,
  "resolutions": [0.2, 0.4, 0.6, 0.8],
  "marker_resolution": 0.4,
  "min_model_source_coverage": 0.8,
  "run_ucell": false
}
```

```bash
Rscript scripts/v7_panel_repair_v1/lineage_workflow.R \
  input_contract.rds config.json NEW_LINEAGE_DIRECTORY
```

Optional `reviewed_keep_file`: TSV with every input `cell_id`, logical TRUE/FALSE `keep`, and an explicit reason for exclusions. Without it all cells remain for initial review. The script does not estimate doublets or automatically remove marker-conflicting cells. Optional `verified_batch_column` enables Harmony only for an independently verified batch variable. Optional `r_library_paths` and `ucell_maxRank` support official UCell when `run_ucell=true`. It refuses unavailable UCell instead of substituting a homemade scoring formula.

Outputs include HVGs/PCA deviations, multi-resolution counts, common-measured and source-specific RNA DE, all-panel-gene coverage-aware expression, and source/cluster same-cell core-positive evidence at raw≥1 and raw≥2. Private cell objects, QC manifests and individual UCell scores must remain private. `complete.json` contains actual local config paths; redact before public release.

No automatic regression of cancer/IFN/stress, fixed target cluster number, inherited final parent, or final subtype assignment. Source-specific DE does not establish patient replication. Final sibling comparisons, doublet assessment and study-held-out stability require lineage-specific review.

## Reproducible verification

```bash
python scripts/v7_panel_repair_v1/test_evidence_engine.py
python scripts/v7_panel_repair_v1/run_smoke.py \
  --rscript Rscript \
  --registry resources/annotation/v7_panel_repair_v1/registry.json \
  --out NEW_SYNTHETIC_TEST_DIRECTORY
python scripts/v7_panel_repair_v1/replay_aggregates.py \
  --repo . --out NEW_PUBLIC_REPLAY_DIRECTORY
```

The synthetic fixture contains360 cells/2100genes/two sources and an explicit four-cell QC exclusion. The independent audit recomputes detection, full-library CP10K, Treg same-cell positive evidence and source-dependent missingness. It tests implementation, not biological annotation performance. No input or output directory is silently overwritten.
