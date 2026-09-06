# PanCancer annotation V7 release

- Version: `pancancer_annotation_v7_cnv_rerun_20260905`
- Run: `20260905_v7_cnv_rerun_freeze_v1`
- Status: `FROZEN_PROJECT_REFERENCE`
- Cell universe: 4,676,787 cells, 1,322 samples, 54 GSE and 43 cancers
- Cell-index SHA-256: `3f802eb8cc52f04beb27fa2987d2ffe3ddf75f4b4f51dec8e6448b26df8d1f16`

The full logic and statistical definitions are in
[`methods/V7_REFERENCE_AND_ANNOTATION_LOGIC_CN.md`](methods/V7_REFERENCE_AND_ANNOTATION_LOGIC_CN.md).
Top-level summary tables are copied into `tables/`; per-sample matrices,
Copykat_python outputs and cell-level objects are intentionally excluded.

## Canonical analysis fields

- Major identity: `major_celltype_v7`
- Fine identity: `final_annotation_v7`
- Identity confidence: `annotation_confidence_v7`
- Primary malignancy call: `malignancy_call_v7`
- Evidence tier: `malignancy_tier_v7`
- Display label: `final_annotation_with_malignancy_v7`

Strict filters:

```text
strict_major_eligible_v7 == True and exclude_from_primary_major_v2 == False
strict_downstream_subtype_eligible_v7 == True
strict_malignant_downstream_eligible_v7 == True
```

The Python files are exact archival copies for methods review, not a standalone
relocatable workflow. Reproduction must use a new append-only run directory and
the original local dependency layout after the local-only preflight passes.
Verify heavy inputs and outputs against `external_artifacts.tsv`.
