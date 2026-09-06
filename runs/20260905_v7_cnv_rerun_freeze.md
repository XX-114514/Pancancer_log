# Run 20260905_v7_cnv_rerun_freeze_v1

## Status

`FROZEN_PROJECT_REFERENCE` as of 2026-09-06.

## Observed inputs and scope

- 4,676,787 cells, 1,322 original samples, 54 GSE and 43 cancer labels.
- Raw non-negative integer counts from the 5M merge Zarr.
- Uncertainty-resolution V2 broad/fine identities.
- Cancer-specific candidate/reference policies.

## Execution

inferCNVpy and fixed-commit Copykat_python were run per original sample using only
local dependencies. inferCNVpy succeeded for 938 samples and explicitly skipped
384 for insufficient candidate/reference counts. Copykat_python succeeded for 668,
skipped 653 for insufficient cells and skipped one for `all_cells_filtered`;
there were no failed terminal samples.

## Main outputs

- Primary malignant: 515,031 cells.
- Uncertain: 484,919.
- Not evaluable candidate: 45,284.
- Two-method high-confidence intersection: 238,418.
- Sensitive union: 513,509.
- Strict broad annotation eligible: 4,057,641.
- Strict subtype eligible: 2,455,424.

## Validation

The cell-index hash is
`3f802eb8cc52f04beb27fa2987d2ffe3ddf75f4b4f51dec8e6448b26df8d1f16`.
The checkpoint and obs-only object hashes are recorded in
[the V7 external artifact manifest](../releases/annotation_v7_20260905/external_artifacts.tsv).
Mapping conflict and unmapped-cell tables are empty apart from headers.

## Interpretation boundary

Skip means not evaluable, not benign. The primary cluster malignancy call, strict
two-method intersection and sensitive union are distinct prespecified definitions.
