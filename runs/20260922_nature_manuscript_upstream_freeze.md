# Run: 20260922_nature_manuscript_upstream_freeze

## Metadata

- Run ID: `20260922_nature_manuscript_upstream_freeze`
- Run type: publication evidence freeze / documentation release
- Date and timezone: 2026-09-22, Asia/Shanghai
- Status: `FROZEN_MANUSCRIPT_UPSTREAM_SNAPSHOT_WITH_OPEN_SCIENTIFIC_GATES`
- Analysis computation: none; all external analysis artifacts were inspected read-only
- Records repository baseline: `main` at `52dcf6f`, six commits ahead of `origin/main`, clean before this task
- Project_v3 code commit: not available because Project_v3 is not a Git repository

## Objective

Freeze the currently defensible annotation upstream for manuscript writing and
prepare a Nature Article submission package without promoting incomplete V8/V9
candidates or relabeling V5-derived downstream analyses.

## Inputs and validation

- V7 final checkpoint, obs-only H5AD and cell-index hashes were matched to the
  existing frozen release.
- V7 identity parents, gene coordinates, actual runner and dynamically loaded CNV
  engines were independently hashed.
- The 14,207,049,627-byte raw-count Zarr contained 13,696 files. Its deterministic
  directory root was computed from `LC_ALL=C` path-sorted `sha256sum` output using
  `./` relative paths: `33440b3325f37ceea03143f794e9c0163f62628e4b5b791d140eeea8af84faa0`.
- Five V8 lineage-freeze rule/taxonomy/sentinel sets were independently hashed.
- The newer V9 report and manifest were inspected; the manifest says
  `FROZEN_CANDIDATE`, not `FROZEN`.
- Nature's current official author pages were checked on 2026-09-22.

## Decision

- V7 remains the sole global manuscript identity/CNV reference.
- V8 is lineage-taxonomy-only evidence.
- V9 is quarantined as a candidate because manual L3 review, stability gates and
  632 suggested CNV reruns remain unresolved.
- The current communication/spatial/clinical release remains explicitly V5-derived.

## Outputs

- `releases/manuscript_upstream_freeze_nature_20260922/`
- `V7_PUBLICATION_LOCK.tsv`
- `EXTERNAL_ARTIFACTS.tsv`
- Nature manuscript, declaration, figure and supplementary templates

## Warnings

- The historical V7 environment was not fully locked. The preflight records local
  dependency success and Copykat_python commit, but not exact versions for all
  Python/R packages; current environments must not be presented as historical fact.
- The records dataset inventory still reflects the earlier discovery scope rather
  than a publication-ready 54-GSE source/licence manifest.
- V7 code commit is not recorded; script-level hashes are the available provenance.
- The raw-count directory root is a deterministic aggregate hash, not a retained
  per-chunk checksum list.

## Next action

Complete the scientific gates in the release README, then create a new append-only
submission-candidate release rather than modifying this freeze in place.
