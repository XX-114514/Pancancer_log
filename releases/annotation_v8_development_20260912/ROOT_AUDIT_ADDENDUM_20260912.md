# V8 root-audit addendum (2026-09-12)

- Addendum ID: `v8_root_audit_addendum_20260912`
- Applies to: `annotation_v8_development_20260912`
- Record type: append-only external-evidence addendum
- External index: [ROOT_AUDIT_EXTERNAL_ARTIFACTS_20260912.tsv](ROOT_AUDIT_EXTERNAL_ARTIFACTS_20260912.tsv)

This addendum records newly inspected external root-audit evidence without changing
the release's existing `README.md`, `VERSION.json`, or `external_artifacts.tsv`.
The external artifacts remain outside Git; this release records only logical paths,
byte sizes, SHA-256 values, and observed scope/status statements.

## Observed root gate

The directly inspected machine-readable root gate records
`PASS_DEVELOPMENT_SOFTWARE_CANDIDATE_NOT_FROZEN`. Its recorded default project
reference is `V7`; automatic reference promotion and formal downstream mutation
remain false. The observed verification record reports `167` V8 pytest tests
passed. The records repository did not re-execute those Project_v3 tests in this
indexing task; it independently recomputed the five artifact sizes and SHA-256
values listed in the external index.

The independent-reference materials remain
`BLOCKED_NO_LOCAL_CELLXGENE_CENSUS_CLIENT`. They are a metadata-only contract and
pending candidate catalog, not an available independent reference, trained scANVI
model, mapping result, or truth set.

## Non-promotion boundary

This evidence supports a development software candidate only. It does not create a
biological V8 freeze, independently establish identity or malignancy accuracy, or
authorize replacement of V7 or mutation of existing downstream results. The V7
frozen reference remains the project default. A biological reference decision still
requires a new, versioned validation run and explicit promotion review.

## Evidence scope

- Root-audit report: human-readable synthesis of candidate integrity, test results,
  gate limitations, and the preserved V7 default.
- Root-release gate: machine-readable status, 167-pytest verification count, and
  the non-promotion controls.
- Independent-reference report, Census contract, and candidate inventory: an active
  metadata-only planning surface whose execution status remains blocked; no
  expression matrix was read or model training/prediction performed.

The source files were inspected read-only from Project_v3 and were not copied into
this repository.
