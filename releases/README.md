# Versioned lightweight releases

This directory contains Git-safe snapshots for writing, Prism figure assembly and
supplementary-material preparation. Releases contain small reports, summary tables,
figures and configurations. Cell-level objects and raw data remain outside Git and
are represented by logical identifiers, sizes and SHA-256 checksums.

Current releases:

- `annotation_v7_20260905`: frozen 5M-cell annotation and CNV reference.
- `annotation_v8_development_20260912`: `NOT_FROZEN` full-metadata software
  candidate index; it preserves V7 as the default reference and is not a biological
  release.
- `annotation_v8_lineage_audit_20260917`: append-only V8 lineage-taxonomy audit
  addendum containing the generated status and five exact freeze manifests; V8
  remains globally `NOT_FROZEN`.
- `annotation_v9_candidate_20260919`: latest safe V9 candidate report and open L3
  review template; it is `FROZEN_CANDIDATE_NOT_PROMOTED` and does not replace V7.
- `communication_evidence_chain_20260906`: manuscript-facing communication,
  spatial and clinical evidence summaries and figures.
- `manuscript_upstream_freeze_nature_20260922`: V7-based publication upstream
  lock, V8/V9 candidate boundaries and Nature Article preparation files.

The annotation release is the current identity/malignancy reference. The
communication release was derived from the earlier V5 downstream CoVarNet/LIANA
branch and must not be described as recomputed from V7.
