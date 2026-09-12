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
- `communication_evidence_chain_20260906`: manuscript-facing communication,
  spatial and clinical evidence summaries and figures.

The annotation release is the current identity/malignancy reference. The
communication release was derived from the earlier V5 downstream CoVarNet/LIANA
branch and must not be described as recomputed from V7.
