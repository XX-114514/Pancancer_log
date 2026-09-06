# Versioned lightweight releases

This directory contains Git-safe snapshots for writing, Prism figure assembly and
supplementary-material preparation. Releases contain small reports, summary tables,
figures and configurations. Cell-level objects and raw data remain outside Git and
are represented by logical identifiers, sizes and SHA-256 checksums.

Current releases:

- `annotation_v7_20260905`: frozen 5M-cell annotation and CNV reference.
- `communication_evidence_chain_20260906`: manuscript-facing communication,
  spatial and clinical evidence summaries and figures.

The annotation release is the current identity/malignancy reference. The
communication release was derived from the earlier V5 downstream CoVarNet/LIANA
branch and must not be described as recomputed from V7.
