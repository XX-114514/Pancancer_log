# Preprint v1 public aggregate tables (2026-11-27)

| File | Content | Derived from |
|---|---|---|
| `Table_S1a_R4_K15_module_members.tsv` | Members of all 15 R4/K15 modules, with high/low unit counts and source coverage | `inventories/analysis/2026-11-19/WEHI_R4_CM_validation` |
| `Table_S1a_R4_K15_module_summary.tsv` | One row per module: member count, sources with high-activity units, largest single-source share, whether Tph/Treg/Tex are all members | as above |
| `Table_S2_source_registry.tsv` | 72 sources: cells in integration v2 and TME v1, malignancy calls, capture, treatment, sampling site, publication metadata, scope decisions | `inventories/tme_v1`, `inventories/source_registry` |
| `PREREG_phaseB.md`, `PREREG_bulk_v2.md` | Pre-specifications written before the corresponding tests (bulk v2 includes Amendment A) | this task |

The tables are source- or module-level aggregates. They contain no cell- or patient-level rows. Curated-state basis weights and edge tables (Supplementary Tables S1b/S1c) and the manuscript remain private until posting. Code: [scripts/preprint_v1_20261127](../../../../scripts/preprint_v1_20261127/README.md). Checksums: `SHA256SUMS`.
