# Preprint v1: module-discovery audit, curated T states, CoVarNet re-discovery and bulk v2

Date range: 2026-11-25 to 2026-11-27 (UTC+8). This is a technical record. It gives scope, method, versions and checksums; scientific estimates, manuscript text and cell-level tables stay in the research environment.

## User requests (in order)

1. Read the repository and server manuscript (v12) and propose a coherent story.
2. Check whether enriched captures violate CoVarNet input assumptions, whether the immune-atlas annotation is ready for downstream analysis, and whether the regulatory/cytotoxic co-presence is biologically meaningful.
3. Check whether CellTypist had been run; record the previous annotation logic; define cell states, re-annotate finely and re-analyse.
4. Re-run CoVarNet on the curated states.
5. Rebuild the bulk readout, check treatment-associated changes, and draft text.
6. Organise a DOI-holding preprint v1 consistent with later versions and with the 4.67M-cell atlas.
7. Draft the full manuscript, figures and supplementary tables, and prepare the code for GitHub.

## Starting state

- Repository main at `0f20c3a`.
- Latest manuscript: private v12 (2026-11-22).
- Immune atlas v1 authority Zarr and 2026-10-14 CM outputs read-only on the WEHI labshare.

## Findings recorded as method facts

- **Annotation and capture.**
  - The immune atlas v1 annotation code (`code_v5/immune_atlas_biology.py`) does not run CellTypist. It assigns majors by panel score, clusters each lineage with Leiden 0.6 after source Harmony, and names each cluster after a single best-contrast marker. Every node carries `fine_review_pending`.
  - The four sources of the `T_whole_TME_primary_tumor` scope are whole-tumour captures. CD45/lymphocyte-enriched collections were in separate scopes.
- **CoVarNet membership.** The original runner selects members by top-N over a global network. In the 15-node T scope this admitted members with near-zero basis weight.
- **Data integrity.**
  - The labshare copy of the SU2C-MARK harmonised TPM file is truncated at exactly 8 MiB (10,503 of 57,525 lines). A verified copy was re-fetched by HTTP range from Zenodo 11179623 (CRC32 2991744730; SHA256 `892f8687…6160`, identical to the 2026-10-18 receipt). The truncated original is unchanged.
  - Raw ICI cohort matrices used in October are not on the labshare. iAtlas-harmonised versions were retrieved through the public cBioPortal API (selected genes only).

## Methods executed (see `scripts/preprint_v1_20261127/README.md`)

- **Curated T states.** 276,787 cells; 4,609 identity HVGs after excluding stress, cycle, TCR/IG and sex genes; Harmony; Leiden 1.5 (21 clusters); CellTypist 1.7.1 Low/High majority voting; manual cluster labelling to 12 CD4/CD8 states plus innate-like T.
- **Pre-specified patient-level test** (`PREREG_phaseB.md`).
- **CoVarNet re-discovery.**
  - Original runner unchanged; CoVarNet R source `e59b314`; R 4.3.3, NMF 0.28.
  - Designs: original units, baseline only, and independent cohorts; plus top-N=3 and fixed K=4 sensitivity runs.
- **Bulk readout v2** (`PREREG_bulk_v2.md` with Amendment A).
  - Arm genes selected by a fixed single-cell rule.
  - A first selection run compared against all non-T lineages, including groups of 1–10 cells. It was corrected to the nine pre-specified lineages before any outcome data were read; both outputs are retained.
  - A first analysis run had a sample-ID type mismatch (IMvigor210) and the truncated SU2C file. That run is retained as `run1_incomplete` and was superseded after the fixes.
  - Reproduction check: the v1 readout recomputed through the new loaders matched the October estimates.
- **Single-cell pre/post comparison** (exploratory).
- **Preprint v1 figures and supplementary tables.** The R4 network substructure was read from the frozen membership table and rendered network figure (Figure S03, v12).

## Files changed in this repository

- `scripts/preprint_v1_20261127/`: code with machine paths replaced by `TFINE_ROOT` / `PROJECT_V5_ROOT`. All Python files pass `py_compile`.
- `inventories/publication/2026-11-27/preprint_v1/`: Table S1a (R4 module members and summary), Table S2 (source registry), both pre-registrations, README and SHA256SUMS. These are aggregates derived from public inventories.
- This log and a CHANGELOG entry.

## Not in this repository

The manuscript, curated basis weights and edges (S1b/S1c), cell-level tables, patient-level fractions, model outputs and figures with scientific content.

## Validation

- `bash scripts/validate_repository.sh` was run before commit (result below).
- Diff and path/privacy grep were reviewed.

## Sync status

- Local commit on branch `claude/preprint-v1-tex-module-20261127`.
- Push was not performed from this session: GitHub write access was not available to the agent here. To be pushed by the user.

## Unverified / open

- Formal citations for the cBioPortal/iAtlas studies.
- The three accessions shared between the immune atlas and R4.
- R4 network edge statistics are not exported publicly (only membership).

## Suggested commit message

`Add preprint v1 audit code, curated T-state CoVarNet rerun and public aggregate tables`
