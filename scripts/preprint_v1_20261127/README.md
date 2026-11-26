# Preprint v1 analysis code (2026-11-25 to 2026-11-27)

Code behind the preprint v1 draft: curated T-cell re-annotation, CoVarNet re-discovery on curated states, bulk readout v2, single-cell pre/post comparison, figures and supplementary tables. Scientific outputs, manuscript text, cell-level tables and learned weights stay in the research environment. This directory holds code only.

## Environment variables

| Variable | Meaning |
|---|---|
| `PROJECT_V5_ROOT` | Project_v5 root holding the 2026-10-14 immune atlas v1 authority Zarr and the 2026-10-18 bulk run |
| `TFINE_ROOT` | New working directory for this analysis (outputs are written here) |

## Order of execution

1. `annotation/obs_check.py <ZARR>`: read-only observation metadata export (`obs_small.parquet`; place it in `TFINE_ROOT`).
2. `annotation/extract_T.py <OUT.h5ad> 1000`: T cells from the four CM-eligible discovery sources plus GSE169246 and lymphocyte-enriched GSE176021; random cap of 1,000 cells per sample.
3. `annotation/annotate_T.py`: identity clustering with stress/cycle/TCR/IG/sex genes removed, Harmony (source, sample), Leiden 1.5, CellTypist Immune_All_Low/High majority voting, marker programmes. Writes cluster review tables; labels are assigned after review.
4. `annotation/phaseB.py`: curated state assignment (cluster map recorded in code) and the pre-specified patient-level test (`PREREG_phaseB.md`).
5. `covarnet_rerun/prepare.py`, then the unchanged original runner `run_immune_CM.R` (immune atlas code_v5; SHA256 `d2abc9db…`) with CoVarNet R source (QiangShiPKU/CoVarNet `e59b314`), then `sensitivity.R` (top-N=3; fixed K=4), then `summarise.py`.
6. `bulk_v2/gene_specificity.py`: fixed single-cell arm-gene rule (`PREREG_bulk_v2.md`).
7. `bulk_v2/fetch_cbioportal.py` (public cBioPortal API, selected genes only) and `bulk_v2/fetch_su2c_member.py` (HTTP-range fetch of one Zenodo member, CRC32 and SHA256 checked against the original receipt).
8. `bulk_v2/bulk_v2_analysis.py` and `bulk_v2/sc_longitudinal.py`.
9. `figures/make_figs.py`, `figures/make_supp_tables.py` and `figures/md2docx.py`.

## Validated environment

- Python 3.11: scanpy 1.11.5, anndata 0.11.4, celltypist 1.7.1, harmonypy, statsmodels, zarr 2.x.
- R 4.3.3 with NMF 0.28.
- `md2docx.py` needs python-docx, installed into an isolated folder.

Pre-registrations and the public aggregate tables are in [inventories/publication/2026-11-27/preprint_v1](../../inventories/publication/2026-11-27/preprint_v1/README.md).
