#!/bin/bash
#SBATCH -J v7_ref
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o %x-%j.out
set -euo pipefail
set -a; source ~/.config/pancancer/pancancer.local.sh; source ~/.config/pancancer/wehi_xenium.local.sh; set +a
PY=${PC_ROOT}/envs/conda/scanpy/bin/python
$PY -I $RUN/scripts/build_v7_spatial_reference.py $RUN ${PC_ROOT}/Project_v3/annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr ${XENIUM_DATA_ROOT}/processed/h5ad
