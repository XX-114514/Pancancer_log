#!/bin/bash -l
#SBATCH --job-name=cmsp-c2l
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:A100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,CODE=<repo>/scripts/spatial_cm_validation,RUN=<run dir> jobs/visium_c2l.sbatch
set -euo pipefail
set -a; source ~/.config/pancancer/pancancer.local.sh; source ~/.config/pancancer/wehi_xenium.local.sh; set +a
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l)"
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
"$PC_ROOT/envs/conda/c2l/bin/python" "$CODE/visium_cell2location_v1.py" --reference "${RUN:?}/reference/atlas_reference_v1.h5ad" \
  --visium-root "$XENIUM_DATA_ROOT/raw/V01_GSE273952_Visium_ccRCC" --output "$RUN/visium_c2l"
