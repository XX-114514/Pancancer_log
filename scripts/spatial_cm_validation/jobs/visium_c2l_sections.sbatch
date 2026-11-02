#!/bin/bash -l
#SBATCH --job-name=cmsp-c2ls
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:1
#SBATCH --constraint=A30|A100
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-18
#SBATCH --output=%x-%A_%a.out
# Per-section full-batch cell2location spatial model (reference signatures from visium_c2l.sbatch).
# Usage: sbatch --export=ALL,CODE=<repo>/scripts/spatial_cm_validation,RUN=<run dir> jobs/visium_c2l_sections.sbatch
set -euo pipefail
set -a; source ~/.config/pancancer/pancancer.local.sh; source ~/.config/pancancer/wehi_xenium.local.sh; set +a
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l)"
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
"$PC_ROOT/envs/conda/c2l/bin/python" "$CODE/visium_cell2location_v1.py" --stage spatial --reference "${RUN:?}/reference/atlas_reference_v1.h5ad" \
  --visium-root "$XENIUM_DATA_ROOT/raw/V01_GSE273952_Visium_ccRCC" --output "$RUN/visium_c2l" --section-index "$SLURM_ARRAY_TASK_ID"
