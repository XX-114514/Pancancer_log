#!/bin/bash -l
#SBATCH --job-name=xenium-panel
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
source "${XENIUM_CODE_ROOT:?export XENIUM_CODE_ROOT}/config.sh"
echo "repo_commit: $(git -C "$REPO_ROOT" rev-parse --short=12 HEAD) job: $SLURM_JOB_ID start: $(date -Is)"
"$ENV_PREFIX/bin/python" -I "$XENIUM_CODE_ROOT/panel_coverage.py" "$OUT_ROOT" "$FROZEN_PLAN"
echo "end: $(date -Is)"
