#!/bin/bash -l
#SBATCH --job-name=xenium-stats
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%A_%a.out
# Array index = 1-based row of samples.tsv; always pass --array explicitly at submit time.
set -euo pipefail
source "${XENIUM_CODE_ROOT:?export XENIUM_CODE_ROOT}/config.sh"
echo "repo_commit: $(git -C "$REPO_ROOT" rev-parse --short=12 HEAD) dirty: $(git -C "$REPO_ROOT" status --porcelain -- scripts/wehi_xenium | wc -l) job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} start: $(date -Is)"
export OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4 MPLCONFIGDIR="${TMPDIR:-/tmp}/mpl-$SLURM_JOB_ID"
"$ENV_PREFIX/bin/python" -I "$XENIUM_CODE_ROOT/section_stats.py" "$SAMPLES_TSV" "$SLURM_ARRAY_TASK_ID" "$OUT_ROOT" "$FROZEN_PLAN"
echo "end: $(date -Is)"
