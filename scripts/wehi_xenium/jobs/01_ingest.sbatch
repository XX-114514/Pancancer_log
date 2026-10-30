#!/bin/bash -l
#SBATCH --job-name=xenium-ingest
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%A_%a.out
# Array index = 1-based row of samples.tsv; always pass --array explicitly at submit time.
set -euo pipefail
source "${XENIUM_CODE_ROOT:?export XENIUM_CODE_ROOT}/config.sh"
echo "repo_commit: $(git -C "$REPO_ROOT" rev-parse --short=12 HEAD) dirty: $(git -C "$REPO_ROOT" status --porcelain -- scripts/wehi_xenium | wc -l) job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} start: $(date -Is)"
"$ENV_PREFIX/bin/python" -I "$XENIUM_CODE_ROOT/ingest_xenium.py" "$SAMPLES_TSV" "$SLURM_ARRAY_TASK_ID" "$OUT_ROOT"
echo "end: $(date -Is)"
