#!/bin/bash -l
#SBATCH --job-name=cmsp
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=8:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,CODE=<repo>/scripts/spatial_cm_validation,ENV=<env>,OUTROOT=<run dir> jobs/run_in_sandbox.sbatch SCRIPT [ARGS...]
# Runs a Python script inside the old-path sandbox (frozen JSONs name old-server paths).
# OLD_V5 (old Project_v5 prefix) comes from the local PATH_MAP.tsv, never from the repo.
set -euo pipefail
set -a; source ~/.config/pancancer/pancancer.local.sh; source ~/.config/pancancer/wehi_xenium.local.sh; set +a
MIG=${CODE:?}/../wehi_migration
echo "repo_commit $(git -C "$CODE" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l)"
OLD_V5=$(awk -F'\t' -v n="$PC_ROOT/Project_v5" '$2==n{print $1; exit}' "$PC_PATH_MAP")
T=${SLURM_CPUS_PER_TASK:-2}
script=$1; shift
mkdir -p "${OUTROOT:?}"
"$MIG/wehi_oldpath_sandbox.sh" --rw "$OUTROOT" -- env OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T NUMBA_NUM_THREADS=$T PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 OLD_V5="$OLD_V5" \
  "$PC_ROOT/envs/conda/${ENV:?}/bin/python" -B "$CODE/$script" "$@"
