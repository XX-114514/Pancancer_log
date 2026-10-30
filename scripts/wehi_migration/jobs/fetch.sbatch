#!/bin/bash -l
#SBATCH --job-name=mig-fetch
#SBATCH --partition=datamover
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,MIG_ITEM=<id>,MIG_DEST=<parent dir> jobs/fetch.sbatch
# Submit with MIG_CODE_ROOT exported.
set -euo pipefail
module load rclone
echo "repo_commit $(git -C "${MIG_CODE_ROOT:?export MIG_CODE_ROOT}" rev-parse HEAD) dirty $(git -C "$MIG_CODE_ROOT" status --porcelain | wc -l)"
MIG_JOBS="${SLURM_CPUS_PER_TASK:-8}" bash "$MIG_CODE_ROOT/wehi_fetch_verify.sh" "${MIG_ITEM:?}" "${MIG_DEST:?}"
