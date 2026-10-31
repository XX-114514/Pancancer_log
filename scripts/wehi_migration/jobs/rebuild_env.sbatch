#!/bin/bash -l
#SBATCH --job-name=env-rebuild
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%A_%a.out
# Usage: sbatch --array=0-3 --export=ALL,MIG_CODE_ROOT=<repo>/scripts/wehi_migration,PC_ROOT=<root> jobs/rebuild_env.sbatch
set -euo pipefail
ENVS=(scanpy scrna_r singler_bioc scvi)
E=${ENVS[${SLURM_ARRAY_TASK_ID:?use --array}]}
echo "repo_commit $(git -C "${MIG_CODE_ROOT:?export MIG_CODE_ROOT}" rev-parse HEAD) dirty $(git -C "$MIG_CODE_ROOT" status --porcelain | wc -l)"
bash "$MIG_CODE_ROOT/wehi_rebuild_env.sh" "${PC_ROOT:?}/envs/specs/$E" "$PC_ROOT/envs/conda/$E"
