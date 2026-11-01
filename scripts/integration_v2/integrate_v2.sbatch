#!/bin/bash -l
#SBATCH --job-name=integrate-v2
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=480G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,CODE=<repo>,PC_ROOT=<root>,D=<decontX dir>,OUT=<dir> integrate_v2.sbatch
set -euo pipefail
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l) host $(hostname) start $(date -u +%FT%TZ)"
export OMP_NUM_THREADS=16 NUMBA_NUM_THREADS=16
"$PC_ROOT/envs/conda/scanpy/bin/python" -I -u -W ignore "$CODE/scripts/integration_v2/integrate_v2.py" "${D:?}" "${OUT:?}" "$CODE" 2026110102
echo "end $(date -u +%FT%TZ)"
