#!/bin/bash -l
#SBATCH --job-name=integrate-tme-v1
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=480G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,CODE=<repo>,PC_ROOT=<root>,S=<v2 shard dir>,M=<malignancy dir>,OUT=<dir>[,TME_HARMONY=source:2] integrate_tme_v1.sbatch
set -euo pipefail
echo "harmony ${TME_HARMONY:-default}"; echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l) host $(hostname) start $(date -u +%FT%TZ)"
export OMP_NUM_THREADS=16 NUMBA_NUM_THREADS=16
"$PC_ROOT/envs/conda/scanpy/bin/python" -I -u -W ignore "$CODE/scripts/malignancy_v1/integrate_tme_v1.py" "${S:?}" "${M:?}" "${OUT:?}" "$CODE" 2026110201
echo "end $(date -u +%FT%TZ)"
