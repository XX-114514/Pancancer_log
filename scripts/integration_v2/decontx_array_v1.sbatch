#!/bin/bash -l
#SBATCH --job-name=decontx
#SBATCH --partition=regular
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --constraint=AVX512
#SBATCH --output=%x-%A_%a.out
# Usage: sbatch --array=0-N --export=ALL,CODE=<repo>,PC_ROOT=<root>,LIST=<file of shard paths>,OUT=<dir> decontx_array_v1.sbatch
set -euo pipefail
SHARD=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${LIST:?}")
SRC=$(basename "$SHARD" .h5ad)
[ -s "$OUT/$SRC.h5ad" ] && { echo "done already $SRC"; exit 0; }
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) source $SRC start $(date -u +%FT%TZ)"
W=$OUT/work/$SRC; mkdir -p "$W"
"$PC_ROOT/envs/conda/scanpy/bin/python" -I -u -W ignore "$CODE/scripts/integration_v2/decontx_source_v1.py" "$SHARD" "$W" "$OUT/$SRC.h5ad.part" \
  "$PC_ROOT/envs/conda/scrna_r_full/bin/Rscript" "$CODE/scripts/integration_v2/decontx_library_worker_celda.R" 2026110101
mv "$OUT/$SRC.h5ad.part" "$OUT/$SRC.h5ad"
echo "end $(date -u +%FT%TZ)"
