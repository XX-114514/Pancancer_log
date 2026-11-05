#!/bin/bash -l
#SBATCH --job-name=celltypist72-v1
#SBATCH --partition=regular
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=AVX512
#SBATCH --array=0-71%4
#SBATCH --output=celltypist72-%A_%a.out
set -euo pipefail
: "${PC_ROOT:?}" "${CODE:?}" "${WORKER:?}" "${SOURCE_LIST:?}" "${OUT:?}"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMBA_NUM_THREADS=4
SHARD=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SOURCE_LIST")
test -n "$SHARD"
SRC=$(basename "$SHARD" .h5ad)
echo "start $(date -u +%FT%TZ) source $SRC code $(git -C "$CODE" rev-parse HEAD)"
sha256sum "$WORKER"
"$PC_ROOT/envs/conda/scanpy/bin/python" -I -u "$WORKER" \
  --shard "$SHARD" --models "$PC_ROOT/resources/celltypist/data/models" \
  --repo "$CODE" --out "$OUT/$SRC"
echo "end $(date -u +%FT%TZ)"
