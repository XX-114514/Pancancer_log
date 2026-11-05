#!/bin/bash -l
#SBATCH --job-name=marker-routed-tme-v1
#SBATCH --partition=regular
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=AVX512
#SBATCH --output=marker-routed-%A_%a.out
# Pass --array=0-N%4 using the number of rows in SOURCE_LIST.
set -euo pipefail
: "${PC_ROOT:?}" "${CODE:?}" "${WORKER:?}" "${SOURCE_LIST:?}" "${OUT:?}"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMBA_NUM_THREADS=4
SHARD=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SOURCE_LIST")
test -n "$SHARD"
SRC=$(basename "$SHARD" .h5ad)
echo "start $(date -u +%FT%TZ) source $SRC code $(git -C "$CODE" rev-parse HEAD)"
sha256sum "$WORKER" "$(dirname "$WORKER")/marker_policy.py" "$(dirname "$WORKER")/singler_hpca.R"
"$PC_ROOT/envs/conda/scanpy/bin/python" -I -u "$WORKER" \
    --shard "$SHARD" --malignancy-dir "$PC_ROOT/malignancy_v1" \
    --models "$PC_ROOT/resources/celltypist/data/models" --repo "$CODE" \
    --rscript "$PC_ROOT/envs/conda/singler_bioc/bin/Rscript" \
    --r-cache "$PC_ROOT/resources/R_ExperimentHub" --out "$OUT/$SRC"
echo "end $(date -u +%FT%TZ)"
