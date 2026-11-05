#!/usr/bin/env bash
set -euo pipefail
TASK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASE_R="${BLISA_BASE_R:-${PC_ROOT}/envs/conda/scrna_r_full/bin/Rscript}"
cd "$TASK_ROOT"
mkdir -p env/R-library logs outputs
export R_LIBS_USER="$TASK_ROOT/env/R-library"
export OMP_NUM_THREADS=2
exec "$BASE_R" --vanilla "$@"
