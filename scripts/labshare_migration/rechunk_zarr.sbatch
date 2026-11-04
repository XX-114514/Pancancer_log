#!/bin/bash -l
#SBATCH --job-name=rechunk-zarr
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%A_%a.out
# Usage: sbatch --array=1-N --export=ALL,CODE=<repo>/scripts/labshare_migration,PY=<python>,LIST=<tsv: name<TAB>src zarr>,OUT=<dir> rechunk_zarr.sbatch
# Writes OUT/<name>.zarr with ~16 MiB chunks and checks it against the source (rechunk_zarr.py).
set -euo pipefail
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l)"
IFS=$'\t' read -r NAME SRC <<< "$(grep -v '^#' "${LIST:?}" | sed -n "${SLURM_ARRAY_TASK_ID}p")"
echo "name $NAME src $SRC start $(date -u +%FT%TZ)"
"${PY:?}" -I -u "$CODE/rechunk_zarr.py" "$SRC" "${OUT:?}/$NAME.zarr" --jobs "$SLURM_CPUS_PER_TASK" --report "$OUT/$NAME.rechunk.tsv"
echo "end $(date -u +%FT%TZ)"
