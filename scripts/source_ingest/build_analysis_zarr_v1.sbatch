#!/bin/bash -l
#SBATCH --job-name=build-zarr
#SBATCH --partition=regular
#SBATCH --cpus-per-task=4
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,CODE=<repo>,PC_ROOT=<root>,OUT=<dir> build_analysis_zarr_v1.sbatch
set -euo pipefail
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l) start $(date -u +%FT%TZ)"
"$PC_ROOT/envs/conda/scanpy/bin/python" -I -u -W ignore "$CODE/scripts/source_ingest/build_analysis_zarr_v1.py" \
  "$PC_ROOT/ingest/existing_set_v2/per_source" "$PC_ROOT/ingest/supplement_v1/per_source_qc" "$CODE" \
  "${OUT:?}/pancancer_analysis_ready_v1.zarr" "$OUT/report"
echo "end $(date -u +%FT%TZ)"
