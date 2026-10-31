#!/bin/bash -l
#SBATCH --job-name=gap3-fingerprint
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,CODE=<repo>,PC_ROOT=<root>,OUT=<dir> cell_fingerprint_duplicates.sbatch
set -euo pipefail
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l) start $(date -u +%FT%TZ)"
V=$PC_ROOT/views/current
mkdir -p "${OUT:?}"
cut -f1 "$CODE/inventories/cancer_type/cancer_type_map_v1.tsv" | tail -n +2 > "$OUT/sources.txt"
"$PC_ROOT/envs/conda/scanpy/bin/python" -I -W ignore "$CODE/scripts/data_gaps/cell_fingerprint_duplicates.py" "$OUT/sources.txt" "$OUT" \
  core=$V/core_fullgene_v1.zarr:gse_id:X \
  immune=$V/immune_atlas_v1.zarr:dataset:layers/counts \
  public=$V/public_preprocessed_union_v1.zarr:dataset:layers/counts
echo "end $(date -u +%FT%TZ)"
