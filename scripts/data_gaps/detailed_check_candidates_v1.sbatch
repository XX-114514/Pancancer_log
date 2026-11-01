#!/bin/bash -l
#SBATCH --job-name=src-detail
#SBATCH --partition=regular
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,CODE=<repo>,PC_ROOT=<root>,OUT=<dir>,PRIO=A detailed_check_candidates_v1.sbatch
set -euo pipefail
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l) start $(date -u +%FT%TZ)"
"$PC_ROOT/envs/conda/scanpy/bin/python" -I -u "$CODE/scripts/data_gaps/detailed_check_candidates_v1.py" \
  "$CODE/inventories/source_search/candidates_curated_v1.tsv" "$PC_ROOT/analysis/source_search_v1/known_accessions.txt" "${OUT:?}" "${PRIO:-A}"
echo "end $(date -u +%FT%TZ)"
