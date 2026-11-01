#!/bin/bash -l
#SBATCH --job-name=src-ingest3
#SBATCH --partition=regular
#SBATCH --cpus-per-task=4
#SBATCH --mem=110G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,CODE=<repo>,PC_ROOT=<root>,D=<ingest v3 dir>[,SOURCES="..."] ingest_supplement_v3.sbatch
set -euo pipefail
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l) start $(date -u +%FT%TZ)"
"$PC_ROOT/envs/conda/scanpy/bin/python" -I -u -W ignore "$CODE/scripts/source_ingest/ingest_supplement_v3.py" "$D/raw" "$D/per_source_raw" \
  "$PC_ROOT/analysis/source_search_v3/detailed_v3" ${SOURCES:-}
echo "end $(date -u +%FT%TZ)"
