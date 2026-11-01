#!/bin/bash -l
#SBATCH --job-name=src-ingest
#SBATCH --partition=regular
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=8:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,CODE=<repo>,PC_ROOT=<root>,D=<ingest dir>[,SOURCES="GSE..."] ingest_supplement_v1.sbatch
set -euo pipefail
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l) start $(date -u +%FT%TZ)"
"$PC_ROOT/envs/conda/scanpy/bin/python" -I -u -W ignore "$CODE/scripts/source_ingest/ingest_supplement_v1.py" "$D/raw" "$D/per_source_raw" ${SOURCES:-}
echo "end $(date -u +%FT%TZ)"
