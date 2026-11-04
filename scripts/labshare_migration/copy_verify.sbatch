#!/bin/bash -l
#SBATCH --job-name=labshare-copy
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%A_%a.out
# Usage: sbatch --array=1-N%4 --export=ALL,CODE=<repo>/scripts/labshare_migration,ITEMS=<items.tsv>,SRC_ROOT=<scratch root>,DST_ROOT=<labshare root>,LOG_DIR=<dir> copy_verify.sbatch
# Copies one item (array index = data line of ITEMS) from scratch to labshare with rsync, then
# compares SHA256 of every file on both sides (verify.py). The source is never modified.
set -euo pipefail
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l)"
line=$(grep -v '^#' "${ITEMS:?}" | sed -n "${SLURM_ARRAY_TASK_ID}p")
IFS=$'\t' read -r ITEM REL EXCL <<< "$line"
SRC=${SRC_ROOT:?}/$REL
DST=${DST_ROOT:?}/$REL
echo "item $ITEM src $SRC dst $DST excludes ${EXCL:-none} start $(date -u +%FT%TZ)"
args=()
IFS=':' read -ra ex <<< "${EXCL:-}"
for e in "${ex[@]}"; do [ -n "$e" ] && args+=(--exclude="/$e"); done
mkdir -p "$(dirname "$DST")"
rsync -aH --no-inc-recursive "${args[@]}" "$SRC/" "$DST/"
echo "rsync done $(date -u +%FT%TZ)"
python3 "$CODE/verify.py" "$SRC" "$DST" --exclude "${EXCL:-}" --jobs "$SLURM_CPUS_PER_TASK" \
  --out "${LOG_DIR:?}/$ITEM"
echo "verify done $(date -u +%FT%TZ)"
