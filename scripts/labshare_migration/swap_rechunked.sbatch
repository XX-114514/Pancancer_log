#!/bin/bash -l
#SBATCH --job-name=swap-rechunked
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,CODE=<repo>/scripts/labshare_migration,NEW=<rechunked zarr>,TARGET=<labshare zarr to replace>,LOG_DIR=<dir> swap_rechunked.sbatch
# Copies the verified rechunked store next to TARGET, checks the copy (verify.py), then deletes the
# old small-chunk TARGET and moves the new store into its place. The scratch original stays as is.
set -euo pipefail
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l)"
TMP=${TARGET:?}.rechunked
[ -e "$TMP" ] && { echo "exists $TMP" >&2; exit 2; }
rsync -a "${NEW:?}/" "$TMP/"
python3 -I "$CODE/verify.py" "$NEW" "$TMP" --jobs "$SLURM_CPUS_PER_TASK" --out "${LOG_DIR:?}/$(basename "$TARGET").rechunked"
echo "copy verified $(date -u +%FT%TZ)"
rm -rf "$TARGET"
mv "$TMP" "$TARGET"
echo "swapped $(date -u +%FT%TZ) files_now $(find "$TARGET" -type f | wc -l)"
