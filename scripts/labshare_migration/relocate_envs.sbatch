#!/bin/bash -l
#SBATCH --job-name=labshare-envs
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,CODE=<repo>/scripts/labshare_migration,SRC_PC=<scratch pancancer>,DST_PC=<labshare pancancer>,LOG_DIR=<dir> relocate_envs.sbatch
# Copies the conda prefixes and the py_magic venv from scratch to DST_PC/envs/*.relocated, checks the
# copy (verify.py), then rewrites SRC_PC -> DST_PC inside them and inside DST_PC/envs/rlib
# (relocate_prefix.py). Nothing is swapped into place here; the scratch envs are not touched.
set -euo pipefail
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) dirty $(git -C "$CODE" status --porcelain | wc -l)"
for d in conda py_magic_venv; do
  [ -e "${DST_PC:?}/envs/$d.relocated" ] && { echo "exists $DST_PC/envs/$d.relocated" >&2; exit 2; }
  rsync -aH "${SRC_PC:?}/envs/$d/" "$DST_PC/envs/$d.relocated/"
  python3 -I "$CODE/verify.py" "$SRC_PC/envs/$d" "$DST_PC/envs/$d.relocated" --jobs "$SLURM_CPUS_PER_TASK" --out "${LOG_DIR:?}/envs_$d"
done
echo "copy verified $(date -u +%FT%TZ)"
for d in conda.relocated py_magic_venv.relocated rlib; do
  python3 -I "$CODE/relocate_prefix.py" "$DST_PC/envs/$d" "$SRC_PC" "$DST_PC" --report "$LOG_DIR/relocate_${d%.relocated}.tsv"
done
echo "relocated $(date -u +%FT%TZ)"
