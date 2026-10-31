#!/bin/bash -l
#SBATCH --job-name=replay-py
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch -c <threads> --export=ALL,MIG_CODE_ROOT=<repo>/scripts/wehi_migration,ORIG=<orig run>,REPLAY=<new dir>,ENTRY=<script>,ENV=<env> jobs/replay_python_stage.sbatch
# Re-executes one accepted Python stage of an old-server run on WEHI: config/, scripts/,
# resources/ and vendor/ are copied unchanged; old absolute paths resolve through
# wehi_oldpath_sandbox.sh; outputs are compared file by file with the original.
set -euo pipefail
set -a; source ~/.config/pancancer/pancancer.local.sh; source ~/.config/pancancer/wehi_xenium.local.sh; set +a
echo "repo_commit $(git -C "${MIG_CODE_ROOT:?}" rev-parse HEAD) dirty $(git -C "$MIG_CODE_ROOT" status --porcelain | wc -l)"
[ -e "${REPLAY:?}" ] && { echo "replay dir exists: $REPLAY" >&2; exit 2; }
mkdir -p "$REPLAY"/{checkpoints,outputs,provenance}
for d in config scripts resources vendor; do [ -d "$ORIG/$d" ] && cp -a "$ORIG/$d" "$REPLAY/"; done
chmod -R u+w "$REPLAY"
PY=$PC_ROOT/envs/conda/${ENV:?}/bin/python
T=${SLURM_CPUS_PER_TASK:-2}
start=$(date -u +%FT%TZ)
set +e
"$MIG_CODE_ROOT/wehi_oldpath_sandbox.sh" --rw "$REPLAY" -- env OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T NUMBA_NUM_THREADS=$T PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 \
  "$PY" -B "$REPLAY/scripts/${ENTRY:?}" --output "$REPLAY" > "$REPLAY/provenance/replay_stdout.txt" 2> "$REPLAY/provenance/replay_stderr.txt"
code=$?
set -e
echo "exit=$code start=$start end=$(date -u +%FT%TZ)"
tail -n 20 "$REPLAY/provenance/replay_stderr.txt"
[ $code -eq 0 ] || exit $code
"$PC_REPO_ENV/bin/python" "$MIG_CODE_ROOT/wehi_compare_outputs.py" "$ORIG/outputs" "$REPLAY/outputs" "$REPLAY/provenance/compare_vs_original.tsv"
