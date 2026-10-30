#!/usr/bin/env bash
# Old server: stage M2 (the current working layer) as a hard-link tree, so packing needs no
# extra disk and the originals are never touched.
# Usage: old_stage_m2.sh PROJECT_V5 EXCLUDE_FILE [EXTRA_DIR ...]
#   Seeds are run directories modified since MIG_M2_SINCE (default 2026-10-29) plus
#   MIG_M2_SEEDS; old_m2_closure.py adds everything they reference under runs/.
#   EXCLUDE_FILE lists run names, or paths relative to runs/, left out (re-downloadable inputs,
#   or stores that travel as their own M1/M3 items). EXTRA_DIRs go under extra/.
# Waits while any process command line names a run in the closure, so running chains are not
# copied half-written. Symlinks are listed in symlinks.tsv and not packed.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
v5="$(readlink -f "$1")"; excl="$(readlink -f "$2")"; shift 2
runs="$v5/runs"; work="${MIG_STAGE}/m2_work"; tree="${MIG_STAGE}/m2_tree"
[[ ! -e "$tree" ]] || { echo "exists: $tree" >&2; exit 1; }
mkdir -p "$work"

mapfile -t seeds < <(cd "$runs" && find . -mindepth 1 -maxdepth 1 -type d -newermt "${MIG_M2_SINCE:-2026-10-29}" -printf '%f\n' | LC_ALL=C sort)
seeds+=(${MIG_M2_SEEDS:-})
mig_log "M2 closure from ${#seeds[@]} seed runs"
python3 "$MIG_CODE_ROOT/old_m2_closure.py" "$v5" "$work/closure" "${seeds[@]}" 2> "$work/closure.err"
tail -n +2 "$work/closure.runs.tsv" | cut -f1 | while read -r r; do [[ -d "$runs/$r" ]] && echo "$r"; done \
  | LC_ALL=C sort > "$work/closure_runs.txt"
grep -v '/' "$excl" | grep . | LC_ALL=C sort > "$work/excluded_runs.txt" || true
LC_ALL=C comm -23 "$work/closure_runs.txt" "$work/excluded_runs.txt" > "$work/selected_runs.txt"
mig_log "M2 selected $(wc -l < "$work/selected_runs.txt") of $(wc -l < "$work/closure_runs.txt") runs"

for ((i = 0; ; i++)); do
  busy="$(pgrep -af "Project_v5/runs/" | grep -v -e old_stage_m2 -e old_pack_upload -e old_m2_closure \
    | grep -oFf "$work/selected_runs.txt" | LC_ALL=C sort -u || true)"
  [[ -z "$busy" ]] && break
  (( i < ${MIG_M2_MAX_WAITS:-48} )) || { echo "runs still in use: $busy" >&2; exit 1; }
  mig_log "M2 waiting for running jobs in: $(echo "$busy" | tr '\n' ' ')"
  sleep 1800
done

mkdir -p "$tree/Project_v5_runs" "$tree/extra"
while read -r r; do cp -al "$runs/$r" "$tree/Project_v5_runs/$r"; done < "$work/selected_runs.txt"
{ grep '/' "$excl" || true; } | while read -r p; do
  [[ -n "$p" && "$p" != *..* ]] && rm -rf -- "$tree/Project_v5_runs/$p"
done
for d in "$@"; do cp -al "$(readlink -f "$d")" "$tree/extra/"; done
(cd "$tree" && find . -type l -printf '%p\t%l\n') > "$work/symlinks.tsv"
find "$tree" -type l -delete
cp "$work/closure.runs.tsv" "$work/selected_runs.txt" "$work/symlinks.tsv" "$excl" "$tree/"
mig_log "M2 staged: $(find "$tree" -type f | wc -l) files, $(du -sh --apparent-size "$tree" | cut -f1)"
