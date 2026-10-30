#!/usr/bin/env bash
# Old server: collect M0 (code, configs, run manifests, small text outputs, private wheel) into a
# local tree with rsync, so it can go through old_pack_upload.sh like any other item.
# Usage: old_stage_m0.sh PROJECT_V5_ROOT [EXTRA_FILE ...]
# Arrays, checkpoints, figures and any file over MIG_M0_MAX_SIZE are left out.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
v5="$(readlink -f "$1")"; shift
tree="${MIG_STAGE}/m0_tree/Project_v5_records"
mkdir -p "$tree" "${MIG_STAGE}/m0_tree/extra"
rsync -a --prune-empty-dirs --max-size="${MIG_M0_MAX_SIZE:-20M}" \
  --include='/scripts/***' --include='/config/***' --include='/frozen/***' --include='/sources/***' \
  --include='/docs/***' --include='/tests/***' --include='/*.md' --include='/*.sh' \
  --include='/runs/' --include='/runs/*/' \
  --include='/runs/*/config/***' --include='/runs/*/provenance/***' --include='/runs/*/scripts/***' \
  --include='/runs/*/tests/***' --include='/runs/*/manifests/***' \
  --include='/runs/*/outputs/' \
  --include='/runs/*/outputs/*.json' --include='/runs/*/outputs/*.tsv' --include='/runs/*/outputs/*.csv' \
  --include='/runs/*/outputs/*.md' --include='/runs/*/outputs/*.txt' --include='/runs/*/outputs/*.yaml' \
  --exclude='*' "$v5/" "$tree/"
for f in "$@"; do cp -p "$f" "${MIG_STAGE}/m0_tree/extra/"; done
mig_log "M0 staged: $(find "${MIG_STAGE}/m0_tree" -type f | wc -l) files, $(du -sh --apparent-size "${MIG_STAGE}/m0_tree" | cut -f1)"
