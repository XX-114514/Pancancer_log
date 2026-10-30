#!/usr/bin/env bash
# Old server: build a staging tree for one item from a list of sources, for old_pack_upload.sh.
# Usage: old_stage_list.sh ITEM LIST.tsv
# Each LIST line: MODE<TAB>DEST<TAB>SOURCE, DEST relative to the tree.
#   link  hard-link copy: no extra disk, originals untouched; only for static data on the
#         staging filesystem.
#   copy  real copy: for data that may still change (a consistent snapshot for hashing) or that
#         lives on another filesystem.
# Symlinks are listed in symlinks.tsv and not packed; the list itself is kept in the tree.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
item="$1"; list="$(readlink -f "$2")"
tree="${MIG_STAGE}/${item}_tree"
[[ ! -e "$tree" ]] || { echo "exists: $tree" >&2; exit 1; }
mkdir -p "$tree"
while IFS=$'\t' read -r mode dest src; do
  [[ -z "${mode}" || "${mode}" == \#* ]] && continue
  if [[ "$dest" == /* || "$dest" == *..* || ! -e "$src" ]]; then
    echo "bad line: $mode $dest $src" >&2; exit 1
  fi
  mkdir -p "$tree/$(dirname "$dest")"
  case "$mode" in
    link) cp -al "$src" "$tree/$dest" ;;
    copy) cp -a "$src" "$tree/$dest" ;;
    *) echo "bad mode: $mode" >&2; exit 1 ;;
  esac
done < "$list"
(cd "$tree" && find . -type l -printf '%p\t%l\n') > "$tree/symlinks.tsv"
find "$tree" -type l -delete
cp "$list" "$tree/source_list.tsv"
mig_log "$item staged: $(find "$tree" -type f | wc -l) files, $(du -sh --apparent-size "$tree" | cut -f1)"
