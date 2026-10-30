#!/usr/bin/env bash
# Old server: hash every file of one item, pack it into fixed-size tar chunks and upload them.
# Usage: old_pack_upload.sh ITEM_ID SOURCE_PATH
# Restartable: finished hashing/packing steps are kept; rclone skips chunks already uploaded.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
item="$1"; src="$(readlink -f "$2")"
[[ -e "$src" ]] || { echo "missing source: $src" >&2; exit 1; }
parent="$(dirname "$src")"; name="$(basename "$src")"
out="${MIG_STAGE}/${item}"; mkdir -p "$out/chunks"
dest="${MIG_REMOTE}:${MIG_BUCKET}/${item}"
exec > >(tee -a "$out/pack.log") 2>&1
mig_log "$item start source=$name"

if [[ ! -s "$out/files.sha256" ]]; then
  # Links would be packed as links and not verified by content; record them and stop.
  (cd "$parent" && find "$name" -type l) > "$out/symlinks.txt"
  if [[ -s "$out/symlinks.txt" ]]; then
    echo "symlinks present, see $out/symlinks.txt; resolve before packing" >&2; exit 1
  fi
  mig_log "$item hashing files"
  (cd "$parent" && find "$name" -type f -print0 | mig_sha256_parallel) > "$out/files.sha256.tmp"
  mv "$out/files.sha256.tmp" "$out/files.sha256"
fi

if [[ ! -f "$out/chunks.done" ]]; then
  mig_log "$item packing"
  rm -f "$out/chunks/"*
  tar -C "$parent" -cf - "$name" | split -b "${MIG_CHUNK:-1G}" -d -a 5 - "$out/chunks/${item}.tar."
  (cd "$out/chunks" && find . -type f -name "${item}.tar.*" -printf '%f\0' | mig_sha256_parallel) > "$out/chunks.sha256"
  {
    printf 'item\t%s\n' "$item"
    printf 'source_name\t%s\n' "$name"
    printf 'n_files\t%s\n' "$(wc -l < "$out/files.sha256")"
    printf 'file_bytes\t%s\n' "$(cd "$parent" && find "$name" -type f -printf '%s\n' | awk '{s+=$1} END{print s}')"
    printf 'n_chunks\t%s\n' "$(wc -l < "$out/chunks.sha256")"
    printf 'chunk_size\t%s\n' "${MIG_CHUNK:-1G}"
    printf 'files_sha256_of_list\t%s\n' "$(sha256sum < "$out/files.sha256" | cut -d' ' -f1)"
    printf 'packed_utc\t%s\n' "$(date -u +%FT%TZ)"
  } > "$out/meta.tsv"
  touch "$out/chunks.done"
fi

# Manifests first, so the WEHI side can verify chunks as soon as they land.
for f in meta.tsv files.sha256 chunks.sha256; do
  "$MIG_RCLONE" copyto "$out/$f" "$dest/$f" --retries 20
done
mig_log "$item uploading $(wc -l < "$out/chunks.sha256") chunks"
"$MIG_RCLONE" copy "$out/chunks" "$dest/chunks" \
  --transfers "${MIG_TRANSFERS:-8}" --s3-upload-concurrency "${MIG_PART_CONCURRENCY:-4}" \
  --s3-chunk-size "${MIG_PART_SIZE:-16M}" --retries 50 --low-level-retries 100 \
  --stats 10m --stats-one-line --stats-log-level NOTICE
date -u +%FT%TZ | "$MIG_RCLONE" rcat "$dest/_UPLOADED"
mig_log "$item upload complete"
