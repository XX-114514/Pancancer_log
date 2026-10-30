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

# Budget guard: chunks go up only while the whole bucket stays under MIG_BUCKET_CAP (default 9G,
# inside the R2 free storage tier). WEHI deletes each chunk after verifying it and leaves a
# zero-byte marker under received/, so a chunk is pending until it is in neither place.
chunk_bytes="$(numfmt --from=iec "${MIG_CHUNK:-1G}")"
cap_bytes="$(numfmt --from=iec "${MIG_BUCKET_CAP:-9G}")"
(( chunk_bytes < cap_bytes )) || { echo "MIG_CHUNK must be smaller than MIG_BUCKET_CAP" >&2; exit 1; }
cut -c67- "$out/chunks.sha256" | LC_ALL=C sort > "$out/all_chunks.txt"
mig_log "$item uploading $(wc -l < "$out/all_chunks.txt") chunks, bucket cap ${MIG_BUCKET_CAP:-9G}"
while :; do
  { "$MIG_RCLONE" lsf "$dest/received/" 2>/dev/null || true
    "$MIG_RCLONE" lsf "$dest/chunks/" 2>/dev/null || true; } | LC_ALL=C sort -u > "$out/present.txt"
  LC_ALL=C comm -23 "$out/all_chunks.txt" "$out/present.txt" > "$out/pending.txt"
  [[ -s "$out/pending.txt" ]] || break
  used="$("$MIG_RCLONE" size "${MIG_REMOTE}:${MIG_BUCKET}" --json | grep -o '"bytes":[0-9]*' | cut -d: -f2)"
  room=$(( (cap_bytes - ${used:-0}) / chunk_bytes ))
  if (( room < 1 )); then sleep "${MIG_POLL:-300}"; continue; fi
  head -n "$(( room < ${MIG_TRANSFERS:-8} ? room : ${MIG_TRANSFERS:-8} ))" "$out/pending.txt" > "$out/batch.txt"
  "$MIG_RCLONE" copy "$out/chunks" "$dest/chunks" --files-from "$out/batch.txt" --no-traverse \
    --transfers "${MIG_TRANSFERS:-8}" --s3-upload-concurrency "${MIG_PART_CONCURRENCY:-4}" \
    --s3-chunk-size "${MIG_PART_SIZE:-16M}" --retries 50 --low-level-retries 100 \
    --stats 10m --stats-one-line --stats-log-level NOTICE
done
date -u +%FT%TZ | "$MIG_RCLONE" rcat "$dest/_UPLOADED"
mig_log "$item upload complete"
