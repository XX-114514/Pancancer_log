#!/usr/bin/env bash
# WEHI: download one item's chunks as they arrive, verify each against the old server's chunk
# SHA256, then delete it from the bucket and leave a received/ marker (keeps the bucket under the
# budget cap). After the last chunk: unpack into DEST_PARENT (never overwrites), verify per-file
# SHA256 against the old server's list, write a receipt and remove the item from the bucket.
# A chunk that fails verification is deleted on both sides without a marker, so it is re-uploaded.
# Usage: wehi_fetch_verify.sh ITEM_ID DEST_PARENT
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
item="$1"; dest_parent="$(readlink -m "$2")"
src="${MIG_REMOTE}:${MIG_BUCKET}/${item}"
st="${MIG_STAGE}/${item}"; mkdir -p "$st/chunks"
exec > >(tee -a "$st/fetch.log") 2>&1
rc() { "$MIG_RCLONE" "$@" --retries 20 --low-level-retries 50; }

mig_log "$item waiting for manifests"
until rc lsf "$src/chunks.sha256" 2>/dev/null | grep -q .; do sleep "${MIG_POLL:-300}"; done
for f in meta.tsv files.sha256 chunks.sha256; do rc copyto "$src/$f" "$st/$f"; done
name="$(awk -F'\t' '$1=="source_name"{print $2}' "$st/meta.tsv")"
want="$(awk -F'\t' '$1=="n_chunks"{print $2}' "$st/meta.tsv")"
[[ -n "$name" ]] || { echo "bad meta.tsv" >&2; exit 1; }
[[ ! -e "$dest_parent/$name" ]] || { echo "destination exists: $dest_parent/$name" >&2; exit 1; }
touch "$st/verified.txt"

while (( $(wc -l < "$st/verified.txt") < want )); do
  rc lsf "$src/chunks/" 2>/dev/null | LC_ALL=C sort > "$st/remote.txt" || true
  if [[ ! -s "$st/remote.txt" ]]; then sleep "${MIG_POLL:-300}"; continue; fi
  rc copy "$src/chunks" "$st/chunks" --files-from "$st/remote.txt" --no-traverse \
    --transfers "${MIG_TRANSFERS:-16}" --stats 10m --stats-one-line --stats-log-level NOTICE
  while read -r c; do
    expect="$(awk -v c="$c" 'substr($0,67)==c{print $1}' "$st/chunks.sha256")"
    got="$(sha256sum "$st/chunks/$c" | cut -d' ' -f1)"
    if [[ -n "$expect" && "$expect" == "$got" ]]; then
      rc touch "$src/received/$c"
      rc deletefile "$src/chunks/$c"
      grep -qxF "$c" "$st/verified.txt" || echo "$c" >> "$st/verified.txt"
    else
      mig_log "$item chunk $c failed SHA256; deleting for re-upload"
      rm -f "$st/chunks/$c"; rc deletefile "$src/chunks/$c" || true
    fi
  done < "$st/remote.txt"
  mig_log "$item verified $(wc -l < "$st/verified.txt")/$want chunks"
done

mig_log "$item all chunks verified; re-checking local chunk set"
(cd "$st/chunks" && sha256sum --quiet -c ../chunks.sha256)
mkdir -p "$dest_parent"
mig_log "$item unpacking to $dest_parent/$name"
find "$st/chunks" -type f -name "${item}.tar.*" -print0 | LC_ALL=C sort -z | xargs -0 cat | tar -C "$dest_parent" -xf -
mig_log "$item verifying files"
(cd "$dest_parent" && find "$name" -type f -print0 | mig_sha256_parallel) > "$st/files.wehi.sha256"
if ! cmp -s "$st/files.sha256" "$st/files.wehi.sha256"; then
  diff "$st/files.sha256" "$st/files.wehi.sha256" | head -20 >&2
  echo "per-file SHA256 mismatch" >&2; exit 1
fi
{
  cat "$st/meta.tsv"
  printf 'verified_utc\t%s\n' "$(date -u +%FT%TZ)"
  printf 'files_sha256_of_list_wehi\t%s\n' "$(sha256sum < "$st/files.wehi.sha256" | cut -d' ' -f1)"
  printf 'status\tVERIFIED\n'
} > "$st/receipt.tsv"
mig_log "$item VERIFIED; receipt $st/receipt.tsv"
# Wait for the uploader's completion marker, then clear the item (markers, manifests) from the bucket.
until rc lsf "$src/_UPLOADED" 2>/dev/null | grep -q .; do sleep "${MIG_POLL:-300}"; done
rc purge "$src" && mig_log "$item purged from bucket"
if [[ "${MIG_KEEP_CHUNKS:-0}" != 1 ]]; then rm -rf "$st/chunks"; fi
