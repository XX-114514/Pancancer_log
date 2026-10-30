#!/usr/bin/env bash
# WEHI: download one item's chunks while they are uploaded, verify chunk and per-file SHA256,
# unpack into DEST_PARENT and write a receipt. Never overwrites an existing destination.
# Usage: wehi_fetch_verify.sh ITEM_ID DEST_PARENT
# MIG_PURGE_REMOTE=1 deletes the item from the bucket after full verification.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
item="$1"; dest_parent="$(readlink -m "$2")"
src="${MIG_REMOTE}:${MIG_BUCKET}/${item}"
st="${MIG_STAGE}/${item}"; mkdir -p "$st/chunks"
exec > >(tee -a "$st/fetch.log") 2>&1
pull() {
  "$MIG_RCLONE" copy "$src" "$st" --exclude '_UPLOADED' --transfers "${MIG_TRANSFERS:-16}" \
    --retries 20 --low-level-retries 50 --stats 10m --stats-one-line --stats-log-level NOTICE
}
mig_log "$item waiting for upload; polling every ${MIG_POLL:-600}s"
until "$MIG_RCLONE" lsf "$src/_UPLOADED" 2>/dev/null | grep -q .; do
  pull || true
  sleep "${MIG_POLL:-600}"
done
pull

name="$(awk -F'\t' '$1=="source_name"{print $2}' "$st/meta.tsv")"
want="$(awk -F'\t' '$1=="n_chunks"{print $2}' "$st/meta.tsv")"
have="$(find "$st/chunks" -type f -name "${item}.tar.*" | wc -l)"
[[ -n "$name" && "$want" == "$have" ]] || { echo "chunk count $have != $want" >&2; exit 1; }
mig_log "$item verifying $have chunks"
(cd "$st/chunks" && sha256sum --quiet -c ../chunks.sha256)

[[ ! -e "$dest_parent/$name" ]] || { echo "destination exists: $dest_parent/$name" >&2; exit 1; }
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
if [[ "${MIG_PURGE_REMOTE:-0}" == 1 ]]; then
  "$MIG_RCLONE" purge "$src" && mig_log "$item purged from bucket"
fi
if [[ "${MIG_KEEP_CHUNKS:-0}" != 1 ]]; then rm -rf "$st/chunks"; fi
