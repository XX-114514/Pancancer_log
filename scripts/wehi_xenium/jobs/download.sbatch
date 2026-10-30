#!/bin/bash
#SBATCH --job-name=spatial-dl
#SBATCH --partition=datamover
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out
# Usage: sbatch download.sbatch MANIFEST.tsv
# Manifest columns: http_code  bytes  relative_dir  url  (bytes may be 0 = unknown size, accept any non-empty file)
# Downloads processed tables only; FASTQ / morphology images / transcripts are listed separately when needed.
set -uo pipefail
source "${XENIUM_CODE_ROOT:?export XENIUM_CODE_ROOT}/config.sh"
MANIFEST=${1:?manifest}
DEST="$RAW_ROOT"
mkdir -p "$DEST"
fail=0
while IFS=$'\t' read -r code bytes dir url; do
  [[ -z "$url" || "$code" == "http_code" ]] && continue
  name=$(basename "${url%%\?*}")
  out="$DEST/$dir/$name"
  mkdir -p "$DEST/$dir"
  if [[ -s "$out" && ( "$bytes" == 0 || $(stat -c%s "$out") == "$bytes" ) ]]; then echo "SKIP $dir/$name"; continue; fi
  for try in 1 2 3 4; do
    curl -sSL --fail -C - --retry 3 -o "$out" "$url" && break
    sleep $((2**try))
  done
  got=$(stat -c%s "$out" 2>/dev/null || echo 0)
  if [[ "$got" -gt 0 && ( "$bytes" == 0 || "$got" == "$bytes" ) ]]; then echo "OK   $dir/$name $got"
  else echo "FAIL $dir/$name expected=$bytes got=$got" >&2; fail=1; fi
done < "$MANIFEST"
# Checksums per source directory so later downloads do not rehash everything.
for d in $(cut -f3 "$MANIFEST" | cut -d/ -f1 | sort -u); do
  [[ -d "$DEST/$d" ]] || continue
  (cd "$DEST/$d" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 -r sha256sum > SHA256SUMS)
done
echo "done fail=$fail"; exit $fail
