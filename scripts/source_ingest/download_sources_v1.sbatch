#!/bin/bash -l
#SBATCH --job-name=src-download
#SBATCH --partition=regular
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=6:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,D=<ingest dir with download_manifest.tsv> download_sources_v1.sbatch
# Downloads each manifest file into D/raw/<source>/ (read-only afterwards), checks the byte count
# against Content-Length when known, and writes SHA256SUMS per source.
set -euo pipefail
cd "${D:?}"
tail -n +2 download_manifest.tsv | while IFS=$'\t' read -r src name url bytes; do
  mkdir -p "raw/$src"; f="raw/$src/$name"
  if [ -s "$f" ] && { [ "$bytes" = 0 ] || [ "$(stat -c %s "$f")" = "$bytes" ]; }; then continue; fi
  # -f makes HTTP errors fail (an error page otherwise counts as success); size is checked per try
  for i in 1 2 3 4 5 6; do
    if curl -fsS -L --retry 3 --max-time 7200 -o "$f.part" "$url" && \
       { [ "$bytes" = 0 ] || [ "$(stat -c %s "$f.part")" = "$bytes" ]; }; then mv "$f.part" "$f"; break; fi
    rm -f "$f.part"; sleep $((2**i))
  done
  got=$(stat -c %s "$f"); [ "$bytes" = 0 ] || [ "$got" = "$bytes" ] || { echo "SIZE_MISMATCH $f $got $bytes"; exit 1; }
  echo "ok $src $name $got"
done
for s in raw/*/; do (cd "$s" && sha256sum * > SHA256SUMS && chmod a-w * ); done
echo "done $(date -u +%FT%TZ)"
