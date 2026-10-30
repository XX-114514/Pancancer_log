# Source from every migration script. Real paths, bucket and remote names live in an untracked
# local file on each server so the public repository holds no host paths; see README.md.
MIG_CODE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MIG_LOCAL_CONFIG="${MIG_LOCAL_CONFIG:-${HOME}/.config/pancancer/migration.local.sh}"
if [[ -f "${MIG_LOCAL_CONFIG}" ]]; then
  # shellcheck source=/dev/null
  source "${MIG_LOCAL_CONFIG}"
fi
: "${MIG_STAGE:?set MIG_STAGE in ${MIG_LOCAL_CONFIG}}"
: "${MIG_BUCKET:?set MIG_BUCKET in ${MIG_LOCAL_CONFIG}}"
MIG_RCLONE="${MIG_RCLONE:-rclone}"
MIG_REMOTE="${MIG_REMOTE:-r2mig}"
MIG_JOBS="${MIG_JOBS:-8}"
# The old server's login shell exports a local proxy that is not always up; R2 is reached directly.
unset http_proxy https_proxy all_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY
export MIG_CODE_ROOT MIG_STAGE MIG_BUCKET MIG_RCLONE MIG_REMOTE MIG_JOBS

mig_log() { echo "$(date -u +%FT%TZ) $*"; }

# Parallel sha256sum of NUL-separated paths on stdin; each worker writes its own file so
# lines from different workers never interleave. Output is sorted by path.
mig_sha256_parallel() {
  local tmp; tmp="$(mktemp -d "${MIG_STAGE}/.sha.XXXXXX")"
  xargs -0 -n 256 -P "${MIG_JOBS}" sh -c 'sha256sum -- "$@" > "$0/part.$$.$RANDOM"' "$tmp"
  cat "$tmp"/part.* 2>/dev/null | LC_ALL=C sort -k2
  rm -rf "$tmp"
}
