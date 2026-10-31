#!/usr/bin/env bash
# Run a command with old-server absolute paths visible at their migrated WEHI locations.
# Usage: wehi_oldpath_sandbox.sh [--gpu] [--rw PATH]... -- CMD [ARGS...]
# Builds a bwrap root from the host's top-level dirs, then binds every PATH_MAP.tsv prefix
# read-only, so frozen configs that name old-server absolute paths resolve unchanged.
# The PATH_MAP row whose target is "(rebuild)" gives the old env prefix; each <prefix>/<env>
# maps to the rebuilt ${PC_ROOT}/envs/conda/<env>. Old paths live only in the local PATH_MAP.tsv.
# PATH gains shims/ so `quota -g <old group>` reports the caller's real WEHI user quota.
# --rw PATH (an output directory or existing file under ${PC_ROOT}/Project_v5) is also made
# writable at its old-server path; at its own path it is writable already.
# --gpu passes the host /dev and /sys through (NVIDIA device nodes, driver info for CUDA).
set -euo pipefail
: "${PC_ROOT:?source ~/.config/pancancer/pancancer.local.sh}"
RW=(); GPU=0
while [ $# -gt 0 ] && [ "$1" != -- ]; do
  case $1 in --rw) RW+=("$(readlink -f "$2")"); shift 2;; --gpu) GPU=1; shift;; *) echo "unknown $1" >&2; exit 2;; esac
done
shift
if [ "$GPU" = 1 ]; then
  args=(--die-with-parent --proc /proc --dev-bind /dev /dev --ro-bind /sys /sys --tmpfs /tmp)
else
  args=(--die-with-parent --proc /proc --dev /dev --tmpfs /tmp)
fi
for d in /usr /etc /opt /var /run /vast /stornext /bin /sbin /lib /lib64 /tmp; do
  [ -e "$d" ] && args+=(--bind "$d" "$d")
done
args+=(--bind "$(readlink -f "$HOME")" "$HOME")
ENV_PREFIX=; V5_ALIASES=()
while IFS=$'\t' read -r old new _; do
  [ "$old" = old_prefix ] && continue
  [ "$new" = "(rebuild)" ] && ENV_PREFIX=$old
  [ "$new" = "$PC_ROOT/Project_v5" ] && V5_ALIASES+=("$old")
  [ "${new:0:1}" = / ] && [ -e "$new" ] && args+=(--ro-bind "$new" "$old")
done < "$PC_ROOT/PATH_MAP.tsv"
if [ -n "$ENV_PREFIX" ]; then
  for e in "$PC_ROOT"/envs/conda/*/; do
    e=${e%/}; args+=(--ro-bind "$e" "$ENV_PREFIX/$(basename "$e")")
  done
fi
for d in "${RW[@]}"; do
  # $d stays writable at its own path through the /vast bind; binding it there again would make
  # it a separate mount, and df would then report a sub-path that no quota line matches.
  case $d in "$PC_ROOT"/Project_v5/*)
    for a in "${V5_ALIASES[@]}"; do args+=(--bind "$d" "$a/${d#"$PC_ROOT"/Project_v5/}"); done;;
  esac
done
# quota shim (see shims/quota): old-server group quota query -> real WEHI user quota
args+=(--setenv PATH "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/shims:$PATH")
exec bwrap "${args[@]}" -- "$@"
