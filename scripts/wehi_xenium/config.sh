# Source this file from jobs and interactive shells. Real server paths live in an untracked
# local file so the public repository holds no host paths; see README.md.
XENIUM_CODE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
XENIUM_LOCAL_CONFIG="${XENIUM_LOCAL_CONFIG:-${HOME}/.config/pancancer/wehi_xenium.local.sh}"
if [[ -f "${XENIUM_LOCAL_CONFIG}" ]]; then
  # shellcheck source=/dev/null
  source "${XENIUM_LOCAL_CONFIG}"
fi
: "${XENIUM_DATA_ROOT:?set XENIUM_DATA_ROOT in ${XENIUM_LOCAL_CONFIG}}"
: "${PC_SPATIAL_ENV:?set PC_SPATIAL_ENV in ${XENIUM_LOCAL_CONFIG}}"
RAW_ROOT="${XENIUM_RAW_ROOT:-${XENIUM_DATA_ROOT}/raw}"
OUT_ROOT="${XENIUM_OUT_ROOT:-${XENIUM_DATA_ROOT}/processed}"
META_ROOT="${XENIUM_META_ROOT:-${XENIUM_DATA_ROOT}/metadata}"
SAMPLES_TSV="${META_ROOT}/samples.tsv"
ENV_PREFIX="${PC_SPATIAL_ENV}"
REPO_ROOT="$(cd "${XENIUM_CODE_ROOT}/../.." && pwd)"
# Frozen panel definitions are read from the repository record, not copied.
FROZEN_PLAN="${REPO_ROOT}/inventories/analysis/2026-10-18/priority_immune_CM_bulk_v1/frozen_plan_before_effects.json"
export XENIUM_CODE_ROOT RAW_ROOT OUT_ROOT META_ROOT SAMPLES_TSV ENV_PREFIX REPO_ROOT FROZEN_PLAN
