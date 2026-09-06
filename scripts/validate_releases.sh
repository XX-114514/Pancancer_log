#!/usr/bin/env bash
set -uo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}" || exit 1
errors=0

check_manifest() {
  local manifest="$1"
  local base
  base="$(dirname "${manifest}")"
  while IFS=$'\t' read -r rel expected state; do
    [[ "${rel}" == "release_path" ]] && continue
    if [[ ! -f "${base}/${rel}" ]]; then
      printf 'FAIL missing: %s\n' "${base}/${rel}"
      errors=$((errors + 1))
      continue
    fi
    actual="$(sha256sum "${base}/${rel}" | awk '{print $1}')"
    if [[ "${actual}" != "${expected}" ]]; then
      printf 'FAIL hash: %s\n' "${base}/${rel}"
      errors=$((errors + 1))
    fi
  done < "${manifest}"
}

check_manifest releases/annotation_v7_20260905/SOURCE_PROVENANCE.tsv
check_manifest releases/communication_evidence_chain_20260906/SOURCE_PROVENANCE.tsv

if rg -n '/data[0-9]|/home/|USER[0-9]|xiongxin' releases \
  -g '!communication_evidence_chain_20260906/scripts/**' \
  -g '!annotation_v7_20260905/scripts/run_v7_pipeline_portable.sh'; then
  printf 'FAIL release contains internal absolute paths\n'
  errors=$((errors + 1))
else
  printf 'PASS no internal absolute paths in releasable files\n'
fi

python -m py_compile releases/annotation_v7_20260905/scripts/*.py || errors=$((errors + 1))

printf 'release validation errors: %d\n' "${errors}"
[[ "${errors}" -eq 0 ]]
