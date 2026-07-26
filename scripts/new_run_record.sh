#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
run_name="${*:-run}"
slug="$(printf '%s' "${run_name}" \
  | tr '[:upper:]_' '[:lower:]-' \
  | sed -E 's/[^a-z0-9-]+/-/g; s/-+/-/g; s/^-//; s/-$//')"

if [[ -z "${slug}" ]]; then
  slug="run"
fi

timestamp="${PROJECT_RECORD_TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
year="${timestamp:0:4}"
run_id="${timestamp}_${slug}"
template="${repo_root}/templates/run_record_template.md"
destination="${repo_root}/runs/${year}/${run_id}.md"

if [[ ! -f "${template}" ]]; then
  printf 'ERROR: template not found: %s\n' "${template}" >&2
  exit 1
fi

mkdir -p "$(dirname "${destination}")"
if [[ -e "${destination}" ]]; then
  printf 'EXISTS: %s\n' "${destination}"
  exit 0
fi

sed \
  -e "s/RUN_ID/${run_id}/g" \
  -e "s/^- Run ID:$/- Run ID: ${run_id}/" \
  "${template}" > "${destination}"
printf '%s\n' "${destination}"
