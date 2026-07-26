#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
today="${PROJECT_RECORD_DATE:-$(date +%F)}"
year="${today%%-*}"
template="${repo_root}/templates/daily_log_template.md"
destination="${repo_root}/logs/daily/${year}/${today}.md"

if [[ ! -f "${template}" ]]; then
  printf 'ERROR: template not found: %s\n' "${template}" >&2
  exit 1
fi

mkdir -p "$(dirname "${destination}")"
if [[ -e "${destination}" ]]; then
  printf 'EXISTS: %s\n' "${destination}"
  exit 0
fi

sed "s/YYYY-MM-DD/${today}/g" "${template}" > "${destination}"
printf '%s\n' "${destination}"
