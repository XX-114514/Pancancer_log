#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  printf 'Usage: %s "task name"\n' "$0" >&2
  exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
task_name="$*"
slug="$(printf '%s' "${task_name}" \
  | tr '[:upper:]_' '[:lower:]-' \
  | sed -E 's/[^a-z0-9-]+/-/g; s/-+/-/g; s/^-//; s/-$//')"

if [[ -z "${slug}" ]]; then
  slug="codex-task"
fi

today="${PROJECT_RECORD_DATE:-$(date +%F)}"
year="${today%%-*}"
template="${repo_root}/templates/codex_task_template.md"
destination="${repo_root}/logs/codex/${year}/${today}_${slug}.md"

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
  -e "s/TASK_NAME/${slug}/g" \
  -e "s/^- Date:$/- Date: ${today}/" \
  "${template}" > "${destination}"
printf '%s\n' "${destination}"
