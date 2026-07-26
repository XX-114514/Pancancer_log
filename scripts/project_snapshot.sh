#!/usr/bin/env bash
set -u

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
summary_target="${1:-}"

printf '# Project Snapshot\n\n'
printf -- '- Time: %s\n' "$(date --iso-8601=seconds)"
printf -- '- Host: %s\n' "$(hostname 2>/dev/null || printf 'unavailable')"
printf -- '- Record directory: ${RECORD_ROOT}\n'

if git -C "${repo_root}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  printf -- '- Git branch: %s\n' "$(git -C "${repo_root}" branch --show-current 2>/dev/null || printf 'unavailable')"
  printf -- '- Git commit: %s\n' "$(git -C "${repo_root}" rev-parse HEAD 2>/dev/null || printf 'unborn')"
  printf '\n## Git status\n\n```text\n'
  git -C "${repo_root}" status --short 2>/dev/null || true
  printf '```\n'
else
  printf -- '- Git: not initialized\n'
fi

printf '\n## Environments\n\n'
printf -- '- Conda environment: %s\n' "${CONDA_DEFAULT_ENV:-not-active}"
printf -- '- Python: %s\n' "$(python --version 2>&1 || printf 'unavailable')"
printf -- '- R: %s\n' "$(R --version 2>/dev/null | sed -n '1p' || printf 'unavailable')"

printf '\n## Current user SLURM jobs\n\n```text\n'
if command -v squeue >/dev/null 2>&1; then
  squeue --me --noheader --format='%.18i %.12P %.24j %.8T %.10M %.6D' 2>/dev/null || printf 'squeue query failed\n'
else
  printf 'squeue unavailable\n'
fi
printf '```\n'

if [[ -n "${summary_target}" ]]; then
  printf '\n## Result directory summary\n\n'
  printf -- '- Logical argument: ${RESULTS_ROOT}\n'
  if [[ -d "${summary_target}" ]]; then
    file_count="$(find "${summary_target}" -type f 2>/dev/null | wc -l | tr -d ' ')"
    size_summary="$(du -sh "${summary_target}" 2>/dev/null | awk '{print $1}')"
    printf -- '- File count: %s\n' "${file_count}"
    printf -- '- Total size: %s\n' "${size_summary:-unavailable}"
  else
    printf -- '- Status: directory not found or not readable\n'
  fi
fi

printf '\nReview and redact snapshot output before saving it in Git.\n'
