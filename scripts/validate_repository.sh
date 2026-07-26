#!/usr/bin/env bash
set -uo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
cd "${repo_root}" || exit 1

errors=0
warnings=0

pass() { printf 'PASS: %s\n' "$1"; }
warn() { printf 'WARN: %s\n' "$1"; warnings=$((warnings + 1)); }
fail() { printf 'FAIL: %s\n' "$1"; errors=$((errors + 1)); }

core_files=(
  README.md AGENTS.md STATUS.md ROADMAP.md TODO.md CHANGELOG.md .gitignore
  docs/PROJECT_OVERVIEW.md docs/DIRECTORY_STRUCTURE.md docs/WORKFLOW.md
  docs/NAMING_CONVENTIONS.md docs/GLOSSARY.md docs/DATA_POLICY.md
  docs/REPRODUCIBILITY.md methods/README.md decisions/README.md
  runs/README.md inventories/datasets.tsv inventories/runs.tsv
  inventories/artifacts.tsv
)

for path in "${core_files[@]}"; do
  if [[ -f "${path}" ]]; then
    pass "core file exists: ${path}"
  else
    fail "missing core file: ${path}"
  fi
done

while IFS= read -r md; do
  [[ -n "${md}" ]] || continue
  if [[ ! -s "${md}" ]]; then
    fail "empty Markdown file: ${md}"
  fi
done < <(find . -type f -name '*.md' -not -path './.git/*' -print)

placeholder_hits="$(grep -RInE 'TASK_NAME|RUN_ID|ADR-XXXX|Decision title|Incident: TITLE|待填写' \
  --include='*.md' \
  --exclude-dir='.git' \
  --exclude-dir='templates' \
  . 2>/dev/null || true)"
if [[ -n "${placeholder_hits}" ]]; then
  fail "unreplaced placeholders outside templates:\n${placeholder_hits}"
else
  pass "no obvious unreplaced placeholders outside templates"
fi

scan_files=()
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  while IFS= read -r file; do
    [[ -n "${file}" ]] && scan_files+=("${file}")
  done < <(git ls-files)
  while IFS= read -r file; do
    [[ -n "${file}" ]] && scan_files+=("${file}")
  done < <(git ls-files --others --exclude-standard)
else
  while IFS= read -r file; do
    scan_files+=("${file#./}")
  done < <(find . -type f -not -path './.git/*' -print)
fi

secret_hits=""
path_hits=""
for file in "${scan_files[@]}"; do
  [[ -f "${file}" ]] || continue
  [[ "${file}" == "scripts/validate_repository.sh" ]] && continue
  if LC_ALL=C grep -InE \
    '(BEGIN[[:space:]].*PRIVATE KEY|api[_-]?key[[:space:]]*[:=][[:space:]]*[^[:space:]]+|access[_-]?token[[:space:]]*[:=][[:space:]]*[^[:space:]]+|password[[:space:]]*[:=][[:space:]]*[^[:space:]]+|gh[pousr]_[A-Za-z0-9_]{20,})' \
    "${file}" >/dev/null 2>&1; then
    secret_hits+="${file}"$'\n'
  fi
  if LC_ALL=C grep -InE \
    '(/home/[A-Za-z0-9._-]+/|/data[0-9]*/[^ ${}`]+|ssh://|[A-Za-z0-9._-]+@[A-Za-z0-9._-]+:)' \
    "${file}" >/dev/null 2>&1; then
    path_hits+="${file}"$'\n'
  fi
done

if [[ -n "${secret_hits}" ]]; then
  fail "possible secrets detected in:\n${secret_hits}"
else
  pass "no obvious secret patterns detected"
fi

if [[ -n "${path_hits}" ]]; then
  warn "possible sensitive absolute/internal paths detected in:\n${path_hits}"
else
  pass "no common sensitive absolute/internal path patterns detected"
fi

large_files=""
disallowed_files=""
for file in "${scan_files[@]}"; do
  [[ -f "${file}" ]] || continue
  size="$(wc -c < "${file}" | tr -d ' ')"
  if [[ "${size}" -gt 5242880 ]]; then
    large_files+="${file} (${size} bytes)"$'\n'
  fi
  case "${file}" in
    *.h5ad|*.h5|*.hdf5|*.loom|*.rds|*.rda|*.pkl|*.pickle|*.parquet|*.feather|*.mtx|*.mtx.gz|*.bam|*.bai|*.cram|*.fastq|*.fastq.gz|*.fq|*.fq.gz)
      disallowed_files+="${file}"$'\n'
      ;;
  esac
done

if [[ -n "${large_files}" ]]; then
  fail "tracked/candidate files larger than 5 MiB:\n${large_files}"
else
  pass "no tracked/candidate file exceeds 5 MiB"
fi

if [[ -n "${disallowed_files}" ]]; then
  fail "disallowed data formats found:\n${disallowed_files}"
else
  pass "no disallowed data formats found"
fi

while IFS= read -r tsv; do
  [[ -n "${tsv}" ]] || continue
  if awk -F '\t' '
    NR == 1 { expected = NF; next }
    NF != expected { bad = 1; printf "%s:%d expected=%d actual=%d\n", FILENAME, NR, expected, NF }
    END { exit bad ? 1 : 0 }
  ' "${tsv}"; then
    pass "consistent TSV columns: ${tsv}"
  else
    fail "inconsistent TSV columns: ${tsv}"
  fi
done < <(find . -type f -name '*.tsv' -not -path './.git/*' -print)

broken_links=""
while IFS= read -r md; do
  [[ -n "${md}" ]] || continue
  while IFS= read -r target; do
    [[ -n "${target}" ]] || continue
    case "${target}" in
      http://*|https://*|mailto:*|\#*|\$\{*|/*) continue ;;
    esac
    clean_target="${target%%#*}"
    [[ -n "${clean_target}" ]] || continue
    if [[ ! -e "$(dirname "${md}")/${clean_target}" ]]; then
      broken_links+="${md} -> ${target}"$'\n'
    fi
  done < <(grep -oE '\]\([^)]+' "${md}" 2>/dev/null | sed 's/^](//' || true)
done < <(find . -type f -name '*.md' -not -path './.git/*' -print)

if [[ -n "${broken_links}" ]]; then
  fail "broken local relative links:\n${broken_links}"
else
  pass "no broken local relative Markdown links detected"
fi

printf '\nSUMMARY: %d error(s), %d warning(s)\n' "${errors}" "${warnings}"
if [[ "${errors}" -gt 0 ]]; then
  exit 1
fi
exit 0
