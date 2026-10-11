#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
task_python="${PYTHON:-python3}"
"${task_python}" -c 'import sys; assert sys.version_info >= (3, 10), "Python >=3.10 required"'
PYTHONDONTWRITEBYTECODE=1 "${task_python}" -m unittest discover -s tests -v
"${task_python}" scripts/validate_runtime_snapshot.py inventories/runtime/2026-10-11/runtime_snapshot_v2.json
bash scripts/validate_repository.sh
bash scripts/validate_releases.sh
git diff --check
