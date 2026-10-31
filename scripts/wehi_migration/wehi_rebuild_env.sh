#!/usr/bin/env bash
# Rebuild one old-server conda env on WEHI from its exported specs, then report differences.
# Usage: wehi_rebuild_env.sh SPEC_DIR PREFIX
#   SPEC_DIR  e.g. ${PC_ROOT}/envs/specs/scrna_r (from old_export_env.sh)
#   PREFIX    target env path, e.g. ${PC_ROOT}/envs/conda/scrna_r
# Steps: conda explicit (no solve) -> PyPI pins missing from the env (--no-deps) ->
# copied non-CRAN R library (scrna_r) -> package diff reports in PREFIX/rebuild_report/.
# Versions are never upgraded; anything that cannot be installed is listed as a difference.
set -euo pipefail
SPEC=$1; PREFIX=$2
[ -f "$SPEC/conda_explicit.txt" ] || { echo "missing $SPEC/conda_explicit.txt" >&2; exit 2; }
[ -e "$PREFIX" ] && { echo "refusing to overwrite existing $PREFIX" >&2; exit 2; }

module load micromamba
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-/vast/scratch/users/$USER/micromamba}
export PYTHONNOUSERSITE=1
echo "start $(date -u +%FT%TZ) spec=$SPEC prefix=$PREFIX host=$(hostname)"

micromamba create -y -p "$PREFIX" --file "$SPEC/conda_explicit.txt"
REP=$PREFIX/rebuild_report; mkdir -p "$REP"
run() { micromamba run -p "$PREFIX" "$@"; }

if [ -f "$SPEC/pip_freeze.txt" ] && [ -x "$PREFIX/bin/python" ]; then
  run python -m pip freeze --all 2>/dev/null | grep '==' | tr 'A-Z_' 'a-z-' | sort -u > "$REP/pip_present.txt"
  # Install every PyPI pin whose exact name==version is not already importable (old envs
  # sometimes pip-overrode a conda build, so a present name is not enough).
  grep -v ' @ ' "$SPEC/pip_freeze.txt" | grep '==' | while IFS= read -r line; do
    grep -qx "$(echo "$line" | tr 'A-Z_' 'a-z-')" "$REP/pip_present.txt" || echo "$line"
  done > "$REP/pip_to_install.txt"
  if [ -s "$REP/pip_to_install.txt" ]; then
    while IFS= read -r req; do
      run python -m pip install --no-deps --no-cache-dir "$req" >> "$REP/pip_install.log" 2>&1 \
        || echo "$req" >> "$REP/pip_failed.txt"
    done < "$REP/pip_to_install.txt"
  fi
  run python -m pip freeze --all > "$REP/pip_freeze_wehi.txt" 2>/dev/null || true
  python3 - "$SPEC/pip_freeze.txt" "$REP/pip_freeze_wehi.txt" > "$REP/pip_diff.tsv" <<'EOF'
import sys
def load(p):
    d = {}
    for l in open(p):
        l = l.strip()
        if not l or l.startswith('#'):
            continue
        n, v = (l.split('==', 1) + [''])[:2] if '==' in l else (l.split(' @ ')[0], '@')
        d[n.lower().replace('_', '-')] = v
    return d
old, new = load(sys.argv[1]), load(sys.argv[2])
print('package\told\twehi')
for n in sorted(set(old) | set(new)):
    if old.get(n) != new.get(n) and not (old.get(n) == '@' and n in new):
        print(f'{n}\t{old.get(n, "")}\t{new.get(n, "")}')
EOF
fi

if [ -d "$SPEC/R_nonstandard_library" ] && [ -d "$PREFIX/lib/R/library" ]; then
  for p in "$SPEC"/R_nonstandard_library/*; do
    b=$(basename "$p")
    [ -e "$PREFIX/lib/R/library/$b" ] && echo "$b" >> "$REP/R_library_kept_conda.txt" || cp -a "$p" "$PREFIX/lib/R/library/"
  done
fi

if [ -x "$PREFIX/bin/Rscript" ] && [ -f "$SPEC/R_packages.tsv" ]; then
  run Rscript -e 'ip <- installed.packages()[, c("Package","Version")]; write.table(ip, "'"$REP"'/R_installed_wehi.tsv", sep="\t", quote=FALSE, row.names=FALSE)'
  python3 - "$SPEC/R_packages.tsv" "$REP/R_installed_wehi.tsv" > "$REP/R_diff.tsv" <<'EOF'
import csv, sys
old = {r['Package']: (r['Version'], r.get('origin', '')) for r in csv.DictReader(open(sys.argv[1]), delimiter='\t')}
new = {r['Package']: r['Version'] for r in csv.DictReader(open(sys.argv[2]), delimiter='\t')}
print('package\told\twehi\torigin')
for n in sorted(set(old) | set(new)):
    ov, org = old.get(n, ('', ''))
    if ov != new.get(n, ''):
        print(f'{n}\t{ov}\t{new.get(n, "")}\t{org}')
EOF
fi

for f in pip_diff.tsv R_diff.tsv; do [ -f "$REP/$f" ] && echo "$f: $(( $(wc -l < "$REP/$f") - 1 )) differences" || true; done
[ -f "$REP/pip_failed.txt" ] && echo "pip_failed.txt: $(wc -l < "$REP/pip_failed.txt") packages" || true
echo "end $(date -u +%FT%TZ)"
