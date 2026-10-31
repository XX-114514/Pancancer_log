#!/bin/bash -l
#SBATCH --job-name=replay-covarnet
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,MIG_CODE_ROOT=<repo>/scripts/wehi_migration,REPLAY=<new dir>,ORIG=<orig run>[,RENV=<env>] jobs/replay_covarnet.sbatch
# Replays one accepted CoVarNet fit on WEHI: the original R driver and frozen author code are
# copied unchanged; the original frozen R request is reused with its input/output paths mapped
# through PATH_MAP.tsv; input and author-code SHA256 are checked against the frozen acceptance;
# then every output file is compared with the original.
set -euo pipefail
set -a; source ~/.config/pancancer/pancancer.local.sh; source ~/.config/pancancer/wehi_xenium.local.sh; set +a
echo "repo_commit $(git -C "${MIG_CODE_ROOT:?}" rev-parse HEAD) dirty $(git -C "$MIG_CODE_ROOT" status --porcelain | wc -l)"
RS=$PC_ROOT/envs/conda/${RENV:-scrna_r_covarnet}/bin/Rscript  # conda-forge R 4.3.3 + exact CRAN versions
[ -e "${REPLAY:?}" ] && { echo "replay dir exists: $REPLAY" >&2; exit 2; }
mkdir -p "$REPLAY"/{config,scripts,provenance,outputs}
cp -p "${ORIG:?}"/scripts/*.R "$REPLAY/scripts/"
cd "$REPLAY"
"$PC_REPO_ENV/bin/python" - "$ORIG" "$REPLAY" "$PC_PATH_MAP" <<'PY'
import csv, json, hashlib, sys
orig, rep, pm = sys.argv[1:4]
maps = sorted(((r['old_prefix'], r['new_prefix']) for r in csv.DictReader(open(pm), delimiter='\t') if r['new_prefix'].startswith('/')), key=lambda x: -len(x[0]))
def m(p):
    for o, n in maps:
        if p.startswith(o): return n + p[len(o):]
    return p
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
spec = json.load(open(f'{orig}/config/spec_v1.json'))
req = json.load(open(f'{orig}/config/frozen_R_request_v1.json'))
for k in ('frequency', 'metadata'):
    req[k] = m(req[k])
req['author_code'] = f'{rep}/scripts/frozen_author_CoVarNet_v1.R'
req['output'] = f'{rep}/outputs'
json.dump(req, open(f'{rep}/config/frozen_R_request_v1.json', 'w'), indent=2)
assert sha(req['author_code']) == spec['executed_author_code_SHA256']
acc = m(spec['input_acceptance']['path'])
assert sha(acc) == spec['input_acceptance']['SHA256']
a = json.load(open(acc))
for k in ('frequency', 'metadata'):
    assert sha(req[k]) == a[k]['SHA256'], k
print('inputs and author code match frozen SHA256')
PY
start=$(date -u +%FT%TZ)
"$RS" --vanilla scripts/run_TME_CoVarNet_scheme_v1.R config/frozen_R_request_v1.json > provenance/R_stdout.txt 2> provenance/R_stderr.txt
echo "R exit 0 start=$start end=$(date -u +%FT%TZ)"
"$PC_REPO_ENV/bin/python" "$MIG_CODE_ROOT/wehi_compare_outputs.py" "$ORIG/outputs" outputs provenance/compare_vs_original.tsv
