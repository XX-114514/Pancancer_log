import os
"""Prepare or submit the marker-routed 72-source workflow from a host terminal."""
import argparse
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import pwd
import re
import shutil
import subprocess
import sys

PARENTS = {'GSE125449': '30645863_0', 'GSE183904': '30645933_31',
           'GSE253013': '30645933_55', 'GSE269936': '30645933_59'}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--submit', action='store_true')
    p.add_argument('--run-dir', type=Path)
    p.add_argument('--output-root', type=Path)
    p.add_argument('--repo', type=Path, default=Path(os.environ['CODE_ROOT']))
    args = p.parse_args()
    task = Path(__file__).resolve().parents[1]
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    run = (args.run_dir or task / 'submissions' / stamp).resolve()
    shards = [Path(s) for s in (task / 'sources.txt').read_text().splitlines() if s]
    if len(shards) != 72 or len({s.stem for s in shards}) != 72:
        raise ValueError('Expected 72 distinct source shards')
    if not all(s.is_file() for s in shards):
        raise ValueError('One or more source shards are missing')
    pc = shards[0].parents[3]; malig = pc / 'malignancy_v1'
    py = pc / 'envs/conda/scanpy/bin/python'
    if args.submit:
        # No jobs or scratch directories are created unless host identity,
        # authentication and scheduler connectivity all pass first.
        try:
            pwd.getpwuid(os.getuid())
        except KeyError:
            raise RuntimeError('Current UID is absent from NSS; use the normal HPC login terminal')
        auth = subprocess.run(['munge', '-n'], capture_output=True, text=True)
        if auth.returncode:
            raise RuntimeError('Cannot access MUNGE authentication: ' + auth.stderr.strip())
        ping = subprocess.run(['scontrol', 'ping'], capture_output=True, text=True, timeout=15)
        if ping.returncode or not re.search(r'\bis UP\b', ping.stdout):
            raise RuntimeError('Controller health could not be established: ' + ping.stderr.strip())
        active = subprocess.run(['squeue', '-h', '-u', str(os.getuid()), '-o', '%i|%j|%T'],
                                capture_output=True, text=True, check=True, timeout=15)
        if any('|marker-routed-tme-v1|' in line or '|tme-route-collect-v1|' in line
               for line in active.stdout.splitlines()):
            raise RuntimeError('An annotation/collection job already exists; inspect it before submitting again')
    run.mkdir(parents=True, exist_ok=False)
    (run / 'code').mkdir(); (run / 'logs').mkdir()
    for script in task.joinpath('scripts').iterdir():
        if script.suffix in ['.py', '.R', '.sbatch']:
            shutil.copy2(script, run / 'code' / script.name)
    outroot = (args.output_root or pc / 'annotation_marker_routed_v1' / stamp).resolve()
    sources, ready, pending, reused = [], [], [], []
    inventory = task / 'input_inventory.tsv'
    expected = {}
    if inventory.is_file():
        with inventory.open() as handle:
            expected = {r['source']: int(r['cells']) for r in csv.DictReader(handle, delimiter='\t')}
    for shard in shards:
        s = shard.stem; sidecar = malig / (s + '.tsv.gz')
        old = task / (s + '_v2'); oldsum = old / 'summary.json'
        use_old = False
        if oldsum.is_file() and sidecar.is_file():
            j = json.loads(oldsum.read_text())
            use_old = (not j['pilot'] and j['source'] == s
                and j['malignancy_sha256'] == sha(sidecar)
                and j['script_sha256'] == sha(run / 'code/predict_routed_source.py')
                and j['policy_sha256'] == sha(run / 'code/marker_policy.py'))
            for name, model in j['models'].items():
                if 'sha256' in model:
                    use_old &= model['sha256'] == sha(pc / 'resources/celltypist/data/models' / (name + '.pkl'))
                if name == 'SingleR_HPCA':
                    use_old &= (model['matrix_sha256'] == sha(pc / 'resources/R_ExperimentHub/1d748702c9eed_3508')
                                and model['labels_sha256'] == sha(pc / 'resources/R_ExperimentHub/1d7483d3c0df_3509')
                                and model['worker_sha256'] == sha(run / 'code/singler_hpca.R'))
            if use_old:
                subprocess.run([str(py), '-I', str(run / 'code/check_output.py'),
                                str(old), '--malignancy', str(sidecar)], check=True)
        entry = {'source': s, 'shard': str(shard),
                 'output': str(old if use_old else outroot / s),
                 'status': 'reused_verified' if use_old else
                           ('ready' if sidecar.is_file() else 'awaiting_malignancy'),
                 'parent_job': PARENTS.get(s)}
        if s in expected:
            entry['expected_cells'] = expected[s]
        sources.append(entry)
        (reused if use_old else ready if sidecar.is_file() else pending).append(entry)
    plan = {'created_utc': stamp, 'submitted': False, 'pc_root': str(pc),
            'repo': str(args.repo), 'output_root': str(outroot),
            'malignancy_dir': str(malig), 'sources': sources,
            'official_wrapper_sha256': sha(args.repo / 'scripts/CellTypist_official_unicode_axis_v1.py'),
            'source_counts': {'ready_to_submit': len(ready), 'reused_verified': len(reused),
                              'awaiting_malignancy': len(pending)},
            'code_hashes': {s.name: sha(s) for s in (run / 'code').iterdir() if s.is_file()}}
    planpath = run / 'plan.json'; planpath.write_text(json.dumps(plan, indent=2))
    for name, entries in [('ready', ready), ('pending', pending)]:
        (run / (name + '.txt')).write_text(''.join(e['shard'] + '\n' for e in entries))
    print(json.dumps({'plan': str(planpath), **plan['source_counts'], 'submitted': False}), flush=True)
    if not args.submit:
        return 0
    env = os.environ.copy()
    env.update(PC_ROOT=str(pc), CODE=str(args.repo),
               WORKER=str(run / 'code/predict_routed_source.py'), OUT=str(outroot))
    receipt = {'jobs': [], 'blocked_sources': [], 'plan': str(planpath)}
    receiptpath = run / 'submission_receipt.json'

    def save():
        receiptpath.write_text(json.dumps(receipt, indent=2))

    def submit(options, script, source_list=None):
        if source_list:
            env['SOURCE_LIST'] = str(source_list)
        result = subprocess.run(['sbatch', '--parsable', '--chdir=' + str(run),
            '--output=' + str(run / 'logs/%x-%A_%a.out'), '--export=ALL',
            *options, str(script)], env=env, capture_output=True, text=True, check=True)
        job = result.stdout.strip().split(';')[0]
        if not job.isdigit():
            raise ValueError('Unexpected sbatch job ID')
        receipt['jobs'].append({'job_id': job, 'options': options, 'script': script.name})
        save()
        return job

    save()
    if ready:
        submit(['--array=0-' + str(len(ready)-1) + '%4'],
               run / 'code/predict_array.sbatch', run / 'ready.txt')
    for i, e in enumerate(pending):
        parent = e['parent_job']
        result = subprocess.run(['scontrol', 'show', 'job', parent, '-o'],
                                capture_output=True, text=True) if parent else None
        state = re.search(r'JobState=(\S+)', result.stdout) if result else None
        if state and state.group(1) in ['RUNNING', 'PENDING', 'CONFIGURING', 'COMPLETING']:
            submit(['--array=' + str(i), '--dependency=afterok:' + parent],
                   run / 'code/predict_array.sbatch', run / 'pending.txt')
        else:
            receipt['blocked_sources'].append({'source': e['source'], 'parent_job': parent,
                'observed_state': state.group(1) if state else 'unavailable',
                'reason': 'Final malignancy output absent and live producing job not established'})
            save()
    producer_ids = [e['job_id'] for e in receipt['jobs']]
    collect = run / 'collect.sbatch'
    import shlex
    command = [str(py), '-I', '-u', str(run / 'code/collect_sources.py'), '--plan', str(planpath)]
    collect.write_text('#!/bin/bash\nset -euo pipefail\n' + shlex.join(command) + '\n')
    options = ['--job-name=tme-route-collect-v1', '--partition=regular', '--cpus-per-task=2',
               '--mem=16G', '--time=04:00:00']
    if producer_ids:
        options.append('--dependency=afterany:' + ':'.join(producer_ids))
    submit(options, collect)
    plan['submitted'] = True; planpath.write_text(json.dumps(plan, indent=2))
    print(json.dumps(receipt), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
