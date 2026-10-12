#!/usr/bin/env python3
"""Rerun all audit checks with textual identifiers; preserve the failed proof."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def save(path,data):
    with path.open('x') as out:json.dump(data,out,indent=2)


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--preprocess-root',required=True)
    p.add_argument('--collection',required=True);p.add_argument('--code-git-sha',required=True)
    args=p.parse_args();root=Path(args.preprocess_root)
    original=json.loads((root/'manifests/independent_preparation_audit.json').read_text())
    if original['all_checks_pass']:raise ValueError('Do not replace an already passing proof')
    if original['n_datasets']!=9 or original['n_libraries']!=792:raise ValueError('Unexpected audit scope')
    if not all(l['all_checks_pass'] for d in original['datasets'] for l in d['libraries']):
        raise ValueError('This repair does not bypass cell-level failures')
    backup=root/'config/identifier_audit_repair_v4';backup.mkdir(mode=0o700)
    for rel in ['manifests/independent_preparation_audit.json','manifests/execution_terminal.json',
                'scripts/audit_ambient_preparation.py']:
        path=root/rel;shutil.copyfile(path,backup/path.name)
    script=root/'scripts/audit_ambient_preparation_identifier_v4.py'
    save(root/'manifests/identifier_audit_repair_start_v4.json',dict(
        started_utc=datetime.now(timezone.utc).isoformat(),code_git_sha=args.code_git_sha,
        repair_script_sha256=digest(Path(__file__)),new_audit_sha256=digest(script),
        old_audit_report_sha256=digest(backup/'independent_preparation_audit.json'),
        reason='textual_IDs_preserve_leading_zeros',scientific_worker_parameters_changed=False))
    report=root/'manifests/independent_preparation_audit_identifier_v4.json'
    with (root/'logs/independent_audit_identifier_v4.log').open('x') as log:
        result=subprocess.run([sys.executable,'-u',str(script),'--output-root',str(root/'objects'),
            '--collection',args.collection,'--report',str(report)],stdout=log,stderr=subprocess.STDOUT)
    proof=json.loads(report.read_text()) if report.exists() else {}
    if result.returncode or not proof.get('all_checks_pass'):
        save(root/'manifests/identifier_audit_repair_terminal_v4.json',dict(status='needs_repair',exit_code=result.returncode))
        return 1
    for key in ['n_datasets','n_libraries','n_output','n_ready_technical_cells']:
        if proof[key]!=original[key]:raise ValueError('Scientific counts/gates changed during ID-only audit repair')
    shutil.copyfile(report,root/'manifests/independent_preparation_audit.json')
    shutil.copyfile(script,root/'scripts/audit_ambient_preparation.py')
    terminal=json.loads((backup/'execution_terminal.json').read_text())
    terminal['original_failed_audit_end_time']=terminal['end_time']
    terminal.update(end_time=datetime.now(timezone.utc).isoformat(),status='completed',audit_exit_code=0,
        audit_repair_version='textual_identifiers_v4',audit_code_git_sha=args.code_git_sha,
        original_failed_audit_retained=True,new_audit_sha256=digest(script))
    # The original terminal and report remain sealed in the private backup.
    (root/'manifests/execution_terminal.json').write_text(json.dumps(terminal,indent=2))
    summary=dict(schema='identifier_audit_repair_terminal_v4',completed_utc=terminal['end_time'],
        status='completed',all_checks_pass=True,n_datasets=proof['n_datasets'],n_libraries=proof['n_libraries'],
        n_cells=proof['n_output'],n_ready_technical_cells=proof['n_ready_technical_cells'],
        audit_script_sha256=digest(script),code_git_sha=args.code_git_sha,
        old_failed_report_retained=True,scientific_worker_parameters_changed=False,source_final_AC=False)
    save(root/'manifests/identifier_audit_repair_terminal_v4.json',summary);print(json.dumps(summary));return 0


if __name__=='__main__':raise SystemExit(main())
