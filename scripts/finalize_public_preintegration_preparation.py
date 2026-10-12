#!/usr/bin/env python3
"""Finish the current detached run's readable log and preparation index."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

from build_public_preintegration_manifest import build
from collect_ambient_preparation_status import collect, render


def save(path, value):
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False)


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-root',required=True)
    p.add_argument('--author-history',required=True);p.add_argument('--code-git-sha',required=True)
    args=p.parse_args();root=Path(args.run_root)
    save(root/'manifests/preintegration_finalizer_start_v1.json',dict(
        started_utc=datetime.now(timezone.utc).isoformat(),code_git_sha=args.code_git_sha,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        component_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ['build_public_preintegration_manifest.py','collect_ambient_preparation_status.py']},
        author_history_sha256=hashlib.sha256(Path(args.author_history).read_bytes()).hexdigest(),
        source_final_AC=False,cross_study_integrated=False))
    while not (root/'manifests/execution_terminal.json').exists():
        status=collect(root)
        if not status['orchestrator_alive']:
            save(root/'manifests/preintegration_finalizer_interrupted_v1.json',dict(
                observed_utc=status['observed_utc'],status='orchestrator_stopped_without_terminal',
                source_final_AC=False,cross_study_integrated=False))
            return 1
        time.sleep(20)
    status=collect(root)
    save(root/'manifests/public_status_terminal_v1.json',status)
    with (root/'manifests/PREPROCESS_STATUS_terminal_v1.md').open('x') as f:f.write(render(status))
    if status['status']!='completed' or status['independent_audit_status']!='passed':return 1
    manifest=build(root/'manifests/independent_preparation_audit.json',Path(args.author_history))
    save(root/'manifests/preintegration_preparation_terminal_v1.json',manifest)
    save(root/'manifests/preintegration_finalizer_terminal_v1.json',dict(
        observed_utc=datetime.now(timezone.utc).isoformat(),status='preparation_index_written_pending_source_and_scope_gates',
        n_audited_prepared_cells=manifest['n_audited_prepared_cells'],
        source_final_AC=False,cross_study_integrated=False,
        git_publication='not_automated_requires_reviewed_ordinary_push'))
    return 0


if __name__=='__main__':raise SystemExit(main())
