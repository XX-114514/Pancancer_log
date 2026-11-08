"""Bounded continuation after worker exit; never stops a running worker."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json
import os
import subprocess
import sys
import time
from resumable_public_zip_v1 import exclusive_lock, BusyError, atomic_json


def supervise(root, hours=48, launches=3, launcher=subprocess.run,
              sleeper=time.sleep, monotonic=time.monotonic):
    root = Path(root).resolve()
    out = root / 'download_recovery_v1'
    out.mkdir(exist_ok=True)
    deadline = monotonic() + hours * 3600
    count = 0

    def completed():
        path = out / 'latest_terminal_PRIVATE.json'
        if not path.exists():
            return False
        item = json.loads(path.read_text())
        return item.get('status') == 'ALL_ACCEPTED' and item.get('accepted_native_sources') == 6

    def access_review():
        path = out / 'events_PRIVATE.jsonl'
        if not path.exists():
            return False
        lines = path.read_text().splitlines(keepends=True)
        for index, line in enumerate(lines):
            try:
                item = json.loads(line)
            except json.JSONDecodeError:
                if index == len(lines) - 1 and not line.endswith('\n'):
                    continue  # A crash or current writer may leave an unfinished event.
                return True  # A malformed committed event requires review.
            if item.get('needs_source_or_access_review'):
                return True
        return False

    def event(**data):
        with (out / 'supervisor_events_PRIVATE.jsonl').open('a') as handle:
            handle.write(json.dumps({'UTC': datetime.now(timezone.utc).isoformat(), **data}) + '\n')

    with exclusive_lock(out / 'supervisor.lock'):
        event(status='FINITE_CONTINUATION_STARTED_OLD_WORKER_PRESERVED', hours=hours, max_launches=launches)
        while monotonic() < deadline:
            if completed():
                event(status='ALL6_ACCEPTED_SUPERVISION_COMPLETE', additional_launches=count)
                return 0
            if access_review():
                event(status='SOURCE_OR_ACCESS_REVIEW_REQUIRED_NO_AUTOMATIC_BYPASS')
                return 3
            try:
                with exclusive_lock(out / 'recovery_writer.lock'):
                    pass
            except BusyError:
                sleeper(min(30, max(0, deadline - monotonic())))
                continue
            if count >= launches:
                break
            remaining = deadline - monotonic()
            if remaining <= 0:
                break
            count += 1
            event(status='START_NEXT_WORKER_WITH_SAME_PERSISTENT_CACHE', launch=count)
            with (out / f'supervised_worker_{count}_PRIVATE.log').open('ab') as log:
                result = launcher([sys.executable, '-u', str(Path(__file__).with_name('recover_xenium_downloads_v1.py')),
                                   str(root), '--hours', str(min(18, remaining / 3600)), '--rounds', '3'],
                                  stdout=log, stderr=subprocess.STDOUT,
                                  env={**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'})
            event(status='WORKER_EXITED_KEEP_ALL_PARTIAL_AND_CACHE', exit_code=result.returncode, launch=count)
            if completed():
                event(status='ALL6_ACCEPTED_SUPERVISION_COMPLETE', additional_launches=count)
                return 0
            sleeper(min(300, max(0, deadline - monotonic())))
        event(status='FINITE_SUPERVISION_ENDED_INCOMPLETE_CACHE_RETAINED', additional_launches=count)
        return 2


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_root', type=Path)
    parser.add_argument('--hours', type=float, default=48)
    parser.add_argument('--launches', type=int, default=3)
    args = parser.parse_args()
    if not 0 < args.hours <= 72 or not 1 <= args.launches <= 4:
        parser.error('Finite bounds: hours in(0,72], launches in[1,4]')
    raise SystemExit(supervise(args.run_root, args.hours, args.launches))
