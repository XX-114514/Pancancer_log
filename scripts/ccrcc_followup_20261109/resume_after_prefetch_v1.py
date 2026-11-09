"""Recover acceptance after safe cache-lock contention, without killing writers.

Run beside the initial acquisition and prefetch workers. It waits for prefetch
completion and the original acquisition terminal before a bounded same-cache
resume. A cached block is never considered a complete input on its own.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import subprocess
import sys
import time
from resumable_public_zip_v1 import atomic_json, exclusive_lock

r = Path(sys.argv[1]).resolve()
deadline = time.monotonic() + 6 * 3600


def ready():
    return (r / 'inputs/complete_PRIVATE.json').exists() and (r / 'inputs/boundary_complete_PRIVATE.json').exists()


with exclusive_lock(r / 'resume_after_prefetch.lock'):
    attempts = 0
    status = 'PREFETCH_WAIT_EXPIRED'
    while time.monotonic() < deadline:
        if ready():
            status = 'COMPLETE_INPUTS_READY'
            break
        done = all((r / p).exists() for p in ['core_prefetch_v1/complete.json',
                                              'inputs/boundary_complete_PRIVATE.json',
                                              'acquisition_v1/terminal_PRIVATE.json'])
        if done:
            if attempts >= 3:
                status = 'RESUME_BUDGET_EXHAUSTED_REQUIRES_REVIEW'
                break
            attempts += 1
            with (r / 'logs/resume_after_prefetch_PRIVATE.log').open('ab', buffering=0) as f:
                result = subprocess.run([sys.executable, '-B', str(r / 'code/acquire_ccrcc_inputs_v1.py'), str(r)],
                                        stdout=f, stderr=f, timeout=max(1, deadline - time.monotonic()))
            # A nonzero return (including an existing writer lock) is retained and
            # retried only after a bounded wait, never by unlocking another writer.
        time.sleep(30)
    atomic_json(r / 'resume_after_prefetch_terminal.json', {'UTC': datetime.now(timezone.utc).isoformat(),
                                                          'status': status, 'resume_attempts': attempts})
