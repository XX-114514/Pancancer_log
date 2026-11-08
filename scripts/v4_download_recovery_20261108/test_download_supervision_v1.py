"""Finite-controller tests with an injected clock and synthetic worker exits."""
from pathlib import Path
from types import SimpleNamespace
import fcntl
import json
import sys
from supervise_download_recovery_v1 import supervise

base = Path(sys.argv[1])
base.mkdir(parents=True, exist_ok=False)


def fixture(name):
    root = base / name
    out = root / 'download_recovery_v1'
    out.mkdir(parents=True)
    return root, out


root, out = fixture('retry_then_success')
clock = [0.]
calls = []


def sleep(seconds):
    clock[0] += seconds


def worker(command, **kwargs):
    calls.append(command)
    success = len(calls) == 2
    (out / 'latest_terminal_PRIVATE.json').write_text(json.dumps({
        'status': 'ALL_ACCEPTED' if success else 'INCOMPLETE_RESUMABLE_REPEAT_SAME_COMMAND',
        'accepted_native_sources': 6 if success else 1}))
    return SimpleNamespace(returncode=0 if success else 2)


assert supervise(root, hours=1, launcher=worker, sleeper=sleep, monotonic=lambda: clock[0]) == 0
assert len(calls) == 2 and calls[0][-1] == '3'
assert supervise(root, hours=1, launcher=worker, sleeper=sleep, monotonic=lambda: clock[0]) == 0
assert len(calls) == 2  # Already complete: no extra worker.

root, out = fixture('active_worker_preserved')
clock = [0.]
handle = (out / 'recovery_writer.lock').open('a')
fcntl.flock(handle, fcntl.LOCK_EX)
waited = []


def release_after_wait(seconds):
    waited.append(seconds)
    clock[0] += seconds
    fcntl.flock(handle, fcntl.LOCK_UN)


def done_worker(command, **kwargs):
    assert waited  # Never launches while the prior writer owns its lock.
    (out / 'latest_terminal_PRIVATE.json').write_text(json.dumps({'status': 'ALL_ACCEPTED', 'accepted_native_sources': 6}))
    return SimpleNamespace(returncode=0)


assert supervise(root, hours=1, launcher=done_worker, sleeper=release_after_wait,
                 monotonic=lambda: clock[0]) == 0
handle.close()

root, out = fixture('access_review')
(out / 'events_PRIVATE.jsonl').write_text(json.dumps({'needs_source_or_access_review': True}) + '\n')
assert supervise(root, hours=1, launcher=lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError('must not launch'))) == 3

root, out = fixture('finite_attempts')
clock = [0.]
calls = []


def fail(command, **kwargs):
    calls.append(command)
    return SimpleNamespace(returncode=1)


assert supervise(root, hours=1, launches=1, launcher=fail, sleeper=sleep,
                 monotonic=lambda: clock[0]) == 2
assert len(calls) == 1

root, out = fixture('unfinished_event')
(out / 'events_PRIVATE.jsonl').write_text('{"status":')


def finish(command, **kwargs):
    (out / 'latest_terminal_PRIVATE.json').write_text(json.dumps({'status': 'ALL_ACCEPTED', 'accepted_native_sources': 6}))
    return SimpleNamespace(returncode=0)


assert supervise(root, hours=1, launcher=finish) == 0
print('PASS_SYNTHETIC_ONLY: incomplete worker resumed, completed reused, active writer waited without termination, access failure stops, finite launches enforced')
