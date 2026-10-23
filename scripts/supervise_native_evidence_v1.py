"""Capture the actual exit status of this task's bounded native-RNA child."""
import argparse
import datetime
import json
import os
from pathlib import Path
import signal
import subprocess
import sys


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def startticks(pid):
    path = Path('/proc') / str(pid) / 'stat'
    return path.read_text().rsplit(')', 1)[1].split()[19] if path.exists() else None


def main(args):
    root = Path(args.output)
    command = [sys.executable, '-B', str(root / 'scripts/run_native_coarse_evidence_v1.py'),
               '--source-run', args.source_run, '--output', str(root),
               '--block-size', str(args.block_size), '--wall-seconds', str(args.wall_seconds),
               '--output-cap-bytes', str(args.output_cap_bytes)]
    with (root / 'provenance/science_stdout.txt').open('x') as stdout, (root / 'provenance/science_stderr.txt').open('x') as stderr:
        child = subprocess.Popen(command, stdout=stdout, stderr=stderr, start_new_session=True)
        ticks = startticks(child.pid)
        save(root / 'provenance/child_launch_v1.json',
             {'started_UTC': now(), 'PID': child.pid, 'startticks': ticks})
        timed_out = False
        try:
            code = child.wait(timeout=args.wall_seconds + 180)
        except subprocess.TimeoutExpired:
            timed_out = True
            if child.poll() is None and ticks is not None and startticks(child.pid) == ticks:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    if child.poll() is None and startticks(child.pid) == ticks:
                        os.killpg(child.pid, signal.SIGKILL)
                code = child.wait(timeout=30)
            else:
                raise RuntimeError('Cannot verify bounded child identity for timeout cleanup')
    save(root / 'provenance/supervisor_terminal_v1.json',
         {'completed_UTC': now(), 'actual_child_exit_code': code,
          'wall_timeout': timed_out, 'PID': child.pid, 'startticks': ticks,
          'science_terminal_exists': (root / 'provenance/terminal_v1.json').exists()})
    return 0 if code == 0 and not timed_out else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    parser.add_argument('--source-run', required=True)
    parser.add_argument('--block-size', type=int, default=2048)
    parser.add_argument('--wall-seconds', type=int, default=10800)
    parser.add_argument('--output-cap-bytes', type=int, default=16 * 1024 ** 2)
    sys.exit(main(parser.parse_args()))
