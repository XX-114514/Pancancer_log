"""Supervise a frozen task-owned Python stage and retain its actual exit code."""
import argparse
import datetime
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys


def ticks(pid):
    path = Path('/proc') / str(pid) / 'stat'
    return path.read_text().rsplit(')', 1)[1].split()[19] if path.exists() else None


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def main(root):
    config = json.loads((root / 'config/supervisor_request_v1.json').read_text())
    name = config['entrypoint']
    if not re.fullmatch(r'[A-Za-z0-9_]+\.py', name):
        raise ValueError('Entrypoint must be a task-owned Python filename')
    arguments = config['arguments']
    if not isinstance(arguments, list) or not all(isinstance(x, str) for x in arguments):
        raise ValueError('Frozen argument vector required')
    wall = int(config['wall_seconds'])
    if wall < 1 or wall > 86400:
        raise ValueError('Invalid bounded stage wall budget')
    with (root / 'provenance/science_stdout.txt').open('x') as stdout, (root / 'provenance/science_stderr.txt').open('x') as stderr:
        child = subprocess.Popen([sys.executable, '-B', str(root / 'scripts' / name)] + arguments,
                                 stdout=stdout, stderr=stderr, start_new_session=True)
        start = ticks(child.pid)
        save(root / 'provenance/child_launch_v1.json',
             {'started_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'PID': child.pid, 'startticks': start})
        timed_out = False
        try:
            code = child.wait(timeout=wall)
        except subprocess.TimeoutExpired:
            timed_out = True
            if child.poll() is None and start is not None and ticks(child.pid) == start:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    if child.poll() is None and ticks(child.pid) == start:
                        os.killpg(child.pid, signal.SIGKILL)
                code = child.wait(timeout=30)
            else:
                raise RuntimeError('Cannot verify task child identity at timeout')
    save(root / 'provenance/supervisor_terminal_v1.json',
         {'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'actual_child_exit_code': code, 'wall_timeout': timed_out,
          'science_terminal_exists': (root / 'provenance/terminal_v1.json').exists(),
          'PID': child.pid, 'startticks': start})
    return 0 if code == 0 and not timed_out else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    sys.exit(main(Path(parser.parse_args().output)))
