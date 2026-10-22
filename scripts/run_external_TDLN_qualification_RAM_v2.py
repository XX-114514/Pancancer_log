"""Dispatch qualification with the separately versioned 2 GiB input guard."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import shlex
import subprocess


def save(path, value):
    assert not path.exists()
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    task = Path(os.environ['PANCANCER_EXTERNAL_TDLN_TASK'])
    code = Path(__file__).parent
    frozen = task/'private_frozen_native_qualification_v1.json'
    plan = json.loads(frozen.read_text())
    for name, sha in plan['code_SHA256'].items():
        assert hashlib.sha256((code/name).read_bytes()).hexdigest() == sha
    reader = (code/'read_external_TDLN_qualification_RAM_v1.py').read_text()
    R = '\n'.join((code/name).read_text() for name in plan['R_components'])
    assert hashlib.sha256(R.encode()).hexdigest() == plan['R_code_SHA256']
    wrapper = 'import sys,json,hashlib\nplan=json.loads(sys.stdin.readline())\n'
    wrapper += 'source='+repr(reader)+'\nR_SOURCE='+repr(R)+'\n'
    for variable, name in [('PIPE_SOURCE', 'pipe_native_RDS_RAM_v1.py'),
                           ('RANGE_SOURCE', 'fetch_public_RDS_ranges_resume_RAM_v1.py'),
                           ('RANGE_CORE_SOURCE', 'fetch_public_RDS_ranges_RAM_v2.py')]:
        wrapper += variable+'='+repr((code/name).read_text())+'\n'
    wrapper += 'assert hashlib.sha256(source.encode()).hexdigest()==plan["code_SHA256"]["read_external_TDLN_qualification_RAM_v1.py"]\n'
    wrapper += 'exec(compile(source,"frozen_external_TDLN_reader","exec"))\nmain(plan)\n'
    args = ['ssh', '-T', '-o', 'ClearAllForwardings=yes', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15',
            '-o', 'ServerAliveInterval=30', '-o', 'ServerAliveCountMax=10',
            '-o', 'ProxyCommand='+os.environ['PANCANCER_READONLY_PROXY_COMMAND'], os.environ['PANCANCER_READONLY_SSH_TARGET'],
            shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -c '+shlex.quote(wrapper)]
    child = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    child.stdin.write(json.dumps(plan)+'\n');child.stdin.close()
    terminal = None
    for line in child.stdout:
        assert len(line) <= 1024**2
        frame = json.loads(line); kind = frame['kind']
        if kind == 'SOURCE_QUALIFICATION_STARTED':
            save(task/'private_actual_source_start_v1.json', frame)
            print(json.dumps(frame), flush=True)
        elif kind == 'PUBLIC_INPUT_RAM_BYTES_PROGRESS':
            (task/'private_live_input_progress_v1.json').write_text(json.dumps(frame)+'\n')
        elif kind == 'COMPLETE_INPUT_IDENTITY_ACCEPTED':
            save(task/'private_complete_input_identity_v1.json', frame)
        elif kind == 'COMPLETE_EXTERNAL_TDLN_QUALIFICATION_RAM_RELEASED':
            terminal = frame
        elif kind == 'EXTERNAL_TDLN_QUALIFICATION_NOT_ACCEPTED':
            save(task/'private_source_not_accepted_v1.json', frame)
        else:
            raise ValueError('UNRECOGNIZED_SOURCE_FRAME')
    exitcode = child.wait(timeout=30)
    if exitcode != 0 or terminal is None:
        save(task/'private_actual_failure_exit_v1.json',dict(status='NATIVE_QUALIFICATION_NOT_ACCEPTED',
            observed_UTC=datetime.now(timezone.utc).isoformat(), actual_exit_code=exitcode))
        raise SystemExit(2)
    assert terminal['source_R_exit_code'] == 0 and terminal['source_payload_RAM_released']
    save(task/'private_complete_source_qualification_v1.json', terminal)
    save(task/'private_actual_exit_v1.json', dict(status='SOURCE_EXIT0_QUALIFICATION_AND_RAM_RELEASE_COMPLETE',
        observed_UTC=datetime.now(timezone.utc).isoformat(),actual_exit_code=exitcode,
        frozen_contract_SHA256=hashlib.sha256(frozen.read_bytes()).hexdigest(),
        complete_source_SHA256=hashlib.sha256((task/'private_complete_source_qualification_v1.json').read_bytes()).hexdigest()))
    print(json.dumps(dict(status='COMPLETE_EXTERNAL_TDLN_QUALIFICATION_EXIT0')), flush=True)


if __name__ == '__main__':
    main()
