"""Prospectively fixed secondary family; generic code in, source aggregates out."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import io
import json
import os
import shlex
import subprocess


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    assert not path.exists()
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    task = Path(os.environ['PANCANCER_LN_STRUCTURE_TASK'])
    code = Path(__file__).parent
    frozen = task/'private_frozen_LN_specificity_execution_v1.json'
    plan = json.loads(frozen.read_text())
    for name, sha in plan['code_SHA256'].items():
        assert digest(code/name) == sha
    assert plan['range_module'] == 'fetch_public_RDS_ranges_resume_RAM_v1.py'
    reader = (code/'read_LN_regulatory_specificity_RAM_v1.py').read_text()
    R = '\n'.join((code/name).read_text() for name in plan['R_components'])
    assert hashlib.sha256(R.encode()).hexdigest() == plan['R_code_SHA256']
    wrapper = 'import sys,json,io,hashlib\nplan=json.loads(sys.stdin.readline())\n'
    wrapper += 'source='+repr(reader)+'\nR_SOURCE='+repr(R)+'\n'
    for variable, name in [('PIPE_SOURCE', 'pipe_native_RDS_RAM_v1.py'),
                           ('RANGE_SOURCE', plan['range_module']),
                           ('RANGE_CORE_SOURCE', 'fetch_public_RDS_ranges_RAM_v1.py')]:
        wrapper += variable+'='+repr((code/name).read_text())+'\n'
    wrapper += 'assert hashlib.sha256(source.encode()).hexdigest()==plan["code_SHA256"]["read_LN_regulatory_specificity_RAM_v1.py"]\n'
    wrapper += 'sys.stdin=io.StringIO(json.dumps(plan))\nns={"__name__":"frozen_secondary_source","R_SOURCE":R_SOURCE,"PIPE_SOURCE":PIPE_SOURCE,"RANGE_SOURCE":RANGE_SOURCE,"RANGE_CORE_SOURCE":RANGE_CORE_SOURCE}\n'
    wrapper += 'exec(compile(source,"frozen_secondary_source","exec"),ns)\nns["main"]()\n'
    args = ['ssh', '-T', '-o', 'ClearAllForwardings=yes', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15',
            '-o', 'ServerAliveInterval=30', '-o', 'ServerAliveCountMax=10',
            '-o', 'ProxyCommand='+os.environ['PANCANCER_READONLY_PROXY_COMMAND'], os.environ['PANCANCER_READONLY_SSH_TARGET'],
            shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -c '+shlex.quote(wrapper)]
    child = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    child.stdin.write(json.dumps(plan)+'\n');child.stdin.close()
    terminal = None
    for line in child.stdout:
        assert len(line) <= 1024**2
        frame = json.loads(line);kind = frame['kind']
        assert kind in ['SOURCE_RAM_CONTROL_ROOT_CREATED', 'SOURCE_RAM_WORKER_STARTED', 'PUBLIC_INPUT_RAM_BYTES_PROGRESS',
                        'COMPLETE_PUBLIC_PROTEIN_SECONDARY_SPECIFICITY_RAM_RELEASED', 'PUBLIC_PROTEIN_RDS_SCOPE_NOT_ACCEPTED', 'SOURCE_SPECIFICITY_FAILED']
        if kind == 'SOURCE_RAM_CONTROL_ROOT_CREATED':
            save(task/'private_LN_specificity_control_root_v1.json', frame)
        elif kind == 'SOURCE_RAM_WORKER_STARTED':
            save(task/'private_LN_specificity_actual_worker_v1.json', frame)
            print(json.dumps(frame), flush=True)
        elif kind == 'PUBLIC_INPUT_RAM_BYTES_PROGRESS':
            (task/'private_LN_specificity_live_checkpoint_v1.json').write_text(json.dumps(frame)+'\n')
        elif kind == 'COMPLETE_PUBLIC_PROTEIN_SECONDARY_SPECIFICITY_RAM_RELEASED':
            terminal = frame
        elif kind == 'SOURCE_SPECIFICITY_FAILED':
            save(task/'private_LN_specificity_R_failed_v1.json', frame)
        else:
            save(task/'private_LN_specificity_input_failed_v1.json', frame)
    exitcode = child.wait(timeout=30)
    if exitcode != 0 or terminal is None:
        save(task/'private_LN_specificity_transport_exit_v1.json', dict(status='FROZEN_SECONDARY_WORKFLOW_NOT_ACCEPTED',
            observed_UTC=datetime.now(timezone.utc).isoformat(), actual_exit_code=exitcode, no_individual_payload_received=True))
        raise RuntimeError('FROZEN_SECONDARY_WORKFLOW_NOT_ACCEPTED')
    assert terminal['source_R_exit_code'] == 0 and terminal['source_payload_RAM_released']
    assert terminal['summary']['audit']['status'] == 'PASS_FULL_SECONDARY_SPECIFICITY_NUMERIC_AUDIT'
    save(task/'private_LN_specificity_complete_source_v1.json', terminal)
    save(task/'private_LN_specificity_actual_exit_v1.json', dict(status='SOURCE_EXIT0_SECONDARY_FAMILY_AND_RAM_RELEASE_COMPLETE',
        observed_UTC=datetime.now(timezone.utc).isoformat(), actual_exit_code=exitcode,
        frozen_contract_SHA256=digest(frozen), source_result_SHA256=digest(task/'private_LN_specificity_complete_source_v1.json'),
        no_individual_arrays_or_identifiers_transferred=True))
    print(json.dumps(dict(status='COMPLETE_SOURCE_SECONDARY_FAMILY_AUDIT_EXIT0')), flush=True)


if __name__ == '__main__':
    main()
