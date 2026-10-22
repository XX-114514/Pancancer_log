"""Dispatch generic code only; receive source aggregates, never cell arrays."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import io
import json
import os
import shlex
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    assert not path.exists()
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    task = Path(os.environ['PANCANCER_LN_STRUCTURE_TASK'])
    code = Path(__file__).parent
    frozen = task/'private_frozen_LN_object_structure_v6.json'
    plan = json.loads(frozen.read_text())
    for name, digest in plan['code_SHA256'].items():
        assert sha(code/name) == digest
    python = (code/'read_LN_spatial_protein_RDS_RAM_v5.py').read_text()
    R = (code/'qualify_LN_spatial_RDS_structure_RAM_v1.R').read_text()
    pipe = (code/'pipe_native_RDS_RAM_v1.py').read_text()
    ranges = (code/'fetch_public_RDS_ranges_retry_RAM_v3.py').read_text()
    range_core = (code/'fetch_public_RDS_ranges_RAM_v1.py').read_text()
    wrapper = 'import sys,json,io,hashlib\nplan=json.loads(sys.stdin.readline())\n'
    wrapper += 'source='+repr(python)+'\nR_SOURCE='+repr(R)+'\nPIPE_SOURCE='+repr(pipe)+'\nRANGE_SOURCE='+repr(ranges)+'\nRANGE_CORE_SOURCE='+repr(range_core)+'\n'
    wrapper += 'assert hashlib.sha256(source.encode()).hexdigest()==plan["code_SHA256"]["read_LN_spatial_protein_RDS_RAM_v5.py"]\n'
    wrapper += 'sys.stdin=io.StringIO(json.dumps(plan))\nns={"__name__":"qualified_public_input","R_SOURCE":R_SOURCE,"PIPE_SOURCE":PIPE_SOURCE,"RANGE_SOURCE":RANGE_SOURCE,"RANGE_CORE_SOURCE":RANGE_CORE_SOURCE}\nexec(compile(source,"qualified_public_input","exec"),ns)\nns["main"]()\n'
    args = ['ssh', '-T', '-o', 'ClearAllForwardings=yes', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15', '-o', 'ServerAliveInterval=30', '-o', 'ServerAliveCountMax=10',
        '-o', 'ProxyCommand='+os.environ['PANCANCER_READONLY_PROXY_COMMAND'], os.environ['PANCANCER_READONLY_SSH_TARGET'],
        shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -c '+shlex.quote(wrapper)]
    child = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    child.stdin.write(json.dumps(plan)+'\n')
    child.stdin.close()
    terminal = None
    progress = task/'private_LN_structure_live_checkpoint_v6.json'
    for line in child.stdout:
        assert len(line) <= 1024**2
        frame = json.loads(line)
        kind = frame['kind']
        assert kind in ['SOURCE_RAM_CONTROL_ROOT_CREATED', 'SOURCE_RAM_WORKER_STARTED', 'PUBLIC_INPUT_RAM_BYTES_PROGRESS',
            'COMPLETE_PUBLIC_PROTEIN_RDS_RAM_STRUCTURE_NO_EFFECTS', 'PUBLIC_PROTEIN_RDS_SCOPE_NOT_ACCEPTED', 'SOURCE_RDS_STRUCTURE_FAILED']
        if kind == 'SOURCE_RAM_CONTROL_ROOT_CREATED':
            save(task/'private_LN_structure_control_root_v6.json', frame)
        elif kind == 'SOURCE_RAM_WORKER_STARTED':
            save(task/'private_LN_structure_actual_worker_v6.json', dict(frame, frozen_contract_SHA256=sha(frozen)))
            print(json.dumps(dict(status='ACTUAL_SOURCE_RAM_STRUCTURE_WORKER_STARTED')), flush=True)
        elif kind == 'PUBLIC_INPUT_RAM_BYTES_PROGRESS':
            progress.write_text(json.dumps(frame, indent=2)+'\n')
            print(json.dumps(frame), flush=True)
        elif kind == 'COMPLETE_PUBLIC_PROTEIN_RDS_RAM_STRUCTURE_NO_EFFECTS':
            assert terminal is None
            terminal = frame
        elif kind == 'SOURCE_RDS_STRUCTURE_FAILED':
            save(task/'private_LN_R_structure_failed_v6.json', frame)
        else:
            save(task/'private_LN_structure_failed_v6.json', frame)
    exitcode = child.wait(timeout=30)
    if exitcode != 0 or terminal is None:
        save(task/'private_LN_structure_transport_failed_v6.json', dict(status='SOURCE_STRUCTURE_INPUT_NOT_ACCEPTED',
            observed_UTC=datetime.now(timezone.utc).isoformat(), source_exit_code=exitcode, no_individual_payload_received=True))
        raise RuntimeError('SOURCE_STRUCTURE_INPUT_NOT_ACCEPTED')
    assert terminal['source_R_exit_code'] == 0 and terminal['source_payload_RAM_released']
    save(task/'private_LN_complete_native_structure_v6.json', terminal)
    save(task/'private_LN_structure_actual_exit_v6.json', dict(status='COMPLETE_SOURCE_EXIT0_RAM_RELEASED_STRUCTURE_ONLY',
        observed_UTC=datetime.now(timezone.utc).isoformat(), source_exit_code=exitcode,
        complete_source_result_SHA256=sha(task/'private_LN_complete_native_structure_v6.json'),
        frozen_contract_SHA256=sha(frozen), no_individual_omics_identifiers_or_clinical_values_transferred=True))
    print(json.dumps(dict(status='COMPLETE_PUBLIC_SOURCE_STRUCTURE_AND_RAM_RELEASE_NO_EFFECTS')), flush=True)


if __name__ == '__main__':
    main()
