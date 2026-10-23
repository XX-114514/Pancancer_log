"""Bounded opaque local HTTPS -> existing SSH -> source RAM -> fixed R work."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import re
import resource
import shlex
import signal
import subprocess
import sys
import stream_public_RDS_ordered_RAM_v1 as streamer


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path,value):
    with path.open('x') as stream:
        json.dump(value,stream,indent=2)
        stream.write('\n')


def available_memory_estimate():
    if sys.platform=='darwin':
        raw = subprocess.check_output(['vm_stat'],text=True)
        page = int(re.search(r'page size of (\d+) bytes',raw).group(1))
        return sum(int(re.search(r'^'+re.escape(k)+r':\s+(\d+)',raw,re.M).group(1))
                   for k in ('Pages free','Pages inactive','Pages speculative'))*page
    return int(next(v.split()[1] for v in Path('/proc/meminfo').read_text().splitlines()
                    if v.startswith('MemAvailable:')))*1024


def main():
    task = Path(os.environ['PANCANCER_EXTERNAL_FIXED_TASK']);code = Path(__file__).parent
    frozen = task/'private_frozen_fixed_native_workflow_v1.json';plan = json.loads(frozen.read_text())
    for name,sha in plan['code_SHA256'].items():
        assert hashlib.sha256((code/name).read_bytes()).hexdigest()==sha
    assert available_memory_estimate()>=plan['minimum_local_available_memory_estimate_bytes']
    reader = (code/'receive_external_TDLN_opaque_relay_RAM_v1.py').read_text()
    R = '\n'.join((code/n).read_text() for n in plan['R_components'])
    pipe = (code/'pipe_native_RDS_RAM_v1.py').read_text()
    assert hashlib.sha256(R.encode()).hexdigest()==plan['R_code_SHA256']
    wrapper = 'import sys,json,hashlib\nline=sys.stdin.buffer.readline(65537)\nassert len(line)<=65536 and line.endswith(b"\\n")\nplan=json.loads(line)\n'
    wrapper += 'source='+repr(reader)+'\nR_SOURCE='+repr(R)+'\nPIPE_SOURCE='+repr(pipe)+'\n'
    wrapper += 'assert hashlib.sha256(source.encode()).hexdigest()==plan["code_SHA256"]["receive_external_TDLN_opaque_relay_RAM_v1.py"]\n'
    wrapper += 'exec(compile(source,"frozen_opaque_source_receiver","exec"))\nmain(plan)\n'
    args = ['ssh','-T','-o','ClearAllForwardings=yes','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes',
            '-o','ConnectTimeout=15','-o','ProxyCommand=none','-o','ProxyJump=none','-o','IPQoS=none',
            '-o','ServerAliveInterval=30','-o','ServerAliveCountMax=10',os.environ['PANCANCER_READONLY_SSH_TARGET'],
            shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -c '+shlex.quote(wrapper)]
    child = subprocess.Popen(args,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
    terminal = None;local_error = None;phase = 'SOURCE_CONTROL_AND_READY_GATES'
    def expire(_sig,_frame):
        raise TimeoutError('FINITE_LOCAL_TRANSFER_OR_SOURCE_RESULT_WAIT')
    for sig in (signal.SIGALRM,signal.SIGTERM,signal.SIGHUP):
        signal.signal(sig,expire)
    signal.alarm(plan['maximum_download_seconds'])
    def frame():
        line = child.stdout.readline(1024**2+1)
        if not line:
            return None
        assert len(line)<=1024**2 and line.endswith(b'\n')
        return json.loads(line)
    def collect(item):
        nonlocal terminal
        kind = item['kind']
        if kind=='COMPLETE_INPUT_IDENTITY_ACCEPTED':
            save(task/'private_complete_input_identity_v1.json',item)
        elif kind=='COMPLETE_EXTERNAL_FIXED_WORKFLOW_RAM_RELEASED':
            terminal = item
        elif kind=='EXTERNAL_FIXED_WORKFLOW_NOT_ACCEPTED':
            save(task/'private_source_not_accepted_v1.json',item)
        elif kind=='EXTERNAL_FIXED_WORKFLOW_NATIVE_GATE_FAILURE':
            save(task/'private_native_gate_failure_v1.json',item)
        else:
            raise ValueError('UNRECOGNIZED_SOURCE_FRAME')
    try:
        child.stdin.write((json.dumps(plan)+'\n').encode());child.stdin.flush()
        ready = frame()
        if ready is None or ready['kind']!='SOURCE_OPAQUE_RELAY_READY':
            if ready is not None:
                collect(ready)
            raise ValueError('SOURCE_NOT_READY_NO_PUBLIC_PAYLOAD_SENT')
        assert ready['R_code_SHA256']==plan['R_code_SHA256'] and ready['source_available_memory_gate_passed'] is True
        save(task/'private_actual_source_start_v1.json',ready)
        print(json.dumps(dict(status='SOURCE_OPAQUE_RELAY_READY',started_UTC=ready['started_UTC'])),flush=True)
        phase = 'BOUNDED_LOCAL_OPAQUE_ORDERED_HTTPS_RANGES_TO_SSH'
        def progress(count):
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            if sys.platform!='darwin':
                rss *= 1024
            assert rss<=plan['local_maximum_RSS_bytes'], 'LOCAL_RSS_RESOURCE_GUARD'
            (task/'private_local_relay_progress_v1.json').write_text(json.dumps(dict(
                status='PUBLIC_COMPRESSED_BYTES_SENT_NOT_YET_WHOLE_SOURCE_VERIFIED',bytes_sent=count,
                observed_UTC=now(),local_max_RSS_bytes=rss,no_local_native_parse_or_payload_files=True))+'\n')
        proof,md5,sha = streamer.stream_ranges(plan['public_content_URL'],child.stdin,
            expected_bytes=plan['expected_file_bytes'],expected_MD5=plan['expected_file_MD5'],
            expected_SHA256=plan['expected_file_SHA256'],chunk_bytes=plan['range_chunk_bytes'],workers=plan['range_workers'],
            request_timeout=plan['range_request_timeout_seconds'],maximum_range_seconds=plan['maximum_range_seconds'],
            maximum_attempts=plan['range_maximum_attempts'],maximum_attempt_seconds=plan['range_maximum_attempt_seconds'],
            retry_backoff_seconds=plan['range_retry_backoff_seconds'],progress=progress)
        footer = dict(kind='COMPLETE_ORDERED_PUBLIC_INPUT_PROOF',MD5=md5,SHA256=sha,range_assembly_proof=proof)
        child.stdin.write((json.dumps(footer)+'\n').encode());child.stdin.close()
        save(task/'private_complete_local_opaque_stream_v1.json',dict(observed_UTC=now(),**footer))
        phase = 'WAIT_SAME_SOURCE_FIXED_R_WORKFLOW_AND_AUDIT'
        signal.alarm(plan['maximum_R_seconds']+90)
    except Exception as error:
        reason = getattr(error,'reason',None)
        local_error = dict(status='LOCAL_OPAQUE_RELAY_NOT_ACCEPTED',observed_UTC=now(),phase=phase,
            error_type=type(error).__name__,stream_error_code=getattr(error,'code_name',None),
            HTTP_code=getattr(error,'code',None),URL_reason_class=type(reason).__name__,
            URL_reason_errno=getattr(reason,'errno',None),SSL_verify_code=getattr(reason,'verify_code',None),
            no_raw_error_or_individual_payload=True)
        save(task/'private_local_relay_failure_v1.json',local_error)
        if not child.stdin.closed:
            try:
                child.stdin.close()
            except BrokenPipeError:
                pass
    try:
        while True:
            item = frame()
            if item is None:
                break
            collect(item)
        exitcode = child.wait(timeout=30)
    except Exception as error:
        if local_error is None:
            local_error = dict(status='LOCAL_SOURCE_RESULT_NOT_ACCEPTED',observed_UTC=now(),
                phase=phase,error_type=type(error).__name__,no_raw_error_or_individual_payload=True)
            save(task/'private_local_observation_failure_v1.json',local_error)
    finally:
        signal.alarm(0)
        if child.poll() is None:
            child.terminate();child.wait(timeout=15)
        exitcode = child.returncode
    if local_error is not None or exitcode!=0 or terminal is None:
        save(task/'private_actual_failure_exit_v1.json',dict(status='FIXED_NATIVE_WORKFLOW_NOT_ACCEPTED',
            observed_UTC=now(),actual_exit_code=exitcode))
        print('{"status":"OPAQUE_RELAY_OR_SOURCE_WORKFLOW_FAILED"}',flush=True)
        raise SystemExit(2)
    assert terminal['source_R_exit_code']==0 and terminal['source_payload_RAM_released'] is True
    save(task/'private_complete_source_fixed_workflow_v1.json',terminal)
    save(task/'private_actual_exit_v1.json',dict(status='SOURCE_EXIT0_FIXED_WORKFLOW_AND_RAM_RELEASE_COMPLETE',
        observed_UTC=now(),actual_exit_code=exitcode,frozen_contract_SHA256=hashlib.sha256(frozen.read_bytes()).hexdigest(),
        complete_source_SHA256=hashlib.sha256((task/'private_complete_source_fixed_workflow_v1.json').read_bytes()).hexdigest()))
    print('{"status":"COMPLETE_OPAQUE_RELAY_FIXED_WORKFLOW_EXIT0"}',flush=True)


if __name__=='__main__':
    main()
