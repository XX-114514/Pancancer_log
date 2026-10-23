"""One declared open IMC object; unit/parent gates and fixed spatial family in RAM.

No images, masks, clinical values, identifiers or native arrays are returned.
Read only the existing user upload; no network or additional native cache.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import resource
import signal
import sys


def now():
    return datetime.now(timezone.utc).isoformat()


def main(plan):
    assert sys.platform.startswith('linux')
    for key, source in [('R_code_SHA256', R_SOURCE), ('pipe_code_SHA256', PIPE_SOURCE),
                        ('upload_reader_code_SHA256', UPLOAD_SOURCE)]:
        assert hashlib.sha256(source.encode()).hexdigest() == plan[key]
    resource.setrlimit(resource.RLIMIT_AS, (plan['parent_address_cap_bytes'], plan['inherited_address_hard_cap_bytes']))
    def expired(_signum, _frame):
        raise TimeoutError('FIXED_SCOPE_ELAPSED_BUDGET')
    signal.signal(signal.SIGALRM, expired)
    signal.signal(signal.SIGTERM, expired)
    signal.signal(signal.SIGHUP, expired)
    signal.alarm(plan['maximum_download_seconds'])
    started = now()
    print(json.dumps(dict(kind='SOURCE_FIXED_WORKFLOW_STARTED', started_UTC=started,
        source_pid=os.getpid(), source_start_ticks=int(Path('/proc/self/stat').read_text().split()[21]))), flush=True)
    phase = 'COMPLETE_USER_UPLOAD_FILE_IDENTITY'
    buffer = None
    try:
        available = int(next(v.split()[1] for v in Path('/proc/meminfo').read_text().splitlines() if v.startswith('MemAvailable:')))*1024
        assert available >= plan['minimum_available_memory_bytes']
        upload_ns = {}
        exec(UPLOAD_SOURCE, upload_ns)
        def progress(count):
            if count == plan['expected_file_bytes'] or count % (32*1024**2) == 0:
                print(json.dumps(dict(kind='PUBLIC_INPUT_RAM_BYTES_PROGRESS',
                    bytes_received=count, observed_UTC=now(), input_mode='EXISTING_USER_UPLOAD')), flush=True)
        buffer, md5, sha, input_proof = upload_ns['read_verified_upload'](plan['uploaded_input_path'],
            expected_bytes=plan['expected_file_bytes'], expected_MD5=plan['expected_file_MD5'],
            expected_SHA256=plan['expected_file_SHA256'], maximum_input_bytes=plan['maximum_input_bytes'],
            maximum_seconds=plan['maximum_download_seconds'], progress=progress)
        print(json.dumps(dict(kind='COMPLETE_INPUT_IDENTITY_ACCEPTED', observed_UTC=now(),
            bytes=buffer.tell(), MD5=md5.hexdigest(), SHA256=sha.hexdigest(),
            uploaded_file_read_proof=input_proof)), flush=True)
        phase = 'GATED_UNIT_PARENT_STATE_FIXED_SPATIAL_AND_AUDIT'
        signal.alarm(plan['maximum_R_seconds']+30)
        R = Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript'
        assert R.is_file()
        bootstrap = Path(plan['source_control_root'])/'generic_fixed_spatial.R'
        assert hashlib.sha256(bootstrap.read_bytes()).hexdigest() == plan['R_code_SHA256']
        env = dict(os.environ, PANCANCER_EXTERNAL_FIXED_MODE='NATIVE', PANCANCER_EXTERNAL_FIXED_CONFIG=json.dumps(plan['analysis']),
                   R_MAX_VSIZE=str(plan['R_address_cap_bytes']),
                   OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        def child_limit():
            resource.setrlimit(resource.RLIMIT_AS, (plan['R_address_cap_bytes'], plan['R_address_cap_bytes']))
        pipe_ns = {}
        exec(PIPE_SOURCE, pipe_ns)
        child, pipe_proof = pipe_ns['pipe_native_RDS'](buffer, [str(R), '--vanilla', str(bootstrap)],
            maximum_expanded_bytes=plan['maximum_expanded_bytes'], timeout=plan['maximum_R_seconds'], env=env, preexec_fn=child_limit)
        buffer.close(); buffer = None
        answer = json.loads(child.stdout)
        if child.returncode != 0 or answer.get('kind') != 'SOURCE_ONLY_EXTERNAL_FIXED_WORKFLOW_COMPLETE':
            print(json.dumps(dict(kind='EXTERNAL_FIXED_WORKFLOW_NATIVE_GATE_FAILURE',
                phase=answer.get('phase'),error_class=answer.get('error_class'),
                source_R_exit_code=child.returncode,no_raw_error_or_individual_payload=True)), flush=True)
        assert child.returncode == 0 and answer['kind'] == 'SOURCE_ONLY_EXTERNAL_FIXED_WORKFLOW_COMPLETE'
        assert answer['native_RAM_released']
        signal.alarm(0)
        print(json.dumps(dict(kind='COMPLETE_EXTERNAL_FIXED_WORKFLOW_RAM_RELEASED',
            started_UTC=started, completed_UTC=now(), input_MD5=md5.hexdigest(), input_SHA256=sha.hexdigest(),
            input_bytes=plan['expected_file_bytes'], record_metadata_SHA256=plan['prior_record_metadata_SHA256'],
            R_code_SHA256=plan['R_code_SHA256'], summary=answer['summary'], pipe_technical_proof=pipe_proof,
            uploaded_file_read_proof=input_proof, source_R_exit_code=child.returncode,
            source_payload_RAM_released=True, no_dynamic_followups_or_additional_native_cache=True,
            existing_user_upload_read_only_and_retained=True,
            no_native_individual_clinical_coordinate_state_or_graph_files=True,
            no_individual_values_or_identifiers_transferred=True)), flush=True)
    except Exception as error:
        signal.alarm(0)
        if buffer is not None:
            buffer.close()
        print(json.dumps(dict(kind='EXTERNAL_FIXED_WORKFLOW_NOT_ACCEPTED', observed_UTC=now(), phase=phase,
            error_type=type(error).__name__, uploaded_input_failure_reason=getattr(error, 'code_name', None),
            no_raw_error_or_individual_payload=True)), flush=True)
        raise SystemExit(2)



if __name__ == '__main__':
    root = Path(sys.argv[1])
    plan = json.loads((root/'source_plan.json').read_text())
    for name, sha in plan['source_files_SHA256'].items():
        assert Path(name).name == name
        assert hashlib.sha256((root/name).read_bytes()).hexdigest() == sha
    R_SOURCE = (root/'generic_fixed_spatial.R').read_text()
    PIPE_SOURCE = (root/'pipe_native_RDS_RAM_v1.py').read_text()
    UPLOAD_SOURCE = (root/'read_uploaded_RDS_RAM_v1.py').read_text()
    assert plan['source_control_root'] == str(root)
    main(plan)
