"""Source-side frozen secondary specificity workflow, RAM only; no dynamic followups.

Fetches directly from the declared open Zenodo record. No patient/omics object
is uploaded from the client or returned from this process. Only generic control code, manifests and aggregate JSON may be files.
The distinct previously denied clinical CSV and TCR transport are untouched.
"""
from pathlib import Path
from datetime import datetime, timezone
from urllib.request import Request, build_opener, ProxyHandler
from urllib.parse import urlparse, quote
import hashlib
import io
import json
import os
import resource
import signal
import subprocess
import sys
import time
import tempfile


def now():
    return datetime.now(timezone.utc).isoformat()


def main():
    plan = json.load(sys.stdin)
    assert sys.platform.startswith('linux')
    assert hashlib.sha256(R_SOURCE.encode()).hexdigest() == plan['R_code_SHA256']
    assert hashlib.sha256(PIPE_SOURCE.encode()).hexdigest() == plan['pipe_code_SHA256']
    assert hashlib.sha256(RANGE_SOURCE.encode()).hexdigest() == plan['range_code_SHA256']
    assert hashlib.sha256(RANGE_CORE_SOURCE.encode()).hexdigest() == plan['range_core_code_SHA256']
    resource.setrlimit(resource.RLIMIT_AS, (plan['parent_address_cap_bytes'], plan['inherited_address_hard_cap_bytes']))
    def expired(_signum, _frame):
        raise TimeoutError('FIXED_ELAPSED_BUDGET')
    signal.signal(signal.SIGALRM, expired)
    signal.alarm(plan['maximum_download_seconds'])
    started = now()
    previous = Path('/proc')/str(plan['previous_source_pid'])
    if previous.exists():
        assert int((previous/'stat').read_text().rsplit(')',1)[1].split()[19]) != plan['previous_source_start_ticks'], 'PREVIOUS_SOURCE_WORKER_STILL_LIVE'
    previous_R = Path('/proc')/str(plan['previous_R_pid'])
    if previous_R.exists():
        assert int((previous_R/'stat').read_text().rsplit(')',1)[1].split()[19]) != plan['previous_R_start_ticks'], 'PREVIOUS_R_STILL_LIVE'
    control = Path(tempfile.mkdtemp(prefix='pancancer_LN_specificity_', dir='/dev/shm'))
    os.chmod(control, 0o700)
    print(json.dumps(dict(kind='SOURCE_RAM_CONTROL_ROOT_CREATED', source_control_root=str(control), started_UTC=started)), flush=True)
    print(json.dumps(dict(kind='SOURCE_RAM_WORKER_STARTED', source_pid=os.getpid(),
        source_start_ticks=int(Path('/proc/self/stat').read_text().split()[21]), started_UTC=started)), flush=True)
    phase = 'PUBLIC_RECORD_IDENTITY'
    try:
        available = int(next(v.split()[1] for v in Path('/proc/meminfo').read_text().splitlines() if v.startswith('MemAvailable:')))*1024
        assert available >= plan['minimum_available_memory_bytes']
        # Per-request direct standard TLS. No global environment, proxy, trust
        # chain, DNS, forwarding or authentication setting is changed.
        opener = build_opener(ProxyHandler({}))
        url = 'https://zenodo.org/api/records/'+plan['record']
        with opener.open(Request(url), timeout=40) as response:
            assert response.status == 200 and response.url == url
            raw = response.read(1024**2+1)
            assert len(raw) <= 1024**2
        record = json.loads(raw)
        assert str(record['id']) == plan['record'] and record['metadata']['access_right'] == 'open'
        entries = [v for v in record['files'] if v['key'].lower().endswith('.rds')]
        assert len(entries) == 1
        entry = entries[0]
        assert entry['size'] == plan['expected_file_bytes']
        assert entry['checksum'] == 'md5:'+plan['expected_file_MD5']
        assert hashlib.sha256(entry['key'].encode()).hexdigest() == plan['file_key_SHA256']
        content = 'https://zenodo.org/records/'+plan['record']+'/files/'+quote(entry['key'],safe='')+'?download=1'
        phase = 'PUBLIC_PROTEIN_OBJECT_RANGES_COMPLETE_RAM_ONLY'
        range_ns = {"RANGE_CORE_SOURCE": RANGE_CORE_SOURCE}
        exec(RANGE_SOURCE, range_ns)
        def progress(count):
            print(json.dumps(dict(kind='PUBLIC_INPUT_RAM_BYTES_PROGRESS', bytes_received=count, observed_UTC=now())), flush=True)
        buffer, md5, sha, range_proof = range_ns['fetch_ranges'](content,
            expected_bytes=plan['expected_file_bytes'], expected_MD5=plan['expected_file_MD5'],
            chunk_bytes=plan['range_chunk_bytes'], workers=plan['range_workers'],
            request_timeout=plan['range_request_timeout_seconds'],
            maximum_range_seconds=plan['maximum_range_seconds'], progress=progress,
            opener_factory=globals().get('RANGE_OPENER_FACTORY'),
            maximum_attempts=plan['range_maximum_attempts'], retry_backoff_seconds=plan['range_retry_backoff_seconds'],
            maximum_attempt_seconds=plan['range_maximum_attempt_seconds'])
        assert sha.hexdigest() == plan['expected_file_SHA256'], 'WHOLE_OBJECT_SHA256_MISMATCH'
        (control/'private_verified_public_input_identity.json').write_text(json.dumps(dict(MD5=md5.hexdigest(), SHA256=sha.hexdigest(), bytes=buffer.tell(), R_code_SHA256=plan['R_code_SHA256'], complete_UTC=now(), no_native_object_file=True), indent=2)+'\n')
        signal.alarm(plan['maximum_R_seconds']+30)
        phase = 'R_FROZEN_SECONDARY_SPECIFICITY_AND_NUMERIC_AUDIT'
        R = Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript'
        assert R.is_file()
        rcode = R_SOURCE
        r_deadline = datetime.fromtimestamp(time.time()+plan['maximum_R_seconds'], timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
        env = dict(os.environ, PANCANCER_LN_CONTROL_ROOT=str(control),
                   PANCANCER_LN_SPEC_CONFIG=json.dumps(plan['analysis_config']),
                   PANCANCER_LN_R_DEADLINE_UTC=r_deadline,
                   R_MAX_VSIZE=str(plan['R_address_cap_bytes']), OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        def child_limit():
            resource.setrlimit(resource.RLIMIT_AS, (plan['R_address_cap_bytes'], plan['R_address_cap_bytes']))
        pipe_ns = {}
        exec(PIPE_SOURCE, pipe_ns)
        bootstrap = control/'native_reader.R'
        bootstrap.write_text(rcode)
        child, pipe_proof = pipe_ns['pipe_native_RDS'](buffer, [str(R), '--vanilla', str(bootstrap)],
            maximum_expanded_bytes=plan['maximum_expanded_bytes'], timeout=plan['maximum_R_seconds'],
            env=env, preexec_fn=child_limit)
        buffer.close()
        assert len(child.stdout) <= 1024**2
        answer = json.loads(child.stdout)
        if child.returncode != 0:
            safe = dict(kind='SOURCE_SPECIFICITY_FAILED', source_R_exit_code=child.returncode,
                        phase=answer.get('phase') if answer.get('phase') in ['READ_SAME_VERIFIED_PUBLIC_PROTEIN_OBJECT','FROZEN_PHENOTYPE_REPLAY','FINITE_SECONDARY_SPECIFICITY_PRODUCTION','FULL_SECONDARY_NUMERIC_REPLAY_AUDIT'] else 'UNREVIEWED_PHASE',
                        error_class=answer.get('error_class') if answer.get('error_class') in ['simpleError', 'error', 'packageNotFoundError'] else 'UNREVIEWED_ERROR_CLASS')
            print(json.dumps(safe), flush=True)
        assert child.returncode == 0 and answer['kind'] == 'SOURCE_ONLY_SPECIFICITY_SUMMARY'
        signal.alarm(0)
        proof = dict(kind='COMPLETE_PUBLIC_PROTEIN_SECONDARY_SPECIFICITY_RAM_RELEASED', started_UTC=started, completed_UTC=now(),
            source_object_MD5=md5.hexdigest(), source_object_SHA256=sha.hexdigest(), source_object_bytes=entry['size'],
            R_code_SHA256=plan['R_code_SHA256'], summary=answer['summary'],
            pipe_code_SHA256=plan['pipe_code_SHA256'], pipe_technical_proof=pipe_proof,
            native_reader_code_SHA256=hashlib.sha256(rcode.encode()).hexdigest(),
            range_code_SHA256=plan['range_code_SHA256'], range_assembly_proof=range_proof,
            source_R_exit_code=child.returncode, source_payload_RAM_released=True,
            source_parent_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            source_child_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*1024,
            no_source_object_individual_ID_protein_coordinate_or_clinical_record_files=True,
            no_source_arrays_or_clinical_outcome_values_returned=True,
            original_project_objects_and_denied_actions_untouched=True,
            no_dynamic_followup_nodes=True, frozen_secondary_family_SHA256=plan['analysis_config']['family_contract_SHA256'])
        print(json.dumps(proof), flush=True)
    except Exception as error:
        signal.alarm(0)
        print(json.dumps(dict(kind='PUBLIC_PROTEIN_RDS_SCOPE_NOT_ACCEPTED', phase=phase, observed_UTC=now(),
            error_type=type(error).__name__, HTTP_code=getattr(error, 'code', None),
            range_failure_reason=getattr(error,'code_name',None),
            range_start=getattr(error,'range_start',None), range_end=getattr(error,'range_end',None),
            expected_range_bytes=getattr(error,'expected_range_bytes',None), received_range_bytes=getattr(error,'received_range_bytes',None),
            range_attempt_number=getattr(error,'range_attempt_number',None),
            range_elapsed_seconds=getattr(error,'range_elapsed_seconds',None), no_raw_error_or_individual_payload=True)), flush=True)
        raise SystemExit(2)


if __name__ == '__main__':
    main()
