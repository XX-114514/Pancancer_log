"""One declared open IMC object; unit/parent gates and fixed spatial family in RAM.

No images, masks, clinical values, identifiers or native arrays are returned.
No followup interface or object cache. Authentication/HTTP denial is terminal.
"""
from pathlib import Path
from datetime import datetime, timezone
from urllib.request import Request, build_opener, ProxyHandler
from urllib.parse import quote
import hashlib
import json
import os
import resource
import signal
import sys
import tempfile


def now():
    return datetime.now(timezone.utc).isoformat()


def main(plan):
    assert sys.platform.startswith('linux')
    for key, source in [('R_code_SHA256', R_SOURCE), ('pipe_code_SHA256', PIPE_SOURCE),
                        ('range_code_SHA256', RANGE_SOURCE), ('range_core_code_SHA256', RANGE_CORE_SOURCE)]:
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
    phase = 'PUBLIC_DATASET_RECORD_IDENTITY'
    buffer = None
    try:
        available = int(next(v.split()[1] for v in Path('/proc/meminfo').read_text().splitlines() if v.startswith('MemAvailable:')))*1024
        assert available >= plan['minimum_available_memory_bytes']
        opener = build_opener(ProxyHandler({}))
        url = 'https://zenodo.org/api/records/'+plan['record']
        with opener.open(Request(url), timeout=40) as response:
            assert response.status == 200 and response.url == url
            raw = response.read(1024**2+1)
            assert len(raw) <= 1024**2
        record = json.loads(raw)
        assert str(record['id']) == plan['record'] and record['metadata']['access_right'] == 'open'
        assert any(v.get('identifier') == plan['article_DOI'] and v.get('relation') == 'isSupplementTo'
                   for v in record['metadata']['related_identifiers'])
        entries = [v for v in record['files'] if hashlib.sha256(v['key'].encode()).hexdigest() == plan['file_key_SHA256']]
        assert len(entries) == 1
        entry = entries[0]
        assert entry['size'] == plan['expected_file_bytes']
        assert entry['checksum'] == 'md5:'+plan['expected_file_MD5']
        content = 'https://zenodo.org/records/'+plan['record']+'/files/'+quote(entry['key'], safe='')+'?download=1'
        phase = 'BOUNDED_COMPLETE_PUBLIC_DATASET_INPUT_RAM_ONLY'
        range_ns = {'RANGE_CORE_SOURCE': RANGE_CORE_SOURCE}
        exec(RANGE_SOURCE, range_ns)
        def progress(count):
            print(json.dumps(dict(kind='PUBLIC_INPUT_RAM_BYTES_PROGRESS', bytes_received=count, observed_UTC=now())), flush=True)
        buffer, md5, sha, range_proof = range_ns['fetch_ranges'](content,
            expected_bytes=plan['expected_file_bytes'], expected_MD5=plan['expected_file_MD5'],
            chunk_bytes=plan['range_chunk_bytes'], workers=plan['range_workers'],
            request_timeout=plan['range_request_timeout_seconds'], maximum_range_seconds=plan['maximum_range_seconds'],
            progress=progress, maximum_attempts=plan['range_maximum_attempts'],
            retry_backoff_seconds=plan['range_retry_backoff_seconds'], maximum_attempt_seconds=plan['range_maximum_attempt_seconds'])
        assert sha.hexdigest() == plan['expected_file_SHA256']
        print(json.dumps(dict(kind='COMPLETE_INPUT_IDENTITY_ACCEPTED', observed_UTC=now(),
            bytes=buffer.tell(), MD5=md5.hexdigest(), SHA256=sha.hexdigest(), range_assembly_proof=range_proof)), flush=True)
        phase = 'GATED_UNIT_PARENT_STATE_FIXED_SPATIAL_AND_AUDIT'
        signal.alarm(plan['maximum_R_seconds']+30)
        R = Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript'
        assert R.is_file()
        control = Path(tempfile.mkdtemp(prefix='pancancer_external_TDLN_fixed_', dir='/dev/shm'))
        os.chmod(control, 0o700)
        bootstrap = control/'generic_fixed_workflow.R'
        bootstrap.write_text(R_SOURCE)
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
            input_bytes=entry['size'], record_metadata_SHA256=hashlib.sha256(raw).hexdigest(),
            R_code_SHA256=plan['R_code_SHA256'], summary=answer['summary'], pipe_technical_proof=pipe_proof,
            range_assembly_proof=range_proof, source_R_exit_code=child.returncode,
            source_payload_RAM_released=True, no_dynamic_followups_or_native_cache=True,
            no_native_individual_clinical_coordinate_state_or_graph_files=True,
            no_individual_values_or_identifiers_transferred=True)), flush=True)
    except Exception as error:
        signal.alarm(0)
        if buffer is not None:
            buffer.close()
        print(json.dumps(dict(kind='EXTERNAL_FIXED_WORKFLOW_NOT_ACCEPTED', observed_UTC=now(), phase=phase,
            error_type=type(error).__name__, HTTP_code=getattr(error, 'code', None),
            range_failure_reason=getattr(error, 'code_name', None),
            range_start=getattr(error, 'range_start', None), range_end=getattr(error, 'range_end', None),
            expected_range_bytes=getattr(error, 'expected_range_bytes', None), received_range_bytes=getattr(error, 'received_range_bytes', None),
            no_raw_error_or_individual_payload=True)), flush=True)
        raise SystemExit(2)

