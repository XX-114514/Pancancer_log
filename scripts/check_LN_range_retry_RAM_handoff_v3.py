"""Source-only synthetic check of the frozen reader and RAM stage handoff.

All HTTP responses are supplied by an in-memory fake opener. Only generic
control code, contracts and synthetic aggregate JSON are written. This never
fetches a real object, reads a real individual record or modifies a live job.
"""
from pathlib import Path
from datetime import datetime, timezone
import base64
import hashlib
import io
import json
import os
import subprocess
import sys
import tempfile
import threading
import time
from contextlib import redirect_stdout
from urllib.error import HTTPError
from types import SimpleNamespace


def main():
    r = Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript'
    make = '''
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(jsonlite))
x<-SingleCellExperiment(assays=list(counts=matrix(seq_len(60),6,10,
 dimnames=list(c("CD3","CD4","CD8a","FOXP3","PD1","Ki67"),NULL))),
 colData=S4Vectors::DataFrame(patient_id=rep(c("FAKE_A","FAKE_B"),each=5)))
cat(toJSON(list(raw=base64_enc(serialize(x,NULL))),auto_unbox=TRUE))
'''
    made = subprocess.run([str(r), '--vanilla', '-e', make], capture_output=True, timeout=90)
    if made.returncode:
        raise RuntimeError('SYNTHETIC_HANDOFF_GENERATOR_FAILED')
    payload = base64.b64decode(json.loads(made.stdout)['raw'])
    key, record_id = 'synthetic.rds', 'SYNTHETIC_NO_NETWORK'
    md5 = hashlib.md5(payload).hexdigest()
    metadata_url = 'https://zenodo.org/api/records/'+record_id
    content_url = 'https://zenodo.org/records/'+record_id+'/files/synthetic.rds?download=1'
    test_mode = ['OK']
    attempt_counts, request_lock = {}, threading.Lock()
    metadata = json.dumps(dict(id=record_id, metadata=dict(access_right='open'),
        files=[dict(key=key, size=len(payload), checksum='md5:'+md5,
                    links=dict(self=content_url))])).encode()
    requested, controls, thread_errors = [], [], []

    class Response(io.BytesIO):
        def __init__(self, value, url, status=200, content_range=None, declared_length=None):
            super().__init__(value)
            self.status, self.url = status, url
            self.headers = {'Content-Length': str(len(value) if declared_length is None else declared_length),
                            'Content-Type': 'application/octet-stream'}
            if content_range is not None:
                self.headers['Content-Range'] = content_range

    class FakeOpener:
        def open(self, request, timeout):
            span = request.get_header('Range')
            requested.append((request.full_url, span))
            assert request.full_url in [metadata_url, content_url]
            if request.full_url == metadata_url:
                return Response(metadata, request.full_url)
            assert span and span.startswith('bytes=')
            start, end = map(int, span[6:].split('-'))
            with request_lock:
                attempt_counts[span] = attempt_counts.get(span,0)+1
                attempt = attempt_counts[span]
            if test_mode[0] == 'HTTP403':
                raise HTTPError(request.full_url,403,'SYNTHETIC_ONLY',{},None)
            body = payload[start:end+1]
            cr = f'bytes {start}-{end}/{len(payload)}'
            if test_mode[0] == 'BAD_HEADER': cr = 'bytes 0-0/1'
            declared = len(body)
            if test_mode[0] == 'TRUNCATED' or (test_mode[0] == 'TRANSIENT_SHORT_FIRST' and attempt == 1): body = body[:-1]
            status = 200 if test_mode[0] == 'STATUS200' else 206
            return Response(body, request.full_url, status=status, content_range=cr, declared_length=declared)


    def stage(root):
        try:
            deadline = time.monotonic()+180
            while not (root/'private_native_schema_ready.json').exists():
                if time.monotonic() > deadline:
                    raise TimeoutError('SYNTHETIC_SCHEMA_DEADLINE')
                time.sleep(0.1)
            schema = json.loads((root/'private_native_schema_ready.json').read_text())
            assert schema['status'] == 'NATIVE_SCHEMA_READY_SOURCE_OBJECT_HELD_RAM_ONLY'
            assert schema['summary']['source_cells'] == 10
            code = '''
stopifnot(attr(attr(native_object,"colData",exact=TRUE),"nrows",exact=TRUE)==10L)
writeLines(jsonlite::toJSON(list(status="SYNTHETIC_RAM_STAGE_VERIFIED"),auto_unbox=TRUE),
 file.path(control_root,"private_synthetic_stage_receipt.json"))
'''
            (root/'phase001.R').write_text(code)
            (root/'phase001.json').write_text(json.dumps(dict(
                code_SHA256=hashlib.sha256(code.encode()).hexdigest())))
            while not (root/'private_phase001_complete.json').exists():
                if time.monotonic() > deadline:
                    raise TimeoutError('SYNTHETIC_STAGE_DEADLINE')
                time.sleep(0.1)
            assert json.loads((root/'private_synthetic_stage_receipt.json').read_text())['status'] == 'SYNTHETIC_RAM_STAGE_VERIFIED'
            (root/'release.json').write_text(json.dumps(dict(synthetic_release=True)))
        except BaseException as error:
            thread_errors.append(type(error).__name__)
            (root/'release.json').write_text(json.dumps(dict(synthetic_release_after_failure=True)))

    def control_dir(**kwargs):
        root = Path(tempfile.mkdtemp(prefix='pancancer_LN_SYNTHETIC_control_', dir='/dev/shm'))
        controls.append(root)
        t = threading.Thread(target=stage, args=(root,), daemon=True)
        t.start()
        return str(root)

    plan = dict(record=record_id, expected_file_bytes=len(payload), expected_file_MD5=md5,
        file_key_SHA256=hashlib.sha256(key.encode()).hexdigest(),
        R_schema_code_SHA256=hashlib.sha256(R_SOURCE.encode()).hexdigest(),
        pipe_code_SHA256=hashlib.sha256(PIPE_SOURCE.encode()).hexdigest(),
        parent_address_cap_bytes=4*1024**3, inherited_address_hard_cap_bytes=16*1024**3,
        R_address_cap_bytes=12*1024**3, minimum_available_memory_bytes=24*1024**3,
        maximum_download_seconds=180, maximum_R_seconds=180, RAM_handoff_seconds=30,
        maximum_expanded_bytes=64*1024**2, previous_source_pid=0, previous_source_start_ticks=0,
        range_code_SHA256=hashlib.sha256(RANGE_SOURCE.encode()).hexdigest(),
        range_chunk_bytes=1024, range_workers=4, range_request_timeout_seconds=10, maximum_range_seconds=15,
        range_core_code_SHA256=hashlib.sha256(RANGE_CORE_SOURCE.encode()).hexdigest(),
        range_maximum_attempts=4, range_retry_backoff_seconds=0)
    ns = {'__name__': 'synthetic_frozen_reader', 'R_SOURCE': R_SOURCE,
          'PIPE_SOURCE': PIPE_SOURCE, 'RANGE_SOURCE': RANGE_SOURCE, 'RANGE_CORE_SOURCE': RANGE_CORE_SOURCE,
          'RANGE_OPENER_FACTORY': lambda: FakeOpener()}
    exec(READER_SOURCE, ns)
    ns['build_opener'] = lambda *args, **kwargs: FakeOpener()
    ns['tempfile'] = SimpleNamespace(mkdtemp=control_dir)
    prior_stdin, captured = sys.stdin, io.StringIO()
    try:
        sys.stdin = io.StringIO(json.dumps(plan))
        with redirect_stdout(captured):
            ns['main']()
    except BaseException as error:
        frames = [json.loads(line) for line in captured.getvalue().splitlines()]
        safe = [{k: v for k, v in f.items() if k in
                 ['kind', 'phase', 'error_type', 'error_class', 'source_R_exit_code', 'range_failure_reason']}
                for f in frames]
        print(json.dumps(dict(status='SYNTHETIC_FROZEN_READER_NOT_ACCEPTED',
            exception_type=type(error).__name__, safe_frames=safe,
            thread_error_classes=thread_errors, no_human_data_read=True)))
        raise SystemExit(2)
    finally:
        sys.stdin = prior_stdin
    frames = [json.loads(line) for line in captured.getvalue().splitlines()]
    assert not thread_errors and requested[0] == (metadata_url, None)
    expected_spans = {f'bytes={start}-{min(len(payload)-1,start+1023)}' for start in range(0,len(payload),1024)}
    assert len(requested) == 1+len(expected_spans)
    assert {span for url,span in requested[1:]} == expected_spans
    assert all(url == content_url for url,span in requested[1:])
    final = frames[-1]
    assert final['kind'] == 'COMPLETE_PUBLIC_PROTEIN_RDS_RAM_STRUCTURE_NO_EFFECTS'
    assert final['source_payload_RAM_released'] and final['source_R_exit_code'] == 0
    assert final['summary']['completed_control_phases'] == 1
    allowed = {'private_verified_public_input_identity.json', 'private_native_schema_ready.json',
        'private_control_heartbeat.json', 'phase001.R', 'phase001.json',
        'private_phase001_complete.json', 'private_synthetic_stage_receipt.json', 'release.json', 'native_reader.R'}
    assert len(controls) == 1 and {p.name for p in controls[0].iterdir()} <= allowed
    assert 'FAKE_A' not in captured.getvalue()
    range_ns = {"RANGE_CORE_SOURCE": RANGE_CORE_SOURCE}
    exec(RANGE_SOURCE, range_ns)
    failures = []
    for mode, expected in [('BAD_HEADER','CONTENT_RANGE_MISMATCH'),('TRUNCATED','RANGE_LENGTH_SHORT'),
                           ('STATUS200','NOT_EXACT_OFFICIAL_206'),('HASH','WHOLE_OBJECT_MD5_MISMATCH')]:
        test_mode[0] = mode
        attempt_counts.clear()
        try:
            range_ns['fetch_ranges'](content_url, expected_bytes=len(payload),
                expected_MD5='0'*32 if mode=='HASH' else md5, chunk_bytes=1024, workers=4,
                request_timeout=10, maximum_range_seconds=15, opener_factory=lambda: FakeOpener(),
                maximum_attempts=4, retry_backoff_seconds=0)
        except range_ns['RangeInputError'] as error:
            assert error.code_name == expected
            failures.append(expected)
            assert max(attempt_counts.values()) <= (4 if mode=='TRUNCATED' else 1)
        else:
            raise AssertionError('NEGATIVE_RANGE_INTEGRITY_CASE_NOT_REJECTED')
    test_mode[0] = 'TRANSIENT_SHORT_FIRST'
    attempt_counts.clear()
    recovered, dm, ds, proof = range_ns['fetch_ranges'](content_url, expected_bytes=len(payload),
        expected_MD5=md5, chunk_bytes=1024, workers=4, request_timeout=10,
        maximum_range_seconds=15, opener_factory=lambda: FakeOpener(), maximum_attempts=4, retry_backoff_seconds=0)
    assert recovered.getvalue() == payload and dm.hexdigest() == md5
    assert proof['ranges_requiring_retry'] == len(expected_spans)
    assert proof['maximum_actual_attempts_per_range'] == 2 and proof['actual_request_attempts'] == 2*len(expected_spans)
    recovered.close()
    test_mode[0] = 'HTTP403'
    attempt_counts.clear()
    try:
        range_ns['fetch_ranges'](content_url, expected_bytes=len(payload), expected_MD5=md5,
            chunk_bytes=1024, workers=4, request_timeout=10, maximum_range_seconds=15,
            opener_factory=lambda: FakeOpener(), maximum_attempts=4, retry_backoff_seconds=0)
    except HTTPError as error:
        assert error.code == 403 and max(attempt_counts.values()) == 1
    else:
        raise AssertionError('AUTHENTICATION_ERROR_WAS_NOT_TERMINAL')
    print(json.dumps(dict(status='PASS_BOUNDED_TRANSPORT_RETRY_NATIVE_RAM_HANDOFF_INTEGRITY_AND_NO_AUTH_RETRY',
        observed_UTC=datetime.now(timezone.utc).isoformat(),
        reader_SHA256=hashlib.sha256(READER_SOURCE.encode()).hexdigest(),
        pipe_SHA256=hashlib.sha256(PIPE_SOURCE.encode()).hexdigest(),
        range_SHA256=hashlib.sha256(RANGE_SOURCE.encode()).hexdigest(), range_core_SHA256=hashlib.sha256(RANGE_CORE_SOURCE.encode()).hexdigest(),
        rejected_integrity_cases=failures, transient_short_recovered_with_identical_whole_MD5=True, HTTP403_not_retried=True,
        no_real_source_request_or_individual_payload_files=True,
        synthetic_stage_SHA_verified=True, original_live_job_untouched=True)))


if __name__ == '__main__':
    main()
