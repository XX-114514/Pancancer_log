"""Receive opaque compressed bytes through existing SSH; verify before R.

Only generic R code is written. Native, individual, coordinate and graph
payloads remain in source RAM. No interactive followup or persistent cache.
"""
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, build_opener, ProxyHandler
from urllib.parse import quote
import hashlib
import io
import json
import os
import resource
import signal
import sys
import tempfile


def now():
    return datetime.now(timezone.utc).isoformat()


def receive_payload(stream,plan):
    n = plan['expected_file_bytes'];assert type(n) is int and 0<n<=2*1024**3
    buffer = io.BytesIO();md5 = hashlib.md5();sha = hashlib.sha256()
    try:
        count = 0
        while count<n:
            block = stream.read(min(65536,n-count))
            if not block:
                raise ValueError('RELAY_INPUT_SHORT')
            count += len(block);buffer.write(block);md5.update(block);sha.update(block)
        line = stream.readline(65537)
        assert 0<len(line)<=65536 and line.endswith(b'\n'), 'INVALID_RELAY_FOOTER'
        footer = json.loads(line)
        assert not stream.read(1), 'EXCESS_RELAY_INPUT'
        assert footer['kind']=='COMPLETE_ORDERED_PUBLIC_INPUT_PROOF'
        assert md5.hexdigest()==footer['MD5']==plan['expected_file_MD5']
        assert sha.hexdigest()==footer['SHA256']==plan['expected_file_SHA256']
        proof = footer['range_assembly_proof']
        assert proof['complete_ranges']==(n+plan['range_chunk_bytes']-1)//plan['range_chunk_bytes']
        assert proof['workers']==plan['range_workers']
        for flag in ('exact_ranges_and_lengths','no_gaps_or_overlaps','whole_official_MD5_verified',
                     'whole_known_SHA256_verified','no_range_payload_files',
                     'no_alternate_endpoint_or_authentication_protocol_hash_retries',
                     'ordered_stream_no_local_file_sized_buffer'):
            assert proof[flag] is True
        proof.update(opaque_compressed_relay_over_existing_SSH=True,
                     whole_bytes_independently_hashed_on_source_before_R=True,
                     range_header_proof_from_frozen_local_streamer=True)
        buffer.seek(n)
        return buffer,md5.hexdigest(),sha.hexdigest(),proof
    except BaseException:
        buffer.close()
        raise


def main(plan):
    assert sys.platform.startswith('linux')
    assert hashlib.sha256(R_SOURCE.encode()).hexdigest()==plan['R_code_SHA256']
    assert hashlib.sha256(PIPE_SOURCE.encode()).hexdigest()==plan['pipe_code_SHA256']
    resource.setrlimit(resource.RLIMIT_AS,(plan['parent_address_cap_bytes'],plan['inherited_address_hard_cap_bytes']))
    def expire(_sig,_frame):
        raise TimeoutError('FINITE_RELAY_OR_ANALYSIS_BUDGET')
    for sig in (signal.SIGALRM,signal.SIGTERM,signal.SIGHUP):
        signal.signal(sig,expire)
    signal.alarm(plan['maximum_download_seconds'])
    started = now();buffer = None;phase = 'PUBLIC_RECORD_AND_SOURCE_RESOURCE_GATES'
    try:
        available = int(next(v.split()[1] for v in Path('/proc/meminfo').read_text().splitlines() if v.startswith('MemAvailable:')))*1024
        assert available>=plan['minimum_available_memory_bytes']
        url = 'https://zenodo.org/api/records/'+plan['record']
        with build_opener(ProxyHandler({})).open(Request(url),timeout=40) as response:
            assert response.status==200 and response.url==url
            raw = response.read(1024**2+1);assert len(raw)<=1024**2
        record = json.loads(raw)
        assert str(record['id'])==plan['record'] and record['metadata']['access_right']=='open'
        assert any(v.get('identifier')==plan['article_DOI'] and v.get('relation')=='isSupplementTo'
                   for v in record['metadata']['related_identifiers'])
        entries = [v for v in record['files'] if hashlib.sha256(v['key'].encode()).hexdigest()==plan['file_key_SHA256']]
        assert len(entries)==1
        entry = entries[0]
        assert entry['size']==plan['expected_file_bytes'] and entry['checksum']=='md5:'+plan['expected_file_MD5']
        expected_url = 'https://zenodo.org/records/'+plan['record']+'/files/'+quote(entry['key'],safe='')+'?download=1'
        assert expected_url==plan['public_content_URL']
        print(json.dumps(dict(kind='SOURCE_OPAQUE_RELAY_READY',started_UTC=started,
            source_pid=os.getpid(),source_start_ticks=int(Path('/proc/self/stat').read_text().split()[21]),
            R_code_SHA256=plan['R_code_SHA256'],source_available_memory_gate_passed=True)),flush=True)
        phase = 'RECEIVE_AND_INDEPENDENTLY_HASH_COMPRESSED_INPUT_RAM_ONLY'
        buffer,md5,sha,proof = receive_payload(sys.stdin.buffer,plan)
        print(json.dumps(dict(kind='COMPLETE_INPUT_IDENTITY_ACCEPTED',observed_UTC=now(),
            bytes=plan['expected_file_bytes'],MD5=md5,SHA256=sha,range_assembly_proof=proof)),flush=True)
        phase = 'GATED_UNIT_PARENT_STATE_FIXED_SPATIAL_AND_AUDIT'
        signal.alarm(plan['maximum_R_seconds']+30)
        R = Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript';assert R.is_file()
        control = Path(tempfile.mkdtemp(prefix='pancancer_opaque_relay_',dir='/dev/shm'));os.chmod(control,0o700)
        script = control/'generic_fixed_workflow.R';script.write_text(R_SOURCE)
        env = dict(os.environ,PANCANCER_EXTERNAL_FIXED_MODE='NATIVE',
                   PANCANCER_EXTERNAL_FIXED_CONFIG=json.dumps(plan['analysis']),
                   R_MAX_VSIZE=str(plan['R_address_cap_bytes']),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        def child_limit():
            resource.setrlimit(resource.RLIMIT_AS,(plan['R_address_cap_bytes'],plan['R_address_cap_bytes']))
        pipe = {};exec(PIPE_SOURCE,pipe)
        child,pipe_proof = pipe['pipe_native_RDS'](buffer,[str(R),'--vanilla',str(script)],
            maximum_expanded_bytes=plan['maximum_expanded_bytes'],timeout=plan['maximum_R_seconds'],
            env=env,preexec_fn=child_limit)
        buffer.close();buffer = None;answer = json.loads(child.stdout)
        if child.returncode!=0 or answer.get('kind')!='SOURCE_ONLY_EXTERNAL_FIXED_WORKFLOW_COMPLETE':
            print(json.dumps(dict(kind='EXTERNAL_FIXED_WORKFLOW_NATIVE_GATE_FAILURE',
                phase=answer.get('phase'),error_class=answer.get('error_class'),source_R_exit_code=child.returncode)),flush=True)
        assert child.returncode==0 and answer['kind']=='SOURCE_ONLY_EXTERNAL_FIXED_WORKFLOW_COMPLETE'
        assert answer['native_RAM_released'] is True
        signal.alarm(0)
        print(json.dumps(dict(kind='COMPLETE_EXTERNAL_FIXED_WORKFLOW_RAM_RELEASED',started_UTC=started,
            completed_UTC=now(),input_bytes=entry['size'],input_MD5=md5,input_SHA256=sha,
            record_metadata_SHA256=hashlib.sha256(raw).hexdigest(),R_code_SHA256=plan['R_code_SHA256'],
            summary=answer['summary'],pipe_technical_proof=pipe_proof,range_assembly_proof=proof,
            source_R_exit_code=child.returncode,source_payload_RAM_released=True,
            no_dynamic_followups_or_native_cache=True,no_native_individual_clinical_coordinate_state_or_graph_files=True,
            no_individual_values_or_identifiers_transferred=True)),flush=True)
    except Exception as error:
        signal.alarm(0)
        if buffer is not None:
            buffer.close()
        print(json.dumps(dict(kind='EXTERNAL_FIXED_WORKFLOW_NOT_ACCEPTED',observed_UTC=now(),phase=phase,
            error_type=type(error).__name__,HTTP_code=getattr(error,'code',None),
            URL_reason_class=type(getattr(error,'reason',None)).__name__,
            URL_reason_errno=getattr(getattr(error,'reason',None),'errno',None),
            no_raw_error_or_individual_payload=True)),flush=True)
        raise SystemExit(2)
