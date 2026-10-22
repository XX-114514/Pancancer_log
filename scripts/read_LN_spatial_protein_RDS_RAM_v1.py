"""Source-side public protein object structure qualification, RAM only.

Fetches directly from the declared open Zenodo record. No patient/omics object
is uploaded from the client or returned from this process. No source files.
The distinct previously denied clinical CSV and TCR transport are untouched.
"""
from pathlib import Path
from datetime import datetime, timezone
from urllib.request import Request, build_opener, ProxyHandler
from urllib.parse import urlparse
import hashlib
import io
import json
import os
import resource
import signal
import subprocess
import sys
import time


def now():
    return datetime.now(timezone.utc).isoformat()


def main():
    plan = json.load(sys.stdin)
    assert sys.platform.startswith('linux')
    assert hashlib.sha256(R_SOURCE.encode()).hexdigest() == plan['R_schema_code_SHA256']
    resource.setrlimit(resource.RLIMIT_AS, (plan['parent_address_cap_bytes'], plan['inherited_address_hard_cap_bytes']))
    def expired(_signum, _frame):
        raise TimeoutError('FIXED_ELAPSED_BUDGET')
    signal.signal(signal.SIGALRM, expired)
    signal.alarm(plan['maximum_elapsed_seconds'])
    started = now()
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
        content = entry['links']['self']
        assert urlparse(content).scheme == 'https' and urlparse(content).hostname == 'zenodo.org'
        phase = 'PUBLIC_PROTEIN_OBJECT_COMPLETE_RAM_ONLY'
        buffer, md5, sha = io.BytesIO(), hashlib.md5(), hashlib.sha256()
        with opener.open(Request(content), timeout=45) as response:
            assert response.status == 200 and urlparse(response.url).hostname == 'zenodo.org'
            assert 'text/html' not in response.headers.get('Content-Type', '').lower()
            length = response.headers.get('Content-Length')
            if length is not None:
                assert int(length) == plan['expected_file_bytes']
            last = 0
            while True:
                value = response.read(1024**2)
                if not value:
                    break
                buffer.write(value)
                md5.update(value)
                sha.update(value)
                assert buffer.tell() <= plan['expected_file_bytes']
                if buffer.tell()-last >= 128*1024**2:
                    last = buffer.tell()
                    print(json.dumps(dict(kind='PUBLIC_INPUT_RAM_BYTES_PROGRESS', bytes_received=last, observed_UTC=now())), flush=True)
        assert buffer.tell() == plan['expected_file_bytes'] and md5.hexdigest() == plan['expected_file_MD5']
        phase = 'R_NATIVE_STRUCTURE_ONLY_NO_EFFECTS'
        R = Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript'
        assert R.is_file()
        rcode = R_SOURCE + r'''
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(jsonlite))
phase<-"READ_ONE_PUBLIC_SPATIAL_OBJECT"
tryCatch({
 con<-gzcon(file("stdin","rb"));x<-readRDS(con);close(con)
 phase<-"SOURCE_ONLY_STRUCTURE_AND_MARKER_QUALIFICATION"
 answer<-qualify_structure(x)
 answer$runtime_R<-as.character(getRversion())
 answer$native_object_declared_bytes<-as.numeric(object.size(x))
 rm(x);gc(verbose=FALSE)
 cat(toJSON(list(kind="SOURCE_ONLY_STRUCTURE_SUMMARY",summary=answer),auto_unbox=TRUE,null="null",digits=16))
},error=function(e){
 cat(toJSON(list(kind="SOURCE_RDS_STRUCTURE_FAILED",phase=phase,error_class=class(e)[1]),auto_unbox=TRUE))
 quit(status=2L)
})
'''
        env = dict(os.environ, R_MAX_VSIZE=str(plan['R_address_cap_bytes']), OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        def child_limit():
            resource.setrlimit(resource.RLIMIT_AS, (plan['R_address_cap_bytes'], plan['R_address_cap_bytes']))
        child = subprocess.run([str(R), '--vanilla', '-e', rcode], input=buffer.getvalue(), stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL, env=env, preexec_fn=child_limit, timeout=plan['maximum_R_seconds'])
        buffer.close()
        assert len(child.stdout) <= 1024**2
        answer = json.loads(child.stdout)
        if child.returncode != 0:
            safe = dict(kind='SOURCE_RDS_STRUCTURE_FAILED', source_R_exit_code=child.returncode,
                        phase=answer.get('phase') if answer.get('phase') in ['READ_ONE_PUBLIC_SPATIAL_OBJECT', 'SOURCE_ONLY_STRUCTURE_AND_MARKER_QUALIFICATION'] else 'UNREVIEWED_PHASE',
                        error_class=answer.get('error_class') if answer.get('error_class') in ['simpleError', 'error', 'packageNotFoundError'] else 'UNREVIEWED_ERROR_CLASS')
            print(json.dumps(safe), flush=True)
        assert child.returncode == 0 and answer['kind'] == 'SOURCE_ONLY_STRUCTURE_SUMMARY'
        signal.alarm(0)
        proof = dict(kind='COMPLETE_PUBLIC_PROTEIN_RDS_RAM_STRUCTURE_NO_EFFECTS', started_UTC=started, completed_UTC=now(),
            source_object_MD5=md5.hexdigest(), source_object_SHA256=sha.hexdigest(), source_object_bytes=entry['size'],
            R_schema_code_SHA256=plan['R_schema_code_SHA256'], summary=answer['summary'],
            source_R_exit_code=child.returncode, source_payload_RAM_released=True,
            source_parent_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            source_child_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*1024,
            no_source_object_individual_ID_protein_coordinate_or_clinical_record_files=True,
            no_source_arrays_or_clinical_outcome_values_returned=True,
            original_project_objects_and_denied_actions_untouched=True)
        print(json.dumps(proof), flush=True)
    except Exception as error:
        signal.alarm(0)
        print(json.dumps(dict(kind='PUBLIC_PROTEIN_RDS_SCOPE_NOT_ACCEPTED', phase=phase, observed_UTC=now(),
            error_type=type(error).__name__, HTTP_code=getattr(error, 'code', None), no_raw_error_or_individual_payload=True)), flush=True)
        raise SystemExit(2)


if __name__ == '__main__':
    main()
