"""DISABLED historical remote transport, rejected before execution.

Use local-only nonsequence review. The rejected frozen source is archived
privately with its original SHA; this entry point cannot perform transport.

Only source aggregates are persisted or emitted. Per-material numeric payload
stays private RAM; original identifiers and sequence fields never reach it.
"""
from pathlib import Path
from datetime import datetime,timezone
from urllib.request import Request,ProxyHandler,build_opener
import hashlib
import importlib
import io
import json
import os
import shlex
import subprocess
import time


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def save(path,obj):
    assert not path.exists();path.write_text(json.dumps(obj,indent=2)+'\n')


def main():
    raise SystemExit('REMOTE_TCR_TRANSPORT_DISABLED_AFTER_APPROVAL_REJECTION_USE_LOCAL_ONLY_REVIEW')
    task=Path(os.environ['PANCANCER_TCR_TASK']);code=Path(__file__).parent
    frozen=task/'private_frozen_measurement_v1.json';plan=json.loads(frozen.read_text());item=plan['source']
    for name,digest in plan['code_SHA256'].items():assert sha(code/name)==digest
    program="__name__='embedded_TCR_qualification'\n"+(code/'NSCLC_author_state_counts_RAM_v1.py').read_text()+'\nread_RNA_source=read_source\n'
    program+=(code/'qualify_NSCLC_reported_TCR_expansion_RAM_v1.py').read_text()+'\n'
    program+='PLAN='+repr(plan)+'\nROOT='+repr(os.environ['PANCANCER_SOURCE_ROOT'])+'\n'
    program+='''
import sys
try:
 synthetic()
 available=int(next(v.split()[1] for v in Path('/proc/meminfo').read_text().splitlines() if v.startswith('MemAvailable:')))*1024
 assert available>=4*1024**3
 print(json.dumps(dict(kind='PURE_TCR_FIXTURE_ACCEPTED')),flush=True)
 manifest=json.loads(sys.stdin.buffer.readline());assert manifest['bytes']==PLAN['source']['bytes']
 buffer=io.BytesIO();digest=hashlib.sha256();remaining=manifest['bytes']
 while remaining:
  value=sys.stdin.buffer.read(min(1024**2,remaining))
  if not value:raise EOFError('INPUT_CHANNEL_CLOSED')
  buffer.write(value);digest.update(value);remaining-=len(value)
 assert digest.hexdigest()==manifest['SHA256']
 payload,summary=qualify(buffer.getvalue(),ROOT,PLAN)
 buffer.close()
 print(json.dumps(dict(kind='PRIVATE_MEASUREMENT_HANDOFF',summary=summary,payload=payload,transport_SHA256=digest.hexdigest())),flush=True)
except Exception as error:
 print(json.dumps(dict(kind='SOURCE_TECHNICAL_FAILURE',error_type=type(error).__name__)),flush=True)
 raise SystemExit(2)
'''
    assert len(program.encode())<110000
    args=['ssh','-T','-o','ClearAllForwardings=yes','-o','BatchMode=yes','-o','ConnectTimeout=15',
        '-o','ProxyCommand='+os.environ['PANCANCER_READONLY_PROXY_COMMAND'],os.environ['PANCANCER_READONLY_SSH_TARGET'],
        'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 '+shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -c '+shlex.quote(program)]
    started=now();phase='SOURCE_SYNTHETIC';proc=None;source_error=None;payload=None
    try:
        proc=subprocess.Popen(args,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        message=json.loads(proc.stdout.readline())
        if message['kind']=='SOURCE_TECHNICAL_FAILURE':source_error=message['error_type']
        assert message['kind']=='PURE_TCR_FIXTURE_ACCEPTED'
        opener=build_opener(ProxyHandler({'https':os.environ['PANCANCER_METADATA_HTTP_PROXY']}));phase='HTTP_HEAD'
        with opener.open(Request(item['url'],method='HEAD'),timeout=35) as response:
            assert response.status==200 and response.url==item['url']
            assert int(response.headers['Content-Length'])==item['bytes']
            assert response.headers.get('ETag')==item['etag'] and response.headers.get('Last-Modified')==item['last_modified']
        phase='HTTP_PUBLIC_PROCESSED_ANNOTATION';buffer=io.BytesIO();digest=hashlib.sha256()
        with opener.open(Request(item['url']),timeout=45) as response:
            assert response.status==200 and response.url==item['url']
            assert int(response.headers['Content-Length'])==item['bytes']
            assert response.headers.get('ETag')==item['etag'] and response.headers.get('Last-Modified')==item['last_modified']
            while True:
                value=response.read(1024**2)
                if not value:break
                buffer.write(value);digest.update(value);assert buffer.tell()<=item['bytes']
        assert buffer.tell()==item['bytes'];phase='SOURCE_RAM_TRANSPORT'
        proc.stdin.write((json.dumps(dict(bytes=item['bytes'],SHA256=digest.hexdigest()))+'\n').encode())
        buffer.seek(0)
        for value in iter(lambda:buffer.read(1024**2),b''):proc.stdin.write(value)
        proc.stdin.flush();proc.stdin.close();buffer.close();phase='SOURCE_MEASUREMENT'
        summary=None
        for line in proc.stdout:
            message=json.loads(line)
            if message['kind']=='SOURCE_TECHNICAL_FAILURE':source_error=message['error_type'];raise RuntimeError('SOURCE_MEASUREMENT_NOT_ACCEPTED')
            assert message['kind']=='PRIVATE_MEASUREMENT_HANDOFF' and message['transport_SHA256']==digest.hexdigest()
            summary=message['summary'];payload=message['payload']
        exitcode=proc.wait(timeout=30);assert exitcode==0 and summary is not None
        save(task/'private_complete_measurement_v1.json',summary)
        save(task/'private_measurement_terminal_v1.json',dict(status='COMPLETE_NONSEQUENCE_TCR_SOURCE_EXIT0_PRIVATE_RAM_HANDOFF',started_UTC=started,
            completed_UTC=now(),source_exit_code=exitcode,source_RAM_released=True,frozen_contract_SHA256=sha(frozen),
            summary_SHA256=sha(task/'private_complete_measurement_v1.json'),retrieved_processed_annotation_SHA256=digest.hexdigest(),
            no_individual_records_identifiers_sequence_or_matrix_files_persisted=True,source_header_size_signature_and_gzip_CRC_accepted=True))
        accepted=bool(summary['semantics']['operational_reported_expansion_qualification'])
        print(json.dumps(dict(status='COMPLETE_TCR_MEASUREMENT',operational_measurement_qualified=accepted,private_numeric_payload_RAM_only=True)),flush=True)
        if accepted:
            phase='FROZEN_STATISTICS_HANDOFF';stat=task/'private_frozen_statistics_v1.json';release=task/'private_release_RAM_v1.json';deadline=time.monotonic()+plan['RAM_handoff_seconds']
            while not stat.exists() and not release.exists() and time.monotonic()<deadline:time.sleep(1)
            if stat.exists():
                q=json.loads(stat.read_text());assert q['measurement_terminal_SHA256']==sha(task/'private_measurement_terminal_v1.json')
                for name,digest in q['code_SHA256'].items():assert sha(code/name)==digest
                module=importlib.import_module(q['statistics_module']);module.run(payload,task,q)
            elif release.exists():
                decision=json.loads(release.read_text());assert decision['measurement_terminal_SHA256']==sha(task/'private_measurement_terminal_v1.json')
                phase='RAM_RELEASE_AFTER_MEASUREMENT_REVIEW'
            else:phase='NO_STAT_CONTRACT_RAM_RELEASE'
        else:phase='MEASUREMENT_UNQUALIFIED_NO_STATISTICS'
        del payload
        save(task/'private_transport_terminal_v1.json',dict(status='COMPLETE_TCR_TRANSPORT_AND_PRIVATE_RAM_RELEASE',completed_UTC=now(),last_phase=phase,
            measurement_qualified=accepted,all_payload_RAM_released=True,source_exit_code=exitcode))
        print(json.dumps(dict(status='COMPLETE_PRIVATE_RAM_RELEASE',last_phase=phase)),flush=True)
    except Exception as error:
        if proc is not None and proc.poll() is None:
            proc.terminate()
            try:proc.wait(timeout=10)
            except subprocess.TimeoutExpired:proc.kill();proc.wait()
        save(task/'private_failed_terminal_v1.json',dict(status='TCR_SCOPE_NOT_ACCEPTED',observed_UTC=now(),phase=phase,error_type=type(error).__name__,
            source_error_type=source_error,HTTP_code=getattr(error,'code',None),source_exit_code=proc.returncode if proc is not None else None,
            no_raw_error_or_individual_partial_payload_saved=True))
        raise RuntimeError('TCR_SCOPE_NOT_ACCEPTED') from None


if __name__=='__main__':main()
