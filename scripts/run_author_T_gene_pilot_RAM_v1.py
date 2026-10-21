"""Same existing public pilot, full hash verification, T-row assay in RAM only."""
from pathlib import Path
from datetime import datetime,timezone
from urllib.request import Request,build_opener,ProxyHandler
import hashlib
import json
import os
import shlex
import subprocess
from GEO_public_processed_RAM_v1 import read_public


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def save(path,obj):
    assert not path.exists();path.write_text(json.dumps(obj,indent=2)+'\n')


def main():
    task=Path(os.environ['PANCANCER_AUTHOR_RNA_TASK']);code=Path(__file__).parent
    frozen=task/'private_frozen_specific_gene_pilot_v1.json';plan=json.loads(frozen.read_text());item=plan['source']
    for name,digest in plan['code_SHA256'].items():assert sha(code/name)==digest
    program="__name__='embedded_RNA_pilot'\n"
    program+=(code/'biological_annotation_label_gate_v2.py').read_text()+'\n'
    program+=(code/'qualify_author_spatial_metadata_RAM_v2.py').read_text()+'\n'
    program+=(code/'measure_author_T_gene_capture_RAM_v1.py').read_text()+'\n'
    program+='\nPILOT_PLAN='+repr(plan['measurement'])+'\nSOURCE_SIZE='+str(item['bytes'])+'\nSOURCE_SHA='+repr(item['SHA256'])+'\n'
    program+='''
import sys,h5py
try:
 synthetic()
 print(json.dumps(dict(kind='PURE_RNA_FIXTURE_ACCEPTED')),flush=True)
 buffer=io.BytesIO();digest=hashlib.sha256();remaining=SOURCE_SIZE
 while remaining:
  value=sys.stdin.buffer.read(min(1024**2,remaining))
  if not value:raise EOFError('INCOMPLETE_SOURCE')
  buffer.write(value);digest.update(value);remaining-=len(value)
 assert digest.hexdigest()==SOURCE_SHA and buffer.tell()==SOURCE_SIZE
 buffer.seek(0)
 with h5py.File(buffer,'r') as h:result=measure(h,PILOT_PLAN)
 buffer.close()
 print(json.dumps(dict(kind='COMPLETE_SOURCE_AGGREGATE',source_SHA256=digest.hexdigest(),result=result)),flush=True)
except Exception as error:
 print(json.dumps(dict(kind='SOURCE_TECHNICAL_FAILURE',error_type=type(error).__name__)),flush=True)
 raise SystemExit(2)
'''
    assert len(program.encode())<100000
    args=['ssh','-T','-o','ClearAllForwardings=yes','-o','BatchMode=yes','-o','ConnectTimeout=15',
        '-o','ProxyCommand='+os.environ['PANCANCER_READONLY_PROXY_COMMAND'],os.environ['PANCANCER_READONLY_SSH_TARGET'],
        'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 '+shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -c '+shlex.quote(program)]
    started=now();proc=None;phase='HTTP_HEAD';transferred=0;source_error=None
    try:
        http=build_opener(ProxyHandler({'https':os.environ['PANCANCER_METADATA_HTTP_PROXY']}))
        with http.open(Request(item['url'],method='HEAD'),timeout=35) as response:
            assert response.url==item['url'] and response.status==200
            assert int(response.headers['Content-Length'])==item['bytes']
            assert response.headers.get('ETag')==item['etag'] and response.headers.get('Last-Modified')==item['last_modified']
        phase='SOURCE_SYNTHETIC';proc=subprocess.Popen(args,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        first=json.loads(proc.stdout.readline());
        if first['kind']=='SOURCE_TECHNICAL_FAILURE':source_error=first['error_type']
        assert first['kind']=='PURE_RNA_FIXTURE_ACCEPTED'
        phase='PUBLIC_RAM_RANGES'
        def progress(i,n,received):
            print(json.dumps(dict(status='ORIGINAL_PILOT_RAM_RANGE_PROGRESS',ranges=i,total_ranges=n,bytes=received)),flush=True)
        buffer,digest=read_public(item,os.environ['PANCANCER_METADATA_HTTP_PROXY'],progress)
        assert digest==item['SHA256'] and len(buffer.getbuffer())==item['bytes']<=plan['maximum_public_bytes']
        phase='SOURCE_RAM_DELIVERY'
        for block in iter(lambda:buffer.read(1024**2),b''):proc.stdin.write(block);transferred+=len(block)
        proc.stdin.flush();proc.stdin.close();buffer.close();phase='SOURCE_MEASUREMENT'
        result=None
        for line in proc.stdout:
            message=json.loads(line)
            if message['kind']=='SOURCE_TECHNICAL_FAILURE':source_error=message['error_type'];raise RuntimeError('SOURCE_MEASUREMENT_REJECTED')
            assert message['kind']=='COMPLETE_SOURCE_AGGREGATE' and message['source_SHA256']==item['SHA256'];result=message['result']
        exitcode=proc.wait(timeout=30);assert exitcode==0 and result is not None and transferred==item['bytes']
        save(task/'private_specific_gene_pilot_summary_v1.json',result)
        save(task/'private_specific_gene_pilot_terminal_v1.json',dict(status=result['status'],started_UTC=started,completed_UTC=now(),
            source_exit_code=exitcode,source_RAM_released=True,source_bytes=transferred,both_ends_original_whole_file_SHA256=item['SHA256'],
            frozen_contract_SHA256=sha(frozen),summary_SHA256=sha(task/'private_specific_gene_pilot_summary_v1.json'),
            no_expression_files_individual_arrays_or_identifiers_saved=True,runtime_synthetic_accepted_before_human_bytes=True))
        print(json.dumps(dict(status='COMPLETE_SAME_PILOT_SOURCE_EXIT0_RAM_RELEASED')),flush=True)
    except Exception as error:
        if proc is not None and proc.poll() is None:
            proc.terminate()
            try:proc.wait(timeout=10)
            except subprocess.TimeoutExpired:proc.kill();proc.wait()
        save(task/'private_specific_gene_pilot_failed_v1.json',dict(status='PILOT_NOT_ACCEPTED',observed_UTC=now(),phase=phase,
            error_type=type(error).__name__,source_error_type=source_error,HTTP_code=getattr(error,'code',None),
            source_exit_code=proc.returncode if proc is not None else None,bytes_delivered=transferred,
            no_raw_stderr_or_partial_human_payload_saved=True))
        raise RuntimeError('PILOT_NOT_ACCEPTED') from None


if __name__=='__main__':main()
