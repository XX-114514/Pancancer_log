"""Fixed original public ranges -> existing h5py runtime; optional RAM handoff.

No raw/individual files, identifiers, payload logging, HTTP fallback or retries.
"""
from pathlib import Path
from datetime import datetime, timezone
from urllib.request import Request, ProxyHandler, build_opener
import hashlib
import json
import os
import shlex
import subprocess
import time


def now(): return datetime.now(timezone.utc).isoformat()
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path, obj):
    assert not path.exists()
    path.write_text(json.dumps(obj, indent=2)+'\n')


def main():
    task=Path(os.environ['PANCANCER_AUTHOR_SPATIAL_TASK']); code=Path(__file__).parent
    freeze=task/'private_frozen_measurement_v1.json'; plan=json.loads(freeze.read_text())
    for name, digest in plan['code_SHA256'].items(): assert sha(code/name)==digest
    program="__name__='embedded_metadata_scope'\n"
    program+=(code/'H5AD_metadata_range_reader_v1.py').read_text()+'\nrange_synthetic=synthetic\n'
    program+=(code/'biological_annotation_label_gate_v2.py').read_text()+'\nlabel_synthetic=synthetic\n'
    program+=(code/'qualify_author_spatial_metadata_RAM_v2.py').read_text()+'\n'
    program+='''
import sys, h5py, resource
try:
 range_synthetic(); label_synthetic(); synthetic()
except Exception as error:
 print(json.dumps(dict(kind='SOURCE_TECHNICAL_FAILURE',error_type=type(error).__name__)),flush=True)
 raise SystemExit(2)
print(json.dumps(dict(kind='RUNTIME_SYNTHETIC_ACCEPTED')),flush=True)
payloads=[]; results=[]; totals=0
try:
 for index, size in enumerate(SIZES):
  print(json.dumps(dict(kind='SET_FILE',index=index)),flush=True)
  if sys.stdin.buffer.readline()!=b'READY\\n':raise ValueError('CLIENT_HEADER_NOT_ACCEPTED')
  def fetch(start,length):
   print(json.dumps(dict(kind='RANGE_REQUEST',index=index,start=start,length=length)),flush=True)
   parts=[]; remaining=length
   while remaining:
    part=sys.stdin.buffer.read(remaining)
    if not part: raise EOFError('CHANNEL_CLOSED')
    parts.append(part); remaining-=len(part)
   return b''.join(parts)
  reader=RangeReader(size,fetch,budget=PER_FILE_CAP)
  with h5py.File(reader,'r') as h: payload,result=qualify(h)
  assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss < RSS_CAP_KIB
  payloads.append(payload); results.append(result); totals+=reader.transferred
  del reader,payload,result
  print(json.dumps(dict(kind='FILE_COMPLETE',index=index,metadata_bytes_total=totals)),flush=True)
 summary=aggregate(results)
 assert summary['rows']==EXPECTED_ROWS
 summary['peak_RSS_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 print(json.dumps(dict(kind='COMPLETE_MEASUREMENT',summary=summary,metadata_bytes_total=totals)),flush=True)
 command=json.loads(sys.stdin.buffer.readline())
 if command['kind']=='RELEASE_RAM':
  print(json.dumps(dict(kind='RAM_RELEASED_NO_EFFECTS')),flush=True)
 elif command['kind']=='RUN_FROZEN_DESCRIPTION':
  scope={}; exec(command['code'],scope)
  result=scope['run_description'](payloads,command['contract'])
  print(json.dumps(dict(kind='COMPLETE_DESCRIPTION',result=result)),flush=True)
 else: raise ValueError('UNKNOWN_HANDOFF')
except Exception as error:
 print(json.dumps(dict(kind='SOURCE_TECHNICAL_FAILURE',error_type=type(error).__name__)),flush=True)
 raise SystemExit(2)
'''
    constants='SIZES='+repr([v['bytes'] for v in plan['sources']])+'\n'
    constants+='PER_FILE_CAP='+str(plan['per_file_network_cap'])+'\nRSS_CAP_KIB='+str(plan['RSS_cap_KiB'])+'\nEXPECTED_ROWS='+str(plan['expected_total_rows'])+'\n'
    program=program.replace('payloads=[]; results=[]; totals=0',constants+'payloads=[]; results=[]; totals=0')
    args=['ssh','-T','-o','ClearAllForwardings=yes','-o','BatchMode=yes','-o','ConnectTimeout=15',
        '-o','ProxyCommand='+os.environ['PANCANCER_READONLY_PROXY_COMMAND'],os.environ['PANCANCER_READONLY_SSH_TARGET'],
        'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 '+shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -c '+shlex.quote(program)]
    http=build_opener(ProxyHandler({'https':os.environ['PANCANCER_METADATA_HTTP_PROXY']}))
    started=now(); proc=subprocess.Popen(args,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    transferred=0; requests=0; synthetic_ok=False; active=None; phase='SOURCE_SYNTHETIC'; source_error=None; completed_files=0
    frozen_description=None; description_result=None; measurement_ok=False; released=False
    try:
        for line in proc.stdout:
            msg=json.loads(line); kind=msg['kind']
            if kind=='SOURCE_TECHNICAL_FAILURE':source_error=msg['error_type'];raise RuntimeError('SOURCE_REJECTED')
            if kind=='RUNTIME_SYNTHETIC_ACCEPTED':synthetic_ok=True;continue
            assert synthetic_ok
            if kind=='SET_FILE':
                assert msg['index']==completed_files; active=msg['index']; item=plan['sources'][active]; phase='HTTP_HEAD'
                with http.open(Request(item['url'],method='HEAD'),timeout=35) as response:
                    assert response.url==item['url'] and response.status==200
                    assert int(response.headers['Content-Length'])==item['bytes']
                    assert response.headers.get('ETag')==item['etag'] and response.headers.get('Last-Modified')==item['last_modified']
                proc.stdin.write(b'READY\n');proc.stdin.flush();continue
            if kind=='RANGE_REQUEST':
                assert msg['index']==active; item=plan['sources'][active];lo=msg['start'];length=msg['length'];hi=lo+length-1
                assert 0<=lo<=hi<item['bytes'] and 0<length<=65536
                assert transferred+length<=plan['total_network_cap'] and requests<plan['max_requests']
                phase='HTTP_RANGE';
                with http.open(Request(item['url'],headers={'Range':f'bytes={lo}-{hi}'}),timeout=35) as response:
                    assert response.url==item['url'] and response.status==206
                    assert response.headers.get('Content-Range')==f"bytes {lo}-{hi}/{item['bytes']}"
                    assert response.headers.get('ETag')==item['etag'] and response.headers.get('Last-Modified')==item['last_modified']
                    value=response.read(length+1);assert len(value)==length
                phase='DELIVER_RANGE';proc.stdin.write(value);proc.stdin.flush();transferred+=length;requests+=1;continue
            if kind=='FILE_COMPLETE':
                assert msg['index']==completed_files and msg['metadata_bytes_total']==transferred; completed_files+=1
                print(json.dumps(dict(status='METADATA_FILE_COMPLETE',completed_files=completed_files,network_bytes=transferred)),flush=True);continue
            if kind=='COMPLETE_MEASUREMENT':
                assert completed_files==len(plan['sources']) and msg['metadata_bytes_total']==transferred
                phase='MEASUREMENT_HANDOFF';summary=msg['summary'];assert summary['rows']==plan['expected_total_rows']
                save(task/'private_complete_measurement_v1.json',summary)
                save(task/'private_measurement_terminal_v1.json',dict(status='COMPLETE_MEASUREMENT_RAM_HELD_NO_EFFECTS',
                    started_UTC=started,completed_UTC=now(),frozen_contract_SHA256=sha(freeze),network_bytes=transferred,requests=requests,
                    summary_SHA256=sha(task/'private_complete_measurement_v1.json'),source_has_not_exited=True))
                measurement_ok=True;print(json.dumps(dict(status='COMPLETE_METADATA_MEASUREMENT_RAM_HELD')),flush=True)
                deadline=time.monotonic()+plan['RAM_handoff_seconds']; path=task/'private_frozen_description_v1.json'
                while not path.exists() and time.monotonic()<deadline and proc.poll() is None:time.sleep(1)
                if path.exists():
                    frozen_description=json.loads(path.read_text());script=code/frozen_description['code_filename']
                    assert script.name==frozen_description['code_filename'] and script.parent==code and sha(script)==frozen_description['code_SHA256']
                    assert frozen_description['measurement_terminal_SHA256']==sha(task/'private_measurement_terminal_v1.json')
                    command=dict(kind='RUN_FROZEN_DESCRIPTION',code=script.read_text(),contract=frozen_description)
                    phase='FROZEN_DESCRIPTION'
                else:command=dict(kind='RELEASE_RAM');phase='RAM_RELEASE'
                proc.stdin.write((json.dumps(command)+'\n').encode());proc.stdin.flush();continue
            if kind=='COMPLETE_DESCRIPTION':
                assert frozen_description is not None;description_result=msg['result'];continue
            if kind=='RAM_RELEASED_NO_EFFECTS':released=True;continue
            raise ValueError('UNEXPECTED_PROTOCOL_MESSAGE')
        proc.stdin.close();exitcode=proc.wait(timeout=30)
        assert exitcode==0 and measurement_ok and (released or description_result is not None)
        if description_result is not None:save(task/'private_complete_description_v1.json',description_result)
        save(task/'private_terminal_v1.json',dict(status='COMPLETE_DESCRIPTION_RAM_RELEASED' if description_result is not None else 'COMPLETE_METADATA_RAM_RELEASED_NO_EFFECTS',
            started_UTC=started,completed_UTC=now(),source_exit_code=exitcode,frozen_contract_SHA256=sha(freeze),
            network_bytes=transferred,requests=requests,all_files=completed_files,source_RAM_released=True,
            description_contract_SHA256=sha(task/'private_frozen_description_v1.json') if frozen_description else None,
            no_raw_or_individual_files=True,whole_file_SHA_not_recomputed=True))
        print(json.dumps(dict(status='COMPLETE_SOURCE_EXIT0_RAM_RELEASED')),flush=True)
    except Exception as error:
        if proc.poll() is None:proc.terminate()
        try:proc.wait(timeout=10)
        except subprocess.TimeoutExpired:proc.kill();proc.wait()
        save(task/'private_failed_terminal_v1.json',dict(status='SCOPE_NOT_ACCEPTED',observed_UTC=now(),phase=phase,
            error_type=type(error).__name__,source_error_type=source_error,HTTP_code=getattr(error,'code',None),
            source_exit_code=proc.returncode,network_bytes=transferred,requests=requests,completed_files=completed_files,
            measurement_accepted_before_failure=measurement_ok,no_raw_stderr_or_partial_payload_saved=True))
        raise RuntimeError('SCOPE_NOT_ACCEPTED') from None


if __name__=='__main__':main()
