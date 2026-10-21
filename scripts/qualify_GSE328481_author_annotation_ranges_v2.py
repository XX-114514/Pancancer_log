"""Metadata-only public range bridge to the existing remote h5py runtime.

Native H5AD fragments and category/cell arrays stay in RAM on both hosts.
Public URL/header boundaries are fixed, with no automatic retries or fallback.
"""
from pathlib import Path
from datetime import datetime,timezone
from urllib.request import Request,ProxyHandler,build_opener
import hashlib
import json
import os
import shlex
import subprocess


def main():
    task=Path(os.environ['PANCANCER_AUTHOR_SPATIAL_TASK']);code=Path(__file__).parent
    freeze=task/'private_frozen_pilot_annotation_ranges_v2.json';plan=json.loads(freeze.read_text());item=plan['source']
    for name,digest in plan['code_SHA256'].items():assert hashlib.sha256((code/name).read_bytes()).hexdigest()==digest
    http=build_opener(ProxyHandler({'https':os.environ['PANCANCER_METADATA_HTTP_PROXY']}))
    with http.open(Request(item['url'],method='HEAD'),timeout=30) as response:
        assert response.url==item['url'] and response.status==200
        assert int(response.headers['Content-Length'])==item['bytes']
        assert response.headers.get('ETag')==item['etag'] and response.headers.get('Last-Modified')==item['last_modified']
    program="__name__='embedded_metadata_scope'\n"+(code/'H5AD_metadata_range_reader_v1.py').read_text()+'\nrange_synthetic=synthetic\n'
    program+=(code/'qualify_GSE328481_existing_author_annotation_v1.py').read_text()+'\n'
    program+='''
import sys,h5py
range_synthetic();synthetic()
# Exercise h5py itself through the seek/read adapter on a pure RAM fixture.
b=io.BytesIO()
with h5py.File(b,'w') as h:
    o=h.create_group('obs');o.create_dataset('_index',data=np.array([b'one',b'two']))
    o.create_dataset('anno',data=np.array([b'Treg',b'unknown_private']))
raw=b.getvalue();fake=RangeReader(len(raw),lambda start,length:raw[start:start+length])
with h5py.File(fake,'r') as h:
    result=summarize(h);assert result['fields']['anno']['recognized_counts']=={'Treg':1}
assert 'unknown_private' not in json.dumps(result)
del raw,b,fake,result
print(json.dumps(dict(kind='RUNTIME_SYNTHETIC_ACCEPTED')),flush=True)
def fetch(start,length):
    print(json.dumps(dict(kind='RANGE_REQUEST',start=start,length=length)),flush=True)
    parts=[];remaining=length
    while remaining:
        chunk=sys.stdin.buffer.read(remaining)
        if not chunk:raise EOFError('RANGE_CHANNEL_CLOSED')
        parts.append(chunk);remaining-=len(chunk)
    return b''.join(parts)
'''
    program+='reader=RangeReader('+str(item['bytes'])+',fetch,budget='+str(plan['network_byte_cap'])+')\n'
    program+='with h5py.File(reader,"r") as h:\n result=summarize(h)\n'
    program+='print(json.dumps(dict(kind="COMPLETE_SOURCE_AGGREGATE",result=result,metadata_bytes=reader.transferred)),flush=True)\n'
    args=['ssh','-T','-o','ClearAllForwardings=yes','-o','BatchMode=yes','-o','ConnectTimeout=15','-o','ProxyCommand='+os.environ['PANCANCER_READONLY_PROXY_COMMAND'],os.environ['PANCANCER_READONLY_SSH_TARGET'],
          'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 '+shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -c '+shlex.quote(program)]
    started=datetime.now(timezone.utc).isoformat();proc=subprocess.Popen(args,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    transferred=0;requests=0;accepted=False;aggregate=None;synthetic_ok=False
    try:
        for line in proc.stdout:
            message=json.loads(line)
            if message['kind']=='RUNTIME_SYNTHETIC_ACCEPTED':synthetic_ok=True;continue
            if message['kind']=='COMPLETE_SOURCE_AGGREGATE':
                assert synthetic_ok and message['metadata_bytes']==transferred;aggregate=message['result'];accepted=True;break
            assert synthetic_ok and message['kind']=='RANGE_REQUEST'
            lo=message['start'];length=message['length'];hi=lo+length-1
            assert 0<=lo<=hi<item['bytes'] and 0<length<=65536 and requests<512
            assert transferred+length<=plan['network_byte_cap']
            with http.open(Request(item['url'],headers={'Range':f'bytes={lo}-{hi}'}),timeout=35) as response:
                assert response.url==item['url'] and response.status==206
                assert response.headers.get('Content-Range')==f"bytes {lo}-{hi}/{item['bytes']}"
                assert response.headers.get('ETag')==item['etag'] and response.headers.get('Last-Modified')==item['last_modified']
                value=response.read(length+1);assert len(value)==length
            proc.stdin.write(value);proc.stdin.flush();transferred+=length;requests+=1
        proc.stdin.close();exitcode=proc.wait(timeout=30)
        assert exitcode==0 and accepted
    except Exception as error:
        proc.terminate()
        try:proc.wait(timeout=10)
        except subprocess.TimeoutExpired:proc.kill();proc.wait()
        (task/'private_failed_terminal_v2.json').write_text(json.dumps(dict(status='METADATA_RANGE_SCOPE_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),error_type=type(error).__name__,network_bytes=transferred,requests=requests,no_raw_stderr_or_partial_human_payload_saved=True),indent=2)+'\n')
        raise RuntimeError('METADATA_RANGE_SCOPE_NOT_ACCEPTED') from None
    (task/'private_source_pilot_annotation_summary_v2.json').write_text(json.dumps(aggregate,indent=2)+'\n')
    terminal=dict(status='COMPLETE_METADATA_RANGES_ONE_PILOT_NO_EFFECTS',started_UTC=started,completed_UTC=datetime.now(timezone.utc).isoformat(),
        frozen_contract_SHA256=hashlib.sha256(freeze.read_bytes()).hexdigest(),source_size=item['bytes'],metadata_network_bytes=transferred,range_requests=requests,
        actual_source_exit_code=exitcode,runtime_synthetic_h5py_adapter_and_nonemission_accepted=synthetic_ok,
        full_file_SHA_not_rechecked_in_this_metadata_only_transport=True,no_new_human_files_or_source_row_or_matrix_output=True)
    (task/'private_terminal_v2.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],metadata_network_bytes=transferred,range_requests=requests)),flush=True)


if __name__=='__main__':main()
