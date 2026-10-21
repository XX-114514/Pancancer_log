"""Stream original public processed files to a read-only RAM source worker.

Bounded range buffers only on the client; no new matrix/cell files anywhere.
The scope, every byte hash, and exact calculation code must be frozen first.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import shlex
import subprocess
from urllib.request import Request, ProxyHandler, build_opener


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    task=Path(os.environ['PANCANCER_NATIVE_SPATIAL_TASK']);code=Path(__file__).parent
    freeze=task/'private_frozen_all11_native_description_v1.json';plan=json.loads(freeze.read_text())
    for name,digest in plan['code_SHA256'].items():assert sha(code/name)==digest
    out=task/'native_description_v1';out.mkdir(exist_ok=False)
    save=lambda name,value:(out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    program='import types\nevidence=types.ModuleType("fixed_evidence")\n'
    program+='exec('+repr((code/'observed_T_marker_evidence_v1.py').read_text())+',evidence.__dict__)\n'
    program+=(code/'CM02_LUAD_native_state_spatial_description_v1.py').read_text()+'\n'
    program+='try:\n remote_main('+repr(dict(sources=plan['sources']))+',evidence)\n'
    program+='except Exception as error:\n print(json.dumps(dict(kind="failed",error_type=type(error).__name__)),flush=True)\n raise SystemExit(1)\n'
    args=['ssh','-T','-o','ClearAllForwardings=yes','-o','BatchMode=yes','-o','ConnectTimeout=15',
          '-o','ServerAliveInterval=30','-o','ServerAliveCountMax=3']
    proxy=os.environ.get('PANCANCER_READONLY_PROXY_COMMAND')
    if proxy:args+=['-o','ProxyCommand='+proxy]
    args += [os.environ['PANCANCER_READONLY_SSH_TARGET'],
             'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 '+shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -c '+shlex.quote(program)]
    proc=subprocess.Popen(args,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    started=datetime.now(timezone.utc).isoformat()
    save('private_checkpoint.json',dict(stage='SOURCE_PIPE_STARTED',started_UTC=started))
    try:
        ready=json.loads(proc.stdout.readline())
        assert ready['kind']=='ready','SOURCE_READY_GATE_NOT_ACCEPTED'
        save('private_native_runtime_synthetic_acceptance.json',ready)
        for ordinal,chosen in enumerate(plan['sources']):
            proc.stdin.write((json.dumps(chosen)+'\n').encode());proc.stdin.flush()
            size=chosen['bytes'];bounds=[(i,min(size-1,i+4*1024**2-1)) for i in range(0,size,4*1024**2)]
            def fetch(bound):
                lo,hi=bound
                opener=build_opener(ProxyHandler({'https':os.environ['PANCANCER_PUBLIC_HTTP_PROXY']}))
                with opener.open(Request(chosen['url'],headers={'Range':f'bytes={lo}-{hi}'}),timeout=40) as response:
                    assert response.status==206 and response.headers.get('Content-Range')==f'bytes {lo}-{hi}/{size}'
                    assert response.headers.get('Last-Modified')==chosen['last_modified']
                    if chosen.get('etag'):assert response.headers.get('ETag')==chosen['etag']
                    value=response.read(hi-lo+2);assert len(value)==hi-lo+1
                    return value
            digest=hashlib.sha256();sent=0
            with ThreadPoolExecutor(max_workers=4) as executor:
                pending={i:executor.submit(fetch,bounds[i]) for i in range(min(4,len(bounds)))}
                for i in range(len(bounds)):
                    block=pending.pop(i).result();digest.update(block);proc.stdin.write(block);sent+=len(block)
                    if i+4<len(bounds):pending[i+4]=executor.submit(fetch,bounds[i+4])
                    if i%16==0:save('private_checkpoint.json',dict(stage='BOUNDED_PUBLIC_BYTES_TO_RAM_SOURCE',image_ordinal=ordinal+1,
                        expected_images=len(plan['sources']),bytes_sent=sent,source_bytes=size,observed_UTC=datetime.now(timezone.utc).isoformat()))
            proc.stdin.flush();assert sent==size and digest.hexdigest()==chosen['SHA256']
        proc.stdin.close();aggregate=None
        for line in proc.stdout:
            record=json.loads(line)
            if record['kind']=='source_aggregate':aggregate=record['value']
            else:save('private_checkpoint.json',record)
        rc=proc.wait();errors=proc.stderr.read()  # Never persist/print possibly identifying stderr.
        assert rc==0 and aggregate is not None,'SOURCE_MEASUREMENT_NOT_COMPLETE'
        save('private_complete_source_descriptive_summary.json',aggregate)
        save('private_terminal.json',dict(status='COMPLETE_FROZEN_ALL11_DESCRIPTIVE_SCOPE',started_UTC=started,
            completed_UTC=datetime.now(timezone.utc).isoformat(),freeze_SHA256=sha(freeze),source_exit_code=rc,
            source_aggregate_SHA256=sha(out/'private_complete_source_descriptive_summary.json'),
            all_individual_source_and_graph_arrays_RAM_only=True,new_P_CI_or_clinical_model=False))
        print(json.dumps(dict(status='COMPLETE_ALL11_NATIVE_STATE_DESCRIPTION_SOURCE_SUMMARY_ONLY')),flush=True)
    except Exception as error:
        if proc.stdin and not proc.stdin.closed:proc.stdin.close()
        # Closing stdin lets a finite reader reject a truncated file; no broad process kill.
        save('private_failed_terminal.json',dict(status='SOURCE_DESCRIPTION_NOT_ACCEPTED',error_type=type(error).__name__,
            observed_UTC=datetime.now(timezone.utc).isoformat(),no_partial_scientific_effects_accepted=True))
        raise RuntimeError('SOURCE_DESCRIPTION_NOT_ACCEPTED') from None


if __name__=='__main__':
    main()
