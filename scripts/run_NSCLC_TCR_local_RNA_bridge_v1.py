"""Local RAM only: two official public metadata objects, exact key bridge.

No SSH, remote execution, record upload, sequence-field parsing or raw files.
"""
from pathlib import Path
from datetime import datetime,timezone
from urllib.request import Request,ProxyHandler,build_opener
from urllib.parse import urlparse
import hashlib
import importlib
import io
import json
import os
import resource
import time
import pandas as pd
import qualify_NSCLC_reported_TCR_expansion_RAM_v1 as original
import NSCLC_TCR_local_RNA_bridge_v1 as bridge


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def save(path,obj):
    assert not path.exists();path.write_text(json.dumps(obj,indent=2)+'\n')


def fetch(item,opener):
    p=urlparse(item['url']);assert p.scheme=='https' and p.hostname=='ftp.ncbi.nlm.nih.gov' and p.path.endswith('.csv.gz')
    with opener.open(Request(item['url'],method='HEAD'),timeout=35) as response:
        assert response.status==200 and response.url==item['url']
        size=int(response.headers['Content-Length']);assert 0<size<=item['maximum_bytes']
        etag=response.headers.get('ETag');modified=response.headers.get('Last-Modified')
    b=io.BytesIO();digest=hashlib.sha256()
    with opener.open(Request(item['url']),timeout=45) as response:
        assert response.status==200 and response.url==item['url'] and int(response.headers['Content-Length'])==size
        assert response.headers.get('ETag')==etag and response.headers.get('Last-Modified')==modified
        while True:
            value=response.read(1024**2)
            if not value:break
            b.write(value);digest.update(value);assert b.tell()<=size
    assert b.tell()==size and digest.hexdigest()==item['expected_SHA256'];b.seek(0)
    return b,dict(bytes=size,SHA256=digest.hexdigest(),current_HEAD_and_GET_headers_identical=True)


def main():
    task=Path(os.environ['PANCANCER_LOCAL_TCR_TASK']);code=Path(__file__).parent;frozen=task/'private_frozen_local_RNA_bridge_v1.json';plan=json.loads(frozen.read_text())
    for name,digest in plan['code_SHA256'].items():assert sha(code/name)==digest
    bridge.synthetic();started=now();phase='OFFICIAL_PUBLIC_TCR_LOCAL_RAM'
    try:
        opener=build_opener(ProxyHandler({'https':os.environ['PANCANCER_METADATA_HTTP_PROXY']}))
        tcr,receipt1=fetch(plan['TCR_source'],opener)
        frame=original.parse_selected(tcr.getvalue(),plan['TCR_parse']);tcr.close();del tcr
        phase='OFFICIAL_PUBLIC_RNA_METADATA_LOCAL_RAM';rna,receipt2=fetch(plan['RNA_source'],opener)
        md=pd.read_csv(rna,compression='gzip',usecols=bridge.RNA_COLUMNS,dtype=str,keep_default_na=False,nrows=plan['expected_RNA_rows']+1)
        rna.close();del rna
        phase='EXACT_LOCAL_KEY_AND_REPORTED_STATE_MEASUREMENT';payload,summary=bridge.measure(frame,md,plan)
        del frame,md
        assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<plan['local_RSS_cap_bytes']
        save(task/'private_complete_local_RNA_bridge_v1.json',summary)
        save(task/'private_local_RNA_bridge_terminal_v1.json',dict(status=summary['status'],started_UTC=started,completed_UTC=now(),
            frozen_contract_SHA256=sha(frozen),summary_SHA256=sha(task/'private_complete_local_RNA_bridge_v1.json'),
            TCR_object=receipt1,RNA_object=receipt2,local_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            original_RNA_source_SHA_exactly_reused=True,no_individual_files_or_remote_transmission=True,
            actual_local_synthetic_and_private_payload_RAM_only=True))
        accepted=bool(summary['semantic_profile']['operational_reported_expansion_qualified'])
        print(json.dumps(dict(status='COMPLETE_LOCAL_TCR_RNA_MEASUREMENT_PRIVATE_RAM_HELD',measurement_semantics_qualified=accepted)),flush=True)
        if accepted:
            phase='PRIVATE_LOCAL_STATISTICS_HANDOFF';stat=task/'private_frozen_local_statistics_v1.json';release=task/'private_release_local_RAM_v1.json'
            deadline=time.monotonic()+plan['RAM_handoff_seconds']
            while not stat.exists() and not release.exists() and time.monotonic()<deadline:time.sleep(1)
            if stat.exists():
                q=json.loads(stat.read_text());assert q['measurement_terminal_SHA256']==sha(task/'private_local_RNA_bridge_terminal_v1.json')
                for name,digest in q['code_SHA256'].items():assert sha(code/name)==digest
                module=importlib.import_module(q['statistics_module']);module.run(payload,task,q)
            elif release.exists():
                decision=json.loads(release.read_text());assert decision['measurement_terminal_SHA256']==sha(task/'private_local_RNA_bridge_terminal_v1.json')
                phase='LOCAL_RAM_RELEASE_AFTER_MEASUREMENT_REVIEW'
            else:phase='LOCAL_RAM_RELEASE_NO_STATISTICS_CONTRACT'
        else:phase='LOCAL_RAM_RELEASE_SEMANTICS_UNQUALIFIED'
        del payload
        save(task/'private_local_RNA_transport_terminal_v1.json',dict(status='COMPLETE_LOCAL_TCR_RNA_RAM_RELEASED',completed_UTC=now(),last_phase=phase,
            no_remote_record_or_sequence_transmission=True,local_payload_released=True))
        print(json.dumps(dict(status='COMPLETE_LOCAL_RAM_RELEASED',last_phase=phase)),flush=True)
    except Exception as error:
        save(task/'private_local_RNA_bridge_failed_v1.json',dict(status='LOCAL_RNA_BRIDGE_NOT_ACCEPTED',observed_UTC=now(),phase=phase,
            error_type=type(error).__name__,HTTP_code=getattr(error,'code',None),no_raw_error_or_record_payload_saved=True))
        raise RuntimeError('LOCAL_RNA_BRIDGE_NOT_ACCEPTED') from None


if __name__=='__main__':main()
