"""Read only official small ROI metadata; records stay in server RAM.

No cell intensities, learned matrices or outcomes association are selected.
Clinical categories and units are reviewed as source provenance, not fitted.
"""
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import csv
import hashlib
import importlib.util
import io
import json
import os
import re
import traceback
import zipfile
import time
import socket
import ssl
import http.client
import urllib.error
import urllib.request

import pandas as pd

MEMBER='data/input/meta/ROI_meta_afterQC.csv'


def sha(b):return hashlib.sha256(b).hexdigest()


def clean(s):return re.sub('[^a-z0-9]','',s.lower())



class BoundedRetryRangeFile:
    """Same official validated ranges; at most three transient transport tries.

    Budget counts every requested byte before a request, including failed tries.
    TLS certificate and HTTP status failures are never retried or bypassed.
    """
    def __init__(self,source,url,size):
        self.source=source;self.url=url;self.size=size;self.position=0
        self.transferred=0;self.requested=0;self.requests=0;self.identity=None;self.failures=[]
        assert url.startswith(source.BASE)
    def tell(self):return self.position
    def seek(self,offset,whence=0):
        self.position=offset if whence==0 else self.position+offset if whence==1 else self.size+offset
        assert 0<=self.position<=self.size
        return self.position
    def seekable(self):return True
    def read(self,n=-1):
        if n<0:n=self.size-self.position
        n=min(n,self.size-self.position)
        if n==0:return b''
        assert n<=16*1024**2
        start=self.position;end=start+n-1
        for attempt in range(3):
            assert self.requested+n<=64*1024**2
            self.requested+=n;self.requests+=1
            try:
                request=urllib.request.Request(self.url,headers={'Range':f'bytes={start}-{end}'})
                with self.source.opener().open(request,timeout=35) as response:
                    assert response.status==206
                    assert response.headers['Content-Range']==f'bytes {start}-{end}/{self.size}'
                    identity=(response.headers.get('ETag'),response.headers.get('Last-Modified'))
                    if self.identity is None:self.identity=identity
                    else:assert self.identity==identity
                    raw=response.read(n+1);assert len(raw)==n
                self.position+=n;self.transferred+=n
                return raw
            except (TimeoutError,socket.timeout,ConnectionResetError,http.client.IncompleteRead) as error:
                self.failures.append(type(error).__name__)
                if attempt==2:raise
                time.sleep(attempt+1)
        raise AssertionError('unreachable')


def main():
    helper=Path(os.environ['PANCANCER_IMC_TRANSPORT_HELPER'])
    assert sha(helper.read_bytes())=='9fb200d98e321584fde8bca6cb92f72151ac33ec00d64e04d93e7e2b33e870e2'
    spec=importlib.util.spec_from_file_location('published_archive_transport',helper)
    source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert not out.exists();out.mkdir(mode=0o700)
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2)+'\n')
    save('plan_before_ROI_registry_rows.json',dict(frozen_UTC=datetime.now(timezone.utc).isoformat(),
        source_SHA256=sha(Path(__file__).read_bytes()),selected_member=MEMBER,
        scope='One complete official metadata member, max64KiB; field names and source aggregates only.',
        clinical_definition='Author RECIST1.1 CR/PR/SD=responder, PD=nonresponder; no alternative recoding.',
        source_unit='Do not infer identities or outcomes from coded IDs. Candidate fields only; no cells or effects until unit semantics and exact joining accepted.',
        no_cell_or_intensity_rows_or_learnt_weights_or_TCR_selected=True,no_clinical_association=True,
        no_individual_records_written=True,no_author_code_executed=True,
        transport_retries='At most3 per range; transient timeout/reset/incomplete only; count all requested bytes; certificate/HTTP failures not bypassed.'))
    raw=source.fetch(source.BASE+'api/records/6838169',2*1024**2);d=json.loads(raw)
    assert d['id']==6838169 and d['metadata']['access_right']=='open'
    f=next(x for x in d['files'] if x['key']=='Melanoma_code_data.zip')
    assert f['size']==668227028 and f['checksum']=='md5:da74b7228763b76dcd2d9830c1a10cd9'
    rf=BoundedRetryRangeFile(source,f['links']['self'],f['size'])
    with zipfile.ZipFile(rf) as archive:
        info=archive.getinfo(MEMBER);assert info.file_size==4392 and info.file_size<=65536
        blob=archive.read(info);assert len(blob)==info.file_size
        rows=list(csv.DictReader(io.StringIO(blob.decode('utf-8-sig'))))
    columns=list(rows[0]);assert len(columns)==len(set(columns))
    assert all(set(row)==set(columns) for row in rows)
    frame=pd.read_csv(io.BytesIO(blob),dtype=str,keep_default_na=False)
    assert len(frame)==len(rows)
    for col in columns:
        other=frame.iloc[:,columns.index(col)].tolist()
        assert other==[row[col] for row in rows]
    fields={clean(c):c for c in columns if c}
    role_candidates={role:[fields[name] for name in names if name in fields] for role,names in {
        'case':['patient','patientid','patientname','sample','samplename'],
        'ROI':['roi','roiid','sampleid','image','imageid'],
        'region':['region','location','position','site'],
        'response':['response','responder','group','efficacy','recist']}.items()}
    profiles=[]
    for role,candidates in role_candidates.items():
        for col in candidates:
            values=[row[col] for row in rows];counter=Counter(values)
            safe_categories=None
            allowed={'R','NR','CR','PR','SD','PD','CT','IM','core','invasive','Core','Invasive'}
            if role in ['region','response'] and set(counter)<=allowed:safe_categories=dict(sorted(counter.items()))
            profiles.append(dict(role=role,field=col,rows=len(values),unique_values=len(counter),
                empty_values=counter.get('',0),source_category_counts=safe_categories,
                unknown_or_identifier_values_not_exported=True))
    relationships=[]
    for case in role_candidates['case']:
        n=len({r[case] for r in rows});sizes=Counter(r[case] for r in rows)
        sizes_summary=[min(sizes.values()),max(sizes.values())] if n>=5 else None
        for response in role_candidates['response']:
            groups=defaultdict(set)
            for row in rows:groups[row[case]].add(row[response])
            pandas=frame.groupby(case,sort=False,dropna=False)[response].nunique()
            assert {k:len(v) for k,v in groups.items()}==pandas.to_dict()
            consistent=all(len(v)==1 for v in groups.values())
            counts=Counter(next(iter(v)) for v in groups.values()) if consistent else Counter()
            allowed={'R','NR','CR','PR','SD','PD'}
            relationships.append(dict(candidate_case_field=case,response_field=response,
                source_units=n,all_units_single_consistent_response=consistent,
                source_units_by_response=dict(sorted(counts.items())) if set(counts)<=allowed else None,
                ROIs_per_unit_min_max=sizes_summary,
                authors_published26_units14R12NR_exact=consistent and n==26 and counts=={'R':14,'NR':12},
                no_patient_unit_acceptance_from_counts_alone=True))
    final=json.loads(source.fetch(source.BASE+'api/records/6838169',2*1024**2))
    other=next(x for x in final['files'] if x['key']==f['key'])
    assert all(f.get(k)==other.get(k) for k in ['id','size','checksum','links'])
    result=dict(status='COMPLETE_OFFICIAL_ROI_METADATA_SOURCE_COUNTS_AND_FIELD_CONSISTENCY',
        completed_UTC=datetime.now(timezone.utc).isoformat(),producer_SHA256=sha(Path(__file__).read_bytes()),
        original_complete_member_SHA256=sha(blob),member_CRC32=info.CRC,complete_member_CRC_checked=True,
        schema_columns=columns,source_rows=len(rows),field_role_candidates=role_candidates,
        field_profiles=profiles,candidate_case_response_consistency=relationships,
        entire_file_csv_and_Pandas_string_value_replay=True,published_archive_metadata_unchanged=True,
        exact_cell_ROI_join_and_coordinate_and_positivity_not_yet_accepted=True,
        no_individual_keys_outcomes_or_records_written=True,no_new_clinical_or_spatial_effects=True,
        source_API_SHA256=sha(raw),range_bytes=rf.transferred,
        requested_bytes_including_failed_tries=rf.requested,range_requests_including_failed=rf.requests,
        transient_transport_failure_classes=rf.failures)
    save('private_complete_ROI_metadata_qualification.json',result)
    save('private_terminal.json',dict(status=result['status'],completed_UTC=datetime.now(timezone.utc).isoformat(),
        producer_SHA256=sha(Path(__file__).read_bytes()),files={p.name:dict(SHA256=sha(p.read_bytes()),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status=result['status'],terminal_SHA256=sha((out/'private_terminal.json').read_bytes()))))


if __name__=='__main__':
    try:main()
    except Exception as error:
        p=Path(os.environ['PANCANCER_RAM_OUTPUT']);p.mkdir(mode=0o700,parents=True,exist_ok=True)
        value=dict(status='FAILED_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],
            no_individual_values_in_error_record=True)
        (p/'private_failed_terminal.json').write_text(json.dumps(value,indent=2)+'\n')
        print(json.dumps(dict(status=value['status'],error_class=value['error_class'])))
        raise SystemExit(1)
