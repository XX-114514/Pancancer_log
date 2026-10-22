"""Synthetic protocol/continuation oracles, zero real network and payload files."""
from pathlib import Path
from datetime import datetime, timezone
from types import SimpleNamespace
from email.message import Message
from urllib.error import HTTPError
from http.client import IncompleteRead
import hashlib
import io
import json
import sys
import threading


def main():
    core = globals().get('RANGE_CORE_SOURCE')
    resume = globals().get('RESUME_SOURCE')
    legacy = globals().get('LEGACY_RANGE_SOURCE')
    if core is None:
        code = Path(__file__).parent
        core = (code/'fetch_public_RDS_ranges_RAM_v1.py').read_text()
        resume = (code/'fetch_public_RDS_ranges_resume_RAM_v1.py').read_text()
        legacy = (code/'fetch_public_RDS_ranges_retry_RAM_v4.py').read_text()
    ns = {'RANGE_CORE_SOURCE': core};exec(resume, ns)
    old = {'RANGE_CORE_SOURCE': core};exec(legacy, old)
    data = bytes(range(251))*49
    url = 'https://zenodo.org/records/SYNTHETIC/files/fixture.rds?download=1'
    requested, lock = [], threading.Lock()
    mode = ['NORMAL'];clock = [0.0]

    class Response(io.BytesIO):
        def __init__(self, body, start, end, special=False):
            super().__init__(body);self.status=206;self.url=url;self.calls=0;self.special=special
            self.headers=Message()
            self.headers['content-range']=f'bytes {start}-{end}/{len(data)}'
            self.headers['content-length']=str(end-start+1)
            self.headers['etag']='"SYNTHETIC_VERSION_1"'
        def read1(self, n=-1):
            self.calls += 1
            if self.special and mode[0] in ['TIMEOUT_PREFIX','SOFT_ATTEMPT']:
                if self.calls==1:
                    if mode[0]=='SOFT_ATTEMPT':clock[0]+=91
                    return super().read1(min(n,500))
                raise TimeoutError('SYNTHETIC_TRANSPORT')
            if self.special and mode[0]=='INCOMPLETE_PREFIX':
                raise IncompleteRead(super().read1(500),500)
            if self.special and mode[0]=='LOGICAL_BUDGET':
                clock[0]+=481;return super().read1(min(n,500))
            return super().read1(n)

    class Opener:
        def open(self, req, timeout):
            assert req.full_url==url
            start,end=map(int,req.get_header('Range')[6:].split('-'))
            assert 0<=start<=end<len(data)
            with lock:requested.append((start,end,req.get_header('If-match')))
            first=start%4096==0
            if mode[0] in ['HTTP401','HTTP403','HTTP429']:
                raise HTTPError(url,int(mode[0][4:]),'SYNTHETIC',{},None)
            body=data[start:end+1]
            if mode[0] in ['SHORT_LARGE','ETAG_CHANGE','BAD_SUFFIX_HEADER'] and len(body)>32:
                body=body[:-31]
            if mode[0]=='ZERO_PROGRESS':body=b''
            if mode[0]=='EXCESS':body+=b'x'
            if mode[0]=='CORRUPT':body=bytes([body[0]^1])+body[1:]
            r=Response(body,start,end,first)
            if mode[0]=='BAD_HEADER' or (mode[0]=='BAD_SUFFIX_HEADER' and not first):r.headers.replace_header('content-range','bytes 0-0/1')
            if mode[0]=='BAD_LENGTH':r.headers.replace_header('content-length',str(end-start+2))
            if mode[0]=='STATUS200':r.status=200
            if mode[0]=='FOREIGN_REDIRECT':r.url='https://example.invalid/fixture.rds'
            if mode[0]=='ETAG_CHANGE' and not first:r.headers.replace_header('etag','"SYNTHETIC_VERSION_2"')
            if mode[0]=='PLAIN_LOWER':r.headers={k.lower():v for k,v in r.headers.items()}
            return r

    def run(module,scenario,expected=None):
        mode[0]=scenario;requested.clear();clock[0]=0
        result=module['fetch_ranges'](url,expected_bytes=len(data),expected_MD5=expected or hashlib.md5(data).hexdigest(),
            chunk_bytes=4096,workers=1,request_timeout=10,maximum_range_seconds=480,
            maximum_attempts=4,retry_backoff_seconds=0,maximum_attempt_seconds=90,opener_factory=Opener)
        b,md,sha,proof=result
        assert b.getvalue()==data and md.hexdigest()==hashlib.md5(data).hexdigest() and sha.hexdigest()==hashlib.sha256(data).hexdigest()
        b.close();return proof,list(requested)

    proof, requests = run(ns,'SHORT_LARGE')
    for end in sorted({r[1] for r in requests}):
        group=[r for r in requests if r[1]==end]
        base=end//4096*4096
        assert group[0][:2]==(base,end)
        if end-base+1>32:
            assert len(group)==2 and group[1][0]==end-30
            assert group[1][2]=='"SYNTHETIC_VERSION_1"'
        else:assert len(group)==1
    assert proof['logical_ranges_requiring_continuation']==sum(end-end//4096*4096+1>32 for end in {r[1] for r in requests})
    # The exact same repeatedly short large responses exhaust the previous algorithm.
    try:run(old,'SHORT_LARGE')
    except old['RangeInputError'] as e:assert e.code_name=='RANGE_LENGTH_SHORT' and e.range_attempt_number==4
    else:raise AssertionError('LEGACY_FAILURE_NOT_REPRODUCED')
    for scenario in ['NORMAL','PLAIN_LOWER','TIMEOUT_PREFIX','INCOMPLETE_PREFIX']:
        run(ns,scenario)
    failures=[]
    for scenario,reason in [('BAD_HEADER','CONTENT_RANGE_MISMATCH'),('BAD_SUFFIX_HEADER','CONTENT_RANGE_MISMATCH'),
                             ('BAD_LENGTH','CONTENT_LENGTH_MISMATCH'),('STATUS200','NOT_EXACT_OFFICIAL_206'),
                             ('FOREIGN_REDIRECT','NOT_EXACT_OFFICIAL_206'),('ETAG_CHANGE','OBJECT_ETAG_CHANGED'),
                             ('EXCESS','RANGE_LENGTH_EXCESS'),('ZERO_PROGRESS','RANGE_LENGTH_SHORT')]:
        try:run(ns,scenario)
        except ns['RangeInputError'] as e:
            assert e.code_name==reason
            if scenario=='ZERO_PROGRESS':assert e.range_attempt_number==4 and e.received_range_bytes==0
            failures.append(reason)
        else:raise AssertionError('NEGATIVE_CASE_ACCEPTED_'+scenario)
    for scenario in ['HTTP401','HTTP403','HTTP429']:
        try:run(ns,scenario)
        except HTTPError as e:assert e.code==int(scenario[4:]) and len(requested)==1
        else:raise AssertionError('HTTP_ERROR_RETRIED')
    try:run(ns,'CORRUPT')
    except ns['RangeInputError'] as e:assert e.code_name=='WHOLE_OBJECT_MD5_MISMATCH'
    else:raise AssertionError('WHOLE_HASH_MISMATCH_ACCEPTED')
    virtual=SimpleNamespace(monotonic=lambda:clock[0]);ns['time']=virtual;ns['_core']['time']=virtual
    run(ns,'SOFT_ATTEMPT')
    try:run(ns,'LOGICAL_BUDGET')
    except ns['RangeInputError'] as e:assert e.code_name=='RANGE_ELAPSED_BUDGET'
    else:raise AssertionError('LOGICAL_BUDGET_NOT_ENFORCED')
    print(json.dumps(dict(status='PASS_VERIFIED_SUFFIX_CONTINUATION_SYNTHETIC',
        completed_UTC=datetime.now(timezone.utc).isoformat(),runtime_Python=sys.version.split()[0],
        logical_ranges_exact_and_whole_MD5_SHA256_oracles=True,legacy_short_body_failure_reproduced=True,
        only_missing_suffix_requested_and_IfMatch_preserved=True,
        timeout_and_incomplete_read_prefix_continuation_checked=True,
        protocol_HTTP_identity_hash_and_no_progress_failures_terminal=True,
        soft_attempt_prefix_preserved_but_logical_budget_terminal=True,
        case_insensitive_headers_checked=True,no_real_network_or_payload_files=True,
        code_SHA256={k:hashlib.sha256(v.encode()).hexdigest() for k,v in [('core',core),('resume',resume),('legacy',legacy)]})))


if __name__=='__main__':main()
