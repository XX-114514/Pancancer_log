"""Artificial bytes only: known temporary gateway failures vs terminal denial."""
from pathlib import Path
from datetime import datetime, timezone
from email.utils import format_datetime
from urllib.error import HTTPError
import hashlib
import io
import json
import re


def check():
    code=Path(__file__).parent
    namespace={'RANGE_CORE_SOURCE':(code/'fetch_public_RDS_ranges_RAM_v2.py').read_text()}
    exec((code/'fetch_public_RDS_ranges_resume_RAM_v2.py').read_text(),namespace)
    payload=bytes(range(256))*32;url='https://zenodo.org/records/1/files/ARTIFICIAL'
    class Response(io.BytesIO):
        def __init__(self,value,start,end):
            super().__init__(value);self.status=206;self.url=url
            self.headers={'Content-Range':f'bytes {start}-{end}/{len(payload)}','Content-Length':str(end-start+1)}
        def read1(self,n):return self.read(n)
    class Origin:
        def __init__(self,status=504,always=False,after='0',short=False):
            self.status,self.always,self.after,self.short=status,always,after,short;self.calls=[]
        def open(self,request,timeout):
            start,end=map(int,re.fullmatch(r'bytes=(\d+)-(\d+)',request.get_header('Range')).groups())
            self.calls.append((start,end));n=len(self.calls)
            fail=self.always or (n==(2 if self.short else 1))
            if fail:raise HTTPError(url,self.status,'ARTIFICIAL',{'Retry-After':self.after},io.BytesIO())
            data=payload[start:end+1]
            if self.short and n==1:data=data[:600]
            return Response(data,start,end)
    def fetch(origin,**extra):
        options=dict(expected_bytes=len(payload),expected_MD5=hashlib.md5(payload).hexdigest(),
            chunk_bytes=1024,workers=1,request_timeout=1,maximum_range_seconds=5,
            maximum_attempts=4,retry_backoff_seconds=0,maximum_attempt_seconds=1,
            opener_factory=lambda:origin)
        options.update(extra)
        return namespace['fetch_ranges'](url,**options)
    for status in (502,503,504):
        origin=Origin(status=status,short=True);buffer,m,s,proof=fetch(origin)
        buffer.seek(0);assert buffer.read()==payload;buffer.close()
        assert origin.calls[:3]==[(0,1023),(600,1023),(600,1023)]
        assert proof['transient_HTTP_status_retries']=={str(status):1}
        assert s.hexdigest()==hashlib.sha256(payload).hexdigest()
    def rejected(origin,**extra):
        try:fetch(origin,**extra)
        except Exception as error:return error
        raise AssertionError('ARTIFICIAL_INVALID_TRANSPORT_ACCEPTED')
    for status in (401,403,404,407,429,500):
        origin=Origin(status=status);error=rejected(origin)
        assert isinstance(error,HTTPError) and len(origin.calls)==1
    repeated=Origin(always=True);error=rejected(repeated)
    assert isinstance(error,HTTPError) and len(repeated.calls)==4
    wait=Origin(status=503,after='1000');error=rejected(wait)
    assert getattr(error,'code_name',None)=='SERVER_RETRY_AFTER_EXCEEDS_RANGE_BUDGET' and len(wait.calls)==1
    malformed=Origin(after='not-a-date');error=rejected(malformed)
    assert getattr(error,'code_name',None)=='INVALID_RETRY_AFTER'
    then=datetime(2020,1,1,tzinfo=timezone.utc)
    assert namespace['_retry_after_seconds']({'Retry-After':'2'},then)==2
    assert namespace['_retry_after_seconds']({'Retry-After':format_datetime(datetime(2020,1,1,0,0,3,tzinfo=timezone.utc))},then)==3
    assert namespace['_retry_after_seconds']({'Retry-After':format_datetime(then)},then)==0
    return dict(status='PASS_ARTIFICIAL_502_503_504_SUFFIX_RECOVERY_AND_TERMINAL_DENIAL_GATES',
        short_prefix_not_redownloaded_after_gateway_timeout=True,
        Retry_After_respected_without_exceeding_original_budget=True,
        denial_certificate_protocol_hash_policy_not_broadened=True,
        no_real_native_individual_or_clinical_data_read=True)


if __name__=='__main__':
    print(json.dumps(check()))
