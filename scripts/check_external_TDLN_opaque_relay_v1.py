"""Small artificial byte streams only; never native or individual files."""
from datetime import datetime, timezone
from urllib.error import HTTPError, URLError
import hashlib
import io
import json
import re
import ssl
import threading
import time
import receive_external_TDLN_opaque_relay_RAM_v1 as receiver
import stream_public_RDS_ordered_RAM_v1 as transport


URL = 'https://zenodo.org/records/1/files/ARTIFICIAL_BYTES?download=1'
PAYLOAD = bytes(range(256))*40+b'ARTIFICIAL_PROTOCOL_ONLY'


class Response(io.BytesIO):
    def __init__(self,content,start,end,total,*,etag='"FIXED"',url=URL):
        super().__init__(content)
        self.status,self.url = 206,url
        self.headers = {'Content-Range':f'bytes {start}-{end}/{total}',
                        'Content-Length':str(end-start+1),'ETag':etag}
    def read1(self,n):
        return self.read(n)


class Origin:
    def __init__(self,mode='SHORT'):
        self.mode = mode;self.calls = [];self.lock = threading.Lock()
    def open(self,request,timeout):
        start,end = map(int,re.fullmatch(r'bytes=(\d+)-(\d+)',request.get_header('Range')).groups())
        with self.lock:
            self.calls.append((start,end));call = len(self.calls)
        if self.mode=='DENIED':
            raise HTTPError(URL,403,'ARTIFICIAL_DENIAL',{},None)
        if self.mode=='CERTIFICATE':
            raise URLError(ssl.SSLCertVerificationError('ARTIFICIAL_CERTIFICATE_ERROR'))
        if self.mode=='CONNECTION':
            raise URLError(ConnectionRefusedError(111,'ARTIFICIAL_TRANSIENT'))
        body = PAYLOAD[start:end+1]
        if self.mode in ('SHORT','ETAG') and start%1024==0 and len(body)>600:
            body = body[:600]
        etag = '"CHANGED"' if self.mode=='ETAG' and call>1 else '"FIXED"'
        response = Response(body,start,end,len(PAYLOAD),etag=etag)
        if self.mode=='HEADER':
            response.headers['Content-Range'] = 'bytes 0-1/2'
        if self.mode=='REDIRECT':
            response.url = 'http://other.invalid/ARTIFICIAL'
        # Earlier logical spans deliberately finish later.
        time.sleep(.001*(3-start//1024%4))
        return response


def rejected(call,code=None):
    try:
        call()
    except (AssertionError,ValueError,transport.StreamInputError,HTTPError,URLError) as error:
        if code is not None:
            assert getattr(error,'code_name',None)==code
        return True
    raise AssertionError('INVALID_ARTIFICIAL_STREAM_ACCEPTED')


def check():
    md5,sha = hashlib.md5(PAYLOAD).hexdigest(),hashlib.sha256(PAYLOAD).hexdigest()
    def send(origin,sink=None,**changes):
        options = dict(expected_bytes=len(PAYLOAD),expected_MD5=md5,expected_SHA256=sha,
                       chunk_bytes=1024,workers=4,request_timeout=1,maximum_range_seconds=5,
                       maximum_attempt_seconds=1,maximum_attempts=4,retry_backoff_seconds=0,
                       opener_factory=lambda:origin)
        options.update(changes)
        return transport.stream_ranges(URL,io.BytesIO() if sink is None else sink,**options)
    origin = Origin();sink = io.BytesIO();proof,m,s = send(origin,sink)
    assert sink.getvalue()==PAYLOAD and m==md5 and s==sha
    assert proof['logical_ranges_requiring_continuation']==10
    assert proof['actual_request_attempts']==21 and proof['maximum_actual_attempts_per_range']==2
    assert proof['maximum_pending_ranges']==4 and proof['ordered_stream_no_local_file_sized_buffer'] is True
    assert (600,1023) in origin.calls and (0,599) not in origin.calls
    footer = dict(kind='COMPLETE_ORDERED_PUBLIC_INPUT_PROOF',MD5=m,SHA256=s,range_assembly_proof=proof)
    wire = PAYLOAD+(json.dumps(footer)+'\n').encode()
    plan = dict(expected_file_bytes=len(PAYLOAD),expected_file_MD5=md5,expected_file_SHA256=sha,
                range_chunk_bytes=1024,range_workers=4)
    buffer,rm,rs,rproof = receiver.receive_payload(io.BytesIO(wire),plan)
    buffer.seek(0);assert buffer.read()==PAYLOAD and rm==m and rs==s
    assert rproof['whole_bytes_independently_hashed_on_source_before_R'] is True
    buffer.close()
    rejected(lambda:receiver.receive_payload(io.BytesIO(wire[:-10]),plan))
    rejected(lambda:receiver.receive_payload(io.BytesIO(wire+b'EXCESS'),plan))
    altered = b'X'+wire[1:]
    rejected(lambda:receiver.receive_payload(io.BytesIO(altered),plan))
    rejected(lambda:receiver.receive_payload(io.BytesIO(),dict(plan,expected_file_bytes=2**31+1)))
    for mode in ('DENIED','CERTIFICATE','HEADER','REDIRECT'):
        bad = Origin(mode);rejected(lambda:send(bad,workers=1))
        assert len(bad.calls)==1
    transient = Origin('CONNECTION');rejected(lambda:send(transient,workers=1))
    assert len(transient.calls)==4
    rejected(lambda:send(Origin('ETAG'),workers=1),code='OBJECT_ETAG_CHANGED')
    rejected(lambda:send(Origin('FULL'),expected_MD5='0'*32),code='WHOLE_OBJECT_MD5_MISMATCH')
    rejected(lambda:send(Origin('FULL'),expected_SHA256='0'*64),code='WHOLE_OBJECT_SHA256_MISMATCH')
    class ShortSink:
        def write(self,block):
            return len(block)-1
    rejected(lambda:send(Origin('FULL'),ShortSink()),code='SINK_SHORT_WRITE')
    return dict(status='PASS_ARTIFICIAL_ORDERED_BOUNDED_RELAY_AND_SOURCE_HASH_GATES',
                reordered_completion_suffix_only_recovery_and_ordered_sink_checked=True,
                authentication_certificate_header_origin_and_hash_errors_terminal=True,
                short_excess_corrupted_input_and_resource_cap_rejected_before_R=True,
                source_independent_full_hash_and_no_local_file_sized_buffer_protocol=True,
                no_real_native_individual_clinical_or_payload_files=True,
                observed_UTC=datetime.now(timezone.utc).isoformat())


if __name__=='__main__':
    print(json.dumps(check()))
