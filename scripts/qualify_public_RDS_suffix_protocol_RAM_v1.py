"""One official metadata GET and two tiny ranges; headers only, no object acceptance."""
from urllib.request import Request, build_opener, ProxyHandler
from urllib.parse import quote, urlparse
from datetime import datetime, timezone
import hashlib
import json
import sys


def main():
    plan=json.load(sys.stdin)
    opener=build_opener(ProxyHandler({}))
    metadata_url='https://zenodo.org/api/records/'+plan['record']
    with opener.open(Request(metadata_url),timeout=40) as r:
        assert r.status==200 and r.url==metadata_url
        payload=r.read(1024**2+1);assert len(payload)<=1024**2
    metadata=json.loads(payload)
    assert str(metadata['id'])==plan['record'] and metadata['metadata']['access_right']=='open'
    files=[x for x in metadata['files'] if x['key'].lower().endswith('.rds')]
    assert len(files)==1
    f=files[0]
    assert f['size']==plan['expected_file_bytes'] and f['checksum']=='md5:'+plan['expected_file_MD5']
    assert hashlib.sha256(f['key'].encode()).hexdigest()==plan['file_key_SHA256']
    url='https://zenodo.org/records/'+plan['record']+'/files/'+quote(f['key'],safe='')+'?download=1'
    etag=None;records=[]
    for start,end in [(0,31),(32,63)]:
        headers={'Range':f'bytes={start}-{end}'}
        if etag is not None:headers['If-Match']=etag
        with opener.open(Request(url,headers=headers),timeout=45) as r:
            assert r.status==206 and urlparse(r.url).scheme=='https' and urlparse(r.url).hostname=='zenodo.org'
            assert r.headers.get('Content-Range')==f'bytes {start}-{end}/{plan["expected_file_bytes"]}'
            declared=r.headers.get('Content-Length');assert declared is None or int(declared)==32
            current=r.headers.get('ETag')
            if etag is not None and current is not None:assert current==etag
            used=etag is not None
            if etag is None and current is not None:
                assert '\n' not in current and '\r' not in current and len(current)<=256
                if current.startswith('"') and current.endswith('"'):etag=current
            part=r.read(33);assert len(part)==32
            records.append(dict(exact_range_header=True,exact_body_length=True,HTTP_status=r.status,
                conditional_IfMatch_used=used,payload_SHA256=hashlib.sha256(part).hexdigest()))
    print(json.dumps(dict(status='OFFICIAL_TINY_RANGE_AND_CONDITIONAL_SUFFIX_PROTOCOL_QUALIFIED',
        completed_UTC=datetime.now(timezone.utc).isoformat(),records=records,
        strong_entity_tag_available=etag is not None,
        entity_tag_SHA256=hashlib.sha256(etag.encode()).hexdigest() if etag else None,
        no_authentication_protocol_or_endpoint_retry=True,
        only_64_compressed_bytes_read_in_RAM_no_individual_values_returned=True,
        complete_object_integrity_or_native_analysis_not_accepted=True)))


if __name__=='__main__':main()
