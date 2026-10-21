"""Same frozen RNA eligibility, buffer exact expected public sources in RAM.

Strict full response/length/SHA checks precede parsing. Cache only the two
declared files; source changes, permission errors and wrong complete hashes
are terminal. Bounded transient failures can retry without changing origin.
No compressed source, individual annotation rows or arrays are saved.
"""
from datetime import datetime,timezone
from pathlib import Path
import hashlib
import importlib.util
import io
import json
import os
import signal
import time
import urllib.request
import urllib.error
import http.client

SOURCES=[('GSE286827_metadata_CD4T_TCR.rds.gz','b33adaaba779a271df11aa6834266133aa41a1952b745770124221909db28210'),
         ('GSE286827_metadata_CD8T_TCR.rds.gz','dc78fc347c8b400e6ebf1add2492b1c7198ab68ccf833abae543db4dfbcf1343')]
PREFIX='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/'
CAP=8*1024**2


class SourceDeadlineExceeded(RuntimeError):
    pass


def fetch_exact(opener,url,expected):
    with opener(url,timeout=20) as response:
        assert response.status==200,'not_full_200_response'
        assert response.geturl()==url,'source_redirect_not_accepted'
        total=response.headers.get('Content-Length')
        if total is not None:assert 0<int(total)<=CAP,'source_size_outside_cap'
        chunks=[];size=0
        while True:
            block=response.read(min(65536,CAP+1-size))
            if not block:break
            chunks.append(block);size+=len(block)
            assert size<=CAP,'source_size_outside_cap'
        if total is not None and size!=int(total):raise http.client.IncompleteRead(b'',int(total)-size)
        raw=b''.join(chunks)
        assert raw.startswith(bytes([31,139])),'unexpected_source_encoding'
        assert hashlib.sha256(raw).hexdigest()==expected,'complete_source_SHA_mismatch'
        return raw,dict(HTTP_status=200,full_response_bytes=len(raw),SHA256=expected,
            content_length_exact_or_absent=True,original_destination_preserved=True)


def acquire(opener,url,expected,sleeper=time.sleep):
    failures=[]
    for attempt in range(3):
        try:
            raw,proof=fetch_exact(opener,url,expected)
            proof.update(attempts=attempt+1,transient_failure_classes=failures)
            return raw,proof
        except urllib.error.HTTPError as error:
            # No permission/rate-limit/redirect error is retried by this reader.
            raise
        except (urllib.error.URLError,TimeoutError,http.client.IncompleteRead,ConnectionError) as error:
            failures.append(type(error).__name__)
            if attempt==2:raise
            sleeper(2*(attempt+1))


def main():
    folder=Path(__file__).parent;original=folder/'qualify_CM02_RNA_programme_eligibility_v1.py'
    assert hashlib.sha256(original.read_bytes()).hexdigest()=='a05c1d06b3c873b8abc0e78da20e420b20a0768e241487a43b17439fbfa2d3bd'
    root=folder.parent;diagnostic=root/'private_source_cache_qualification_v2.json'
    assert not diagnostic.exists()
    old_urlopen=urllib.request.urlopen;old_handler=signal.getsignal(signal.SIGALRM)
    def deadline(*args):raise SourceDeadlineExceeded('bounded_public_source_deadline')
    signal.signal(signal.SIGALRM,deadline)
    cache={};proofs=[]
    try:
        for name,expected in SOURCES:
            signal.alarm(100)
            raw,proof=acquire(old_urlopen,PREFIX+name,expected)
            signal.alarm(0);proof['public_file']=name;cache[PREFIX+name]=raw;proofs.append(proof)
        diagnostic.write_text(json.dumps(dict(status='COMPLETE_EXPECTED_SOURCE_BYTES_BUFFERED_RAM_ONLY',
            completed_UTC=datetime.now(timezone.utc).isoformat(),sources=proofs,no_source_payloads_or_rows_saved=True,
            wrapper_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+chr(10))
        def cached(url,*args,**kwargs):
            key=url.full_url if hasattr(url,'full_url') else str(url)
            assert key in cache,'undeclared_network_request_in_frozen_input_loader'
            return io.BytesIO(cache[key])
        urllib.request.urlopen=cached
        spec=importlib.util.spec_from_file_location('same_RNA_input_qualification',original)
        core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core);core.main()
    except BaseException as error:
        # The frozen core already records its detailed sanitized input failure.
        if not diagnostic.exists():
            diagnostic.write_text(json.dumps(dict(status='EXPECTED_SOURCE_BYTES_NOT_ACCEPTED',
                observed_UTC=datetime.now(timezone.utc).isoformat(),completed_sources=proofs,
                error_class=type(error).__name__,HTTP_status=getattr(error,'code',None),
                no_full_error_or_individual_values_returned=True),indent=2)+chr(10))
        raise
    finally:
        urllib.request.urlopen=old_urlopen;signal.alarm(0);signal.signal(signal.SIGALRM,old_handler);cache.clear()


if __name__=='__main__':main()
