"""Consume exact public source bytes from stdin, all source payloads RAM only.

The original NCBI files are independently SHA checked at both ends. Nothing
is extracted to a filesystem. The accepted frozen loader selects its original
four annotation fields; no local private/derived records are supplied.
This is source transport, not a changed annotation or a new clinical family.
"""
from datetime import datetime,timezone
from pathlib import Path
import hashlib
import importlib.util
import io
import json
import sys
import tarfile
import urllib.request

SOURCES={'GSE286827_metadata_CD4T_TCR.rds.gz':'b33adaaba779a271df11aa6834266133aa41a1952b745770124221909db28210',
         'GSE286827_metadata_CD8T_TCR.rds.gz':'dc78fc347c8b400e6ebf1add2492b1c7198ab68ccf833abae543db4dfbcf1343'}
PREFIX='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/'
CAP=8*1024**2


def receive_exact_sources(stream):
    raw=stream.read(20*1024**2+1);assert len(raw)<=20*1024**2
    cache={};proof=[]
    with tarfile.open(fileobj=io.BytesIO(raw),mode='r:') as archive:
        members=archive.getmembers()
        assert len(members)==len(SOURCES) and {x.name for x in members}==set(SOURCES)
        for member in members:
            assert member.isfile() and 0<member.size<=CAP
            assert member.name==Path(member.name).name
            body=archive.extractfile(member).read()
            assert len(body)==member.size and body.startswith(bytes([31,139]))
            digest=hashlib.sha256(body).hexdigest();assert digest==SOURCES[member.name]
            cache[PREFIX+member.name]=body
            proof.append(dict(public_file=member.name,bytes=len(body),SHA256=digest))
    return cache,proof


def with_cached_sources(cache,callback):
    original=urllib.request.urlopen
    def selected(url,*args,**kwargs):
        key=url.full_url if hasattr(url,'full_url') else str(url)
        assert key in cache,'undeclared_source_request'
        return io.BytesIO(cache[key])
    try:
        urllib.request.urlopen=selected
        return callback()
    finally:urllib.request.urlopen=original


def main():
    root=Path(__file__).parent.parent;cache={}
    try:
        cache,proof=receive_exact_sources(sys.stdin.buffer)
        (root/'private_public_source_RAM_receipt_v4.json').write_text(json.dumps(dict(
            status='COMPLETE_ORIGINAL_PUBLIC_BYTES_VERIFIED_IN_RECEIVER_RAM',
            completed_UTC=datetime.now(timezone.utc).isoformat(),sources=proof,
            source_payloads_not_extracted_or_saved=True,annotation_and_full_SHA_unchanged=True,
            no_local_private_or_derived_records_received=True,
            receiver_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+chr(10))
        file=Path(__file__).with_name('qualify_CM02_RNA_programme_eligibility_v1.py')
        assert hashlib.sha256(file.read_bytes()).hexdigest()=='a05c1d06b3c873b8abc0e78da20e420b20a0768e241487a43b17439fbfa2d3bd'
        spec=importlib.util.spec_from_file_location('unchanged_RNA_input_qualification',file)
        core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
        with_cached_sources(cache,core.main)
    finally:cache.clear()


if __name__=='__main__':main()
