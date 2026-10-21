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

SOURCES={'GSE286827_metadata_immunecells_broad.rds.gz':'7835c18873d51f8bd571f6bf21ade5d0f18c4eac78d942d26a752a03f5d419a6',
         'GSE286827_metadata_CD4T_TCR.rds.gz':'b33adaaba779a271df11aa6834266133aa41a1952b745770124221909db28210',
         'GSE286827_metadata_CD8T_TCR.rds.gz':'dc78fc347c8b400e6ebf1add2492b1c7198ab68ccf833abae543db4dfbcf1343'}
PREFIX='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/'
CAP=8*1024**2


def receive_exact_sources(stream):
    raw=stream.read(26*1024**2+1);assert len(raw)<=26*1024**2
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

