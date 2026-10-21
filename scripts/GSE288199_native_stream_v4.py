"""Bounded official Matrix Market stream with R Matrix block replay; RAM only."""
import csv
import gzip
import hashlib
import io
import re
import urllib.request
import urllib.parse
import zlib

import numpy as np
from GSE288199_R_block_reference_v2 import RReference

BASE='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288199/'
PREFIX=65536


def opener():return urllib.request.build_opener(urllib.request.ProxyHandler({}))


def feature_contract(url,genes):
    assert url.startswith(BASE+'suppl/')
    with opener().open(url,timeout=30) as r:
        assert r.status==200
        assert urllib.parse.urlparse(r.url).hostname=='ftp.ncbi.nlm.nih.gov'
        raw=r.read(2*1024**2+1);assert len(raw)<=2*1024**2
    data=gzip.decompress(raw);assert len(data)<=16*1024**2
    rows=list(csv.reader(io.StringIO(data.decode()),delimiter='\t'))
    assert rows and all(len(r)==3 for r in rows)
    kind=np.full(len(rows),-1,np.int8);mito=np.zeros(len(rows),bool);target=[]
    for i,(accession,symbol,typ) in enumerate(rows):
        if typ=='Gene Expression':
            human=re.fullmatch(r'ENSG[0-9]+(?:\.[0-9]+)?',accession) is not None
            kind[i]=0 if human else 1
            mito[i]=human and symbol.startswith('MT-')
        elif typ=='Antibody Capture':kind[i]=2
    assert (kind>=0).all() and (kind==0).any() and mito.any()
    for gene in genes:
        hits=[i for i,r in enumerate(rows) if kind[i]==0 and r[1]==gene]
        assert len(hits)==1;target.append(hits[0])
    return dict(kind=kind,mito=mito,target=np.asarray(target,np.int64),features=len(rows),
        feature_compressed_SHA256=hashlib.sha256(raw).hexdigest(),
        feature_uncompressed_SHA256=hashlib.sha256(data).hexdigest())


def header(stream):
    first=stream.readline().strip()
    assert first in [b'%%MatrixMarket matrix coordinate integer general',b'%%MatrixMarket matrix coordinate real general']
    line=stream.readline()
    while line.startswith(b'%'):line=stream.readline()
    dims=tuple(map(int,line.split()))
    assert len(dims)==3 and all(x>0 for x in dims)
    return first.decode(),dims


def census(url,contract):
    assert url.startswith(BASE+'suppl/')
    req=urllib.request.Request(url,headers={'Range':f'bytes=0-{PREFIX-1}'})
    with opener().open(req,timeout=30) as r:
        assert r.status in [200,206]
        assert urllib.parse.urlparse(r.url).hostname=='ftp.ncbi.nlm.nih.gov'
        size=int(r.headers['Content-Length'])
        if r.status==206:
            match=re.fullmatch(r'bytes 0-([0-9]+)/([0-9]+)',r.headers['Content-Range'])
            assert match;size=int(match[2]);expected=min(PREFIX,size)
            assert int(match[1])+1==expected
        else:expected=min(PREFIX,size)
        raw=r.read(expected);assert len(raw)==expected
        modified=r.headers.get('Last-Modified');assert modified
    decoded=zlib.decompressobj(16+zlib.MAX_WBITS).decompress(raw,1024**2)
    fmt,dims=header(io.BytesIO(decoded))
    assert dims[0]==contract['features'] and dims[1]<=1000000
    return dict(url=url,bytes=size,modified=modified,prefix_SHA256=hashlib.sha256(raw).hexdigest(),
                format=fmt,dims=dims,**contract)


def chunks(stream,size=4*1024**2):
    tail=b''
    while True:
        part=stream.read(size)
        if not part:
            if tail.strip():yield np.fromstring(tail.decode('ascii'),sep=' ').reshape(-1,3)
            return
        block=tail+part;last=block.rfind(b'\n')
        assert last>=0 or len(block)<1024
        if last<0:tail=block;continue
        tail=block[last+1:]
        content=block[:last+1]
        if content.strip():yield np.fromstring(content.decode('ascii'),sep=' ').reshape(-1,3)


class Tracked:
    def __init__(self,raw,limit):
        self.raw=raw;self.limit=limit;self.bytes=0;self.digest=hashlib.sha256();self.prefix=bytearray()
    def read(self,n=-1):
        assert n>=0
        data=self.raw.read(n)
        self.bytes+=len(data);assert self.bytes<=self.limit
        self.digest.update(data)
        self.prefix.extend(data[:max(0,PREFIX-len(self.prefix))])
        return data


def _consume(stream,dims,kind,mito,target,block_size=4*1024**2,reference=None):
    nf,nc,nnz=dims;k=len(target)
    counts=np.zeros((nc,k),np.int64);totals=np.zeros((nc,3),np.int64)
    detected=np.zeros(nc,np.int64);mito_total=np.zeros(nc,np.int64)
    slot=np.full(nf,-1,np.int64);slot[target]=np.arange(k)
    entries=0;zeros=0;cell_possible=True;gene_possible=True
    previous_cell=-1;previous_gene=-1
    for a in chunks(stream,block_size):
        assert len(a) and np.isfinite(a).all() and np.array_equal(a,np.floor(a))
        assert a.max()<2**53 and a.min()>=0
        gene=a[:,0].astype(np.int64)-1;cell=a[:,1].astype(np.int64)-1;value=a[:,2].astype(np.int64)
        assert gene.min()>=0 and gene.max()<nf and cell.min()>=0 and cell.max()<nc
        by_cell=cell*nf+gene;by_gene=gene*nc+cell
        cell_possible=cell_possible and by_cell[0]>previous_cell and np.all(np.diff(by_cell)>0)
        gene_possible=gene_possible and by_gene[0]>previous_gene and np.all(np.diff(by_gene)>0)
        assert cell_possible or gene_possible,'Duplicate_or_unsorted_coordinates'
        previous_cell=int(by_cell[-1]);previous_gene=int(by_gene[-1]);positive=value>0
        for j in range(3):
            use=kind[gene]==j
            sums=np.bincount(cell[use],weights=value[use],minlength=nc)
            assert sums.max(initial=0)<2**53
            totals[:,j]+=sums.astype(np.int64)
        human_positive=(kind[gene]==0)&positive
        detected+=np.bincount(cell[human_positive],minlength=nc)
        use=mito[gene]
        mito_total+=np.bincount(cell[use],weights=value[use],minlength=nc).astype(np.int64)
        use=slot[gene]>=0
        np.add.at(counts.reshape(-1),cell[use]*k+slot[gene[use]],value[use])
        raw=np.bincount(cell[use]*k+slot[gene[use]],weights=value[use],minlength=nc*k).reshape(nc,k)
        block_totals=np.column_stack([np.bincount(cell[kind[gene]==j],weights=value[kind[gene]==j],minlength=nc).astype(np.int64) for j in range(3)])
        block_detected=np.bincount(cell[human_positive],minlength=nc)
        block_mito=np.bincount(cell[mito[gene]],weights=value[mito[gene]],minlength=nc).astype(np.int64)
        reference.check(a,raw.astype(np.int64),block_totals,block_detected,block_mito)
        entries+=len(a);zeros+=int((~positive).sum())
    assert entries==nnz and np.all(mito_total<=totals[:,0])
    assert np.all(counts<=totals[:,0,None])
    assert totals.max(initial=0)<2**53
    ordering='cell_then_feature' if cell_possible else 'feature_then_cell'
    return counts,totals,detected,mito_total,dict(entries=entries,explicit_zero_entries=zeros,ordering=ordering)


def consume(stream,dims,kind,mito,target,block_size=4*1024**2):
    with RReference(dims,kind,mito,target) as reference:
        return _consume(stream,dims,kind,mito,target,block_size,reference)


def read_full(contract):
    with opener().open(contract['url'],timeout=45) as response:
        assert response.status==200 and int(response.headers['Content-Length'])==contract['bytes']
        assert urllib.parse.urlparse(response.url).hostname=='ftp.ncbi.nlm.nih.gov'
        assert response.headers.get('Last-Modified')==contract['modified']
        tracked=Tracked(response,contract['bytes'])
        with gzip.GzipFile(fileobj=tracked,mode='rb') as stream:
            fmt,dims=header(stream)
            assert fmt==contract['format'] and dims==contract['dims']
            result=consume(stream,dims,contract['kind'],contract['mito'],contract['target'])
        assert tracked.bytes==contract['bytes']
        assert hashlib.sha256(tracked.prefix).hexdigest()==contract['prefix_SHA256']
        return result,tracked.digest.hexdigest()


def synthetic():
    entries=[(1,1,2),(2,1,0),(3,1,5),(2,2,3),(1,3,4),(3,3,7)]
    proofs=[]
    for order in ['cell_then_feature','feature_then_cell']:
        rows=sorted(entries,key=(lambda r:(r[1],r[0])) if order=='cell_then_feature' else (lambda r:(r[0],r[1])))
        text='%%MatrixMarket matrix coordinate integer general\n% fixture\n3 3 6\n'+'\n'.join(' '.join(map(str,r)) for r in rows)
        raw=text.encode();f=io.BytesIO(raw);fmt,dims=header(f)
        got=consume(f,dims,np.array([0,0,2]),np.array([False,True,False]),np.array([0,1]),block_size=9)
        expected=np.zeros((3,3),np.int64)
        for gene,cell,value in entries:expected[cell-1,gene-1]+=value
        assert np.array_equal(got[0],expected[:,:2])
        assert np.array_equal(got[1][:,0],expected[:,:2].sum(axis=1))
        assert np.array_equal(got[2],(expected[:,:2]>0).sum(axis=1))
        assert got[-1]['explicit_zero_entries']==1
        proofs.append(order)
    bad=io.BytesIO(b'1 1 2\n1 1 3\n')
    try:consume(bad,(3,3,2),np.array([0,0,2]),np.array([False,True,False]),np.array([0,1]))
    except (AssertionError,ValueError):duplicate_rejected=True
    else:duplicate_rejected=False
    assert duplicate_rejected
    return dict(both_coordinate_orders_and_tiny_chunk_boundaries=proofs,
                zero_entries_human_antibody_mito_scalar_fixture_and_R_Matrix_parity=True,duplicate_coordinates_rejected=True)
