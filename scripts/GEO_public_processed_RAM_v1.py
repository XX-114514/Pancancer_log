"""Bounded four-range public H5AD reader; no disk persistence or project data I/O."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
from urllib.parse import urlparse
from urllib.request import build_opener,ProxyHandler,Request


def read_public(chosen,proxy,progress=None):
    url=chosen['url'];size=int(chosen['bytes']);p=urlparse(url)
    assert p.scheme=='https' and p.hostname=='ftp.ncbi.nlm.nih.gov' and p.path.endswith('.h5ad')
    assert 0<size<2*1024**3
    bounds=[(i,min(size-1,i+4*1024**2-1)) for i in range(0,size,4*1024**2)]
    def fetch(bound):
        lo,hi=bound;o=build_opener(ProxyHandler({'https':proxy} if proxy else {}))
        with o.open(Request(url,headers={'Range':f'bytes={lo}-{hi}'}),timeout=40) as r:
            assert r.status==206 and r.headers.get('Content-Range')==f'bytes {lo}-{hi}/{size}'
            assert r.headers.get('Last-Modified')==chosen['last_modified']
            if chosen.get('etag'):assert r.headers.get('ETag')==chosen['etag']
            b=r.read(hi-lo+2);assert len(b)==hi-lo+1;return b
    buffer=io.BytesIO();digest=hashlib.sha256()
    with ThreadPoolExecutor(max_workers=4) as ex:
        pending={i:ex.submit(fetch,bounds[i]) for i in range(min(4,len(bounds)))}
        for i in range(len(bounds)):
            b=pending.pop(i).result();buffer.write(b);digest.update(b)
            if i+4<len(bounds):pending[i+4]=ex.submit(fetch,bounds[i+4])
            if progress:progress(i+1,len(bounds),buffer.tell())
    assert buffer.tell()==size and hashlib.sha256(buffer.getbuffer()).hexdigest()==digest.hexdigest()
    buffer.seek(0);assert buffer.read(8)==b'\x89HDF\r\n\x1a\n';buffer.seek(0)
    return buffer,digest.hexdigest()
