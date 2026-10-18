"""Read the same effect-blind public pilot in RAM; persist only small schema receipts.

No existing project matrix is exported, no obs values or clonotype strings read.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import time
from urllib.request import Request,ProxyHandler,build_opener

import h5py
import numpy as np
from H5AD_public_schema_math_v1 import string,matrix_structure,numerical

URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE328nnn/GSE328481/soft/GSE328481_family.soft.gz'
EXPECTED_SOFT='1dd821744732f4d8b371873bd26082bfda3f75d033aacfcb1310b772a98cb418'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def opener(proxy):return build_opener(ProxyHandler({'https':proxy} if proxy else {}))


def main(args):
    out=args.output;assert not out.exists();out.mkdir(mode=0o700)
    with opener(args.proxy).open(URL,timeout=30) as r:source=r.read(1024**2+1)
    assert len(source)<=1024**2 and hashlib.sha256(source).hexdigest()==EXPECTED_SOFT
    text=gzip.decompress(source).decode();links=[];gsm=None
    for line in text.splitlines():
        if line.startswith('^SAMPLE = '):gsm=line.split(' = ',1)[1]
        if re.match(r'^!Sample_supplementary_file(?:_\d+)? = ',line):
            u=line.split(' = ',1)[1]
            if u!='NONE':links.append(dict(GSM=gsm,url=re.sub(r'^ftp:','https:',u)))
    assert len(links)==11
    def head(row):
        with opener(args.proxy).open(Request(row['url'],method='HEAD'),timeout=30) as r:
            assert r.status==200
            return dict(**row,bytes=int(r.headers['Content-Length']),last_modified=r.headers.get('Last-Modified'),etag=r.headers.get('ETag'))
    with ThreadPoolExecutor(max_workers=2) as ex:inventory=list(ex.map(head,links))
    chosen=min(inventory,key=lambda r:(r['bytes'],r['url']));size=chosen['bytes'];assert size==108991639 and size<128*1024**2
    (out/'public_official_all11_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    panels=json.loads(args.panel_plan.read_text())['panels']
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),shared_math_sha256=sha(Path(__file__).with_name('H5AD_public_schema_math_v1.py')),
        source_SOFT_SHA256=EXPECTED_SOFT,selection='All11 official links independently reconstructed. Same smallest byte count/URL tie rule, fixed before expression. Same pilot as pending remote attempts.',
        chosen=chosen,maximum_public_input_bytes=128*1024**2,threads=4,range_bytes=4*1024**2,local_expression_file_written=False,
        panel_plan_sha256=sha(args.panel_plan),scope='Full public file in RAM, every X/raw/layer numerical value and sparse structure, every available gene-symbol axis, all spatial coordinates; no obs metadata values, clones, outcomes or existing server matrices read.',
        limitations='Same numerical parser kernel as remote schema code; independent byte retrieval, not independent algorithm or biological truth. A single pilot does not qualify all11 tissues, raw scale, cell states or patient independence.')
    (out/'plan_before_in_memory_public_expression_read.json').write_text(json.dumps(plan,indent=2)+'\n')
    bounds=[(i,min(size-1,i+4*1024**2-1)) for i in range(0,size,4*1024**2)]
    def fetch(bound):
        lo,hi=bound
        with opener(args.proxy).open(Request(chosen['url'],headers={'Range':f'bytes={lo}-{hi}'}),timeout=40) as r:
            assert r.status==206 and r.headers.get('Content-Range')==f'bytes {lo}-{hi}/{size}'
            assert r.headers.get('Last-Modified')==chosen['last_modified']
            if chosen['etag']:assert r.headers.get('ETag')==chosen['etag']
            a=r.read(hi-lo+2);assert len(a)==hi-lo+1;return a
    buffer=io.BytesIO();digest=hashlib.sha256();started=time.monotonic()
    with ThreadPoolExecutor(max_workers=4) as ex:
        pending={i:ex.submit(fetch,bounds[i]) for i in range(min(4,len(bounds)))}
        for i in range(len(bounds)):
            b=pending.pop(i).result();buffer.write(b);digest.update(b)
            if i+4<len(bounds):pending[i+4]=ex.submit(fetch,bounds[i+4])
            print(json.dumps(dict(stage='public_RAM_only_download',ranges=i+1,total=len(bounds),bytes=buffer.tell(),elapsed_seconds=time.monotonic()-started)),flush=True)
    assert buffer.tell()==size and hashlib.sha256(buffer.getbuffer()).hexdigest()==digest.hexdigest();buffer.seek(0)
    assert buffer.read(8)==b'\x89HDF\r\n\x1a\n';buffer.seek(0)
    structures=[];matrices={};coords={};axes={}
    with h5py.File(buffer,'r') as h:
        def visit(name,obj):
            if isinstance(obj,h5py.Dataset):structures.append(dict(path=name,shape=list(obj.shape),dtype=str(obj.dtype),compression=obj.compression))
        h.visititems(visit)
        paths=['X']+(['raw/X'] if 'raw/X' in h else [])+(['layers/'+n for n in h['layers']] if 'layers' in h else [])
        for name in paths:matrices[name]=dict(structure=matrix_structure(h[name]),numerical=numerical(h[name]))
        obs=h['obs'];var=h['var'];oi=string(obs.attrs.get('_index','_index'));vi=string(var.attrs.get('_index','_index'))
        nr,nc=len(obs[oi]),len(var[vi]);assert matrices['X']['structure']['shape']==[nr,nc]
        for prefix in ['var','raw/var']:
            if prefix not in h:continue
            group=h[prefix];index=string(group.attrs.get('_index','_index'))
            for name in sorted(set([index])|{'gene_symbol','gene_symbols','gene_name','gene_names','SYMBOL','symbol'}):
                if name not in group:continue
                ds=group[name]
                if not isinstance(ds,h5py.Dataset) or ds.dtype.kind not in ['O','S','U']:continue
                values=np.array([string(v) for v in ds[:]]);counter=Counter(values);unique={g for g,n in counter.items() if n==1}
                axes[prefix+'/'+name]=dict(n=len(values),unique=len(counter),ensembl_like=int(np.char.startswith(values,'ENSG').sum()),sha256=hashlib.sha256(json.dumps(values.tolist(),ensure_ascii=False,separators=(',',':')).encode()).hexdigest(),fixed_panels={p:dict(required=len(gs),uniquely_measured=len(set(gs)&unique),complete=set(gs)<=unique) for p,gs in panels.items()})
        if 'obsm' in h:
            for name in ['spatial','X_spatial']:
                if name not in h['obsm']:continue
                a=np.asarray(h['obsm'][name]);assert len(a)==nr and a.ndim==2 and a.shape[1]>=2 and np.isfinite(a).all()
                coords[name]=dict(shape=list(a.shape),unique_xy=int(len(np.unique(a[:,:2],axis=0))),minimum=a.min(axis=0).tolist(),maximum=a.max(axis=0).tolist(),physical_unit='UNVERIFIED')
        summary=dict(observations=nr,variables=nc,obs_columns=list(obs.keys()),var_columns=list(var.keys()),layers=list(h['layers'].keys()) if 'layers' in h else [],uns_keys=list(h['uns'].keys()) if 'uns' in h else [],obsm_keys=list(h['obsm'].keys()) if 'obsm' in h else [],matrices=matrices,gene_axis_candidates=axes,coordinates=coords,raw_scale='NOT_PROVEN_BY_INTEGERS_OR_SCHEMA',cell_vs_bin='UNVERIFIED',obs_values_read=False,independent_patients_reviewed=False)
    buffer.close();(out/'private_full_HDF5_dataset_schema.json').write_text(json.dumps(structures,indent=2)+'\n')
    (out/'private_complete_pilot_schema_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    files=[p for p in out.iterdir() if p.is_file()];assert sum(p.stat().st_size for p in files)<8*1024**2
    report=dict(status='COMPLETE_FULL_PUBLIC_PILOT_RAM_SCHEMA_ALL_VALUES_SCANNED_NO_LOCAL_MATRIX',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),shared_math_sha256=plan['shared_math_sha256'],public_file_sha256=digest.hexdigest(),public_file_bytes=size,h5py_version=h5py.__version__,numpy_version=np.__version__,local_expression_persisted=False,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},new_effects_or_all_source_qualification=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(status=report['status'],public_file_sha256=report['public_file_sha256'],terminal_sha256=sha(out/'private_terminal.json'))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--panel-plan',type=Path,required=True);p.add_argument('--proxy',default='');main(p.parse_args())
