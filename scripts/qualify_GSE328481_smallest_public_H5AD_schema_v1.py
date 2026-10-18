"""Effect-blind smallest public processed H5AD schema and numerical qualification.

Reads matrix/gene/coordinate structure, not cell metadata values or clonotypes.
This pilot cannot qualify all tissues, raw counts, cell states or biology.
"""
import gzip
import hashlib
import json
import os
import re
from datetime import datetime,timezone
from pathlib import Path
from urllib.request import Request,urlopen

import h5py
import numpy as np

from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_native_MFA_interpretation_inputs_v1 import sha

OUT='GSE328481_smallest_public_H5AD_schema_20261018_v1'


def string(x):
    return x.decode('utf-8') if isinstance(x,bytes) else str(x)


def matrix_structure(group):
    if isinstance(group,h5py.Dataset):return dict(storage='dense',shape=list(group.shape),dtype=str(group.dtype))
    encoding=string(group.attrs.get('encoding-type','unknown'))
    assert encoding in ['csr_matrix','csc_matrix'], 'Unsupported matrix encoding: '+encoding
    shape=list(map(int,group.attrs['shape']));assert len(shape)==2
    ptr=np.asarray(group['indptr'],np.int64);axis=0 if encoding=='csr_matrix' else 1
    assert len(ptr)==shape[axis]+1 and ptr[0]==0 and np.all(np.diff(ptr)>=0) and ptr[-1]==len(group['data'])==len(group['indices'])
    limit=shape[1-axis]
    for start in range(0,len(group['indices']),1000000):
        a=np.asarray(group['indices'][start:start+1000000]);assert (a>=0).all() and (a<limit).all()
    return dict(storage=encoding,shape=shape,dtype=str(group['data'].dtype),nonzero_stored=int(ptr[-1]),empty_axis_entries=int((np.diff(ptr)==0).sum()))


def numerical(group):
    d=group if isinstance(group,h5py.Dataset) else group['data']
    minimum=np.inf;maximum=-np.inf;negative=zeros=noninteger=total=0;sum_values=0.
    if d.ndim==1:
        chunks=(np.asarray(d[start:start+1000000]) for start in range(0,len(d),1000000))
    else:chunks=(np.asarray(d[start:start+128,:]) for start in range(0,d.shape[0],128))
    for a in chunks:
        assert np.isfinite(a).all(),'Nonfinite matrix value'
        total+=a.size
        if a.size:minimum=min(minimum,float(a.min()));maximum=max(maximum,float(a.max()))
        negative+=int((a<0).sum());zeros+=int((a==0).sum());noninteger+=int((a!=np.floor(a)).sum());sum_values+=float(a.astype(float).sum())
    return dict(all_finite=True,stored_values=total,min=None if total==0 else minimum,max=None if total==0 else maximum,negative_stored_values=negative,zero_stored_values=zeros,noninteger_stored_values=noninteger,sum_stored_values=sum_values,
                integer_values_alone_do_not_prove_raw=True)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,512*1024**2)
    src=root/'GSE328481_public_processed_spatial_metadata_20261018_v2'
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    text=gzip.open(src/'public_series_family.soft.gz','rt').read();assert hashlib.sha256((src/'public_series_family.soft.gz').read_bytes()).hexdigest()=='1dd821744732f4d8b371873bd26082bfda3f75d033aacfcb1310b772a98cb418'
    # Independent regex reconstruction of every supplementary link (numbered or not).
    expected=[];gsm=None
    for line in text.splitlines():
        if line.startswith('^SAMPLE = '):gsm=line.split(' = ',1)[1]
        if re.match(r'^!Sample_supplementary_file(?:_\d+)? = ',line):
            url=line.split(' = ',1)[1]
            if url!='NONE':expected.append((gsm,re.sub(r'^ftp:','https:',url)))
    files=json.loads((src/'public_all_processed_file_inventory.json').read_text())
    assert expected==[(r['GSM'],r['url']) for r in files] and len(files)==11
    assert all(r['status']==200 and r['filename'].endswith('.h5ad') and int(r['content_length'])>0 for r in files)
    chosen=min(files,key=lambda r:(int(r['content_length']),r['url']));size=int(chosen['content_length'])
    assert size+8*1024**2<512*1024**2
    out=root/OUT;assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        metadata_terminal_sha256=sha(src/'private_terminal.json'),all11_links_independently_reparsed=True,
        selection_rule='Smallest exact Content-Length H5AD; ties URL lexical. Fixed before expression, schema pilot only, no favorable biology selection.',
        all_processed_file_bytes=sum(int(r['content_length']) for r in files),chosen=chosen,
        download='Exactly one official public processed H5AD, no controlled FASTQ or clonotype sequence values. Expected byte count/Last-Modified, SHA256 tracked, HDF5 signature and full numerical layer scan.',
        metadata_privacy='Visit dataset names/shapes only. Read gene axes and spatial numeric coordinates; never read obs cell values, clone strings, patient rows or clinical outcomes.',
        gate='No raw/normalization/cell-state/independent-patient assertion from a filename, integer values or one pilot. Full-source transfer qualification and finite spatial hypotheses remain separate.',
        new_effect_tests_or_annotation_changes=False)
    (out/'plan_before_pilot_expression_download_or_read.json').write_text(json.dumps(plan,indent=2)+'\n')
    target=out/chosen['filename'];digest=hashlib.sha256();read=0;checkpoint=0
    with urlopen(Request(chosen['url'],headers={'User-Agent':'PanCancer-public-processed-schema/1.0'}),timeout=45) as response, target.open('xb') as dest:
        assert response.status==200 and int(response.headers['Content-Length'])==size
        assert response.headers.get('Last-Modified')==chosen['last_modified']
        while True:
            block=response.read(1024**2)
            if not block:break
            dest.write(block);digest.update(block);read+=len(block);assert read<=size
            if read//(8*1024**2)>checkpoint:
                checkpoint=read//(8*1024**2)
                (out/'private_download_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),bytes_received=read,expected_bytes=size))+'\n')
    assert read==size and target.stat().st_size==size and sha(target)==digest.hexdigest()
    with target.open('rb') as f:assert f.read(8)==b'\x89HDF\r\n\x1a\n'
    structures=[];matrices={};coords={};genes={}
    with h5py.File(target,'r') as h:
        def visit(name,obj):
            if isinstance(obj,h5py.Dataset):structures.append(dict(path=name,shape=list(obj.shape),dtype=str(obj.dtype),compression=obj.compression))
        h.visititems(visit)
        paths=['X']+(['raw/X'] if 'raw/X' in h else [])
        if 'layers' in h:paths+=['layers/'+n for n in h['layers']]
        for name in paths:
            matrices[name]=dict(structure=matrix_structure(h[name]),numerical=numerical(h[name]))
        obs=h['obs'];var=h['var'];oi=string(obs.attrs.get('_index','_index'));vi=string(var.attrs.get('_index','_index'))
        nr=len(obs[oi]);nc=len(var[vi]);assert matrices['X']['structure']['shape']==[nr,nc]
        for name in sorted(set([vi])|{'gene_symbol','gene_symbols','gene_name','gene_names','SYMBOL','symbol'}):
            if name not in var or not isinstance(var[name],h5py.Dataset):continue
            ds=var[name]
            if ds.dtype.kind not in ['O','S','U']:continue
            values=np.array([string(v) for v in ds[:]])
            genes[name]=dict(n=len(values),unique=int(len(set(values))),empty=int((values=='').sum()),ensembl_like=int(np.char.startswith(values,'ENSG').sum()),sha256=hashlib.sha256(json.dumps(values.tolist(),ensure_ascii=False,separators=(',',':')).encode()).hexdigest())
        if 'obsm' in h:
            for name in ['spatial','X_spatial']:
                if name not in h['obsm']:continue
                ds=h['obsm'][name];assert isinstance(ds,h5py.Dataset)
                a=np.asarray(ds);assert len(a)==nr and a.ndim==2 and a.shape[1]>=2 and np.isfinite(a).all()
                coords[name]=dict(shape=list(a.shape),unique_xy=int(len(np.unique(a[:,:2],axis=0))),min=a.min(axis=0).tolist(),max=a.max(axis=0).tolist(),physical_unit='UNVERIFIED')
        summary=dict(observations=nr,variables=nc,obs_columns=list(obs.keys()),var_columns=list(var.keys()),layers=list(h['layers'].keys()) if 'layers' in h else [],uns_keys=list(h['uns'].keys()) if 'uns' in h else [],obsm_keys=list(h['obsm'].keys()) if 'obsm' in h else [],matrices=matrices,gene_axis_candidates=genes,coordinates=coords,raw_RNA_scale='NOT_CONFIRMED_BY_SCHEMA',single_cell_vs_bin='UNVERIFIED',cell_state_annotation_values_read=False,patient_identity_reviewed=False)
    (out/'private_full_HDF5_dataset_schema.json').write_text(json.dumps(structures,indent=2)+'\n')
    (out/'private_pilot_matrix_axes_coordinate_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    resource(root,8*1024**2)
    outputs=[p for p in out.iterdir() if p.is_file()]
    report=dict(status='COMPLETE_ONE_EFFECT_BLIND_PUBLIC_H5AD_SCHEMA_RAW_SCALE_AND_CELL_STATE_UNVERIFIED',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),download_sha256=digest.hexdigest(),download_bytes=size,
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in outputs},no_controlled_raw_clonotype_or_outcome_values_read=True,no_spatial_effect_or_full_source_qualification=True)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],completed_UTC=report['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))

if __name__=='__main__':main()
