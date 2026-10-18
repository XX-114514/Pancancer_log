"""All/pilot source qualification: declared cell-bin metadata and fixed RNA capture.

Public file only in RAM; no clones, cell IDs, patients, spatial effects or refits.
"""
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from GEO_public_processed_RAM_v1 import read_public
from H5AD_public_schema_math_v1 import string,matrix_structure,numerical

FIELDS=['anno','anno2','bin50 structural annotation']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def annotation(obs,name,n):
    if name not in obs:return None
    d=obs[name]
    if isinstance(d,h5py.Group):
        assert string(d.attrs.get('encoding-type',''))=='categorical'
        cats=[string(s) for s in d['categories'][:]];codes=np.asarray(d['codes'],int)
    else:
        values=np.array([string(s) for s in d[:]]);cats,codes=np.unique(values,return_inverse=True);cats=cats.tolist()
    assert len(codes)==n and (codes>=-1).all() and (codes<len(cats)).all()
    assert len(cats)<=256 and all(len(s)<256 for s in cats),'Unexpected high-cardinality annotation; do not infer identity'
    return cats+['UNANNOTATED'],np.where(codes>=0,codes,len(cats))


def scalar(h,path):
    if path not in h:return dict(status='MISSING')
    d=h[path]
    if not isinstance(d,h5py.Dataset) or d.size!=1:return dict(status='NONSCALAR_NOT_READ')
    x=np.asarray(d)[()]
    if isinstance(x,np.ndarray):x=x.ravel()[0]
    if isinstance(x,(bytes,str,np.bytes_)):x=string(x)
    elif isinstance(x,np.generic):x=x.item()
    return dict(status='AUTHOR_SCALAR_READ',value=x)


def measure(h,panels,axes):
    x=h['X'];s=matrix_structure(x);assert s['storage']=='csr_matrix'
    n,nf=s['shape'];obs=h['obs'];var=h['var'];index=string(var.attrs.get('_index','_index'))
    names=[string(v) for v in var[index][:]];assert len(names)==nf
    unique=Counter(names);lookup={v:i for i,v in enumerate(names) if unique[v]==1}
    requested=sorted(set(g for gs in panels.values() for g in gs));gi={v:i for i,v in enumerate(requested)}
    target=np.full(nf,-1,int)
    for g,j in gi.items():
        if g in lookup:target[lookup[g]]=j
    columns={p:np.array([gi[g] for g in gs],int) for p,gs in panels.items()}
    complete={p:all(g in lookup for g in gs) for p,gs in panels.items()}
    anno={f:annotation(obs,f,n) for f in FIELDS};anno={f:v for f,v in anno.items() if v is not None}
    total=np.zeros(n);detected=np.zeros(n,int)
    positive1={p:np.zeros(n,bool) for p in panels};positive2={p:np.zeros(n,bool) for p in panels}
    positive_strict2={p:np.zeros(n,bool) for p in panels}
    raw_identical=None
    if 'raw/X' in h:
        r=h['raw/X'];raw_identical=True
        for key in ['data','indices','indptr']:
            if x[key].shape!=r[key].shape:raw_identical=False;continue
            for start in range(0,len(x[key]),1000000):
                if not np.array_equal(x[key][start:start+1000000],r[key][start:start+1000000]):raw_identical=False
    for start in range(0,n,4096):
        stop=min(start+4096,n);ptr=np.asarray(x['indptr'][start:stop+1],np.int64);lo,hi=int(ptr[0]),int(ptr[-1]);ptr-=lo
        val=np.asarray(x['data'][lo:hi]);idx=np.asarray(x['indices'][lo:hi],int)
        assert np.isfinite(val).all() and (val>0).all() and (val==np.floor(val)).all(),'Deposited count-like positive integer assay required; normalized/negative data need separate qualification'
        row=np.repeat(np.arange(stop-start),np.diff(ptr))
        # Do not silently count duplicate stored coordinates as separate genes.
        if len(idx)>1:assert ((idx[1:]>idx[:-1])|(row[1:]!=row[:-1])).all(),'Unsorted/duplicate source coordinates require explicit canonicalization'
        prefix=np.r_[0.,np.cumsum(val,dtype=np.float64)]
        total[start:stop]=prefix[ptr[1:]]-prefix[ptr[:-1]];detected[start:stop]=np.diff(ptr)
        values=np.zeros((stop-start,len(requested)))
        take=target[idx]>=0;values[row[take],target[idx[take]]]=val[take]
        for p,js in columns.items():
            if not complete[p]:continue
            a=values[:,js]
            positive1[p][start:stop]=(a>0).sum(axis=1)>=1
            positive2[p][start:stop]=(a>0).sum(axis=1)>=2
            positive_strict2[p][start:stop]=(a>=2).sum(axis=1)>=2
    summaries=[];pairs=[]
    # All-observation group is a measurement denominator, not a cell identity.
    anno['ALL_OBSERVATIONS']=(['ALL'],np.zeros(n,int))
    for field,(labels,codes) in anno.items():
        for j,label in enumerate(labels):
            use=codes==j;nn=int(use.sum())
            if not nn:continue
            base=dict(annotation_field=field,label=label,n=nn,median_deposited_total=float(np.median(total[use])),median_detected_genes=float(np.median(detected[use])),zero_deposited_total=int((total[use]==0).sum()))
            for p in panels:summaries.append(dict(**base,panel=p,required_genes=len(panels[p]),unique_measured_genes=sum(g in lookup for g in panels[p]),axis_complete=complete[p],cells_any_gene=int(positive1[p][use].sum()) if complete[p] else None,cells_ge2_genes=int(positive2[p][use].sum()) if complete[p] else None,cells_ge2_genes_each_value_ge2=int(positive_strict2[p][use].sum()) if complete[p] else None))
            for axis,(a,b) in axes.items():pairs.append(dict(annotation_field=field,label=label,n=nn,axis=axis,both_panels_complete=bool(complete[a] and complete[b]),any_gene_both_panels=int((positive1[a][use]&positive1[b][use]).sum()) if complete[a] and complete[b] else None,ge2_genes_both_panels=int((positive2[a][use]&positive2[b][use]).sum()) if complete[a] and complete[b] else None))
    return dict(observations=n,genes=nf,author_scalars={p:scalar(h,'uns/'+p) for p in ['bin_type','bin_size','resolution']},cell_border_shape=list(h['obsm']['cell_border'].shape) if 'obsm/cell_border' in h else None,X_raw_X_exactly_identical=raw_identical,annotation_fields_read=list(anno),fixed_panel_capture=summaries,fixed_axis_same_cell_detection=pairs,physical_unit_and_patient_identity='NOT_INDEPENDENTLY_CONFIRMED',deposited_integer_assay_is_not_automatically_verified_original_UMI=True)


def main(args):
    out=args.output;assert not out.exists();out.mkdir(mode=0o700)
    inventory=json.loads(args.inventory.read_text());assert len(inventory)==11
    fixed=json.loads(args.panel_plan.read_text());todo=inventory if args.scope=='all' else [min(inventory,key=lambda r:(r['bytes'],r['url']))]
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),inventory_sha256=sha(args.inventory),fixed_panel_plan_sha256=sha(args.panel_plan),scope=args.scope,sources=[r['GSM'] for r in todo],
        declared_metadata_fields=FIELDS,author_scalar_fields=['bin_type','bin_size','resolution'],
        measurements='All original7 panels and original2 axes, gene-axis completeness, within-cell panel detection>=1/>=2 genes and >=2genes each deposited value>=2; per-annotation total/detected RNA medians. Complete X/raw equality and positive sorted integer assay gate.',
        privacy='Public matrices in RAM only; no existing project matrices exported or raw cell/clone/patient IDs read. Only low-cardinality annotation fields above; no local per-cell arrays persisted.',
        limits='Count-like deposited integer values do not prove original UMI. No new state annotations, spatial associations, P/CI, clinical refits or independent patient claims. All sources/statuses retained; no expression-based source selection.')
    (out/'plan_before_declared_cellbin_capture_read.json').write_text(json.dumps(plan,indent=2)+'\n')
    statuses=[]
    for chosen in todo:
        name=chosen['GSM'];print(json.dumps(dict(stage='source_started',GSM=name,expected_bytes=chosen['bytes'])),flush=True)
        buffer,digest=read_public(chosen,args.proxy)
        if args.scope=='pilot':assert digest=='fbf4c6f5955851e7dad41a39e70729a4c859be61fd3e09e05f79dc395050711f'
        try:
            with h5py.File(buffer,'r') as h:result=measure(h,fixed['panels'],fixed['axes'])
            result.update(status='COMPLETE_DECLARED_METADATA_AND_FIXED_PANEL_CAPTURE_NO_EFFECTS',GSM=name,public_file_sha256=digest)
        except (ValueError,AssertionError,KeyError) as e:result=dict(status='UNQUALIFIED_REQUIRES_REVIEW',GSM=name,public_file_sha256=digest,reason=type(e).__name__+': '+str(e))
        buffer.close();p=out/(name+'_private_qualification.json');p.write_text(json.dumps(result,indent=2)+'\n');statuses.append(dict(GSM=name,status=result['status'],sha256=sha(p)))
        print(json.dumps(dict(stage='source_finished',GSM=name,status=result['status'])),flush=True)
    files=[p for p in out.iterdir() if p.is_file()];assert sum(p.stat().st_size for p in files)<8*1024**2
    terminal=dict(status='COMPLETE_DECLARED_CELLBIN_CAPTURE_SCOPE_NO_SPATIAL_EFFECTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),scope=args.scope,source_statuses=statuses,local_matrix_or_cell_arrays_persisted=False,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},patient_or_raw_scale_or_cell_state_validation=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(dict(status=terminal['status'],completed_UTC=terminal['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inventory',type=Path,required=True);p.add_argument('--panel-plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--scope',choices=['pilot','all'],required=True);p.add_argument('--proxy',default='');main(p.parse_args())
