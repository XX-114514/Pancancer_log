"""Full native gene units preserve original identities; prior coverage is not inference.

PANCANCER_PREFLIGHT_ROOT=<run> python fullgene_receiver_target_readiness_v3.py
Native cell/library arrays stay on the server. No new statistical p values.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np
import pandas as pd
import h5py,pyreadr
from scipy.sparse import csr_matrix

def now():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def rds(p):
    obj=pyreadr.read_r(str(p));assert len(obj)==1
    d=next(iter(obj.values()));assert isinstance(d,pd.DataFrame)
    return d

def prior_fields(row,lr,targets,tf):
    lg=row.ligand.split('_');rg=row.receptor.split('_');pairs=set(zip(lr['from'],lr['to']))
    exact=(row.ligand,row.receptor) in pairs
    component_hits=[a+'>'+b for a in lg for b in rg if (a,b) in pairs]
    complex_pair=len(lg)>1 or len(rg)>1
    status='EXACT_SINGLETON_LR_PRIOR' if exact and not complex_pair else 'COMPLEX_NOT_ASSEMBLED_BY_FLAT_PRIOR' if complex_pair else 'NO_EXACT_LR_PRIOR'
    has=row.ligand in targets.columns
    all_t=[]
    if has:
        v=targets[row.ligand]
        positive=v[v>0].sort_values(ascending=False,kind='stable')
        # Stable alphabetical ties, independent of current expression or effect.
        all_t=sorted(positive.index,key=lambda g:(-float(positive[g]),g))
    excluded=set(lg+rg+['CCL5','IL7R','MS4A1','FCER1G'])
    top=[g for g in all_t if g not in excluded][:200]
    tf_genes=sorted(tf.index[tf[row.ligand]>0]) if row.ligand in tf.columns else []
    origin=lr[(lr['from']==row.ligand)&(lr['to']==row.receptor)]
    return dict(ligand=row.ligand,receptor=row.receptor,LR_status=status,exact_singleton_LR_prior=bool(exact and not complex_pair),component_edges=';'.join(component_hits),component_edges_do_not_prove_assembled_complex=complex_pair,prior_databases=';'.join(sorted(set(origin.database.astype(str)))) if exact else '',prior_sources=';'.join(sorted(set(origin.source.astype(str)))) if exact else '',exact_ligand_target_column=has,nonzero_target_genes=len(all_t),top200_after_exclusions=';'.join(top),top200_count=len(top),exact_ligand_TF_column=row.ligand in tf.columns,nonzero_ligand_tf_prior_row_symbols=len(tf_genes),receptor_specific_targets_available=False,prior_is_unsigned=True)

def aggregate_source(rec,old,out):
    cache=Path(rec['selected_cache']);source=Path(rec['native_H5'])
    assert sha(cache)==rec['cache_sha256'] and sha(source)==rec['native_H5_sha256']
    z=dict(np.load(cache,allow_pickle=False));nodes=sorted(set(z['subCluster']),key=lambda v:int(v.rsplit('_',1)[1]));samples=z['samples'].astype(str)
    ni={v:j for j,v in enumerate(nodes)};nf={v:str(z['families'][int(z['family_code'][np.flatnonzero(z['subCluster']==v)[0]])]) for v in nodes}
    g=np.array([ni[v] for v in z['subCluster']],int)*len(samples)+z['sample_code'];groups=len(nodes)*len(samples)
    nc=np.bincount(g,minlength=groups).astype(np.int64);rna=np.zeros(groups,np.int64);ng=np.zeros(groups,np.int64)
    with h5py.File(source,'r') as h:
        genes=np.asarray(h['var'][h['var'].attrs['_index']].asstr()[:],dtype=str);assert len(set(genes))==len(genes)
        shape=tuple(map(int,h['X'].attrs['shape']));assert shape[0]==rec['native_source_cells'] and h['X'].attrs['encoding-type']=='csr_matrix'
        old_cols=pd.Index(genes).get_indexer(z['genes']);assert (old_cols>=0).all()
        rawsum=np.zeros((groups,len(genes)),np.int64);positive=np.zeros((groups,len(genes)),np.int32)
        cps=np.zeros((groups,len(genes)),np.float64);logs=np.zeros_like(cps)
        source_rows=z['source_row'];assert len(source_rows)==rec['selected_cells'] and (np.diff(source_rows)>0).all()
        lookup=np.full(shape[0],-1,np.int64);lookup[source_rows]=np.arange(len(source_rows));ptr=h['X/indptr'][:]
        for lo in range(0,shape[0],4096):
            hi=min(lo+4096,shape[0]);take=np.flatnonzero(lookup[lo:hi]>=0)
            if not len(take):continue
            zz=lookup[lo+take];p0,p1=int(ptr[lo]),int(ptr[hi]);data=h['X/data'][p0:p1]
            assert np.isfinite(data).all() and (data>=0).all() and (data==np.floor(data)).all()
            x=csr_matrix((data,h['X/indices'][p0:p1],ptr[lo:hi+1]-p0),shape=(hi-lo,len(genes)))[take].astype(np.int64)
            full=np.asarray(x.sum(axis=1)).ravel();det=np.asarray((x>0).sum(axis=1)).ravel()
            assert np.array_equal(full,z['full_RNA_total'][zz]) and np.array_equal(det,z['native_detected_genes'][zz]) and (full>0).all()
            assert np.array_equal(x[:,old_cols].toarray(),z['raw_counts'][zz])
            unique,inverse=np.unique(g[zz],return_inverse=True);inc=csr_matrix((np.ones(len(zz),np.int64),(inverse,np.arange(len(zz)))),shape=(len(unique),len(zz)))
            rawsum[unique]+=(inc@x).toarray();positive[unique]+=(inc@(x>0).astype(np.int32)).toarray()
            cp=x.astype(np.float64);cp.data*=np.repeat(10000/full,np.diff(cp.indptr));lp=cp.copy();lp.data=np.log1p(lp.data)
            cps[unique]+=(inc@cp).toarray();logs[unique]+=(inc@lp).toarray()
            np.add.at(rna,g[zz],full);np.add.at(ng,g[zz],det)
        assert np.array_equal(rawsum.sum(axis=1),rna)
    mean=np.divide(cps,nc[:,None],out=np.full_like(cps,np.nan),where=nc[:,None]>0);logmean=np.divide(logs,nc[:,None],out=np.full_like(logs,np.nan),where=nc[:,None]>0)
    cat=pd.DataFrame([dict(unit_index=j*len(samples)+k,subCluster=node,family=nf[node],sampleID=sample,n_cells=int(nc[j*len(samples)+k])) for j,node in enumerate(nodes) for k,sample in enumerate(samples)])
    cat['background_diagnostics_eligible']=nc>=20
    cat['detected_genes1pct_n20']=np.where(nc>=20,(positive>=.01*nc[:,None]).sum(axis=1),np.nan)
    cat['detected_genes10pct_n20']=np.where(nc>=20,(positive>=.1*nc[:,None]).sum(axis=1),np.nan)
    sid=rec['gse_id'];cat.to_csv(out/'private_source_unit_catalogs'/(sid+'.tsv'),sep='\t',index=False)
    np.savez_compressed(out/'private_source_unit_arrays'/(sid+'.npz'),genes=genes,nodes=np.asarray(nodes),node_family=np.asarray([nf[v] for v in nodes]),samples=samples,n_cells=nc,sum_raw_counts=rawsum,positive_cell_counts=positive,mean_CP10K=mean,mean_log1p_CP10K=logmean,full_RNA_total_sum=rna,full_detected_gene_sum=ng)
    summary=dict(gse_id=sid,source_cells=shape[0],selected_cells=int(nc.sum()),libraries=len(samples),original_nodes=len(nodes),unit_rows=groups,nonempty_units=int((nc>0).sum()),units_ge20=int((nc>=20).sum()),native_genes=len(genes),gene_axis_unique=True,native_input_sha256=rec['native_H5_sha256'],original_selected_cache_sha256=rec['cache_sha256'])
    family_expr={}
    for fam in sorted(set(nf.values())):
        mask=cat.family.eq(fam).to_numpy();n=int(nc[mask].sum());hits=positive[mask].sum(axis=0)
        family_expr[fam]=dict(n_cells=n,genes1=set(genes[hits>=.01*n]) if n else set(),genes10=set(genes[hits>=.1*n]) if n else set())
    del rawsum,positive,cps,logs,mean,logmean;gc.collect()
    return summary,set(genes),family_expr

def coverage_fields(row,sid,genes,families):
    proposed=row.top200_after_exclusions.split(';') if row.top200_after_exclusions else []
    expr=families.get(row.receiver_family);observed=bool(expr and expr['n_cells']>0);has=bool(row.exact_ligand_target_column)
    measured=set(proposed)&genes
    status='MISSING_EXACT_LIGAND_TARGET_PRIOR' if not has else 'RECEIVER_FAMILY_NOT_CAPTURED' if not observed else 'NO_TOP_TARGET_IN_NATIVE_AXIS' if not measured else 'DESCRIPTIVE_COVERAGE_AVAILABLE'
    return dict(model_key=row.model_key,native_gse_id=sid,discovery_gse_id=row.gse_id,ligand=row.ligand,receptor=row.receptor,receiver_family=row.receiver_family,LR_status=row.LR_status,exact_ligand_target_column=has,target_coverage_status=status,receiver_family_observed=observed,n_top200_targets=len(proposed),targets_in_native_axis=len(measured) if has else None,targets_receiver_detected1pct=len(set(proposed)&expr['genes1']) if observed and has else None,targets_receiver_detected10pct=len(set(proposed)&expr['genes10']) if observed and has else None,receiver_family_cells=expr['n_cells'] if observed else 0,receiver_background_genes1pct=len(expr['genes1']) if observed else None,receiver_background_genes10pct=len(expr['genes10']) if observed else None,independent_validation_source=False)

def main():
    start=now();root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);p=root/'immune_same_name_family_source_review_v1';out=p/'fullgene_receiver_target_readiness_v1';old=p/'supported_LR_native_cell_panel_v1';plan=read(out/'predeclared_plan.json')
    assert not (out/'inventory.json').exists() and read(root/'CM_program_screen_evidence_register_v1.json')['schema']=='native_LR_program_evidence_register_v16'
    prior_files={}
    for item in plan['prior']['files']:
        f=out/'verified_prior_inputs_v2'/item['key'];assert f.stat().st_size==item['size'] and hashlib.md5(f.read_bytes()).hexdigest()==item['checksum'].split(':')[1]
        prior_files[item['key']]=dict(sha256=sha(f),official_md5=item['checksum'],bytes=item['size'])
    lr=rds(out/'verified_prior_inputs_v2/lr_network_human_21122021.rds');target=rds(out/'verified_prior_inputs_v2/ligand_target_matrix_nsga2r_final.rds');tf=rds(out/'verified_prior_inputs_v2/ligand_tf_matrix_nsga2r_final.rds')
    assert set(lr.columns)=={'from','to','database','source'}
    for matrix in [target,tf]:
        matrix.index=matrix.index.astype(str);matrix.columns=matrix.columns.astype(str)
        assert matrix.index.is_unique and matrix.columns.is_unique and np.isfinite(matrix.to_numpy()).all() and (matrix.to_numpy()>=0).all()
    d=pd.read_csv(old/'fixed336_reference_supported_candidates.tsv',sep='\t');d['model_key']=d[['context_index','program','sender_family','receiver_family','ligand','receptor']].astype(str).agg('|'.join,axis=1);assert len(d)==336 and d.model_key.is_unique
    pairs=d[['ligand','receptor']].drop_duplicates().sort_values(['ligand','receptor']);pr=pd.DataFrame([prior_fields(r,lr,target,tf) for r in pairs.itertuples()]);assert len(pr)==209
    public_cols=['model_key','gse_id','cancerType','program','sender_family','receiver_family','ligand','receptor']
    all336=d[public_cols].merge(pr,on=['ligand','receptor'],validate='many_to_one');all336.to_csv(out/'all336_prior_coverage.tsv',sep='\t',index=False,na_rep='NA')
    f=pd.read_csv(p/'axis_validation_feasibility_v1/automatic_next_validation_records.tsv',sep='\t');all15=f[['model_key']].merge(all336,on='model_key',validate='one_to_one');assert len(all15)==15;all15.to_csv(out/'all15_followup_prior_coverage.tsv',sep='\t',index=False,na_rep='NA')
    inv=read(old/'inventory.json');assert len(inv['source_records'])==3
    for folder in ['private_source_unit_arrays','private_source_unit_catalogs']: (out/folder).mkdir(exist_ok=False)
    summaries=[];comparisons=[];source_records=[];observed_original_nodes=set()
    for rec in inv['source_records']:
        sm,genes,families=aggregate_source(rec,old,out);summaries.append(sm);sid=rec['gse_id']
        with np.load(rec['selected_cache'],allow_pickle=False) as accepted:observed_original_nodes.update(accepted['subCluster'].tolist())
        for row in all336.itertuples():comparisons.append(coverage_fields(row,sid,genes,families))
        source_records.append(dict(gse_id=sid,source=str(rec['native_H5']),source_sha256=rec['native_H5_sha256'],accepted_selected_cache=str(rec['selected_cache']),cache_sha256=rec['cache_sha256']))
        write(out/'progress.json',dict(status='PREPARING_FULL_NATIVE_ORIGINAL_NODE_UNITS',completed_sources=len(summaries),selected_cells=sum(s['selected_cells'] for s in summaries)));print(json.dumps(sm),flush=True)
    assert sum(s['selected_cells'] for s in summaries)==285360 and sum(s['libraries'] for s in summaries)==164 and len(comparisons)==1008
    coverage=pd.DataFrame(comparisons)
    for col in ['n_top200_targets','targets_in_native_axis','targets_receiver_detected1pct','targets_receiver_detected10pct','receiver_family_cells','receiver_background_genes1pct','receiver_background_genes10pct']:coverage[col]=coverage[col].astype('Int64')
    coverage.to_csv(out/'all1008_source_target_coverage.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(summaries).to_csv(out/'fullgene_source_summary.tsv',sep='\t',index=False)
    outputs={str(f.relative_to(out)):sha(f) for f in out.rglob('*') if f.is_file()}
    inv=dict(schema='fullgene_receiver_target_inventory_v1',status='PRODUCER_COMPLETE_PENDING_FULLGENE_AND_PRIOR_AUDIT',started_UTC=start,finished_UTC=now(),producer_sha256=sha(__file__),plan_sha256=sha(out/'predeclared_plan.json'),accepted_native_inventory_sha256=sha(old/'inventory.json'),prior_files=prior_files,prior_dimensions=dict(LR_rows=len(lr),target_genes=len(target),target_ligands=len(target.columns),ligand_tf_prior_row_symbols=len(tf),ligand_tf_prior_ligands=len(tf.columns)),scope=dict(records=336,unique_LR=209,followup_records=15,source_comparisons=1008,selected_cells=285360,libraries=164,source_cells=809784,original_catalogued_nodes=13,observed_original_node_union=len(observed_original_nodes)),source_records=source_records,outputs=outputs,ligand_tf_prior_is_not_a_curated_TF_inventory=True,formal_states=False,statistical_models_fitted=False,scientific_release=False,axis_validation_complete=False)
    write(out/'inventory.json',inv);print(json.dumps({k:v for k,v in inv.items() if k not in ['outputs','source_records','prior_files']}),flush=True)

if __name__=='__main__':main()
