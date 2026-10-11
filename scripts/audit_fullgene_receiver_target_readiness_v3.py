"""Independent complete source CSR/unit sums and literal prior coverage audit.

Does not import producer. Computational preparation is not biological validation.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd,h5py,pyreadr
from scipy.sparse import csr_matrix
def now():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def sha(p):
    d=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):d.update(b)
    return d.hexdigest()
def single(p):
    r=pyreadr.read_r(str(p));assert len(r)==1
    return next(iter(r.values()))
def table(p):return pd.read_csv(p,sep='\t',keep_default_na=False)

def main():
    start=now();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';O=P/'fullgene_receiver_target_readiness_v1';N=P/'supported_LR_native_cell_panel_v1'
    assert not (O/'independent_audit.json').exists();v=read(O/'inventory.json');plan=read(O/'predeclared_plan.json');assert sha(O/'predeclared_plan.json')==v['plan_sha256']
    for n,h in v['outputs'].items():assert sha(O/n)==h
    old=read(N/'inventory.json');assert sha(N/'inventory.json')==v['accepted_native_inventory_sha256']
    for item in plan['prior']['files']:
        q=O/'verified_prior_inputs_v2'/item['key'];assert q.stat().st_size==item['size'] and hashlib.md5(q.read_bytes()).hexdigest()==item['checksum'].split(':')[1]
    lr=single(O/'verified_prior_inputs_v2/lr_network_human_21122021.rds');tg=single(O/'verified_prior_inputs_v2/ligand_target_matrix_nsga2r_final.rds');tf=single(O/'verified_prior_inputs_v2/ligand_tf_matrix_nsga2r_final.rds')
    for m in [tg,tf]:m.index=m.index.astype(str);m.columns=m.columns.astype(str);assert m.index.is_unique and m.columns.is_unique and np.isfinite(m.to_numpy()).all() and (m.to_numpy()>=0).all()
    assert v['prior_dimensions']==dict(LR_rows=len(lr),target_genes=len(tg),target_ligands=len(tg.columns),ligand_tf_prior_row_symbols=len(tf),ligand_tf_prior_ligands=len(tf.columns))
    fixed=table(N/'fixed336_reference_supported_candidates.tsv');fixed['model_key']=fixed[['context_index','program','sender_family','receiver_family','ligand','receptor']].astype(str).agg('|'.join,axis=1)
    d=table(O/'all336_prior_coverage.tsv');f=table(O/'all15_followup_prior_coverage.tsv');sc=table(O/'all1008_source_target_coverage.tsv');ss=table(O/'fullgene_source_summary.tsv');original15=table(P/'axis_validation_feasibility_v1/automatic_next_validation_records.tsv')
    assert len(d)==336 and len(f)==15 and len(sc)==1008 and len(ss)==3 and set(d.model_key)==set(fixed.model_key) and set(f.model_key)==set(original15.model_key)
    assert d.model_key.is_unique and f.model_key.is_unique and not sc[['model_key','native_gse_id']].duplicated().any()
    cols=['model_key','gse_id','cancerType','program','sender_family','receiver_family','ligand','receptor']
    pd.testing.assert_frame_equal(d[cols].sort_values('model_key').reset_index(drop=True),fixed[cols].sort_values('model_key').reset_index(drop=True),check_dtype=False)
    assert set(zip(sc.model_key,sc.native_gse_id))=={(key,rec['gse_id']) for key in fixed.model_key for rec in old['source_records']}
    ref=d.set_index('model_key')
    for row in sc.itertuples():
        rr=ref.loc[row.model_key]
        assert row.discovery_gse_id==rr.gse_id
        for col in ['ligand','receptor','receiver_family','LR_status','exact_ligand_target_column']:assert getattr(row,col)==rr[col]
    priors={}
    for row in fixed.itertuples():
        key=(row.ligand,row.receptor)
        if key in priors:continue
        lg=row.ligand.split('_');rg=row.receptor.split('_');exact=lr[(lr['from']==row.ligand)&(lr['to']==row.receptor)]
        hit=[x+'>'+y for x in lg for y in rg if ((lr['from']==x)&(lr['to']==y)).any()]
        complex_pair=len(lg)>1 or len(rg)>1
        status='COMPLEX_NOT_ASSEMBLED_BY_FLAT_PRIOR' if complex_pair else 'EXACT_SINGLETON_LR_PRIOR' if len(exact)>0 else 'NO_EXACT_LR_PRIOR'
        positives=[]
        if row.ligand in tg:
            for gene,value in tg[row.ligand].items():
                if value>0:positives.append((gene,float(value)))
        ranked=[g for g,w in sorted(positives,key=lambda t:(-t[1],t[0]))];exclude=set(lg+rg+['CCL5','IL7R','MS4A1','FCER1G']);top=[x for x in ranked if x not in exclude][:200]
        priors[key]=dict(status=status,top=top,positives=len(positives),tf_count=int((tf[row.ligand]>0).sum()) if row.ligand in tf else 0,complex=complex_pair,component=';'.join(hit),databases=';'.join(sorted(set(exact.database.astype(str)))),sources=';'.join(sorted(set(exact.source.astype(str)))))
    for row in d.itertuples():
        p=priors[(row.ligand,row.receptor)];assert row.LR_status==p['status'] and row.top200_after_exclusions==';'.join(p['top']) and int(row.top200_count)==len(p['top']) and int(row.nonzero_target_genes)==p['positives'] and int(row.nonzero_ligand_tf_prior_row_symbols)==p['tf_count']
        assert row.exact_ligand_target_column==(row.ligand in tg) and row.exact_ligand_TF_column==(row.ligand in tf) and not row.receptor_specific_targets_available
        assert row.exact_singleton_LR_prior==(p['status']=='EXACT_SINGLETON_LR_PRIOR') and row.component_edges_do_not_prove_assembled_complex==p['complex'] and row.component_edges==p['component'] and row.prior_databases==p['databases'] and row.prior_sources==p['sources'] and row.prior_is_unsigned
    pd.testing.assert_frame_equal(f.sort_values('model_key').reset_index(drop=True),d[d.model_key.isin(original15.model_key)].sort_values('model_key').reset_index(drop=True),check_dtype=False)
    all_checks=[];total_cells=0;total_units=0;total_values=0;node_union=set()
    member=table(P/'membership_review_v2/source_node_analysis_membership.tsv')
    for rec in old['source_records']:
        sid=rec['gse_id'];cache=Path(rec['selected_cache']);assert sha(cache)==rec['cache_sha256'] and sha(rec['native_H5'])==rec['native_H5_sha256']
        z=dict(np.load(cache,allow_pickle=False));b=dict(np.load(O/'private_source_unit_arrays'/(sid+'.npz'),allow_pickle=False));cat=table(O/'private_source_unit_catalogs'/(sid+'.tsv'))
        node_union.update(z['subCluster'].tolist())
        assert b['sum_raw_counts'].dtype==np.dtype('int64') and b['positive_cell_counts'].dtype==np.dtype('int32') and b['mean_CP10K'].dtype==np.dtype('float64') and b['mean_log1p_CP10K'].dtype==np.dtype('float64')
        assert np.array_equal(z['samples'],b['samples']) and set(b['nodes'])==set(z['subCluster']) and len(cat)==len(b['n_cells'])
        eligible=member[(member.gse_id==sid)&member.analysis_family_main_membership]
        assert set(z['subCluster']).issubset(set(eligible.subCluster))
        with h5py.File(rec['native_H5'],'r') as h:
            genes=h['var'][h['var'].attrs['_index']].asstr()[:];assert np.array_equal(genes,b['genes'])
            x=csr_matrix((h['X/data'][:],h['X/indices'][:],h['X/indptr'][:]),shape=tuple(map(int,h['X'].attrs['shape'])))[z['source_row']].astype(np.float64)
        assert np.isfinite(x.data).all() and (x.data>=0).all() and np.array_equal(x.data,np.floor(x.data))
        full=np.asarray(x.sum(axis=1)).ravel();detected=np.asarray((x>0).sum(axis=1)).ravel();assert np.array_equal(full,z['full_RNA_total']) and np.array_equal(detected,z['native_detected_genes'])
        original_index=pd.Index(genes).get_indexer(z['genes']);assert (original_index>=0).all() and np.array_equal(x[:,original_index].toarray(),z['raw_counts'])
        cp=x.multiply((10000/full)[:,None]).tocsr();ln=cp.copy();ln.data=np.log1p(ln.data)
        per_node_family={n:str(z['families'][int(z['family_code'][np.flatnonzero(z['subCluster']==n)[0]])]) for n in b['nodes']}
        sample_ids=z['samples'][z['sample_code']]
        for row in cat.itertuples():
            i=int(row.unit_index);pos=np.flatnonzero((z['subCluster']==row.subCluster)&(sample_ids==row.sampleID));n=len(pos)
            assert (row.subCluster,row.sampleID)==(b['nodes'][i//len(b['samples'])],b['samples'][i%len(b['samples'])])
            assert n==int(row.n_cells)==b['n_cells'][i] and row.family==per_node_family[row.subCluster]
            assert np.array_equal(np.asarray(x[pos].sum(axis=0)).ravel().astype(np.int64),b['sum_raw_counts'][i])
            assert np.array_equal(np.asarray((x[pos]>0).sum(axis=0)).ravel().astype(np.int32),b['positive_cell_counts'][i])
            assert full[pos].sum()==b['full_RNA_total_sum'][i] and detected[pos].sum()==b['full_detected_gene_sum'][i]
            assert bool(row.background_diagnostics_eligible)==(n>=20)
            if n>=20:
                dc=np.asarray((x[pos]>0).sum(axis=0)).ravel()
                assert int(float(row.detected_genes1pct_n20))==int((dc>=.01*n).sum()) and int(float(row.detected_genes10pct_n20))==int((dc>=.1*n).sum())
            else:assert row.detected_genes1pct_n20==row.detected_genes10pct_n20==''
            if n:
                assert np.allclose(np.asarray(cp[pos].mean(axis=0)).ravel(),b['mean_CP10K'][i],rtol=1e-7,atol=1e-9)
                assert np.allclose(np.asarray(ln[pos].mean(axis=0)).ravel(),b['mean_log1p_CP10K'][i],rtol=1e-7,atol=1e-9)
            else:assert np.isnan(b['mean_CP10K'][i]).all() and np.isnan(b['mean_log1p_CP10K'][i]).all()
        fam={}
        for family in set(per_node_family.values()):
            mask=cat.family.eq(family).to_numpy();n=int(b['n_cells'][mask].sum());counts=b['positive_cell_counts'][mask].sum(axis=0)
            fam[family]=(n,set(genes[counts>=.01*n]) if n else set(),set(genes[counts>=.1*n]) if n else set())
        for row in sc[sc.native_gse_id==sid].itertuples():
            top=priors[(row.ligand,row.receptor)]['top'];n,e1,e10=fam.get(row.receiver_family,(0,set(),set()))
            has=row.ligand in tg;observed=n>0
            def count_or_na(value,expected):assert value=='NA' if expected is None else int(value)==expected
            assert bool(row.receiver_family_observed)==observed
            count_or_na(row.n_top200_targets,len(top));count_or_na(row.targets_in_native_axis,len(set(top)&set(genes)) if has else None)
            count_or_na(row.targets_receiver_detected1pct,len(set(top)&e1) if has and observed else None);count_or_na(row.targets_receiver_detected10pct,len(set(top)&e10) if has and observed else None)
            assert int(row.receiver_family_cells)==n and not row.independent_validation_source
            count_or_na(row.receiver_background_genes1pct,len(e1) if observed else None);count_or_na(row.receiver_background_genes10pct,len(e10) if observed else None)
            expected_status='MISSING_EXACT_LIGAND_TARGET_PRIOR' if not has else 'RECEIVER_FAMILY_NOT_CAPTURED' if not observed else 'NO_TOP_TARGET_IN_NATIVE_AXIS' if not set(top)&set(genes) else 'DESCRIPTIVE_COVERAGE_AVAILABLE'
            assert row.target_coverage_status==expected_status
        s=ss[ss.gse_id==sid].iloc[0];assert int(s.selected_cells)==len(z['source_row']) and int(s.native_genes)==len(genes) and int(s.libraries)==len(b['samples']) and int(s.original_nodes)==len(b['nodes']) and int(s.unit_rows)==len(cat)
        total_cells+=len(z['source_row']);total_units+=len(cat);total_values+=len(cat)*len(genes)*4
        all_checks.append(dict(gse_id=sid,selected_cells=len(z['source_row']),unit_rows=len(cat),native_genes=len(genes),aggregate_gene_values_checked=len(cat)*len(genes)*4,status='PASS_COMPLETE_NATIVE_SOURCE_UNIT_RAW_DETECTION_CP_LOG_AND_PRIOR_COVERAGE'))
        print(json.dumps(all_checks[-1]|{'at_UTC':now()}),flush=True);del x,cp,ln,b,z;gc.collect()
    assert total_cells==285360 and v['scope']['observed_original_node_union']==len(node_union) and v['scope']['original_catalogued_nodes']==13
    pd.DataFrame(all_checks).to_csv(O/'independent_all_source_checks.tsv',sep='\t',index=False)
    audit=dict(schema='fullgene_receiver_target_audit_v1',status='PASS_ALL336_15_1008_PRIOR_AND_FULLGENE_ORIGINAL_NODE_UNIT_AGGREGATES',started_UTC=start,finished_UTC=now(),selected_cells=total_cells,unit_rows=total_units,aggregate_gene_values_checked=total_values,candidate_records=336,followup_records=15,source_comparisons=1008,producer_inventory_sha256=sha(O/'inventory.json'),auditor_sha256=sha(__file__),method='No producer import; full source CSR sliced by accepted native rows and independently looped original node/library units versus producer block-incidence aggregation; exact unsigned literal prior edges and independently sorted target weights',statistical_models_fitted=False,formal_states=False,scientific_release=False,axis_validation_complete=False)
    (O/'independent_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit),flush=True)
if __name__=='__main__':main()
