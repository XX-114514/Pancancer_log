"""Read all saved cells independently; accept data retained after count-gate bug.
Never overwrites v1, changes labels or imports its producer.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd,anndata as ad,zarr
from scipy.sparse import load_npz,csr_matrix


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def main():
    start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);A=Path(os.environ['PANCANCER_ANNOTATION_ROOT'])
    O=R/'state_singlecell_inputs_acceptance_20261014_v2';D=R/'state_singlecell_inputs_20261014_v1';C=R/'state_singlecell_inputs_20261014_v1_code'
    plan=json.loads((C/'state_input_plan_v1.json').read_text());v2=json.loads(Path(os.environ['PANCANCER_STATE_ACCEPT_PLAN']).read_text())
    assert sha(__file__)==v2['auditor_sha256'] and not O.exists()
    e=json.loads((R/'state_singlecell_inputs_20261014_v1.execution.json').read_text())
    assert e['status']=='TERMINAL_FAILED_INPUT_ASSEMBLY' and e['exit_code']==1
    for k in ['worker_pid','child_pid']:assert not Path('/proc',str(e[k])).exists()
    log=(R/'state_singlecell_inputs_20261014_v1.stdout.log').read_text()
    assert "sum(s['libraries'] for s in summaries)==424" in log and 'AssertionError' in log
    assert not (D/'inventory.json').exists()
    U=A/'20261010_uniform_expression_annotation_v1';N=R/'immune_same_name_family_source_review_v1/mac_T_IFNG_receptor_single_cell_v1'
    for n,h in plan['frozen_input_sha256'].items():assert sha(N/n)==h
    expected=pd.read_csv(N/'private_all1832_sample_node_RNA_counts.tsv',sep='\t').set_index('unit_index')
    z=np.load(N/'private_selected_cell_gene_bits.npz',allow_pickle=False)
    core_contract=json.loads((U/'input_contract.json').read_text());assert sha(U/'input_contract.json')==plan['input_contract_sha256']
    core=zarr.open_group(core_contract['core_zarr'],mode='r');core_index=core['obs'][core['obs'].attrs['_index']]
    block_map={};review_map={}
    for cid,digest in plan['context_sha256'].items():
        file=R/'immune_same_name_family_source_review_v1/major_receiver_TF_preflight_v3/private_context_arrays'/(cid+'.npz');assert sha(file)==digest
        q=np.load(file,allow_pickle=False)
        for sample,block,review in zip(q['samples'],q['blocks'],q['identity_reviewed']):
            sample=str(sample);block_map[sample]=str(block);review_map[sample]=bool(review)
    original=pd.read_csv(U/'all_node_number_catalog.tsv',sep='\t');train=ad.read_h5ad(D/'private_training_common_counts.h5ad')
    assert train.n_vars==16625 and train.n_obs==12000 and train.obs_names.is_unique
    assert not train.obs.gse_id.eq(plan['heldout_source']).any()
    arrays=[];samples=set();blocks=set();totals=[];checked=0;source_files={};source_genes=[]
    for sid in plan['sources']:
        p=D/'private_native_sources'/(sid+'.h5ad');a=ad.read_h5ad(p);X=a.X.tocsr();meta=a.obs
        gt=pd.read_csv(U/'by_gse'/sid/'genes.tsv',sep='\t');assert np.array_equal(a.var_names,gt.gene)
        assert a.obs_names.is_unique and np.isfinite(X.data).all() and (X.data>=0).all() and np.array_equal(X.data,np.floor(X.data))
        positions=meta.core_position.to_numpy();selected=z['global_rows'][expected.loc[z['unit_index'],'gse_id'].to_numpy()==sid]
        assert np.array_equal(positions,selected) and (np.diff(positions)>0).all()
        assert np.array_equal(a.obs_names,core_index.get_orthogonal_selection((positions,)).astype(str))
        prior_indices=z['unit_index'][np.isin(z['global_rows'],positions)];assert np.array_equal(meta.unit_index,prior_indices)
        assert np.array_equal(meta.analysis_block.astype(str),meta.sample_uid.map(block_map).astype(str))
        assert np.array_equal(meta.patient_identity_reviewed,meta.sample_uid.map(review_map))
        for field in ['sample_uid','node_id','subCluster','major','gse_id']:
            assert np.array_equal(meta[field].astype(str),expected.loc[prior_indices,field].astype(str))
        for begin in range(0,a.n_obs,512):
            end=min(begin+512,a.n_obs);rr=positions[begin:end];lo=int(rr[0]);hi=int(rr[-1])+1
            pointer=core['X/indptr'][lo:hi+1].astype(np.int64);left,right=int(pointer[0]),int(pointer[-1])
            fresh=csr_matrix((core['X/data'][left:right],core['X/indices'][left:right],pointer-left),shape=(hi-lo,core_contract['raw_shape'][1]))
            fresh=fresh[rr-lo][:,gt.core_column.to_numpy()]
            diff=fresh-X[begin:end];diff.eliminate_zeros();assert diff.nnz==0
        full=np.asarray(X.sum(axis=1)).ravel();det=np.asarray((X>0).sum(axis=1)).ravel()
        assert np.array_equal(full,meta.full_RNA_total) and np.array_equal(det,meta.full_detected_genes) and (full>0).all()
        cp=X.astype(np.float64).multiply((10000/full)[:,None]).tocsr();cp.data=np.log1p(cp.data)
        native=pd.read_csv(U/'by_gse'/sid/'units.tsv',sep='\t');native['node_id']=original.node_id.to_numpy()[native.node_index]
        native=native[native.candidate_eligible].set_index(['sample_uid','node_id']);assert native.index.is_unique
        old_positive=load_npz(U/'by_gse'/sid/'positive_counts.npz').tocsr();old_log=load_npz(U/'by_gse'/sid/'sum_log1p_CP10K.npz').tocsr()
        for unit,group in meta.reset_index().groupby('unit_index',sort=True):
            row=expected.loc[unit];ii=group.index.to_numpy();assert len(ii)==row.n_cells
            old=int(native.loc[(row.sample_uid,row.node_id)].unit_row)
            assert np.array_equal(np.asarray((X[ii]>0).sum(axis=0)).ravel(),old_positive[old].toarray().ravel())
            assert np.allclose(np.asarray(cp[ii].sum(axis=0)).ravel(),old_log[old].toarray().ravel(),rtol=1e-9,atol=1e-8)
            checked+=2*a.n_vars
        chosen=meta.state_training_selected.to_numpy(bool);tt=train.obs.gse_id.eq(sid).to_numpy()
        assert np.array_equal(a.obs_names[chosen],train.obs_names[tt])
        if chosen.any():
            assert not sid==plan['heldout_source']
            diff=X[chosen][:,pd.Index(a.var_names).get_indexer(train.var_names)]-train.X[tt]
            assert diff.nnz==0
            assert np.array_equal(meta.full_RNA_total.to_numpy()[chosen],train.obs.full_RNA_total.to_numpy()[tt])
            for major in plan['majors']:assert int((chosen&meta.major.eq(major).to_numpy()).sum())<=1500
        elif sid!=plan['heldout_source']:raise AssertionError('Discovery source unexpectedly has no selected cells')
        samples.update(meta.sample_uid);blocks.update(meta.analysis_block);arrays.extend(positions.tolist());source_genes.append(set(a.var_names))
        source_files[sid]=dict(path=str(p),sha256=sha(p))
        totals.append(dict(gse_id=sid,cells=a.n_obs,genes=a.n_vars,libraries=meta.sample_uid.nunique(),sample_major_units=len(meta[['sample_uid','major']].drop_duplicates()),analysis_blocks=meta.analysis_block.nunique(),training_cells=int(chosen.sum()),status='PASS_ALL_CELLS_NATIVE_GENES_RNA_KEYS_AND_UNIT_LOOP_AGGREGATES'))
        print(json.dumps(totals[-1]),flush=True);del a,X,cp,old_positive,old_log;gc.collect()
    assert len(samples)==212 and sum(t['sample_major_units'] for t in totals)==424
    assert len(arrays)==240533 and np.array_equal(np.sort(arrays),z['global_rows'])
    assert set.intersection(*source_genes)==set(train.var_names)
    O.mkdir();pd.DataFrame(totals).to_csv(O/'all_source_readback_checks.tsv',sep='\t',index=False)
    audit=dict(status='PASS_ALL240533_CELL_NATIVE_COUNTS_FULLRNA_AND1832_UNIT_AGGREGATES',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),parent_execution_status=e['status'],parent_failure='Library count212 mistaken for sample-major units424 in final summary assert; no data transform changed',cells=240533,libraries=212,sample_major_units=424,sample_node_units=1832,analysis_blocks=len(blocks),training_cells=12000,common_genes=16625,all240533_raw_rows_directly_compared_to_canonical=True,aggregate_values_checked=checked,source_files=source_files,training_file=dict(path=str(D/'private_training_common_counts.h5ad'),sha256=sha(D/'private_training_common_counts.h5ad')),all_source_checks=totals,state_model_trained=False,state_labels_assigned=False,scientific_release=False)
    (O/'independent_readback_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    (O/'accepted_input_inventory.json').write_text(json.dumps(dict(schema='accepted_state_singlecell_input_inventory_v2',status='ACCEPTED_AUDITED_INPUTS_NOT_STATE_FIT',audit_sha256=sha(O/'independent_readback_audit.json'),verified_dataset_files=source_files,training_file=audit['training_file'],cells=240533,libraries=212,sample_major_units=424,training_cells=12000,common_genes=16625,heldout_source=plan['heldout_source'],heldout_is_internal_only=True,state_labels_assigned=False,scientific_release=False),indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k not in ['source_files','training_file','all_source_checks']}),flush=True)


if __name__=='__main__':main()
