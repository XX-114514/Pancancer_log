"""Official pinned Spectra GPU, prespecified multiseed identity-aware candidates.
Training is not state acceptance, communication validation or release.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc,random
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,torch
from Spectra.Spectra_gpu import est_spectra


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def main():
    started=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    plan_file=Path(os.environ['PANCANCER_SPECTRA_PLAN']);plan=json.loads(plan_file.read_text())
    assert sha(__file__)==plan['producer_sha256'];I=R/'state_singlecell_inputs_acceptance_20261014_v2'
    assert sha(I/'accepted_input_inventory.json')==plan['accepted_input_inventory_sha256']
    inv=json.loads((I/'accepted_input_inventory.json').read_text());assert inv['status']=='ACCEPTED_AUDITED_INPUTS_NOT_STATE_FIT'
    audit=json.loads((I/'independent_readback_audit.json').read_text());assert audit['all240533_raw_rows_directly_compared_to_canonical']
    for name in ['state_singlecell_readback_20261014_v2.execution.json','state_model_runtime_20261014_v1.execution.json']:
        job=json.loads((R/name).read_text());assert job['exit_code']==0 and job['status'].startswith('TERMINAL_SUCCESS')
        for key in ['pid','worker_pid','child_pid']:
            if key in job:assert not Path('/proc',str(job[key])).exists()
    runtime=R/'state_model_runtime_20261014_v1';resources=runtime/'resources';receipt=json.loads((resources/'runtime_resource_receipts.json').read_text())
    assert receipt['Spectra_repository_commit']==plan['Spectra_commit']
    for item in receipt['source_files']:assert sha(resources/'spectra_snapshot'/item['path'])==item['sha256']
    assert torch.cuda.is_available();torch.set_num_threads(2)
    path=Path(inv['training_file']['path']);assert sha(path)==inv['training_file']['sha256']
    raw=ad.read_h5ad(path);assert raw.n_obs==12000 and raw.n_vars==16625 and raw.obs_names.is_unique
    assert not raw.obs.gse_id.eq('GSE183904').any()
    log=raw.X.astype(np.float32).multiply((10000/raw.obs.full_RNA_total.to_numpy())[:,None]).tocsr();log.data=np.log1p(log.data)
    base=ad.AnnData(X=log,obs=raw.obs.copy(),var=raw.var.copy())
    reference=json.loads((resources/'spectra_snapshot/src/Spectra/data/default_gene_sets.json').read_text())
    dictionary={'global':reference['global'],'T cell':{scope+'__'+name:genes for scope in ['CD4_T','CD8_T','Treg'] for name,genes in reference[scope].items()},'Macrophage':reference['MDC']}
    common=set(base.var_names);dictionary={k:{name:sorted(set(gs)&common) for name,gs in value.items() if len(set(gs)&common)>=plan['min_reference_genes']} for k,value in dictionary.items()}
    hvg=set()
    for major in ['T cell','Macrophage']:
        subset=base[base.obs.major.eq(major)].copy()
        sc.pp.highly_variable_genes(subset,flavor='seurat',batch_key='gse_id',n_top_genes=plan['HVG_per_major'])
        hvg.update(subset.var_names[subset.var.highly_variable]);del subset
    prior=set(g for value in dictionary.values() for gs in value.values() for g in gs)
    features=sorted(hvg|prior);data=base[:,features].copy();data.X=data.X.toarray().astype(np.float32)
    L={k:len(value)+(plan['global_free_factors'] if k=='global' else plan['lineage_free_factors']) for k,value in dictionary.items()}
    assert {k:len(v) for k,v in dictionary.items()}==plan['expected_prior_counts']
    O=R/'Spectra_state_discovery_20261014_v1';assert not O.exists();O.mkdir();(O/'private_seed_models').mkdir()
    (O/'fixed_feature_dictionary.json').write_text(json.dumps(dict(features=features,HVG=sorted(hvg),reference=dictionary,L=L,full_RNA_denominator=True),indent=2)+'\n')
    summaries=[]
    for seed in plan['seeds']:
        random.seed(seed);np.random.seed(seed);torch.manual_seed(seed);torch.cuda.manual_seed_all(seed)
        a=data.copy();start=datetime.now(timezone.utc).isoformat();print(json.dumps(dict(seed=seed,started_UTC=start,features=len(features),L=L,training_cells=a.n_obs)),flush=True)
        model=est_spectra(a,dictionary,L=L,use_highly_variable=False,cell_type_key='major',filter_sets=False,n_top_vals=plan['top_genes'],lam=plan['lam'],delta=plan['delta'],rho=plan['rho'],num_epochs=plan['num_epochs'],batch_size=plan['batch_size'],skip=plan['lr_skip'],verbose=False)
        scores=np.asarray(a.obsm['SPECTRA_cell_scores']);factors=np.asarray(a.uns['SPECTRA_factors'])
        assert scores.shape==(12000,sum(L.values())) and factors.shape==(sum(L.values()),len(features))
        assert np.isfinite(scores).all() and np.isfinite(factors).all() and (scores>=0).all() and (factors>=0).all()
        mask=model.internal_model.alpha_mask.detach().cpu().numpy();assert mask.shape==scores.shape
        types=[]
        for j in range(scores.shape[1]):
            permitted=set(a.obs.major.to_numpy()[mask[:,j]>0]);types.append('global' if len(permitted)==2 else next(iter(permitted)))
        root=O/'private_seed_models'/('seed'+str(seed));root.mkdir()
        np.savez_compressed(root/'factors_and_discovery_scores.npz',factors=factors,cell_scores=scores,cell_ids=np.asarray(a.obs_names,dtype=str),genes=np.asarray(features,dtype=str),factor_scope=np.asarray(types,dtype=str))
        torch.save(model.internal_model.state_dict(),root/'model_state_dict.pt')
        pd.DataFrame(a.uns['SPECTRA_markers']).to_csv(root/'top_gene_markers.tsv',sep='\t',index=True)
        a.obs.to_csv(root/'private_discovery_metadata.tsv',sep='\t',index=True)
        summary=dict(seed=seed,started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),factors=scores.shape[1],features=len(features),score_finite=True,trained_candidate_not_state_release=True)
        summaries.append(summary);(O/'progress.json').write_text(json.dumps(dict(status='MULTISEED_TRAINING_PENDING_STABILITY',completed_seeds=len(summaries),summary=summaries),indent=2)+'\n');print(json.dumps(summary),flush=True)
        del model,a,scores,factors,mask;gc.collect();torch.cuda.empty_cache()
    inventory=dict(status='PRODUCED_MULTISEED_SPECTRA_CANDIDATES_PENDING_STABILITY_AND_PROJECTION',started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(plan_file),accepted_input_inventory_sha256=sha(I/'accepted_input_inventory.json'),Spectra_commit=plan['Spectra_commit'],seeds=summaries,training_cells=12000,features=len(features),L=L,reference_counts={k:len(v) for k,v in dictionary.items()},internal_heldout_source='GSE183904',cNMF_contrast_completed=False,UCell_projection_completed=False,all9CM_mapping_completed=False,DIALOGUE_completed=False,spatial_validation_completed=False,formal_state_labels=False,scientific_release=False,outputs={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()})
    (O/'inventory.json').write_text(json.dumps(inventory,indent=2)+'\n');print(json.dumps({k:v for k,v in inventory.items() if k!='outputs'}),flush=True)


if __name__=='__main__':main()
