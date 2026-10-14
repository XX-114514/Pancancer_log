"""Official pyUCell on complete native axes, fixed candidates, no type changes."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc,importlib.metadata,subprocess,re
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import spearmanr
import pyucell as uc

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def rho(x,y):return float(spearmanr(x,y).statistic) if np.std(x)>0 and np.std(y)>0 else np.nan
def main():
    started=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=Path(os.environ['PANCANCER_UCELL_PROJECTION_PLAN']);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256'];assert importlib.metadata.version('pyucell')=='0.8.0'
    D=R/'state_programme_dictionary_20261014_v1';di=json.loads((D/'inventory.json').read_text());assert sha(D/'frozen_programme_dictionary.json')==di['dictionary_sha256'];dictionary=json.loads((D/'frozen_programme_dictionary.json').read_text());records=dictionary['programmes'];ids=[r['programme_id'] for r in records];assert len(set(ids))==98
    I=R/'state_singlecell_inputs_acceptance_20261014_v2/accepted_input_inventory.json';assert sha(I)==plan['input_inventory_sha256'];inv=json.loads(I.read_text());O=R/'state_ucell_native_projection_20261014_v1';assert not O.exists();O.mkdir();(O/'private_source_scores').mkdir();coverage=[];stats=[];node_stats=[];files=[];total=0
    for sid,item in inv['verified_dataset_files'].items():
        quota=subprocess.check_output(['/usr/bin/quota','-g'],text=True)
        q=re.findall(r'^\s*(\d+)\s+0\s+(\d+)\s+\d+\s+0\s+0\s*$',quota,re.M);assert len(q)==1
        free_GiB=(int(q[0][1])-int(q[0][0]))/1048576
        assert free_GiB>=plan['reserved_free_GiB']+plan['projection_disk_budget_GiB'],('reserved_capacity_gate',free_GiB)
        p=Path(item['path']);assert sha(p)==item['sha256'];a=ad.read_h5ad(p);assert a.var_names.is_unique and a.obs_names.is_unique;assert np.isfinite(a.X.data).all() and (a.X.data>=0).all();assert np.array_equal(np.asarray(a.X.sum(axis=1)).ravel(),a.obs.full_RNA_total);total+=a.n_obs;native_obs=a.obs.copy();primary={};tpm_sens={};det=np.zeros((a.n_obs,98),np.uint8);fraction=np.zeros((a.n_obs,98),np.float32);mask=np.zeros((a.n_obs,98),bool)
        for j,r in enumerate(records):
            ix=a.var_names.get_indexer(r['primary_genes']);present=ix>=0;complete=present.mean()>=plan['min_measured_fraction'] and present.sum()>=plan['min_measured_genes'];allowed=a.obs.major.isin(r['allowed_majors']).to_numpy();mask[:,j]=allowed&complete;x=a.X[:,ix[present]].toarray();det[:,j]=(x>0).sum(axis=1);fraction[:,j]=x.sum(axis=1)/a.obs.full_RNA_total.to_numpy();coverage.append(dict(gse_id=sid,programme_id=r['programme_id'],required_genes=50,measured_genes=int(present.sum()),scorable_gene_coverage=bool(complete),missing_genes=';'.join(np.array(r['primary_genes'])[~present])))
            if complete:primary[r['programme_id']]=np.array(r['primary_genes'])[present].tolist()
            if r['method']=='cNMF':
                ix2=a.var_names.get_indexer(r['TPM_sensitivity_genes']);present2=ix2>=0
                if present2.mean()>=plan['min_measured_fraction'] and present2.sum()>=plan['min_measured_genes']:tpm_sens[r['programme_id']]=np.array(r['TPM_sensitivity_genes'])[present2].tolist()
        matrices={}
        for maxrank in plan['max_rank_values']:
            uc.compute_ucell_scores(a,signatures=primary,max_rank=maxrank,ties_method='average',missing_genes='skip',chunk_size=500,n_jobs=1,device='cpu',suffix='_UCell')
            y=np.full((a.n_obs,98),np.nan,np.float32)
            for j,pid in enumerate(ids):
                if pid in primary:y[:,j]=a.obs[pid+'_UCell'].to_numpy();y[~mask[:,j],j]=np.nan
            assert np.all((y[np.isfinite(y)]>=-1e-6)&(y[np.isfinite(y)]<=1+1e-6));matrices['scores_maxrank'+str(maxrank)]=y
            print(json.dumps(dict(gse_id=sid,cells=a.n_obs,maxrank=maxrank,phase='NATIVE_AXIS_PRIMARY_SIGNATURES_SCORED',UTC=datetime.now(timezone.utc).isoformat())),flush=True)
        uc.compute_ucell_scores(a,signatures=tpm_sens,max_rank=1500,ties_method='average',missing_genes='skip',chunk_size=500,n_jobs=1,device='cpu',suffix='_TPM_UCell');t=np.full((a.n_obs,98),np.nan,np.float32)
        for j,pid in enumerate(ids):
            if pid in tpm_sens:t[:,j]=a.obs[pid+'_TPM_UCell'].to_numpy();t[~mask[:,j],j]=np.nan
        matrices['TPM_signature_scores_maxrank1500']=t;root=O/'private_source_scores'/sid;root.mkdir();np.savez_compressed(root/'cell_programme_scores.npz',cell_ids=a.obs_names.to_numpy(str),programme_ids=np.array(ids),detected_programme_genes=det,programme_RNA_fraction=fraction,model_scope_and_coverage_mask=mask,**matrices);native_obs.to_csv(root/'native_cell_identity_metadata.tsv',sep='\t',index=True)
        # Frozen raw identities are preserved, all phenotype annotations are sidecars.
        assert native_obs[['major','subCluster','node_id']].equals(a.obs[['major','subCluster','node_id']]);y=matrices['scores_maxrank1500']
        for major in ['T cell','Macrophage']:
            keep=native_obs.major.eq(major).to_numpy()
            for j,r in enumerate(records):
                valid=keep&mask[:,j];values=y[valid,j]
                if not len(values):continue
                informative=valid&(det[:,j]>=plan['min_detected_genes_for_support']);stats.append(dict(gse_id=sid,major=major,programme_id=r['programme_id'],cells=int(keep.sum()),valid_scores=len(values),informative_gene_support_cells=int(informative.sum()),gene_support_threshold=plan['min_detected_genes_for_support'],score_median=float(np.median(values)),score_q10=float(np.quantile(values,.1)),score_q90=float(np.quantile(values,.9)),score_RNA_spearman=rho(values,native_obs.full_RNA_total.to_numpy()[valid]),score_detected_spearman=rho(values,native_obs.full_detected_genes.to_numpy()[valid]),maxrank1000_rho=rho(values,matrices['scores_maxrank1000'][valid,j]),maxrank3000_rho=rho(values,matrices['scores_maxrank3000'][valid,j]),TPM_signature_rho=rho(values,t[valid,j]) if r['method']=='cNMF' else np.nan,stability_main_candidate=r['stability_main_candidate'],foreign_RNA_sensitive=r['foreign_RNA_sensitive'],basic_expression_review_flag=r['basic_expression_review_flag']))
            for node,meta in native_obs.loc[keep].groupby('subCluster',observed=True):
                ii=native_obs.index.get_indexer(meta.index)
                for j,r in enumerate(records):
                    if major not in r['allowed_majors']:continue
                    vals=y[ii,j];good=np.isfinite(vals)
                    if not good.any():continue
                    node_stats.append(dict(gse_id=sid,major=major,subCluster=node,programme_id=r['programme_id'],cells=len(ii),valid_scores=int(good.sum()),score_median=float(np.median(vals[good])),score_q10=float(np.quantile(vals[good],.1)),score_q90=float(np.quantile(vals[good],.9)),score_RNA_spearman=rho(vals[good],native_obs.full_RNA_total.to_numpy()[ii][good]),informative_gene_support_cells=int((det[ii,j]>=plan['min_detected_genes_for_support']).sum())))
        # Private per sample-node means with exact inherited block and no independent-patient claim.
        aggregates=[]
        for (sample,node),meta in native_obs.groupby(['sample_uid','node_id'],observed=True):
            ii=native_obs.index.get_indexer(meta.index);major=meta.major.iloc[0];assert meta.major.nunique()==1 and meta.analysis_block.nunique()==1
            for j,r in enumerate(records):
                if major not in r['allowed_majors']:continue
                v=y[ii,j];ok=np.isfinite(v)
                aggregates.append(dict(sample_uid=sample,node_id=node,major=major,analysis_block=meta.analysis_block.iloc[0],patient_identity_reviewed=bool(meta.patient_identity_reviewed.iloc[0]),programme_id=r['programme_id'],cells=len(ii),valid_score_cells=int(ok.sum()),programme_mean=float(v[ok].mean()) if ok.any() else np.nan,programme_median=float(np.median(v[ok])) if ok.any() else np.nan,informative_gene_support_cells=int(((det[ii,j]>=plan['min_detected_genes_for_support'])&ok).sum()),eligible_cell_count20=len(ii)>=20))
        pd.DataFrame(aggregates).to_csv(root/'private_sample_node_programme_means.tsv.gz',sep='\t',index=False,compression='gzip');files.append(dict(gse_id=sid,cells=a.n_obs,genes=a.n_vars,score_file_sha256=sha(root/'cell_programme_scores.npz'),metadata_sha256=sha(root/'native_cell_identity_metadata.tsv'),sample_node_means_sha256=sha(root/'private_sample_node_programme_means.tsv.gz')));(O/'progress.json').write_text(json.dumps(dict(completed_sources=len(files),summary=files),indent=2)+'\n');del a,matrices,y,t,det,fraction,mask;gc.collect()
    assert total==240533 and len(coverage)==490
    for n,x in [('all490_native_gene_coverage.tsv',pd.DataFrame(coverage)),('all_source_major_programme_distribution.tsv',pd.DataFrame(stats)),('all_source_node_programme_distribution.tsv',pd.DataFrame(node_stats))]:x.to_csv(O/n,sep='\t',index=False)
    out=dict(status='COMPLETE_NATIVE_UCELL_CANDIDATE_SCORES_PENDING_INDEPENDENT_AUDIT',started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),dictionary_sha256=di['dictionary_sha256'],input_inventory_sha256=sha(I),cells=total,programmes=98,main_stability_candidates=65,version='pyucell0.8.0',primary_maxrank=1500,ties_method='average',missing_genes='skip_with_explicit_measurement_gate',native_gene_axes_preserved=True,original_identities_changed=False,formal_state_labels=False,high_state_percentile_forced=False,internal_source_projection='GSE183904 previously_used_not_independent_validation',source_files=files,public_output_sha256={p.name:sha(p) for p in O.glob('*.tsv')});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
