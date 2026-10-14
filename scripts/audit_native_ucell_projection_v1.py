"""Fresh native genes: full quantities/masks/means and manual tied-rank scores."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd,anndata as ad
from scipy.sparse import csr_matrix

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def manual_ranks(row,cap):
    ranks=np.full(len(row),cap,dtype=float);nonzero=np.flatnonzero(row>0);order=nonzero[np.argsort(-row[nonzero],kind='stable')];start=0
    while start<len(order):
        end=start+1
        while end<len(order) and row[order[end]]==row[order[start]]:end+=1
        rank=(start+1+end)/2; ranks[order[start:end]]=min(rank,cap);start=end
    return ranks
def manual_score(ranks,idx,cap):
    n=len(idx);baseline=n*(n+1)/2;return 1-(sum(float(ranks[j]) for j in idx)-baseline)/(n*cap-baseline)
def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'state_ucell_native_projection_20261014_v2';I=json.loads((O/'inventory.json').read_text());D=R/'state_programme_dictionary_20261014_v2/frozen_programme_dictionary.json';assert sha(D)==I['dictionary_sha256'];records=json.loads(D.read_text())['programmes'];ids=[r['programme_id'] for r in records];rawinv=json.loads((R/'state_singlecell_inputs_acceptance_20261014_v2/accepted_input_inventory.json').read_text());job=json.loads((R/'native_ucell_projection_20261014_v2.execution.json').read_text());assert job['exit_code']==0
    for n,h in I['public_output_sha256'].items():assert sha(O/n)==h
    scalar_checks=0;quantity_checks=0;mean_checks=0;max_error=0.;means_error=0.;rng=np.random.default_rng(20261014);summaries=[]
    for item in I['source_files']:
        sid=item['gse_id'];root=O/'private_source_scores'/sid;p=Path(rawinv['verified_dataset_files'][sid]['path']);assert sha(p)==rawinv['verified_dataset_files'][sid]['sha256'];assert sha(root/'cell_programme_scores.npz')==item['score_file_sha256'];assert sha(root/'native_cell_identity_metadata.tsv')==item['metadata_sha256'];assert sha(root/'private_sample_node_programme_means.tsv.gz')==item['sample_node_means_sha256'];a=ad.read_h5ad(p);z=dict(np.load(root/'cell_programme_scores.npz',allow_pickle=False));assert np.array_equal(z['cell_ids'],a.obs_names.to_numpy(str)) and np.array_equal(z['programme_ids'],np.array(ids));meta=pd.read_csv(root/'native_cell_identity_metadata.tsv',sep='\t',index_col=0,keep_default_na=False)
        for key in ['major','subCluster','node_id','sample_uid','core_position','analysis_block','patient_identity_reviewed']:assert np.array_equal(meta[key].to_numpy(str),a.obs[key].to_numpy(str)),key
        union=sorted({g for r in records for g in r['primary_genes'] if g in a.var_names});x=a.X[:,a.var_names.get_indexer(union)].tocsr();x.sum_duplicates();lookup={g:i for i,g in enumerate(union)};rows=[];cols=[]
        for j,r in enumerate(records):
            for g in r['primary_genes']:
                if g in lookup:rows.append(lookup[g]);cols.append(j)
        A=csr_matrix((np.ones(len(rows),np.int64),(rows,cols)),shape=(len(union),98));count=(x@A).toarray();det=((x>0).astype(np.int64)@A).toarray();assert np.array_equal(det,z['detected_programme_genes']);fraction=count/a.obs.full_RNA_total.to_numpy()[:,None];assert np.allclose(fraction,z['programme_RNA_fraction'],atol=1e-7);quantity_checks+=2*count.size
        mask=np.zeros((a.n_obs,98),bool)
        for j,r in enumerate(records):
            present=a.var_names.isin(r['primary_genes']).sum();mask[:,j]=a.obs.major.isin(r['allowed_majors']).to_numpy()&(present/50>=.8)&(present>=10)
        assert np.array_equal(mask,z['model_scope_and_coverage_mask'])
        for cap in [1000,1500,3000]:
            y=z['scores_maxrank'+str(cap)];assert y.shape==(a.n_obs,98) and np.array_equal(np.isfinite(y),mask);assert np.all((y[mask]>=-1e-6)&(y[mask]<=1+1e-6))
        selected=[]
        for major in ['T cell','Macrophage']:
            ii=np.flatnonzero(a.obs.major.eq(major).to_numpy());selected.extend(rng.choice(ii,min(20,len(ii)),replace=False).tolist())
        primary_indices=[a.var_names.get_indexer([g for g in r['primary_genes'] if g in a.var_names]) for r in records];tpm_indices=[a.var_names.get_indexer([g for g in r['TPM_sensitivity_genes'] if g in a.var_names]) for r in records]
        for cell in selected:
            row=a.X[cell].toarray().ravel()
            for cap in [1000,1500,3000]:
                ranks=manual_ranks(row,cap)
                for j in range(98):
                    if not mask[cell,j]:continue
                    expected=manual_score(ranks,primary_indices[j],cap);reported=z['scores_maxrank'+str(cap)][cell,j];assert np.isclose(expected,reported,atol=1e-6);scalar_checks+=1;max_error=max(max_error,abs(expected-reported))
            ranks=manual_ranks(row,1500)
            for j,r in enumerate(records):
                if r['method']!='cNMF' or not mask[cell,j]:continue
                expected=manual_score(ranks,tpm_indices[j],1500);reported=z['TPM_signature_scores_maxrank1500'][cell,j];assert np.isclose(expected,reported,atol=1e-6);scalar_checks+=1;max_error=max(max_error,abs(expected-reported))
        groups=pd.MultiIndex.from_frame(a.obs[['sample_uid','node_id']]);codes,unique=pd.factorize(groups,sort=True);inc=csr_matrix((np.ones(a.n_obs),(codes,np.arange(a.n_obs))),shape=(len(unique),a.n_obs));y=z['scores_maxrank1500'];sums=inc@np.nan_to_num(y);num=inc@mask.astype(float);means=np.divide(sums,num,out=np.full(sums.shape,np.nan),where=num>0);genecount=inc@((det>=5)&mask).astype(float);tab=pd.read_csv(root/'private_sample_node_programme_means.tsv.gz',sep='\t',keep_default_na=False);unit_map={v:i for i,v in enumerate(unique)};prog_map={v:i for i,v in enumerate(ids)}
        for rr in tab.itertuples():
            i=unit_map[rr.sample_uid,rr.node_id];j=prog_map[rr.programme_id];reported=float(rr.programme_mean) if rr.programme_mean!='' else np.nan;assert np.isclose(means[i,j],reported,atol=2e-6,equal_nan=True);assert rr.valid_score_cells==num[i,j] and rr.informative_gene_support_cells==genecount[i,j];mean_checks+=1
            if np.isfinite(reported):means_error=max(means_error,abs(means[i,j]-reported))
        summaries.append(dict(gse_id=sid,cells=a.n_obs,native_genes=a.n_vars,manual_cells=len(selected),sample_node_means_checked=len(tab),all_masks_quantities_identity_exact=True));print(json.dumps(summaries[-1]),flush=True);del a,z,x,count,det,fraction,mask,y;gc.collect()
    E=R/'native_ucell_projection_audit_20261014_v1';assert not E.exists();E.mkdir();result=dict(status='PASS_NATIVE_UCELL_ALL_CELL_SUPPORT_MASKS_AND_MANUAL_REAL_CELL_SCORES',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),cells=240533,programmes=98,quantity_values_checked=quantity_checks,manual_scalar_scores_checked=scalar_checks,max_manual_score_difference=max_error,sample_node_means_checked=mean_checks,max_sample_mean_difference=means_error,summary=summaries,all_cell_score_values_independently_recomputed=False,all_public_aggregate_correlations_independently_recomputed=False,formal_state_labels=False,biological_or_mechanistic_validation=False);(E/'independent_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
