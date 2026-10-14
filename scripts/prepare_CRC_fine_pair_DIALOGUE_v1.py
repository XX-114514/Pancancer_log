"""Patient-first frozen CRC CM01/03 pair; discovery/holdout cells stay private.
No outcome or expression is used to choose patient/library or partition.
"""
import os,json,hashlib,gzip
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmwrite
from sklearn.decomposition import PCA
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);C=Path(os.environ['PANCANCER_CRC_DIALOGUE_PLAN']);plan=json.loads(C.read_text());assert sha(__file__)==plan['producer_sha256'];I=R/'state_singlecell_inputs_20261014_v1';A=R/'state_singlecell_inputs_acceptance_20261014_v2';assert A.exists();proof_files=list(A.glob('*.json'));proofs=[json.loads(f.read_text()) for f in proof_files];assert any(j.get('status')=='PASS_ALL240533_CELL_NATIVE_COUNTS_FULLRNA_AND1832_UNIT_AGGREGATES' for j in proofs), 'Full raw input readback acceptance required'
 src=I/'private_native_sources/GSE178341.h5ad';O=R/'CRC_fine_pair_DIALOGUE_20261014_v1';assert not O.exists();a=ad.read_h5ad(src);assert a.n_obs==54566 and a.n_vars==26895;obs=a.obs.copy();native=R/'immune_same_name_family_source_review_v1/major_pair_native_assembly_v1/private_group_arrays/MAJ158.npz';n=np.load(native);assert (n['corrected_identity_reviewed_for_next_inference']).all();samples=n['samples'];labels=plan['labels'];take=obs[obs.sample_uid.isin(samples)&obs.subCluster.isin(labels)];tab=take.groupby(['analysis_block','sample_uid','subCluster'],observed=True).size().unstack(fill_value=0);tab=tab[(tab[labels]>=20).all(axis=1)];chosen=tab.reset_index().sort_values(['analysis_block','sample_uid']).drop_duplicates('analysis_block');assert len(tab)==42 and len(chosen)==38
 rng=np.random.default_rng(plan['seed']);blocks=np.sort(chosen.analysis_block.astype(str).to_numpy());discovery=set(rng.choice(blocks,26,replace=False));chosen['partition']=chosen.analysis_block.map(lambda b:'discovery' if b in discovery else 'holdout');selected=set(chosen.sample_uid);obs=obs[obs.sample_uid.isin(selected)&obs.subCluster.isin(labels)].copy();obs['partition']=obs.analysis_block.map(lambda b:'discovery' if b in discovery else 'holdout');O.mkdir();(O/'private_inputs').mkdir();chosen.to_csv(O/'private_inputs/patient_library_partition.tsv',sep='\t',index=False);summary=[]
 for type_name,label in zip(plan['cell_type_names'],labels):
  part=obs[obs.subCluster.eq(label)];positions=[]
  for _,v in part.groupby('analysis_block',observed=True):positions.extend(rng.choice(a.obs_names.get_indexer(v.index),min(plan['max_cells_per_patient_type'],len(v)),replace=False))
  positions=np.sort(positions);md=obs.loc[a.obs_names[positions]].copy();assert md.patient_identity_reviewed.all();raw=a.X[positions].tocsr();total=np.asarray(raw.sum(1)).ravel();assert np.array_equal(total,md.full_RNA_total.to_numpy());cp=raw.astype(float);cp.data*=np.repeat(10000/total,np.diff(cp.indptr));cp.data=np.log1p(cp.data);train=md.partition.eq('discovery').to_numpy();genes=a.var_names.astype(str).to_numpy();active=np.asarray((raw[train]>0).mean(0)).ravel()>=.01;active&=~pd.Index(genes).str.match(r'^(MT-|RPL|RPS)');ix=np.flatnonzero(active);X=cp[:,ix].tocsr();mean=np.asarray(X[train].mean(0)).ravel();var=np.asarray(X[train].power(2).mean(0)).ravel()-mean**2;order=np.argsort(-var,kind='stable');hvg=order[:min(2000,len(order))];pca=PCA(n_components=20,svd_solver='randomized',random_state=plan['seed']);pca.fit(X[train][:,hvg].toarray());pc=pca.transform(X[:,hvg].toarray());cells=md.index.to_numpy(dtype=str)
  quality=np.column_stack((np.log1p(total),np.log1p(md.full_detected_genes.to_numpy())));mu=quality[train].mean(0);sd=quality[train].std(0);assert (sd>0).all();z=(quality-mu)/sd;metadata=pd.DataFrame(dict(cell=cells,sample=md.analysis_block.astype(str).to_numpy(),cellQ=z[:,0],geneQ=z[:,1],partition=md.partition.to_numpy()))
  for partition in ['discovery','holdout']:
   m=metadata.partition.eq(partition).to_numpy();prefix=O/'private_inputs'/(type_name+'_'+partition)
   with gzip.open(str(prefix)+'_tpm.mtx.gz','wb',compresslevel=3) as f:mmwrite(f,X[m].T)
   metadata.loc[m].to_csv(str(prefix)+'_metadata.tsv',sep='\t',index=False);pd.DataFrame(pc[m],index=cells[m],columns=['PC'+str(i+1) for i in range(20)]).to_csv(str(prefix)+'_PC.tsv',sep='\t');pd.Series(genes[ix]).to_csv(str(prefix)+'_genes.tsv',sep='\t',header=False,index=False);summary.append(dict(cell_type=type_name,partition=partition,n_cells=int(m.sum()),n_patient_blocks=metadata.loc[m,'sample'].nunique(),n_candidate_genes=len(ix),n_HVG=len(hvg),n_PC=20))
  np.savez_compressed(O/'private_inputs'/(type_name+'_frozen_projection.npz'),hvg_gene_names=genes[ix[hvg]],pca_mean=pca.mean_,pca_components=pca.components_,quality_mean=mu,quality_sd=sd)
 pd.DataFrame(summary).to_csv(O/'input_scope_summary.tsv',sep='\t',index=False);j=dict(status='PATIENT_SPLIT_AND_REAL_INPUTS_PREPARED_NO_MCP_FIT_YET',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(C),source_input_sha256=sha(src),acceptance_json_sha256={f.name:sha(f) for f in proof_files},patient_blocks=38,discovery_patient_blocks=26,holdout_patient_blocks=12,original_labels=labels,summary=summary,clinical_outcomes_used=False,already_used_discovery_source=True,not_project_wide_independent_validation=True,outputs={str(f.relative_to(O)):sha(f) for f in O.rglob('*') if f.is_file()});(O/'input_inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps({k:v for k,v in j.items() if k!='outputs'}))
if __name__=='__main__':main()
