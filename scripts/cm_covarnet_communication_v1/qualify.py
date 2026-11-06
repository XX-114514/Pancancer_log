from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime,zarr
from anndata.io import read_elem
r=Path(__import__('os').environ['PC_ROOT']);v=r/'Project_v3/runs/20260818_pancancer_5M_merge_prep/10_final_v5_malignancy_covarnet_liana/20260825_v5_final_v1';f=r/'Project_v5/frozen/covarnet';o=r/'Project_v3/runs/20261106_V5_K9_CoVarNet_communication_v1';o.mkdir(exist_ok=True)
for d in ['audit','private','inputs','results','logs','config','scripts']: (o/d).mkdir(exist_ok=True)
W=pd.read_csv(f/'nmf_W_subcluster_by_cm.tsv',sep='\t',index_col=0);H=pd.read_csv(f/'nmf_H_cm_by_sample.tsv',sep='\t',index_col=0);N=pd.read_csv(f/'covarnet_each_nodes.tsv',sep='\t');N.to_csv(o/'audit/network_members.tsv',sep='\t',index=False)
assert list(H.index)==list(W.columns) and H.columns.is_unique and W.index.is_unique
valid=np.isfinite(H).all(axis=0)&(H>=0).all(axis=0)&(H.sum(axis=0)>0);ties=np.isclose(H,H.max(axis=0).values[None,:],rtol=1e-10,atol=1e-12).sum(axis=0)>1;valid &= ~ties
s=pd.DataFrame({'CMT':H.idxmax().where(valid),'eligible_H':valid,'H_tie':ties});s.index.name='sampleID'
rows=[]
for cm in W:
 top=W[cm].sort_values(ascending=False).head(10).index.tolist();nodes=N.loc[N.cm==cm,'subCluster'].tolist()
 for rank,n in enumerate(top,1):rows.append({'cm':cm,'subCluster':n,'top10_rank':rank,'network_member':n in nodes,'W':W.loc[n,cm]})
 for n in set(nodes)-set(top):rows.append({'cm':cm,'subCluster':n,'top10_rank':np.nan,'network_member':True,'W':W.loc[n,cm]})
pd.DataFrame(rows).to_csv(o/'audit/member_comparison.tsv',sep='\t',index=False)
manifest=[]
for p in [f/'nmf_W_subcluster_by_cm.tsv',f/'nmf_H_cm_by_sample.tsv',f/'covarnet_each_nodes.tsv',f/'frequency_normalized.tsv',f/'frequency_raw.tsv',v/'config/run_config.json',v/'checkpoints/final_annotation_v6_identity_v5_malignancy_v3.pkl.gz']:
 h=hashlib.sha256()
 with p.open('rb') as handle:
  for b in iter(lambda:handle.read(8*1024*1024),b''):h.update(b)
 manifest.append({'path':str(p),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
 print('HASHED',p.name,flush=True)
pd.DataFrame(manifest).to_csv(o/'private/input_manifest.tsv',sep='\t',index=False)
q=pd.read_pickle(v/'checkpoints/final_annotation_v6_identity_v5_malignancy_v3.pkl.gz');print('LOADED',q.shape,flush=True)
assert q.index.is_unique and np.array_equal(q.global_position.values,np.arange(len(q)))
z=zarr.open(str(r/'Project_v3/runs/20260818_pancancer_5M_merge_prep/02_merge_new_cohort/objects/5M_merge_concat.zarr'),mode='r')
idx=read_elem(z['obs']['_index']);assert np.array_equal(q.index.astype(str),idx.astype(str))
cols=['sampleID','gse_id','cancer_type','final_annotation_with_malignancy_v6','global_position']
q=q[cols].copy();q['CMT']=q.sampleID.astype(str).map(s.CMT);q.to_pickle(o/'private/cell_metadata.pkl.gz')
meta=q.groupby('sampleID',observed=True).agg(gse_id=('gse_id','first'),cancer_type=('cancer_type','first'),n_cells=('global_position','size'));s=s.join(meta)
s.to_csv(o/'private/sample_eligibility.tsv',sep='\t');q.groupby(['CMT','final_annotation_with_malignancy_v6'],observed=True).size().reset_index(name='n_cells').to_csv(o/'audit/CMT_identity_cell_counts.tsv',sep='\t',index=False)
pd.crosstab(s.CMT,s.cancer_type).to_csv(o/'audit/CMT_cancer_counts.tsv',sep='\t');pd.crosstab(s.CMT,s.gse_id).to_csv(o/'audit/CMT_source_counts.tsv',sep='\t')
res={'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'W_shape':W.shape,'H_shape':H.shape,'H_nonfinite':int((~np.isfinite(H)).sum().sum()),'H_sum_range':[float(H.sum().min()),float(H.sum().max())],'ties':int(ties.sum()),'CMT_counts':s.CMT.value_counts().to_dict(),'network_members':N.groupby('cm').size().to_dict(),'cells':len(q),'identity_labels':int(q.iloc[:,3].nunique()),'index_equal':True,'patient_mapping':'not_yet_verified'}
(o/'audit/qualification.json').write_text(json.dumps(res,indent=2));print(json.dumps(res),flush=True)
