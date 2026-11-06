import json
from datetime import datetime,timezone
import pandas as pd,numpy as np,anndata as ad
r=PROJECT_ROOT/'runs/20261106_v7_other_unresolved_review_v1';o=r/'stage_g_ILC_candidate_followup';o.mkdir(exist_ok=False)
q=pd.read_pickle(r/'stage_c_QC_views/QC_decisions_PRIVATE.pkl.gz');t=pd.read_pickle(r/'stage_f_working_annotations/ILC_working_annotation_PRIVATE.pkl.gz');t=t[t.working_parent.eq('ILC_like_candidate')];assert len(t)==166
a=ad.read_h5ad(r/'stage_a_candidates/candidate_marker_raw_PRIVATE.h5ad');a=a[t.index].copy();x=pd.DataFrame(a.X.toarray(),index=a.obs_names,columns=a.var_names);meta=q.loc[t.index]
type2=x.GATA3.ge(2)&x[['PTGDR2','IL1RL1']].ge(2).any(axis=1)
type3=x.RORC.ge(2)&x[['IL23R','CCR6','IL22']].ge(2).any(axis=1)
assert (type2|type3).all() and x.PTPRC.ge(2).all() and not meta.T_positive_source_aware.any()
program=np.select([type2&type3,type2,type3],['type2_and_type3_shared_programmes','type2_programme_candidate','type3_programme_candidate'],default='unexpected');assert not (program=='unexpected').any()
v=t.copy();v['candidate_programme']=program;v.to_pickle(o/'ILC_166_candidate_programmes_PRIVATE.pkl.gz');pd.testing.assert_frame_equal(v,pd.read_pickle(o/'ILC_166_candidate_programmes_PRIVATE.pkl.gz'))
rows=[];markers=[]
genes=['PTPRC','IL7R','GATA3','PTGDR2','IL1RL1','IL5','IL13','RORC','IL23R','CCR6','IL22','KIT','KLRB1','CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','NKG7','GNLY','KLRD1']
for label in sorted(set(program)):
 take=program==label;ids=v.index[take];m=meta.loc[ids];rows.append({'candidate_programme':label,'n_cells':len(ids),'n_sources':m.gse_id.nunique(),'n_res0_4_partitions':v.loc[ids,'diagnostic_cluster'].nunique(),'max_source_fraction':float(m.gse_id.value_counts().iloc[0]/len(ids)),'raw_library_median':float(m.full_raw_library_total.median()),'decontx_median':float(pd.to_numeric(m.decontx_contamination,errors='coerce').median()),'n_T_apparatus_positive':int(m.T_positive_source_aware.sum()),'final_ILC_subtype_accepted':False})
 for gene in genes:markers.append({'candidate_programme':label,'gene':gene,'n_cells':len(ids),'n_raw_ge1':int(x.loc[ids,gene].ge(1).sum()),'n_raw_ge2':int(x.loc[ids,gene].ge(2).sum()),'interpretation':'Positive RNA evidence only; absent counts are not lineage-negative proof. No source-measured detection rate claimed in this table.'})
pd.DataFrame(rows).to_csv(o/'candidate_programme_summary.tsv',sep='\t',index=False)
pd.DataFrame(markers).to_csv(o/'candidate_positive_RNA_counts.tsv',sep='\t',index=False)
v.groupby(['diagnostic_cluster','candidate_programme'],observed=True).size().rename('n_cells').reset_index().to_csv(o/'candidate_programmes_by_partition.tsv',sep='\t',index=False)
out={'status':'candidate_positive_programmes_readback_verified_not_final_ILC_annotation','UTC':datetime.now(timezone.utc).isoformat(),'n_candidates':166,'n_sources':meta.gse_id.nunique(),'n_type2_programme':int(type2.sum()),'n_type3_programme':int(type3.sum()),'n_both':int((type2&type3).sum()),'all_have_PTPRC_raw_ge2':True,'n_current_T_positive':int(meta.T_positive_source_aware.sum()),'no_Lin_negative_proof_from_dropout':True,'final_ILC_subtype_accepted':False};(o/'complete.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
