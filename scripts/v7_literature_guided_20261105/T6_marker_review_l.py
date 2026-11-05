import os,json,zlib,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad,zarr,h5py
from anndata.io import read_elem
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))));r=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';o=r/'stage_l_T6_marker_review_v1';o.mkdir(exist_ok=False)
k=r/'stage_k_T6_source_robust_diagnostic_v1/fit';terminal=json.loads((r/'launch/T6_diagnostic_k_terminal.json').read_text());assert terminal['exit_code']==0
a=ad.read_h5ad(r/'stage_j_T6_fullgene/T6_v7_fullgene_raw.h5ad');md=pd.read_csv(k/'metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id');assert a.obs_names.equals(md.index) and a.shape==(156631,95038)
spec=json.loads((r/'config/T6_marker_panels_v1.json').read_text());genes=list(dict.fromkeys(g for p in spec['panels']+spec['reference_panels'] for g in p['genes']));ids=a.var_names.get_indexer(genes);raw=np.zeros((a.n_obs,len(genes)),np.int32);found=ids>=0;raw[:,found]=a.X[:,ids[found]].toarray();tot=np.asarray(a.X.sum(axis=1)).ravel();cp=raw*(10000/np.maximum(tot,1)[:,None]);ix={g:j for j,g in enumerate(genes)}
obs=a.obs.copy();cl=md['source_snn_res.0.4'].astype(int).to_numpy();nc=cl.max()+1;assert set(cl)==set(range(nc));obs['T6_review_cluster_source_res0.4']=cl
gses=sorted(obs.gse_id.astype(str).unique());av=[]
for gs in gses:
 with h5py.File(ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad'),'r') as f:var=pd.Index(read_elem(f['var']).index.astype(str))
 av.append(pd.Index(genes).isin(var))
av=np.asarray(av).T;gidx=pd.Index(gses).get_indexer(obs.gse_id.astype(str));measured=av[:,gidx].T;assert not ((raw>0)&~measured).any()
H=sp.csr_matrix((np.ones(a.n_obs,np.int32),(cl,np.arange(a.n_obs))),shape=(nc,a.n_obs));S=sp.csr_matrix((np.ones(a.n_obs,np.int32),(cl,gidx)),shape=(nc,len(gses)));den=np.asarray(S@av.T.astype(np.int32));nd=np.asarray(H@(raw>0).astype(np.int32));sums=np.asarray(H@cp);n=np.bincount(cl,minlength=nc)
e=pd.DataFrame([{'cluster':c,'gene':g,'n_cells':int(n[c]),'n_measured':int(den[c,j]),'coverage_fraction':float(den[c,j]/n[c]),'n_expressing':int(nd[c,j]),'pct_expressing_measured':float(nd[c,j]/den[c,j]) if den[c,j] else np.nan,'mean_cp10k_measured':float(sums[c,j]/den[c,j]) if den[c,j] else np.nan} for c in range(nc) for j,g in enumerate(genes)]);e.to_csv(o/'T6_measured_marker_expression.tsv',sep='\t',index=False,na_rep='NA')
# Positive marker programs are independent evidence, never a hard annotation rule.
def hits(gs,cut):return (raw[:,[ix[g] for g in gs]]>=cut).sum(axis=1)
parent=(hits(['CD3D','CD3E','CD3G'],2)>=1)&(hits(['TRAC','TRBC1','TRBC2'],2)>=1)
obs['T_parent_cross_component_raw2']=parent
for cut in [1,2]:obs['Treg_core_positive_raw'+str(cut)]=parent&(raw[:,ix['FOXP3']]>=cut)&(hits(['IL2RA','IKZF2','IKZF4'],cut)>=1)
comp={p['id']:p['genes'] for p in spec['panels'] if p['level']=='competing_parent'};anycomp=np.zeros(a.n_obs,bool)
for name,gs in comp.items():
 flag=(hits(gs,1)>=3)|(hits(gs,2)>=2);obs[name+'_program_review']=flag;anycomp|=flag
core=zarr.open_group(str(ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r');pos=obs.global_position.to_numpy().astype(int);assert pd.Index(read_elem(core['obs/_index']))[pos].equals(obs.index)
cal=pd.DataFrame({key:read_elem(core['obs/'+key]) for key in ['gse_id','sampleID','doublet_method','doublet_score']});cal['doublet_score']=pd.to_numeric(cal.doublet_score,errors='coerce');grp=cal.groupby(['gse_id','sampleID','doublet_method'],dropna=False,observed=True)['doublet_score'];pct=grp.rank(pct=True,method='average').to_numpy()[pos];size=grp.transform('count').to_numpy()[pos]
obs['doublet_within_sample_method_percentile']=pct;obs['doublet_calibration_evaluable']=(size>=100)&np.isfinite(pct);obs['joint_competitor_doublet_candidate']=parent&anycomp&obs.doublet_calibration_evaluable.to_numpy()&(pct>=.95)
obs['parent_identity_evidence']=np.where(parent,'T_marker_supported_candidate','T_candidate_low_signal');obs['identity_subtype']='not_finally_assigned';cf=obs[[name+'_program_review' for name in comp]].to_numpy();names=list(comp);obs['parent_conflicts']=[','.join(name for name,yes in zip(names,row) if yes) for row in cf]
obs.to_pickle(o/'T6_marker_review_PRIVATE.pkl.gz')
rows=[];refs=[]
for c in range(nc):
 t=obs.iloc[np.flatnonzero(cl==c)];z={'cluster':int(c),'n_cells':len(t),'n_sources':t.gse_id.nunique(),'max_source_fraction':float(t.gse_id.astype(str).value_counts().iloc[0]/len(t)),'n_T_parent_supported_raw2':int(t.T_parent_cross_component_raw2.sum()),'n_Treg_core_raw1':int(t.Treg_core_positive_raw1.sum()),'n_Treg_core_raw2':int(t.Treg_core_positive_raw2.sum()),'n_joint_doublet_candidates':int(t.joint_competitor_doublet_candidate.sum()),'doublet_percentile_median':float(t.doublet_within_sample_method_percentile.median()),'decontx_median':float(pd.to_numeric(t.decontx_contamination,errors='coerce').median())}
 for name in comp:z['n_'+name]=int(t[name+'_program_review'].sum())
 rows.append(z)
 for name in ['ev_celltypist_label','scatomic_full_label']:
  if name in t:
   for label,number in t[name].fillna('missing').astype(str).value_counts().head(5).items():refs.append({'cluster':int(c),'reference_field':name,'label':label,'n_cells':int(number),'total_cluster_cells':len(t)})
pd.DataFrame(rows).to_csv(o/'T6_cluster_marker_QC_summary.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(refs).to_csv(o/'T6_cluster_reference_top5.tsv',sep='\t',index=False)
d=pd.read_csv(k/'RNA_FindAllMarkers_source_res0.4.tsv',sep='\t');d=d[(d.avg_log2FC>0)&(d.p_val_adj<.05)].sort_values(['cluster','avg_log2FC'],ascending=[True,False]).groupby('cluster').head(10);d.to_csv(o/'T6_top10_RNA_DEG.tsv',sep='\t',index=False,na_rep='NA')
manifest={'status':'completed_T6_marker_and_QC_evidence_review_not_final_subtypes','UTC':datetime.now(timezone.utc).isoformat(),'n_cells_retained':a.n_obs,'n_sources':len(gses),'n_clusters':int(nc),'n_marker_genes':len(genes),'n_marker_records':len(e),'genes_absent_from_full_native_axis':[g for g,j in zip(genes,ids) if j<0],'n_T_parent_supported_raw2':int(parent.sum()),'n_Treg_core_positive_raw1':int(obs.Treg_core_positive_raw1.sum()),'n_Treg_core_positive_raw2':int(obs.Treg_core_positive_raw2.sum()),'n_joint_doublet_candidates':int(obs.joint_competitor_doublet_candidate.sum()),'tiny_clusters_under3':[{'cluster':int(c),'n_cells':int(n[c])} for c in range(nc) if n[c]<3],'old_fields_preserved':True,'fine_final_label_changes':0,'limits':['All outputs scoped to old V7 T6 only, not a whole-T or Treg census.','FOXP3 plus another core marker is positive RNA evidence, not a calibrated exclusive Treg classifier.','Competitor programs and doublet quantiles create review candidates only.','Quick-TRANSfer warnings in the fit remain for source/resolution stability review.']}
(o/'complete.json').write_text(json.dumps(manifest,indent=2));packet={'manifest':manifest,'K_complete':json.loads((k/'complete.json').read_text()),'K_terminal':terminal,'tables':{name:pd.read_csv(o/name,sep='\t').to_json(orient='records') for name in ['T6_measured_marker_expression.tsv','T6_cluster_marker_QC_summary.tsv','T6_cluster_reference_top5.tsv','T6_top10_RNA_DEG.tsv']}}
b=zlib.compress(json.dumps(packet).encode(),9);(o/'aggregate_L_transfer_v1.zlib').write_bytes(b);print(json.dumps({'manifest':manifest,'packet_bytes':len(b),'packet_sha256':hashlib.sha256(b).hexdigest()}))
