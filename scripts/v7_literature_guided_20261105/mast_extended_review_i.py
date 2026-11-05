import os,json,hashlib,zlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp,h5py
from anndata.io import read_elem
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))));r=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';o=r/'stage_i_mast_extended_review_v1';o.mkdir(exist_ok=False)
spec=json.loads((r/'config/mast_marker_panels_v4.json').read_text());h=r/'stage_h_mast_source_robust_refit_v1';assert json.loads((r/'launch/mast_refit_h_terminal.json').read_text())['exit_code']==0
obs=pd.read_pickle(h/'mast_cell_review_v3_PRIVATE.pkl.gz');a=ad.read_h5ad(r/'stage_b_mast_fullgene/mast_v7_fullgene_raw.h5ad');assert a.obs_names.equals(obs.index)
md=pd.read_csv(h/'refit/metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id');fit=obs.Mast_source_robust_refit_candidate_v2.to_numpy();assert md.index.equals(obs.index[fit]);cl=md['source_snn_res.0.4'].to_numpy().astype(int);assert set(cl)==set(range(15))
genes=list(dict.fromkeys([g for p in spec['panels']+spec['panMC_reference_panels'] for g in p['genes']]+spec['extra_review_genes']));jj=a.var_names.get_indexer(genes);assert (jj>=0).all();raw=a.X[:,jj].toarray();ix={g:i for i,g in enumerate(genes)};tot=np.asarray(a.X.sum(axis=1)).ravel();cp=raw*(10000/np.maximum(tot,1)[:,None])
gses=sorted(obs.gse_id.astype(str).unique());av=[]
for gs in gses:
 with h5py.File(ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad'),'r') as f:var=pd.Index(read_elem(f['var']).index.astype(str))
 av.append(pd.Index(genes).isin(var))
av=np.asarray(av).T;gidx=pd.Index(gses).get_indexer(obs.gse_id.astype(str))
comp={'T':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],'B':['MS4A1','CD79A','CD79B','CD19','BANK1'],'Plasma':['MZB1','JCHAIN','SDC1','DERL3'],'Myeloid':['LST1','CSF1R','CD14','CD163','MS4A7','C1QA'],'Stromal':['COL1A1','COL1A2','DCN','LUM'],'Epithelial':['EPCAM','KRT8','KRT18','KRT19','KRT14','KRT17'],'Endothelial':['PECAM1','VWF','CLDN5','KDR'],'Platelet':['PF4','PPBP','GP9','ITGA2B']}
anycomp=np.zeros(len(obs),bool);trigger={}
for name,panel in comp.items():
 ids=[ix[g] for g in panel];n1=(raw[:,ids]>=1).sum(axis=1);n2=(raw[:,ids]>=2).sum(axis=1);flag=(n1>=3)|(n2>=2);trigger[name]=flag;anycomp|=flag
 obs['extended_'+name+'_n_detected_genes']=n1;obs['extended_'+name+'_n_genes_ge2UMI']=n2;obs['extended_'+name+'_program_review']=flag
parent=obs.Mast_tryptase_plus_support_raw2.to_numpy();cal=obs.doublet_calibration_evaluable.to_numpy();pct=obs.doublet_within_sample_method_percentile.to_numpy();joint=parent&anycomp&cal&(pct>=.95)
obs['review_route_v4']=np.select([~parent&anycomp,~parent,joint,parent&anycomp],['weak_Mast_with_competing_program_review','low_Mast_marker_review','joint_doublet_candidate_review','Mast_with_competing_program_review'],default='marker_supported_Mast_coarse')
obs['reviewed_major_identity_v1']=np.where(obs.review_route_v4.eq('marker_supported_Mast_coarse'),'Mast','review_required')
obs['reviewed_annotation_level_v1']=np.where(obs.review_route_v4.eq('marker_supported_Mast_coarse'),'major_identity_supported_subtype_not_assigned','identity_or_QC_review')
obs['reviewed_fine_subtype_v1']='not_assigned'
for name in ['MCTC','MCT','cycling','IFN','HLAII']:
 p=next(p for p in spec['panels'] if p['id']==name);ids=[ix[g] for g in p['genes']];obs['positive_'+name+'_marker_program']=parent&((raw[:,ids]>=2).sum(axis=1)>=2)
obs.to_pickle(o/'Mast_annotation_review_v1_PRIVATE.pkl.gz')
trans=obs.groupby(['final_annotation_v7','review_route_v4','reviewed_major_identity_v1'],observed=True).size().reset_index(name='n_cells');assert trans.n_cells.sum()==45036;trans.to_csv(o/'old_V7_to_Mast_reviewed_transition.tsv',sep='\t',index=False)
routes=obs.review_route_v4.value_counts().rename_axis('route').reset_index(name='n_cells');routes.to_csv(o/'Mast_reviewed_route_counts.tsv',sep='\t',index=False)
# Marker table for the 15 new source-refit groups, with actual source availability.
nc=int(cl.max()+1);H=sp.csr_matrix((np.ones(fit.sum(),np.int32),(cl,np.arange(fit.sum()))),shape=(nc,fit.sum()));S=sp.csr_matrix((np.ones(fit.sum(),np.int32),(cl,gidx[fit])),shape=(nc,len(gses)))
den=np.asarray(S@av.T.astype(np.int32));nd=np.asarray(H@(raw[fit]>0).astype(np.int32));sums=np.asarray(H@cp[fit]);n=np.bincount(cl,minlength=nc)
rows=[]
for c in range(nc):
 for j,g in enumerate(genes):
  dd=int(den[c,j]);rows.append({'cluster':c,'gene':g,'n_cells':int(n[c]),'n_measured':dd,'coverage_fraction':float(dd/n[c]),'n_expressing':int(nd[c,j]),'pct_expressing_measured':float(nd[c,j]/dd) if dd else np.nan,'mean_cp10k_measured':float(sums[c,j]/dd) if dd else np.nan})
e=pd.DataFrame(rows);e.to_csv(o/'Mast_H15_expanded_measured_markers.tsv',sep='\t',index=False,na_rep='NA')
# Hand review based on full RNA DEGs, retained source structure and marker programs.
labels={0:'Mast (coarse identity)',1:'Mast; NFKB1/VEGFA program [state]',2:'Mast; heat-shock/IEG program [state]',3:'Mast; ITGA2B/RGS16 profile [review]',4:'Mast (coarse identity)',5:'Mast + T program [review]',6:'IFN-associated Mast [state]',7:'Mast; PHGR1 profile [review]',8:'Mast; Ig-variable RNA [review]',9:'Mast + myeloid program [review]',10:'Cycling Mast [state]',11:'Mast + stromal/endo program [review]',12:'Mast; Ig-constant RNA [review]',13:'Mast + B-cell program [review]',14:'Mast; KLF2/IGHG4 profile [review]'}
pd.DataFrame([{'cluster':k,'display_annotation':v,'annotation_level':'working_major_or_state_not_final_subtype'} for k,v in labels.items()]).to_csv(o/'Mast_H15_display_annotations.tsv',sep='\t',index=False)
cs=[]
for c in range(nc):
 t=obs.iloc[np.flatnonzero(fit)[cl==c]];z={'cluster':c,'display_annotation':labels[c],'n_cells':len(t),'n_sources':t.gse_id.nunique(),'n_joint_doublet_candidates':int((t.review_route_v4=='joint_doublet_candidate_review').sum()),'doublet_percentile_median':float(t.doublet_within_sample_method_percentile.median()),'decontx_contamination_median':float(t.decontx_contamination.median())}
 for name in comp:z['n_'+name+'_program_review']=int(t['extended_'+name+'_program_review'].sum())
 cs.append(z)
pd.DataFrame(cs).to_csv(o/'Mast_H15_extended_competitor_summary.tsv',sep='\t',index=False)
programs=[]
for name in ['MCTC','MCT','cycling','IFN','HLAII']:
 take=obs['positive_'+name+'_marker_program']&obs.review_route_v4.eq('marker_supported_Mast_coarse');t=obs[take];programs.append({'program':name,'n_marker_supported_Mast_cells':len(t),'n_sources':t.gse_id.nunique(),'n_sample_units':t.groupby(['gse_id','sampleID'],observed=True).ngroups,'status':'positive RNA marker program; not a mature subtype assignment'})
pd.DataFrame(programs).to_csv(o/'Mast_positive_state_programs_v1.tsv',sep='\t',index=False)
manifest={'status':'completed_extended_Mast_marker_review','UTC':datetime.now(timezone.utc).isoformat(),'n_all_cells_retained':len(obs),'n_source_refit_cells':int(fit.sum()),'n_refit_clusters':nc,'n_expanded_genes':len(genes),'n_marker_cluster_records':len(e),'route_counts':dict(zip(routes.route,routes.n_cells.astype(int))),'n_old_Mast_unresolved_now_marker_supported_coarse':int((obs.final_annotation_v7.eq('Mast_unresolved')&obs.review_route_v4.eq('marker_supported_Mast_coarse')).sum()),'old_V7_fields_preserved':True,'final_fine_subtypes_assigned':0,'limits':['Expanded competitor review corrects the earlier T/myeloid/stromal/epithelial-only panel scope.','Three detected genes >=1 UMI or two genes >=2 UMI flags an extended competitor program; it is review evidence, not proof of a doublet.','Mast major identity is marker-supported within the old V7 Mast candidate pool; outside-parent recruitment and final atlas release are not complete.','No corrected full-gene decontX expression exists; raw competitor RNA may be ambient.']}
(o/'complete.json').write_text(json.dumps(manifest,indent=2))
packet={'manifest':manifest,'H_complete':json.loads((h/'refit/complete.json').read_text()),'tables':{}}
for name in ['old_V7_to_Mast_reviewed_transition.tsv','Mast_reviewed_route_counts.tsv','Mast_H15_expanded_measured_markers.tsv','Mast_H15_display_annotations.tsv','Mast_H15_extended_competitor_summary.tsv','Mast_positive_state_programs_v1.tsv']:packet['tables'][name]=pd.read_csv(o/name,sep='\t').to_json(orient='records')
b=zlib.compress(json.dumps(packet).encode(),9);(o/'aggregate_I_transfer_v1.zlib').write_bytes(b)
print(json.dumps({'manifest':manifest,'packet_bytes':len(b),'packet_sha256':hashlib.sha256(b).hexdigest()}))
