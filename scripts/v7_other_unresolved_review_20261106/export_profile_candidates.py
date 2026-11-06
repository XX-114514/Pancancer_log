import os,json,hashlib,time
from datetime import datetime,timezone
import numpy as np,pandas as pd,scipy.sparse as sp,zarr,h5py,anndata as ad
from anndata.io import read_elem
r=PROJECT_ROOT;run=r/'runs/20261106_v7_other_unresolved_review_v1';o=run/'stage_a_candidates';o.mkdir(exist_ok=False);policy=json.loads((run/'review_policy_v1.json').read_text());start=datetime.now(timezone.utc).isoformat()
v=pd.read_pickle(r/'runs/20260818_pancancer_5M_merge_prep/11_v7_cnv_freeze/20260905_v7_cnv_rerun_freeze_v1/checkpoints/final_annotation_v7.pkl.gz');g=zarr.open_group(str(r/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r');assert pd.Index(read_elem(g['obs/_index'])).equals(v.index);m=np.load(r/policy['prior_QC_mask']);scope=v.major_celltype_v7.astype(str).isin(policy['scope_majors']).to_numpy();rows=np.flatnonzero(scope&m['primary_keep']);assert len(v)==4676787
var=pd.DataFrame(index=pd.Index(read_elem(g['var/_index']),name='gene_symbol'));assert len(var)==95038 and var.index.is_unique;genes=policy['genes'];ids=var.index.get_indexer(genes);found=ids>=0
obs=v.iloc[rows][[k for k in ['gse_id','sampleID','cancer_type','cohort','cluster_uid','major_celltype_v7','final_annotation_v7','annotation_confidence_v7','major_resolution_reason_v2','celltypist_major','celltypist_confidence'] if k in v]].copy();obs['global_position']=rows
for key in ['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_hb','doublet_score','doublet_method','doublet_status','decontx_contamination','decontx_status','ev_celltypist_label','ev_celltypist_conf','scatomic_full_label','scatomic_full_confidence']:
 if key in g['obs']:obs[key]=np.asarray(read_elem(g['obs/'+key]))[rows]
# Keep old, full frozen QC-unit calibration comparable to previous Mast/T6.
cal=v[['gse_id','sampleID']].copy();cal['doublet_method']=read_elem(g['obs/doublet_method']);cal['doublet_score']=pd.to_numeric(read_elem(g['obs/doublet_score']),errors='coerce');gg=cal.groupby(['gse_id','sampleID','doublet_method'],observed=True,dropna=False)['doublet_score'];obs['doublet_within_sample_method_percentile']=gg.rank(pct=True,method='average').to_numpy()[rows];obs['doublet_calibration_n']=gg.transform('count').to_numpy()[rows];obs['doublet_calibration_evaluable']=(obs.doublet_calibration_n>=100)&obs.doublet_within_sample_method_percentile.notna();del cal,gg
ptr=g['X/indptr'][:];lo=ptr[rows];hi=ptr[rows+1];length=hi-lo;outptr=np.r_[0,np.cumsum(length)].astype(np.int64);data=np.empty(int(length.sum()),np.int32);indices=np.empty_like(data);blocks=[];a=0
for j in range(1,len(rows)):
 if hi[j]-lo[a]>4000000 or lo[j]-hi[j-1]>226095:blocks.append((a,j));a=j
blocks.append((a,len(rows)))
for bi,(a,b) in enumerate(blocks):
 left=int(lo[a]);right=int(hi[b-1]);dd=g['X/data'][left:right];ii=g['X/indices'][left:right]
 for k in range(a,b):data[outptr[k]:outptr[k+1]]=dd[int(lo[k]-left):int(hi[k]-left)];indices[outptr[k]:outptr[k+1]]=ii[int(lo[k]-left):int(hi[k]-left)]
 if bi%50==0:(o/'progress.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'blocks_done':bi+1,'blocks_total':len(blocks),'n_cells':len(rows)}))
assert (data>=0).all();X=sp.csr_matrix((data,indices,outptr),shape=(len(rows),95038));X.check_format(full_check=True);total=np.asarray(X.sum(axis=1)).ravel();obs['full_raw_library_total']=total;raw=np.zeros((len(rows),len(genes)),np.int32);raw[:,found]=X[:,ids[found]].toarray();gx={x:i for i,x in enumerate(genes)}
gses=sorted(obs.gse_id.astype(str).unique());available=[];per=r/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'
for gs in gses:
 with h5py.File(per/(gs+'_symbol_counts.h5ad'),'r') as f:source_var=pd.Index(read_elem(f['var']).index.astype(str))
 available.append(var.index.isin(source_var))
av=np.asarray(available);gi=pd.Index(gses).get_indexer(obs.gse_id.astype(str));marker_av=np.zeros((len(gses),len(genes)),bool);marker_av[:,found]=av[:,ids[found]];assert not ((raw>0)&~marker_av[gi]).any()
flags={}
for p in policy['panels']:
 for cut in [1,2]:flags[p['id']+'_core_positive_raw'+str(cut)]=np.logical_and.reduce([(raw[:,[gx[k] for k in group]]>=cut).any(axis=1) for group in p['core_groups']])
 flags[p['id']+'_core_evaluable']=np.logical_and.reduce([marker_av[:,[gx[k] for k in group]].any(axis=1)[gi] for group in p['core_groups']])
def hit(gs,cut=2):return (raw[:,[gx[k] for k in gs]]>=cut).any(axis=1)
def nhit(gs,cut=2):return (raw[:,[gx[k] for k in gs]]>=cut).sum(axis=1)
T=(hit(['CD3D','CD3E','CD3G'])&hit(['TRAC','TRBC1','TRBC2']))|(nhit(['CD3D','CD3E','CD3G'])>=2)
NK=hit(['KLRD1','NCR1'])&hit(['NKG7','GNLY','PRF1']);B=hit(['MS4A1','CD19','BANK1'])&hit(['CD79A','CD79B']);ASC=hit(['MZB1','DERL3'])&hit(['JCHAIN','SDC1']);MY=hit(['CSF1R','CD14','CD163','C1QA'])&hit(['LST1','MS4A7']);EN=hit(['VWF','CLDN5'])&hit(['PECAM1','KDR']);EP=hit(['EPCAM'])&hit(['KRT8','KRT18','KRT19']);FB=hit(['COL1A1','COL1A2'])&hit(['DCN','LUM'])
ery_struct=hit(['GYPA','SLC4A1','RHAG','GYPB','EPB42'])&hit(['AHSP','KLF1','ALAS2']);ery_nonHB=nhit(['AHSP','KLF1','ALAS2','GYPA','SLC4A1'])>=2;plate=hit(['PF4','PPBP'])&hit(['ITGA2B','GP9','GP1BA','GP1BB','TUBB1']);hspc=hit(['SPINK2','HLF','MEIS1','DNTT'])&hit(['CD34','PROM1','GATA2']);hspc_heme=hspc&hit(['PTPRC']);ilc_specific=flags['ILC2_core_positive_raw2']|flags['ILC3_core_positive_raw2']
for key,value in {'T_positive_source_aware':T,'NK_receptor_effector_positive':NK,'B_distinct_components_positive':B,'ASC_distinct_components_positive':ASC,'Myeloid_mature_components_positive':MY,'Endothelial_distinct_components_positive':EN,'Epithelial_distinct_components_positive':EP,'Fibroblast_distinct_components_positive':FB,'Ery_structural_synthesis_nonHB_positive':ery_struct,'Ery_nonHB_multigene_positive':ery_nonHB,'Platelet_secreted_structural_positive':plate,'HSPC_primitive_components_positive':hspc,'HSPC_hematopoietic_components_positive':hspc_heme,'ILC2_or_ILC3_specific_program_positive':ilc_specific}.items():flags[key]=value
obs=pd.concat([obs,pd.DataFrame(flags,index=obs.index)],axis=1);high=obs.doublet_calibration_evaluable.to_numpy()&(obs.doublet_within_sample_method_percentile.to_numpy()>=.95)
rare=np.stack([ery_nonHB,plate,hspc_heme,ilc_specific],axis=1);focus=policy['focus'];external=~obs.major_celltype_v7.astype(str).isin(focus).to_numpy();external_ok=external&(rare.sum(axis=1)==1)&~(T|B|ASC|MY);obs['model_pool']='not_in_focus_models';obs['recruitment_basis']='old_V7_candidate_parent'
for j,name in enumerate(focus):
 take=obs.major_celltype_v7.astype(str).eq(name).to_numpy()|(external_ok&rare[:,j]);obs.loc[take,'model_pool']=name;obs.loc[take&external,'recruitment_basis']='external_unique_positive_components_candidate_not_final_identity'
joint_ery=ery_struct&(T|B|ASC|MY|EN|EP|FB);joint_plate=plate&(T|B|ASC|MY|EN|EP|FB);joint_hspc=hspc_heme&(EN|EP|FB);joint_ilc=ilc_specific&~(T|NK)&(B|ASC|EN|EP|FB)
obs['joint_disjoint_compartment_QC_candidate']=high&np.select([obs.model_pool.eq('Erythroid'),obs.model_pool.eq('Platelet'),obs.model_pool.eq('HSPC'),obs.model_pool.eq('ILC')],[joint_ery,joint_plate,joint_hspc,joint_ilc],default=False)
obs.to_pickle(o/'candidate_marker_QC_evidence_PRIVATE.pkl.gz');marker=ad.AnnData(sp.csr_matrix(raw),obs=obs[['gse_id','major_celltype_v7','final_annotation_v7','global_position','full_raw_library_total']].copy(),var=pd.DataFrame(index=pd.Index(genes,name='gene_symbol')));marker.write_h5ad(o/'candidate_marker_raw_PRIVATE.h5ad',compression='gzip')
agg=[];me=[];cp=raw*(10000/np.maximum(total,1)[:,None])
for label,idx in obs.groupby('final_annotation_v7',observed=True).indices.items():
 t=obs.iloc[idx];row={'label_v7':str(label),'major':str(t.major_celltype_v7.iloc[0]),'n_cells':len(t),'n_sources':t.gse_id.nunique(),'n_joint_disjoint_QC_candidates':int(t.joint_disjoint_compartment_QC_candidate.sum())}
 for key in ['T_positive_source_aware','NK_receptor_effector_positive','Ery_structural_synthesis_nonHB_positive','Ery_nonHB_multigene_positive','Platelet_secreted_structural_positive','HSPC_hematopoietic_components_positive','ILC2_or_ILC3_specific_program_positive']:row['n_'+key]=int(t[key].sum())
 agg.append(row);den=marker_av[gi[idx]].sum(axis=0);nz=(raw[idx]>0).sum(axis=0);ss=cp[idx].sum(axis=0)
 for j,gene in enumerate(genes):me.append({'label_v7':str(label),'gene':gene,'n_cells':len(idx),'n_measured':int(den[j]),'coverage_fraction':float(den[j]/len(idx)),'n_expressing':int(nz[j]),'pct_expressing_measured':float(nz[j]/den[j]) if den[j] else np.nan,'mean_cp10k_measured':float(ss[j]/den[j]) if den[j] else np.nan})
pd.DataFrame(agg).to_csv(o/'old_label_positive_evidence_summary.tsv',sep='\t',index=False);pd.DataFrame(me).to_csv(o/'old_label_measured_marker_expression.tsv',sep='\t',index=False,na_rep='NA')
exports=[]
for name in focus:
 take=obs.model_pool.eq(name).to_numpy();folder=o/name;folder.mkdir();obs_sub=obs.iloc[np.flatnonzero(take)].copy();counts=X[take].tocsr();a=ad.AnnData(counts,obs=obs_sub,var=var.copy());a.uns['scope']='Old parent plus external unique marker-nominated candidates; diagnostic input before additional QC adjudication.';a.write_h5ad(folder/(name+'_candidate_full_raw.h5ad'),compression='gzip');sources=sorted(obs_sub.gse_id.astype(str).unique());coverage=av[pd.Index(gses).get_indexer(sources)].mean(axis=0);pd.DataFrame({'gene_symbol':var.index[coverage>=.8]}).to_csv(folder/'model_genes.tsv',sep='\t',index=False);pd.DataFrame({'gene_symbol':var.index[coverage==1]}).to_csv(folder/'common_measured_DE_genes.tsv',sep='\t',index=False);exports.append({'lineage':name,'n_cells':int(take.sum()),'n_old_parent':int((take&~external).sum()),'n_external_candidates':int((take&external).sum()),'n_sources':len(sources),'n_joint_disjoint_QC_candidates':int(obs_sub.joint_disjoint_compartment_QC_candidate.sum())})
out={'status':'completed_candidate_export_and_marker_QC_profile_not_identity_acceptance','UTC_start':start,'UTC_end':datetime.now(timezone.utc).isoformat(),'n_candidates':len(obs),'n_marker_genes':len(genes),'n_panels':len(policy['panels']),'n_read_blocks':len(blocks),'n_QC_joint_candidates':int(obs.joint_disjoint_compartment_QC_candidate.sum()),'exports':exports,'gene_axis':95038,'canonical_V7_changed':False};(o/'complete.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
