"""Full V7 positive identity audit, keeping missing measurement separate from zero."""
import os,json,hashlib,time
from datetime import datetime,timezone
import numpy as np,pandas as pd,scipy.sparse as sp,zarr,h5py
from anndata.io import read_elem

r=PROJECT_ROOT;run=r/'runs/20261106_v7_publication_marker_annotation_v1';o=run/'stage_a_full_axis';o.mkdir(exist_ok=False);start=datetime.now(timezone.utc).isoformat()
policy=json.loads((run/'identity_policy.json').read_text());resource=json.loads((run/'marker_registry.json').read_text());rules=policy['rules']
genes=list(dict.fromkeys([g for p in resource['panels'] for g in p['positive_genes']]+[g for p in rules for group in p['groups'] for g in group]+['DNAJB1','HSPA1A','HSPA1B','FOS','JUN','ISG15','IFIT1','IFIT3','MX1','MKI67','TOP2A','UBE2C','CENPF','IGHD','IGHM','TCL1A','FCER2','CD27','AIM2','TNFRSF13B','FCRL2','BCL6','AICDA','RGS13','MEF2B','TBX21','ITGAX','FCRL5','LAMP3','FSCN1','CCR7','CST3']))
vpath=r/'runs/20260818_pancancer_5M_merge_prep/11_v7_cnv_freeze/20260905_v7_cnv_rerun_freeze_v1/checkpoints/final_annotation_v7.pkl.gz'
v=pd.read_pickle(vpath);cur=pd.read_pickle(r/'runs/20261106_v7_B_unresolved_review_v1/stage_g_consolidated_annotation/current_annotation_all_cells_PRIVATE.pkl.gz');assert cur.index.equals(v.index)
g=zarr.open_group(str(r/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r');assert pd.Index(read_elem(g['obs/_index'])).equals(v.index)
var=pd.Index(read_elem(g['var/_index']));ids=var.get_indexer(genes);look=np.full(len(var),-1,np.int32);look[ids[ids>=0]]=np.flatnonzero(ids>=0);N=len(v);G=len(genes);gx={s:i for i,s in enumerate(genes)}
np.save(o/'gene_symbols.npy',np.asarray(genes,dtype='U32'));raw=np.lib.format.open_memmap(o/'marker_counts_PRIVATE.npy',mode='w+',dtype=np.int32,shape=(N,G));total=np.zeros(N,np.int64);ptr=g['X/indptr'][:];block=20000
for left in range(0,N,block):
 right=min(N,left+block);lo=int(ptr[left]);hi=int(ptr[right]);ind=g['X/indices'][lo:hi];data=g['X/data'][lo:hi];ip=ptr[left:right+1]-lo;X=sp.csr_matrix((data,ind,ip),shape=(right-left,len(var)));assert (data>=0).all();total[left:right]=np.asarray(X.sum(axis=1)).ravel();raw[left:right]=X[:,ids[ids>=0]].toarray() if (ids>=0).all() else 0
 if not (ids>=0).all():raw[left:right,np.flatnonzero(ids>=0)]=X[:,ids[ids>=0]].toarray()
 if left//block%10==0:(o/'progress.json').write_text(json.dumps({'stage':'full_RNA_marker_scan','UTC':datetime.now(timezone.utc).isoformat(),'n_cells_done':right,'n_cells_total':N,'n_genes':G}))
raw.flush();np.save(o/'full_raw_library_totals_PRIVATE.npy',total)
source=v.gse_id.astype(str);sources=sorted(source.unique());si=pd.Index(sources).get_indexer(source);av=[]
for gs in sources:
 with h5py.File(r/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad'),'r') as f:available=pd.Index(read_elem(f['var']).index.astype(str))
 av.append(pd.Index(genes).isin(available))
av=np.asarray(av);np.save(o/'source_gene_measured.npy',av);pd.DataFrame({'source':sources}).to_csv(o/'source_axis.tsv',sep='\t',index=False)
for left in range(0,N,block):
 right=min(left+block,N);assert not ((raw[left:right]>0)&~av[si[left:right]]).any()
positive=np.zeros((N,len(rules)),bool);evaluable=np.zeros_like(positive)
for j,p in enumerate(rules):
 hits=[];ev=[]
 for group in p['groups']:
  ii=[gx[x] for x in group];hits.append((raw[:,ii]>=policy['raw_UMI_per_component_gene']).any(axis=1));ev.append(av[:,ii].any(axis=1)[si])
 positive[:,j]=np.logical_and.reduce(hits);evaluable[:,j]=np.logical_and.reduce(ev)
# Cell-level corroboration rules; shared T/NK effectors do not override positive T apparatus.
t_rule=next(j for j,x in enumerate(rules) if x['parent']=='Tcell');cd3=[gx[x] for x in ['CD3D','CD3E','CD3G']];positive[:,t_rule]|=(raw[:,cd3]>=2).sum(axis=1)>=2
q=v[['gse_id','cluster_uid','major_celltype_v7','final_annotation_v7','annotation_confidence_v7']].copy();q.index.name='cell_id';q['current_parent']=cur.current_parent.astype(object);q['current_QC_keep']=cur.current_primary_QC_keep.to_numpy();q['already_reviewed_broad']=cur.refit_identity_reviewed_this_cycle.to_numpy()&cur.identity_review_status.eq('broad_working').to_numpy();q['old_current_annotation_layer']=cur.annotation_layer.astype(object);q['diagnostic_model']=cur.diagnostic_model.astype(object);q['diagnostic_cluster']=cur.diagnostic_cluster
for key in ['doublet_score','doublet_method','doublet_status','decontx_contamination','ev_celltypist_label','scatomic_full_label']:
 if key in g['obs']:q[key]=read_elem(g['obs/'+key])
for j,p in enumerate(rules):q['positive_'+str(j)]=positive[:,j];q['evaluable_'+str(j)]=evaluable[:,j]
group=q.final_annotation_v7.astype(str)+'|L1_'+q.cluster_uid.astype(str)
modeled=cur.refit_identity_reviewed_this_cycle.to_numpy();group.loc[modeled]=q.loc[modeled,'diagnostic_model'].astype(str)+'|C_'+q.loc[modeled,'diagnostic_cluster'].astype(str)
q['audit_group']=group;cats=pd.Categorical(group);ci=cats.codes;names=cats.categories;K=len(names);take=q.current_QC_keep.to_numpy();counts=np.bincount(ci[take],minlength=K);keepids=np.flatnonzero(take)
H=sp.csr_matrix((np.ones(len(keepids),np.int32),(ci[take],np.arange(len(keepids)))),shape=(K,len(keepids)));S=sp.csr_matrix((np.ones(len(keepids),np.int32),(ci[take],si[take])),shape=(K,len(sources)));den=np.asarray(S@av.astype(np.int32));nz=np.zeros((K,G),np.int64);sums=np.zeros((K,G));sumraw=np.zeros((K,G),np.int64)
for left in range(0,N,block):
 right=min(left+block,N);use=take[left:right];positions=np.arange(left,right)[use];h=sp.csr_matrix((np.ones(len(positions),np.int32),(ci[positions],np.arange(len(positions)))),shape=(K,len(positions)));x=np.asarray(raw[positions]);nz+=np.asarray(h@(x>0).astype(np.int32));sumraw+=np.asarray(h@x.astype(np.int64));sums+=np.asarray(h@(x*(10000.0/np.maximum(total[positions],1)[:,None])))
marker=[];groups=[];refs=[];origin=[]
for k,name in enumerate(names):
 positions=np.flatnonzero((ci==k)&take);n=len(positions)
 if not n:continue
 t=q.iloc[positions];src=t.gse_id.astype(str).value_counts();z={'audit_group':name,'n_cells':n,'n_sources':len(src),'max_source_fraction':float(src.iloc[0]/n),'old_major_mode':str(t.major_celltype_v7.value_counts().index[0]),'current_parent_mode':str(t.current_parent.value_counts().index[0]),'n_already_reviewed_broad':int(t.already_reviewed_broad.sum()),'raw_library_median':float(np.median(total[positions]))}
 for j,p in enumerate(rules):z['n_positive_'+str(j)]=int(positive[positions,j].sum());z['n_evaluable_'+str(j)]=int(evaluable[positions,j].sum())
 groups.append(z)
 for j,gene in enumerate(genes):marker.append({'audit_group':name,'gene':gene,'n_cells':n,'n_measured':int(den[k,j]),'n_expressing':int(nz[k,j]),'coverage_fraction':float(den[k,j]/n),'pct_expressing_measured':float(nz[k,j]/den[k,j]) if den[k,j] else np.nan,'mean_cp10k_measured':float(sums[k,j]/den[k,j]) if den[k,j] else np.nan,'raw_UMI_sum':int(sumraw[k,j])})
 for field in ['ev_celltypist_label','scatomic_full_label']:
  if field in t:
   for label,number in t[field].astype(object).fillna('missing').astype(str).value_counts().head(5).items():refs.append({'audit_group':name,'reference_field':field,'label':label,'n_cells':int(number),'n_group':n})
 for label,number in t.final_annotation_v7.astype(str).value_counts().items():origin.append({'audit_group':name,'old_label':label,'n_cells':int(number)})
q['full_raw_library_total']=total;q.to_pickle(o/'all_cells_identity_evidence_PRIVATE.pkl.gz');pd.DataFrame(groups).to_csv(o/'all_group_component_summary.tsv',sep='\t',index=False);pd.DataFrame(marker).to_csv(o/'all_group_marker_expression.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(refs).to_csv(o/'all_group_reference_top5.tsv',sep='\t',index=False);pd.DataFrame(origin).to_csv(o/'all_group_old_label_counts.tsv',sep='\t',index=False)
out={'status':'completed_full_axis_positive_identity_profile_not_freeze_acceptance','UTC_start':start,'UTC_end':datetime.now(timezone.utc).isoformat(),'n_cells':N,'n_QC_keep':int(take.sum()),'n_genes':G,'n_rules':len(rules),'n_groups_with_QC_keep':len(groups),'n_marker_records':len(marker),'source_measured_zero_violations':True,'original_current_axis_preserved':True,'policy_sha256':hashlib.sha256((run/'identity_policy.json').read_bytes()).hexdigest(),'main_input':'20261106_v7_B_unresolved_review_v1/stage_g_consolidated_annotation/current_annotation_all_cells_PRIVATE.pkl.gz'}
(o/'complete.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
