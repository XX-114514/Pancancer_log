"""Independent all-cell raw marker, component and source-denominator audit."""
import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd,zarr,scipy.sparse as sp,h5py
from anndata.io import read_elem
r=PROJECT_ROOT;run=r/'runs/20261106_v7_publication_marker_annotation_v1';a=run/'stage_a_full_axis';assert json.loads((run/'launch/full_axis_profile_terminal.json').read_text())['exit_code']==0
q=pd.read_pickle(a/'all_cells_identity_evidence_PRIVATE.pkl.gz');cur=pd.read_pickle(r/'runs/20261106_v7_B_unresolved_review_v1/stage_g_consolidated_annotation/current_annotation_all_cells_PRIVATE.pkl.gz');assert q.index.equals(cur.index) and len(q)==4676787
assert np.array_equal(q.current_QC_keep.to_numpy(),cur.current_primary_QC_keep.to_numpy())
genes=np.load(a/'gene_symbols.npy').tolist();gx={g:i for i,g in enumerate(genes)};x=np.load(a/'marker_counts_PRIVATE.npy',mmap_mode='r');tot=np.load(a/'full_raw_library_totals_PRIVATE.npy');policy=json.loads((run/'identity_policy.json').read_text());source=pd.read_csv(a/'source_axis.tsv',sep='\t').source.astype(str);av=np.load(a/'source_gene_measured.npy');si=pd.Index(source).get_indexer(q.gse_id.astype(str));assert (si>=0).all()
g=zarr.open_group(str(r/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r');var=pd.Index(read_elem(g['var/_index']));ii=var.get_indexer(genes);found=ii>=0;ip=g['X/indptr'][:];block=40000
for left in range(0,len(q),block):
 right=min(left+block,len(q));lo=int(ip[left]);hi=int(ip[right]);data=g['X/data'][lo:hi];indices=g['X/indices'][lo:hi];p=ip[left:right+1]-lo;raw=sp.csr_matrix((data,indices,p),shape=(right-left,len(var)));assert np.array_equal(np.asarray(raw.sum(axis=1)).ravel(),tot[left:right]);assert np.array_equal(raw[:,ii[found]].toarray(),x[left:right][:,found]);assert not x[left:right][:,~found].any();assert not ((x[left:right]>0)&~av[si[left:right]]).any()
 for j,rule in enumerate(policy['rules']):
  flags=np.ones(right-left,bool);evaluable=np.ones(right-left,bool)
  for comp in rule['groups']:
   jj=[gx[t] for t in comp];flags&=(x[left:right][:,jj]>=2).any(axis=1);evaluable&=av[si[left:right]][:,jj].any(axis=1)
  if rule['parent']=='Tcell':flags|=(x[left:right][:,[gx[t] for t in ['CD3D','CD3E','CD3G']]]>=2).sum(axis=1)>=2
  assert np.array_equal(flags,q['positive_'+str(j)].iloc[left:right].to_numpy());assert np.array_equal(evaluable,q['evaluable_'+str(j)].iloc[left:right].to_numpy())
 if left//block%10==0:(a/'audit_progress.json').write_text(json.dumps({'n_all_cells_raw_verified':right,'n_total':len(q),'UTC':datetime.now(timezone.utc).isoformat()}))
for i,gs in enumerate(source):
 with h5py.File(r/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad'),'r') as f:actual=pd.Index(read_elem(f['var']).index.astype(str))
 assert np.array_equal(pd.Index(genes).isin(actual),av[i])
m=pd.read_csv(a/'all_group_marker_expression.tsv',sep='\t');s=pd.read_csv(a/'all_group_component_summary.tsv',sep='\t');maxerr=0.;K=len(s);ids=pd.Index(s.audit_group);ci=ids.get_indexer(q.audit_group);assert (ci[q.current_QC_keep]>=0).all();count=np.zeros(K,np.int64);den=np.zeros((K,len(genes)),np.int64);nz=np.zeros_like(den);sumraw=np.zeros_like(den);sumcp=np.zeros_like(den,dtype=float)
for left in range(0,len(q),block):
 right=min(left+block,len(q));pos=np.arange(left,right)[q.current_QC_keep.iloc[left:right].to_numpy()];c=ci[pos];H=sp.csr_matrix((np.ones(len(pos),np.int64),(c,np.arange(len(pos)))),shape=(K,len(pos)));xx=x[pos];count+=np.bincount(c,minlength=K);den+=np.asarray(H@av[si[pos]].astype(np.int64));nz+=np.asarray(H@(xx>0).astype(np.int64));sumraw+=np.asarray(H@xx.astype(np.int64));sumcp+=np.asarray(H@(xx*(10000.0/np.maximum(tot[pos],1)[:,None])))
for z in m.to_dict('records'):
 i=ids.get_loc(z['audit_group']);j=gx[z['gene']];assert count[i]==z['n_cells'] and den[i,j]==z['n_measured'] and nz[i,j]==z['n_expressing'] and sumraw[i,j]==z['raw_UMI_sum'];assert abs(den[i,j]/count[i]-z['coverage_fraction'])<1e-12
 if den[i,j]:err=abs(sumcp[i,j]/den[i,j]-z['mean_cp10k_measured']);maxerr=max(maxerr,err);assert err<1e-7 and abs(nz[i,j]/den[i,j]-z['pct_expressing_measured'])<1e-12
assert count.sum()==4668312 and len(m)==K*len(genes)
out={'status':'PASS','UTC':datetime.now(timezone.utc).isoformat(),'n_all_cells_raw_verified':len(q),'n_QC_keep':int(count.sum()),'n_genes':len(genes),'n_groups':K,'n_marker_records':len(m),'all_original_counts_and_raw_library_totals_exact':True,'all_positive_components_and_evaluable_flags_recomputed':True,'all_source_measured_gene_axes_rechecked':True,'all_marker_statistics_independently_recomputed':True,'max_absolute_CP10K_error':maxerr};(a/'independent_profile_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
