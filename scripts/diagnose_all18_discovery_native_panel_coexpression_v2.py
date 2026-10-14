"""Fixed diagnostic panels on all accepted discovery cells; not population prevalence."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd,anndata as ad,zarr
from scipy.sparse import csr_matrix,vstack
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);P=Path('/tmp/all18_discovery_panel_coexpression_plan_v2.json');plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256'];A=R/'all_lineage_state_native_input_audit_20261014_v2/accepted_input_inventory.json';assert sha(A)==plan['accepted_input_inventory_sha256'];iv=json.loads(A.read_text());panels=plan['panels'];AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);U=AR/'20261010_uniform_expression_annotation_v1';assert sha(U/'input_contract.json')==plan['canonical_contract_sha256'];contract=json.loads((U/'input_contract.json').read_text());g=zarr.open_group(contract['core_zarr'],mode='r');xx=g['X'];shape=contract['raw_shape'];O=V/'all18_discovery_native_panel_coexpression_20261014_v2';assert not O.exists();O.mkdir();rows=[];joint=[];coverage=[];total=0;maxdiff=0
for item in iv['files']:
 assert sha(item['path'])==item['sha256'];a=ad.read_h5ad(item['path']);major=item['major'];assert a.obs.major.eq(major).all();total+=a.n_obs;tot=a.obs.full_RNA_total.to_numpy();values={};ok={}
 rr=a.obs.core_position.to_numpy(int);assert np.all(np.diff(rr)>0);parts=[];start=0
 while start<len(rr):
  end=min(start+500,len(rr),int(np.searchsorted(rr,rr[start]+5000,side='left')));take=rr[start:end];lo=int(take[0]);hi=int(take[-1])+1;ip=np.asarray(xx['indptr'][lo:hi+1],np.int64);aa,bb=int(ip[0]),int(ip[-1]);z=csr_matrix((xx['data'][aa:bb],xx['indices'][aa:bb],ip-aa),shape=(hi-lo,shape[1]))[take-lo].tocsr();z.sort_indices();assert z.has_canonical_format and np.isfinite(z.data).all() and (z.data>=0).all() and np.array_equal(z.data,np.floor(z.data));parts.append(z.astype(np.int64));start=end
 full=vstack(parts,format='csr');del parts;assert np.array_equal(np.asarray(full.sum(axis=1)).ravel(),tot)
 values={name:(np.zeros(a.n_obs,np.int64),np.zeros(a.n_obs),np.zeros(a.n_obs,bool)) for name in panels}
 for sid in sorted(a.obs.gse_id.unique()):
  inds=np.flatnonzero(a.obs.gse_id.eq(sid).to_numpy());sub=full[inds];gp=U/'by_gse'/sid/'genes.tsv';manifest=json.loads((U/'by_gse'/sid/'manifest.json').read_text());assert sha(gp)==manifest['output_sha256']['genes.tsv'];gt=pd.read_csv(gp,sep='\t').set_index('gene');allowed=np.zeros(shape[1],bool);allowed[gt.core_column.to_numpy(int)]=True;assert allowed[sub.indices].all()
  for name,genes in panels.items():
   mapped=gt.core_column.reindex(genes);present=mapped.notna().to_numpy();cols=mapped.dropna().to_numpy(int);coverage.append(dict(gse_id=sid,major=major,panel=name,required_genes=len(genes),measured_genes=int(present.sum()),complete_gene_axis=bool(present.all()),missing_genes=';'.join(np.array(genes)[~present])));x=sub[:,cols];det=np.asarray((x>0).sum(axis=1)).ravel();fraction=np.asarray(x.sum(axis=1)).ravel()/tot[inds];ok[sid,name]=bool(present.all());values[name][0][inds]=det;values[name][1][inds]=fraction;values[name][2][inds]=(det>=3)&(fraction>=.01)
   dd=np.zeros(len(inds),np.int64);ss=np.zeros(len(inds))
   for col in cols:
    v=sub[:,int(col)].toarray().ravel();dd+=(v>0);ss+=v
   assert np.array_equal(det,dd) and np.allclose(fraction,ss/tot[inds],atol=1e-12);maxdiff=max(maxdiff,float(np.max(np.abs(fraction-ss/tot[inds]),initial=0)))
 for keys,positions in a.obs.groupby(['gse_id','node_id','subCluster'],observed=True,sort=True).indices.items():
  idx=np.asarray(positions);sid,node,label=keys
  for name,(det,frac,strong) in values.items():rows.append(dict(gse_id=sid,major=major,node_id=node,subCluster=label,cells=len(idx),panel=name,all_panel_genes_measured=ok[sid,name],cells_detecting_atleast3_panel_genes=int((det[idx]>=3).sum()) if ok[sid,name] else np.nan,cells_atleast3genes_and1pct_native_RNA=int(strong[idx].sum()) if ok[sid,name] else np.nan,panel_native_RNA_fraction_median=float(np.median(frac[idx])) if ok[sid,name] else np.nan,diagnostic_marker_support_not_state_positive_rate=True,balanced_discovery_not_cohort_prevalence=True))
  for left,right in plan['joint_panels']:
   complete=ok[sid,left] and ok[sid,right];n=int((values[left][2][idx]&values[right][2][idx]).sum()) if complete else np.nan;joint.append(dict(gse_id=sid,major=major,node_id=node,subCluster=label,cells=len(idx),panel_A=left,panel_B=right,complete_gene_axes=complete,joint3genes_and1pct_each_cells=n,joint_fraction=n/len(idx) if complete else np.nan,not_contamination_or_doublet_diagnosis=True))
 print(json.dumps(dict(major=major,cells=a.n_obs,phase='PANELS_CHECKED',UTC=datetime.now(timezone.utc).isoformat())),flush=True);del a,full,sub;gc.collect()
assert total==100766
for name,tab in [('all_major_panel_gene_coverage.tsv',coverage),('all_source_node_panel_support.tsv',rows),('all_source_node_joint_panel_support.tsv',joint)]:pd.DataFrame(tab).to_csv(O/name,sep='\t',index=False)
out=dict(status='COMPLETE_ALL100766_NATIVE_PANEL_DIAGNOSTICS_NOT_STATE_LABELS',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),accepted_input_inventory_sha256=sha(A),cells=total,majors=18,panels=len(panels),source_node_panel_rows=len(rows),joint_panel_rows=len(joint),all_cell_panel_quantities_independently_accumulated=True,max_fraction_difference=maxdiff,full_cohort_prevalence=False,source_native_gene_axes=True,native_denominators_independently_recomputed=True,formal_state=False,original_labels_changed=False,output_sha256={p.name:sha(p) for p in O.glob('*.tsv')});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
