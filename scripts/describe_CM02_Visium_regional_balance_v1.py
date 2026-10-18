"""Fixed CM02 mixed-region RNA balance, complete libraries; descriptive only."""
import argparse,json
from pathlib import Path
from datetime import datetime,timezone
import h5py,numpy as np,pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import rankdata
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

def dummy(v):
 levels=np.unique(v);return np.column_stack([v==k for k in levels[1:]]).astype(float) if len(levels)>1 else np.empty((len(v),0))

def describe(scores,z,rna,det,g,legacy_names):
 degree=np.asarray(g.sum(1)).ravel();r=scores[:,0];exposure=np.divide(g@r,degree,out=np.full(len(r),np.nan),where=degree>0)
 rows=[]
 for lib in np.unique(z['library_codes']):
  ix=np.flatnonzero((z['library_codes']==lib)&~z['folded_excluded']&(degree>0));n=len(ix)
  region=np.column_stack([dummy(z['TLS_primary'][ix]),dummy(z['tumor_zone'][ix])]);xy=np.column_stack([z['pixel_x'][ix],z['pixel_y'][ix]]).astype(float)
  xy=(xy-xy.mean(0))/np.where(xy.std(0)>0,xy.std(0),1);poly=np.column_stack([xy,xy[:,0]**2,xy[:,1]**2,xy[:,0]*xy[:,1]])
  continuous=np.column_stack([scores[ix,2:5],np.log1p(rna[ix]),np.log1p(det[ix]),poly])
  legacy=z['legacy6_major_fractions'][ix][:,[j for j,k in enumerate(legacy_names) if k!='Epithelial']]
  for model in ['quality_T_IFN_cycle_regions_xy','plus_incomplete_legacy_composition']:
   cc=np.column_stack([continuous,legacy]) if model.endswith('composition') else continuous
   for contrast in ['same_spot_R_C','neighbor_R_C_given_own_R']:
    vals=cc if contrast=='same_spot_R_C' else np.column_stack([cc,r[ix]])
    X=np.column_stack([np.ones(n),rankdata(vals,axis=0,method='average'),region]);valid=np.isfinite(X).all() and np.isfinite(scores[ix,:5]).all() and np.isfinite(exposure[ix]).all()
    row=dict(library_code=int(lib),model=model,contrast=contrast,n_spots=n,n_columns=X.shape[1],status='unavailable',reason='',partial_Spearman_rho=np.nan)
    if n<200:row['reason']='fewer_than200_nonfolded_connected_spots'
    elif not valid:row['reason']='missing_fixed_panel_or_covariate'
    else:
     # Standardize continuous design, preserve binary author categories.
     sd=X[:,1:].std(0);X[:,1:]=(X[:,1:]-X[:,1:].mean(0))/np.where(sd>0,sd,1)
     u,s,v=np.linalg.svd(X,full_matrices=False);rank=int((s>s[0]*max(X.shape)*np.finfo(float).eps).sum());cond=float(s[0]/s[-1]) if s[-1]>0 else np.inf
     row.update(rank=rank,condition_number=cond)
     if rank<X.shape[1] or cond>=1e8:row['reason']='full_required_design_not_identifiable'
     elif n-rank<100:row['reason']='fewer_than100_residual_df'
     else:
      x=rankdata(r[ix] if contrast=='same_spot_R_C' else exposure[ix],method='average');y=rankdata(scores[ix,1],method='average');rx=x-u@(u.T@x);ry=y-u@(u.T@y);den=float(np.linalg.norm(rx)*np.linalg.norm(ry))
      if den<=1e-12:row['reason']='constant_residual'
      else:row.update(status='estimable_descriptive_only',reason='',partial_Spearman_rho=float(rx@ry/den))
    rows.append(row)
 return pd.DataFrame(rows)

def main(a):
 budget=resource(a.root,64*1024**2);assert not a.output.exists();a.output.mkdir(mode=0o700)
 plan=json.loads(a.plan.read_text());assert plan['n_slots']==76 and not plan['P_values']
 V=Path(json.loads((a.root/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
 I=V/'existing_GSE273952_Visium_raw_preflight_20261016_v2';A=V/'existing_GSE273952_Visium_input_audit_20261016_v1';G=V/'GSE273952_spatial_neighbourhood_preflight_20261016_v1';GA=V/'GSE273952_spatial_geometry_audit_20261016_v1'
 iv=json.loads((I/'inventory.json').read_text());ai=json.loads((A/'inventory.json').read_text());ga=json.loads((GA/'inventory.json').read_text());assert ai['status'].startswith('PASS_ALL71398') and ga['status'].startswith('PASS_ALL71398')
 assert sha(A/'private_all_spot_derived_QC_and_grid.npz')==ai['private_QC_sha256'];q=np.load(A/'private_all_spot_derived_QC_and_grid.npz',allow_pickle=False)
 assert sha(G/'private_spatial_graph_and_partial_composition.npz')==ga['private_geometry_sha256'];z=np.load(G/'private_spatial_graph_and_partial_composition.npz',allow_pickle=False)
 lookup=json.loads((I/'private_native_gene_and_library_schema.json').read_text())['feature_lookup'];schema=json.loads((G/'private_region_and_partial_reference_schema.json').read_text())
 fp=a.root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';panels=json.loads(fp.read_text())['panels'];genes=sorted(set(sum(panels.values(),[])));assert len(genes)==34 and all(g in lookup for g in genes)
 old=a.root.parents[2]/'Project_v2/Downstream/0716_Data/immunotherapy_4cohort_analysis/20260803_v1';path=Path(json.loads((old/'outputs/spatial/GSE273952_spatial_contract_audit.json').read_text())['input_h5ad']);assert sha(path)==iv['input_file_sha256']
 bound=dict(plan=plan,started_utc=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),plan_sha256=sha(a.plan),input_sha256=sha(path),fixed_panel_sha256=sha(fp),raw_input_audit_sha256=sha(A/'inventory.json'),geometry_audit_sha256=sha(GA/'inventory.json'),budget=budget)
 (a.output/'bound_plan_before_regional_effects.json').write_text(json.dumps(bound,indent=2)+'\n')
 n=len(q['RNA_total']);assert n==71398;raw=np.zeros((n,34),np.int64)
 with h5py.File(path,'r') as h:
  group=h['layers/counts'];nf=int(group.attrs['shape'][1])
  for start in range(0,n,2048):
   stop=min(n,start+2048);p=group['indptr'][start:stop+1];lo,hi=int(p[0]),int(p[-1]);block=csr_matrix((group['data'][lo:hi],group['indices'][lo:hi],p-p[0]),shape=(stop-start,nf));raw[start:stop]=block[:,[lookup[g] for g in genes]].toarray()
 score_names=['regulatory_associated','cytotoxic','pan_T','interferon','proliferation','memory_associated','GZMK_effector'];score=np.full((n,7),np.nan);coverage=[];rna=q['RNA_total'];det=q['detected_genes']
 g=csr_matrix((np.ones(len(z['graph_indices'])),z['graph_indices'],z['graph_indptr']),shape=(n,n));degree=np.asarray(g.sum(1)).ravel()
 for lib in np.unique(z['library_codes']):
  ix=np.flatnonzero((z['library_codes']==lib)&~z['folded_excluded']&(degree>0));expr=raw[ix]/rna[ix,None]*10000;ranked=rankdata(expr,axis=0,method='average')/len(ix) if len(ix) else expr
  for j,k in enumerate(score_names):
   ci=[genes.index(s) for s in panels[k]];variable=bool(len(ix)>0 and np.all(np.ptp(expr[:,ci],axis=0)>0));coverage.append(dict(library_code=int(lib),panel=k,all_genes_native_measured=True,all_genes_variable=variable,n_spots=len(ix)))
   if variable:score[ix,j]=ranked[:,ci].mean(1)
 result=describe(score,z,rna,det,g,schema['legacy_majors']);assert len(result)==76
 result.to_csv(a.output/'private_all76_regional_descriptors.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(a.output/'private_all133_library_panel_measurement.tsv',sep='\t',index=False)
 np.savez_compressed(a.output/'private_fixed34_counts_and7_scores.npz',raw_counts=raw,genes=np.array(genes,str),scores=score,score_names=np.array(score_names,str))
 terminal=dict(status='ALL76_FIXED_REGIONAL_DESCRIPTORS_COMPLETE_INDEPENDENT_REPLAY_PENDING',finished_utc=datetime.now(timezone.utc).isoformat(),n_slots=len(result),n_spots=n,no_P_or_FDR=True,not_single_cell_contact=True,files={p.name:sha(p) for p in a.output.iterdir() if p.is_file()});(a.output/'terminal_public.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps({k:v for k,v in terminal.items() if k!='files'}))

if __name__=='__main__':
 p=argparse.ArgumentParser()
 for name in ['root','plan','output']:p.add_argument('--'+name,type=Path,required=True)
 main(p.parse_args())
