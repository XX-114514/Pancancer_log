"""Original fixed-gene, alternate percentile, edge-sum and least-squares replay."""
import argparse,json
from pathlib import Path
from datetime import datetime,timezone
import h5py,numpy as np,pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha

def main(a):
 assert not a.output.exists();a.output.mkdir(mode=0o700)
 term=json.loads((a.producer/'terminal_public.json').read_text());assert term['n_slots']==76
 for n,s in term['files'].items():assert sha(a.producer/n)==s
 bound=json.loads((a.producer/'bound_plan_before_regional_effects.json').read_text());V=Path(json.loads((a.root/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
 I=V/'existing_GSE273952_Visium_raw_preflight_20261016_v2';A=V/'existing_GSE273952_Visium_input_audit_20261016_v1';G=V/'GSE273952_spatial_neighbourhood_preflight_20261016_v1';GA=V/'GSE273952_spatial_geometry_audit_20261016_v1'
 assert sha(A/'inventory.json')==bound['raw_input_audit_sha256'] and sha(GA/'inventory.json')==bound['geometry_audit_sha256']
 q=np.load(A/'private_all_spot_derived_QC_and_grid.npz',allow_pickle=False);z=np.load(G/'private_spatial_graph_and_partial_composition.npz',allow_pickle=False)
 lookup=json.loads((I/'private_native_gene_and_library_schema.json').read_text())['feature_lookup'];legacy_names=json.loads((G/'private_region_and_partial_reference_schema.json').read_text())['legacy_majors']
 fp=a.root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';assert sha(fp)==bound['fixed_panel_sha256'];panels=json.loads(fp.read_text())['panels'];genes=sorted(set(sum(panels.values(),[])))
 old=a.root.parents[2]/'Project_v2/Downstream/0716_Data/immunotherapy_4cohort_analysis/20260803_v1';path=Path(json.loads((old/'outputs/spatial/GSE273952_spatial_contract_audit.json').read_text())['input_h5ad']);assert sha(path)==bound['input_sha256']
 n=len(q['RNA_total']);raw=np.zeros((n,34),np.int64)
 with h5py.File(path,'r') as h:
  x=h['layers/counts'];nf=int(x.attrs['shape'][1]);physical=np.full(nf,-1,np.int32)
  for j,gene in enumerate(genes):assert physical[lookup[gene]]==-1;physical[lookup[gene]]=j
  for start in range(0,n,4096):
   stop=min(n,start+4096);p=x['indptr'][start:stop+1];lo,hi=int(p[0]),int(p[-1]);c=x['indices'][lo:hi];v=x['data'][lo:hi];rr=np.repeat(np.arange(stop-start),np.diff(p));mapped=physical[c];keep=mapped>=0;block=np.zeros((stop-start,34),np.int64);np.add.at(block,(rr[keep],mapped[keep]),v[keep]);raw[start:stop]=block
 saved=np.load(a.producer/'private_fixed34_counts_and7_scores.npz',allow_pickle=False);assert np.array_equal(saved['genes'],genes) and np.array_equal(saved['raw_counts'],raw)
 degree=np.diff(z['graph_indptr']);score=np.full((n,7),np.nan);names=['regulatory_associated','cytotoxic','pan_T','interferon','proliferation','memory_associated','GZMK_effector']
 for lib in sorted(set(z['library_codes'])):
  ix=np.where((z['library_codes']==lib)&~z['folded_excluded']&(degree>0))[0];expr=raw[ix]/q['RNA_total'][ix,None]*10000;ranked=pd.DataFrame(expr).rank(method='average',pct=True).to_numpy()
  for j,k in enumerate(names):
   cols=[genes.index(g) for g in panels[k]]
   if len(ix) and (np.max(expr[:,cols],axis=0)>np.min(expr[:,cols],axis=0)).all():score[ix,j]=ranked[:,cols].mean(1)
 assert np.array_equal(np.isnan(score),np.isnan(saved['scores'])) and np.allclose(score,saved['scores'],rtol=0,atol=1e-14,equal_nan=True)
 edge_row=np.repeat(np.arange(n),degree);edge_col=z['graph_indices'];sums=np.bincount(edge_row,weights=score[edge_col,0],minlength=n);exposure=np.divide(sums,degree,out=np.full(n,np.nan),where=degree>0)
 old_rows=pd.read_csv(a.producer/'private_all76_regional_descriptors.tsv',sep='\t',keep_default_na=False).set_index(['library_code','model','contrast']);checked=0;errors=[]
 for lib in sorted(set(z['library_codes'])):
  ix=np.where((z['library_codes']==lib)&~z['folded_excluded']&(degree>0))[0];nn=len(ix)
  region=[]
  for field in ['TLS_primary','tumor_zone']:
   cats=sorted(set(z[field][ix]));region.append(pd.get_dummies(pd.Categorical(z[field][ix],categories=cats),drop_first=True,dtype=float).to_numpy())
  region=np.column_stack(region);xy=np.column_stack([z['pixel_x'][ix],z['pixel_y'][ix]]);sd=xy.std(0);xy=(xy-xy.mean(0))/np.where(sd>0,sd,1)
  base=np.column_stack([score[ix,2:5],np.log1p(q['RNA_total'][ix]),np.log1p(q['detected_genes'][ix]),xy,xy[:,0]**2,xy[:,1]**2,xy[:,0]*xy[:,1]])
  legacy=z['legacy6_major_fractions'][ix][:,[i for i,k in enumerate(legacy_names) if k!='Epithelial']]
  for model in ['quality_T_IFN_cycle_regions_xy','plus_incomplete_legacy_composition']:
   cont=np.column_stack([base,legacy]) if model.endswith('composition') else base
   for contrast in ['same_spot_R_C','neighbor_R_C_given_own_R']:
    values=cont if contrast=='same_spot_R_C' else np.column_stack([cont,score[ix,0]])
    X=np.column_stack([np.ones(nn),pd.DataFrame(values).rank(method='average').to_numpy(),region]);saved_row=old_rows.loc[(lib,model,contrast)];reason=''
    if nn<200:reason='fewer_than200_nonfolded_connected_spots'
    elif not(np.isfinite(X).all() and np.isfinite(score[ix,:5]).all() and np.isfinite(exposure[ix]).all()):reason='missing_fixed_panel_or_covariate'
    else:
     sd=X[:,1:].std(0);X[:,1:]=(X[:,1:]-X[:,1:].mean(0))/np.where(sd>0,sd,1);rank=np.linalg.matrix_rank(X);cond=np.linalg.cond(X)
     if rank<X.shape[1] or cond>=1e8:reason='full_required_design_not_identifiable'
     elif nn-rank<100:reason='fewer_than100_residual_df'
     else:
      xx=pd.Series(score[ix,0] if contrast=='same_spot_R_C' else exposure[ix]).rank(method='average').to_numpy();yy=pd.Series(score[ix,1]).rank(method='average').to_numpy();rx=xx-X@np.linalg.lstsq(X,xx,rcond=None)[0];ry=yy-X@np.linalg.lstsq(X,yy,rcond=None)[0];den=np.linalg.norm(rx)*np.linalg.norm(ry)
      if den<=1e-12:reason='constant_residual'
      else:
       rho=float(np.dot(rx,ry)/den);err=abs(rho-float(saved_row.partial_Spearman_rho));errors.append(err);assert err<1e-8
    assert saved_row.reason==reason and saved_row.status==('unavailable' if reason else 'estimable_descriptive_only');assert int(saved_row.n_spots)==nn;checked+=1
 assert checked==76
 receipt=dict(status='ALL76_REGIONAL_DESCRIPTORS_ORIGINAL_FIXED_RNA_AND_ALTERNATE_NUMERIC_REPLAY_ACCEPTED',finished_utc=datetime.now(timezone.utc).isoformat(),n_slots=76,raw_fixed_gene_values_re_read=int(raw.size),all7_panel_scores_compared=int(score.size),maximum_rho_difference=max(errors,default=0),producer_terminal_sha256=sha(a.producer/'terminal_public.json'),auditor_sha256=sha(__file__),no_P_or_FDR=True,not_patient_or_single_cell_validation=True)
 (a.output/'independent_audit_public.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))

if __name__=='__main__':
 p=argparse.ArgumentParser()
 for name in ['root','producer','output']:p.add_argument('--'+name,type=Path,required=True)
 main(p.parse_args())
