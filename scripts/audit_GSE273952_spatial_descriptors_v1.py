"""Independent least-squares replay of all fixed spatial descriptors."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,time
import numpy as np,pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import rankdata
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);D=V/'GSE273952_spatial_RNA_descriptors_20261016_v1';P=V/'GSE273952_frozen_programme_projection_20261016_v1';G=V/'GSE273952_spatial_neighbourhood_preflight_20261016_v1';sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for tick in range(480):
 if (D/'inventory.json').exists():break
 e=V/(D.name+'.execution.json')
 if e.exists() and json.loads(e.read_text()).get('exit_code') not in [None,0]:raise RuntimeError('Spatial producer failed')
 time.sleep(30)
else:raise TimeoutError('Four-hour spatial descriptor wait')
di=json.loads((D/'inventory.json').read_text());assert di['pvalues_computed']==False;pi=json.loads((P/'inventory.json').read_text());z=np.load(G/'private_spatial_graph_and_partial_composition.npz',allow_pickle=False);schema=json.loads((G/'private_region_and_partial_reference_schema.json').read_text());n=71398;score=np.full((n,817),np.nan);RNA=np.zeros(n);det=np.zeros(n);counts=np.zeros((n,4));ids=None
for q in pi['private_outputs']:
 p=P/q['relative_file'];assert sha(p)==q['sha256'];s=np.load(p,allow_pickle=False);ii=s['spot_rows'];score[ii]=s['scores'];RNA[ii]=s['RNA_total'];det[ii]=s['detected_genes'];counts[ii]=s['raw_NAMPT_ITGA5_ITGB1_DDIT3_counts'];ids=s['programme_ids'] if ids is None else ids
adj=csr_matrix((np.ones(len(z['graph_indices'])),z['graph_indices'],z['graph_indptr']),shape=(n,n));degree=np.diff(z['graph_indptr']);norm=np.log1p(counts/RNA[:,None]*10000);exp=np.divide(adj@norm[:,0],degree,out=np.full(n,np.nan),where=degree>0);mt=pd.read_csv(D/'private_all_programme_library_Moran.tsv',sep='\t');lt=pd.read_csv(D/'private_all_library_neighbour_LR_RNA.tsv',sep='\t');dt=pd.read_csv(D/'private_all_model_designs.tsv',sep='\t');max_m=0.;max_r=0.;checked_m=0;checked_r=0;designs=0
models=['quality_author_regions','plus_incomplete_legacy6_composition','plus_spatial_polynomial_gradient']
def scale(a):
 a=np.asarray(a,float);s=a.std(axis=0);return (a-a.mean(axis=0))/np.where(s>0,s,1)
def dummies(a):
 cats=sorted(set(a));return np.column_stack([a==v for v in cats[1:]]).astype(float) if len(cats)>1 else np.empty((len(a),0))
def model_status(X):
 if len(X)<200:return 'TOO_FEW_SPOTS'
 if np.linalg.matrix_rank(X)<X.shape[1] or np.linalg.cond(X)>=1e8:return 'DESIGN_UNIDENTIFIABLE'
 if len(X)-X.shape[1]<100:return 'INSUFFICIENT_RESIDUAL_DF'
 return 'ESTIMABLE_DESCRIPTIVE_ONLY'
for lib in sorted(set(z['library_codes'])):
 ix=np.flatnonzero((z['library_codes']==lib)&~z['folded_excluded']&(degree>0));nn=len(ix);sub=adj[ix][:,ix];reg=np.column_stack([dummies(z['TLS_primary'][ix]),dummies(z['tumor_zone'][ix])]);base=np.column_stack([np.log1p(RNA[ix]),np.log1p(det[ix])]);fra=z['legacy6_major_fractions'][ix][:,[j for j,v in enumerate(schema['legacy_majors']) if v!='Epithelial']];xy=scale(np.column_stack([z['pixel_x'][ix],z['pixel_y'][ix]]));pol=np.column_stack([xy,xy[:,0]**2,xy[:,1]**2,xy[:,0]*xy[:,1]])
 for name in models:
  cont=base if name==models[0] else np.column_stack([base,fra if name==models[1] else pol]);X=np.column_stack([np.ones(nn),scale(cont),reg]);st=model_status(X);dr=dt[(dt.library_code==lib)&(dt.model==name)&(dt.kind=='Moran')];assert len(dr)==1 and dr.iloc[0].status==st;designs+=1;rows=mt[(mt.library_code==lib)&(mt.model==name)].set_index('programme_id');assert len(rows)==817;finite=np.isfinite(score[ix]).all(0)
  if st=='ESTIMABLE_DESCRIPTIVE_ONLY':
   y=score[ix][:,finite];r=y-X@np.linalg.lstsq(X,y,rcond=None)[0];ss=np.einsum('ij,ij->j',r,r);expected=nn/sub.sum()*np.einsum('ij,ij->j',r,sub@r)/np.where(ss>1e-16,ss,np.nan)
   for key,val,dd in zip(ids[finite],expected,ss):
    actual=rows.loc[key];assert actual.status==('ESTIMABLE_DESCRIPTIVE_ONLY' if dd>1e-16 else 'CONSTANT_RESPONSE');assert np.isclose(val,actual.Moran_I,atol=1e-8,equal_nan=True)
    if np.isfinite(val):max_m=max(max_m,abs(float(val)-float(actual.Moran_I)))
  else:assert (rows.loc[ids[finite]].status==st).all()
  assert (rows.loc[ids[~finite]].status=='UNMEASURABLE').all();checked_m+=817
  Xr=np.column_stack([np.ones(nn),scale(rankdata(np.column_stack([cont,norm[ix,:3]]),axis=0,method='average')),reg]);st=model_status(Xr);dr=dt[(dt.library_code==lib)&(dt.model==name)&(dt.kind=='neighbour_LR_RNA')];assert len(dr)==1 and dr.iloc[0].status==st;designs+=1
  for j in [814,815,816]:
   actual=lt[(lt.library_code==lib)&(lt.model==name)&(lt.endpoint==ids[j])];assert len(actual)==1;actual=actual.iloc[0];status=st;val=np.nan
   if st=='ESTIMABLE_DESCRIPTIVE_ONLY':
    yy=np.column_stack([rankdata(exp[ix]),rankdata(score[ix,j])]);res=yy-Xr@np.linalg.lstsq(Xr,yy,rcond=None)[0];den=float(np.linalg.norm(res[:,0])*np.linalg.norm(res[:,1]))
    if den>1e-12:val=float(res[:,0]@res[:,1]/den)
    else:status='CONSTANT_RESIDUAL_RESPONSE_OR_EXPOSURE'
   assert actual.status==status and np.isclose(val,actual.partial_Spearman_rho,atol=1e-8,equal_nan=True)
   if np.isfinite(val):max_r=max(max_r,abs(val-float(actual.partial_Spearman_rho)))
   checked_r+=1
 print('audited_library',int(lib),flush=True)
# Independently replay every cohort summary, including NA and all LOO medians.
for table,keys,value,public in [(mt,['model','programme_id'],'Moran_I','public_all2451_programme_spatial_descriptors.tsv'),(lt,['model','endpoint'],'partial_Spearman_rho','public_all9_neighbour_NAMPT_target_RNA_descriptors.tsv')]:
 tab=pd.read_csv(D/public,sep='\t').set_index(keys)
 for key,group in table.groupby(keys):
  row=tab.loc[key];v=group[value].dropna().to_numpy();assert row.libraries_total==len(group) and row.libraries_estimable==len(v) and row.positive_libraries==sum(v>0) and row.negative_libraries==sum(v<0)
  vals=[np.median(v),np.quantile(v,.25),np.quantile(v,.75),min(v),max(v)] if len(v) else [np.nan]*5
  loo=[np.median(np.delete(v,k)) for k in range(len(v))] if len(v)>1 else []
  vals.extend([min(loo),max(loo)] if loo else [np.nan,np.nan])
  for name,val in zip(['median','q25','q75','minimum','maximum','leave_one_library_median_min','leave_one_library_median_max'],vals):assert np.isclose(val,row[name],atol=1e-8,equal_nan=True)
O=V/'GSE273952_spatial_RNA_descriptors_audit_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);out=dict(status='PASS_ALL_FROZEN_SPATIAL_MORAN_LR_RNA_DESIGNS_EFFECTS_AND_COHORT_SUMMARIES_INDEPENDENT_LEAST_SQUARES',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(D/'inventory.json'),Moran_records_checked=checked_m,LR_RNA_records_checked=checked_r,model_designs_checked=designs,max_Moran_error=max_m,max_partial_rho_error=max_r,pvalues_or_patient_independence_claim=False,biological_validation=False);(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
