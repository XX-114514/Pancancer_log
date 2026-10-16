"""Prespecified mixed-spot spatial RNA descriptors; no independence-dependent significance."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,time,gc
import numpy as np,pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import rankdata
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);P=Path('/tmp/GSE273952_spatial_RNA_descriptive_plan_v1.json');plan=json.loads(P.read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();assert sha(__file__)==plan['producer_sha256'];D=V/'GSE273952_frozen_programme_projection_20261016_v1';A=V/'GSE273952_frozen_programme_independent_audit_20261016_v1';G=V/'GSE273952_spatial_neighbourhood_preflight_20261016_v1';GA=V/'GSE273952_spatial_geometry_audit_20261016_v1'
for tick in range(480):
 if (A/'inventory.json').exists() and (GA/'inventory.json').exists():break
 for name in [A.name,GA.name,D.name]:
  e=V/(name+'.execution.json')
  if e.exists() and json.loads(e.read_text()).get('exit_code') not in [None,0]:raise RuntimeError('Upstream task failed '+name)
 time.sleep(30)
else:raise TimeoutError('Four-hour dependency wait')
ai=json.loads((A/'inventory.json').read_text());di=json.loads((D/'inventory.json').read_text());ga=json.loads((GA/'inventory.json').read_text());assert ai['status'].startswith('PASS_ALL71398') and ga['status'].startswith('PASS_ALL71398');z=np.load(G/'private_spatial_graph_and_partial_composition.npz',allow_pickle=False);assert sha(G/'private_spatial_graph_and_partial_composition.npz')==ga['private_geometry_sha256'];schema=json.loads((G/'private_region_and_partial_reference_schema.json').read_text());n=71398;Y=np.full((n,817),np.nan,float);raw=np.zeros((n,4),float);rna=np.zeros(n);det=np.zeros(n);seen=np.zeros(n,bool);ids=None
for info in di['private_outputs']:
 p=D/info['relative_file'];assert sha(p)==info['sha256'];a=np.load(p,allow_pickle=False);ix=a['spot_rows'];assert not seen[ix].any();seen[ix]=True;Y[ix]=a['scores'];raw[ix]=a['raw_NAMPT_ITGA5_ITGB1_DDIT3_counts'];rna[ix]=a['RNA_total'];det[ix]=a['detected_genes'];ids=a['programme_ids'] if ids is None else ids;assert np.array_equal(ids,a['programme_ids'])
assert seen.all();g=csr_matrix((np.ones(len(z['graph_indices'])),z['graph_indices'],z['graph_indptr']),shape=(n,n));degree=np.asarray(g.sum(1)).ravel();logRNA=np.log1p(rna);logdet=np.log1p(det);native=np.log1p(raw/rna[:,None]*10000);exposure=np.divide(g@native[:,0],degree,out=np.full(n,np.nan),where=degree>0);oldtotal=None
# Historical total is intentionally not present in normalised input: use only five explicitly partial fractions.
# Amend frozen sensitivity definition before any outcome: total abundance cannot be recovered from fractions.
assert plan['partial_composition_covariates']=='five legacy fractions, Epithelial fixed reference; total abundance unavailable and not fabricated'
models=['quality_author_regions','plus_incomplete_legacy6_composition','plus_spatial_polynomial_gradient'];morans=[];effects=[];designrows=[]
def standard(x):
 x=np.asarray(x,float);sd=x.std(0);return (x-x.mean(0))/np.where(sd>0,sd,1)
def dummy(values):
 levels=np.unique(values);return np.column_stack([(values==k).astype(float) for k in levels[1:]]) if len(levels)>1 else np.empty((len(values),0))
def geometry(X):
 if len(X)<200:return None,'TOO_FEW_SPOTS',np.nan
 u,s,v=np.linalg.svd(X,full_matrices=False);rank=int(np.sum(s>s[0]*max(X.shape)*np.finfo(float).eps));cond=float(s[0]/s[-1]) if s[-1]>0 else np.inf
 if rank<X.shape[1] or cond>=1e8:return None,'DESIGN_UNIDENTIFIABLE',cond
 if len(X)-rank<100:return None,'INSUFFICIENT_RESIDUAL_DF',cond
 return u,'ESTIMABLE_DESCRIPTIVE_ONLY',cond
for lib in np.unique(z['library_codes']):
 ix=np.flatnonzero((z['library_codes']==lib)&~z['folded_excluded']&(degree>0));nn=len(ix);sg=g[ix][:,ix];s0=float(sg.sum());region=np.column_stack([dummy(z['TLS_primary'][ix]),dummy(z['tumor_zone'][ix])]);quality=np.column_stack([logRNA[ix],logdet[ix]]);legacy=z['legacy6_major_fractions'][ix][:,[k for k,v in enumerate(schema['legacy_majors']) if v!='Epithelial']];xy=standard(np.column_stack([z['pixel_x'][ix],z['pixel_y'][ix]]));poly=np.column_stack([xy,xy[:,0]**2,xy[:,1]**2,xy[:,0]*xy[:,1]])
 for model in models:
  cont=quality.copy()
  if model==models[1]:cont=np.column_stack([cont,legacy])
  elif model==models[2]:cont=np.column_stack([cont,poly])
  X=np.column_stack([np.ones(nn),standard(cont),region]);u,status,cond=geometry(X);designrows.append(dict(library_code=int(lib),model=model,kind='Moran',spots=nn,columns=X.shape[1],status=status,condition_number=cond));good=np.all(np.isfinite(Y[ix]),axis=0);mi=np.full(817,np.nan);ystat=np.full(817,'UNMEASURABLE',dtype='U40')
  if u is not None:
   yy=Y[ix][:,good];res=yy-u@(u.T@yy);den=(res*res).sum(0);val=nn/s0*(res*(sg@res)).sum(0)/np.where(den>1e-16,den,np.nan);mi[good]=val;ystat[good]=np.where(den>1e-16,'ESTIMABLE_DESCRIPTIVE_ONLY','CONSTANT_RESPONSE')
  else:ystat[good]=status
  for j,key in enumerate(ids):morans.append(dict(library_code=int(lib),model=model,programme_id=str(key),spots=nn,Moran_I=float(mi[j]),status=str(ystat[j])))
  # Rank continuous quantities; categorical author labels remain dummy variables.
  lc=np.column_stack([cont,native[ix,:3]]);rx=np.column_stack([np.ones(nn),standard(rankdata(lc,axis=0,method='average')),region]);qu,st,cn=geometry(rx);designrows.append(dict(library_code=int(lib),model=model,kind='neighbour_LR_RNA',spots=nn,columns=rx.shape[1],status=st,condition_number=cn));x=rankdata(exposure[ix],method='average');x=x-qu@(qu.T@x) if qu is not None else x
  for j in [814,815,816]:
   val=np.nan;ss=st
   if qu is not None:
    y=rankdata(Y[ix,j],method='average');y=y-qu@(qu.T@y);den=float(np.linalg.norm(x)*np.linalg.norm(y))
    if den>1e-12:val=float(x@y/den)
    else:ss='CONSTANT_RESIDUAL_RESPONSE_OR_EXPOSURE'
   effects.append(dict(library_code=int(lib),model=model,endpoint=str(ids[j]),spots=nn,partial_Spearman_rho=val,status=ss))
 print(json.dumps(dict(library_ordinal=int(lib),processed_spots=nn)),flush=True);gc.collect()
O=V/'GSE273952_spatial_RNA_descriptors_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);pd.DataFrame(morans).to_csv(O/'private_all_programme_library_Moran.tsv',sep='\t',index=False);pd.DataFrame(effects).to_csv(O/'private_all_library_neighbour_LR_RNA.tsv',sep='\t',index=False);pd.DataFrame(designrows).to_csv(O/'private_all_model_designs.tsv',sep='\t',index=False)
def summarize(df,keys,value):
 out=[]
 for key,grp in df.groupby(keys,sort=True):
  key=key if isinstance(key,tuple) else (key,);v=grp[value].dropna().to_numpy();d=dict(zip(keys,key));d.update(libraries_total=len(grp),libraries_estimable=len(v),positive_libraries=int((v>0).sum()),negative_libraries=int((v<0).sum()),median=float(np.median(v)) if len(v) else np.nan,q25=float(np.quantile(v,.25)) if len(v) else np.nan,q75=float(np.quantile(v,.75)) if len(v) else np.nan,minimum=float(v.min()) if len(v) else np.nan,maximum=float(v.max()) if len(v) else np.nan,leave_one_library_median_min=float(min(np.median(np.delete(v,k)) for k in range(len(v)))) if len(v)>1 else np.nan,leave_one_library_median_max=float(max(np.median(np.delete(v,k)) for k in range(len(v)))) if len(v)>1 else np.nan,inference_pvalue_computed=False,cell_intrinsic_or_mechanism_validation=False);out.append(d)
 return pd.DataFrame(out)
summarize(pd.DataFrame(morans),['model','programme_id'],'Moran_I').to_csv(O/'public_all2451_programme_spatial_descriptors.tsv',sep='\t',index=False);summarize(pd.DataFrame(effects),['model','endpoint'],'partial_Spearman_rho').to_csv(O/'public_all9_neighbour_NAMPT_target_RNA_descriptors.tsv',sep='\t',index=False);out=dict(status='ALL19_LIBRARIES_FROZEN_SPATIAL_RNA_DESCRIPTORS_COMPLETE_PENDING_INDEPENDENT_AUDIT_NOT_MECHANISM',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),projection_inventory_sha256=sha(D/'inventory.json'),projection_audit_sha256=sha(A/'inventory.json'),geometry_audit_sha256=sha(GA/'inventory.json'),programmes=817,library_model_programme_records=len(morans),fixed_LR_RNA_records=len(effects),model_designs=len(designrows),pvalues_computed=False,patient_independence_reviewed=False,legacy6_composition_incomplete=True,formal_spatial_cellstate_or_mechanism_validation=False,public_outputs={p.name:sha(p) for p in O.glob('public_*.tsv')});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
