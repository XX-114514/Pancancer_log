"""Independent exact SOFT-derived subset full-X coefficients and resampling."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
from scipy.stats import t
from scipy.linalg import qr
from statsmodels.stats.multitest import multipletests
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def exact_int(value):
 v=float(value);assert np.isfinite(v) and v.is_integer();return int(v)
def main():
 start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';L=P/'CM_edge_LR_extraction_v2';O=P/'GSE161529_primary_patient_sensitivity_v1';M=P/'GSE161529_official_metadata_recon_v2';planpath=Path(os.environ['PANCANCER_BREAST_PLAN']);plan=json.loads(planpath.read_text());assert sha(__file__)==plan['auditor_sha256'];inv=json.loads((O/'inventory.json').read_text());assert inv['plan_sha256']==sha(planpath) and not (O/'independent_audit.json').exists()
 for n,h in inv['outputs'].items():assert sha(O/n)==h
 d=pd.read_csv(O/'all4712_primary_author_patient_models.tsv.gz',sep='\t',keep_default_na=False);assert len(d)==4712;meta=pd.read_csv(M/'private_all69_official_metadata_exact_mapping.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('sampleID');checked=0
 base=pd.read_csv(P/'CM_LR_full_control_association_v1/all26994_full_control_CM_LR_associations.tsv.gz',sep='\t');base=base[(base.gse_id=='GSE161529')&(base.variant=='endpoint_excluded_library')];assert len(base)==1178
 for variant,g in d.groupby('variant'):
  assert set(zip(g.case_id,g.cm))==set(zip(base.case_id,base.cm)) and not g.duplicated(['case_id','cm']).any()
 for gid,g in d.groupby('group_id'):
  b=dict(np.load(L/'private_group_arrays'/(gid+'.npz'),allow_pickle=False));mm=meta.loc[b['samples']];primary=(mm.official_material_role=='tumor_mammary_unresolved_repeat_specimens').to_numpy();lr=list(zip(b['ligands'],b['receptors']))
  for row in g.itertuples():
   j=lr.index((row.ligand,row.receptor));mask=primary&np.isfinite(b['scores'][:,j,0]);N=int(mask.sum());nendpoint=int(primary.sum());blocks=('GSE161529:tumor_author:'+mm.author_patient_token[mask]).to_numpy(dtype=str);G=len(set(blocks));assert exact_int(row.n_observed_primary_libraries)==N and exact_int(row.n_author_patient_blocks)==G
   if row.status=='INSUFFICIENT_PRIMARY_PATIENT_COVERAGE':assert nendpoint<20 or N<20 or G<10 or N/nendpoint<.8;continue
   if row.status=='NATIVE_SUBUNIT_NOT_MEASURED':assert not b['native_measured'][j];continue
   y0=np.log1p(b['scores'][mask,j,0])
   if row.status=='NO_RESPONSE_VARIATION':assert y0.std(ddof=1)<=1e-12;continue
   ci=list(b['CMs']).index(row.cm);activity=b['endpoint_excluded_activity' if row.variant.startswith('endpoint_excluded_library') else 'original_H'][mask,ci];raw=np.column_stack([activity,b['source_within_major'][mask],b['target_within_major'][mask],np.log1p(b['source_n'][mask]),np.log1p(b['target_n'][mask])]);cat=mm.author_cancer_type[mask];levels=sorted(set(cat))
   if row.variant.endswith('subtype') and len(levels)>1:raw=np.column_stack([raw]+[cat.eq(k).to_numpy(dtype=float) for k in levels[1:]])
   sd=raw.std(axis=0,ddof=1);X=np.column_stack([np.ones(N),np.divide(raw-raw.mean(axis=0),sd,out=np.zeros_like(raw),where=sd>1e-12)]);pinv=np.linalg.pinv(X);rank=np.linalg.matrix_rank(X);ident=np.diag(pinv@X)[1]
   if row.status=='CM_EXPOSURE_NOT_IDENTIFIABLE':assert ident<1-1e-6 or sd[0]<=1e-12;continue
   if row.status=='ILL_CONDITIONED_OR_HIGH_INFLUENCE_DESIGN':assert np.diag(X@pinv).max()>.95-1e-9;continue
   assert row.status=='ESTIMABLE_FULL_CONTROL_SPAN_ASSOCIATION' and ident>1-1e-7 and rank==exact_int(row.design_rank);y=np.log1p(b['scores'][mask,j,0]);y=(y-y.mean())/y.std(ddof=1);coef=pinv@y;err=y-X@coef;a=pinv[1];se=np.sqrt(sum((a[blocks==k]@err[blocks==k])**2 for k in set(blocks))*G/(G-1)*(N-1)/(N-rank));p=2*t.sf(abs(coef[1]/se),G-1);assert np.isclose(coef[1],float(row.beta),atol=1e-8,rtol=1e-6) and np.isclose(se,float(row.CR1_se),atol=1e-8,rtol=1e-6) and np.isclose(p,float(row.analytic_p),atol=1e-10,rtol=1e-6);checked+=1
 for variant,g in d.groupby('variant'):
  assert len(g)==1178;pp=pd.to_numeric(g.residual_block_permutation_p,errors='coerce');q=multipletests(pp.fillna(1),method='fdr_bh')[1];valid=pp.notna();assert np.allclose(pd.to_numeric(g.loc[valid,'full1178_BH_q']),q[valid])
 # Fresh explicit block-resampled full-X fit for every saved example.
 examples=[]
 for file in sorted((O/'private_resampling_examples').glob('*.npz')):
     b=dict(np.load(file,allow_pickle=False));raw=b['raw_predictors'];Y=b['raw_response'];blocks=b['patient_blocks'];N=len(raw);sd=raw.std(axis=0,ddof=1);X=np.column_stack([np.ones(N),np.divide(raw-raw.mean(axis=0),sd,out=np.zeros_like(raw),where=sd>1e-12)]);yy=(Y-Y.mean())/Y.std(ddof=1);Z=X[:,[0]+list(range(2,X.shape[1]))];fit0=Z@np.linalg.lstsq(Z,yy,rcond=None)[0];e0=yy-fit0;co=np.linalg.lstsq(X,yy,rcond=None)[0][1];order=list(dict.fromkeys(blocks));ids=[np.flatnonzero(blocks==k) for k in order];G=len(ids);bins={}
     for i,index in enumerate(ids):bins.setdefault(len(index),[]).append(i)
     rng=np.random.default_rng(int(b['resampling_seed']));perm=[]
     for draw in range(2000):
         vec=e0.copy()
         for members in bins.values():
             shuffle=rng.permutation(members)
             for target,source in zip(members,shuffle):vec[ids[target]]=e0[ids[source]]
         perm.append(np.linalg.lstsq(X,fit0+vec,rcond=None)[0][1])
     Qz,Rz,piv=qr(Z,mode='economic',pivoting=True);rz=np.linalg.matrix_rank(Z);basis=np.column_stack([Qz[:,:rz],X[:,1]])
     pp=(1+np.count_nonzero(abs(np.array(perm))>=abs(co)-1e-12))/2001;draws=rng.integers(0,G,size=(2000,G));bootstrap=[];rank=np.linalg.matrix_rank(X)
     for draw in draws:
         take=np.concatenate([ids[j] for j in draw]);xx=X[take];r=np.linalg.matrix_rank(xx)
         if r!=rank:continue
         pinv=np.linalg.pinv(xx);H=xx@pinv
         if np.diag(H).max()>.95:continue
         # Effective nonzero singular-value condition is a numerical guard.
         ev=np.linalg.eigvalsh(basis[take].T@basis[take]);
         if ev[0]<=1e-12*ev[-1] or ev[-1]/ev[0]>=1e16:continue
         bootstrap.append((pinv@yy[take])[1])
     row=d[(d.case_id==str(b['case_id']))&(d.cm==str(b['cm']))&(d.variant==str(b['variant']))].iloc[0];assert np.isclose(pp,float(row.residual_block_permutation_p));assert len(bootstrap)==exact_int(row.boot_valid),(file.name,len(bootstrap),row.boot_valid)
     if len(bootstrap)>=1800:assert np.allclose(np.quantile(bootstrap,[.025,.975]),[float(row.bootstrap_ci_low),float(row.bootstrap_ci_high)],rtol=1e-6,atol=1e-8)
     examples.append(dict(case_id=str(b['case_id']),cm=str(b['cm']),variant=str(b['variant']),status='PASS_EXPLICIT_FULLX2000_BLOCK_PERMUTATIONS_AND_CASE_BOOTSTRAPS',valid_bootstraps=len(bootstrap)))
 report=dict(status='PASS_ALL4712_PRIMARY_AUTHOR_PATIENT_SUBSET_MODELS_AND_EXPLICIT_RESAMPLING',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),estimated_models_independently_recomputed=checked,resampling_examples=examples,global_registry_changed=False,scientific_validation=False);(O/'independent_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
