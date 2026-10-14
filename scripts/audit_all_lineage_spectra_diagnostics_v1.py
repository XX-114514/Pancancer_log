"""Independent factor matching, scalar ranks, dummy regression and sparse incidence."""
from pathlib import Path
from datetime import datetime,timezone
from itertools import combinations
import os,json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr
from scipy.sparse import csr_matrix
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);block=os.environ['PANCANCER_STATE_AUDIT_BLOCK'];D=R/'all_lineage_Spectra_diagnostics_20261014_v1'/block;S=R/'all_lineage_Spectra_state_discovery_20261014_v1'/block;sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();iv=json.loads((D/'inventory.json').read_text());assert iv['status'].startswith('COMPLETE_BLOCK_MULTISEED');plan=json.loads(Path(os.environ['PANCANCER_STATE_AUDIT_PLAN']).read_text());assert sha(os.environ['PANCANCER_STATE_AUDIT_PLAN'])==iv['plan_sha256']
for name,h in iv['output_sha256'].items():assert sha(D/name)==h
for name,h in iv['input_seed_sha256'].items():assert sha(S.parent/name)==h
models=[]
for seed in plan['seeds']:
 with np.load(S/('seed'+str(seed))/'factors_and_discovery_scores.npz') as z:models.append({k:z[k].copy() for k in z.files})
meta=pd.read_csv(S/('seed'+str(plan['seeds'][0]))/'private_discovery_metadata.tsv.gz',sep='\t',index_col=0,keep_default_na=False);ref=models[0];assert np.array_equal(meta.index.to_numpy(str),ref['cell_ids']);scope=ref['factor_scope'];norm=[m['factors'].astype(float)/np.sqrt(np.einsum('ij,ij->i',m['factors'].astype(float),m['factors'].astype(float)))[:,None] for m in models];pairs=pd.read_csv(D/'all_factor_matches.tsv',sep='\t');score=pd.read_csv(D/'all_major_score_reproducibility.tsv',sep='\t');diag=pd.read_csv(D/'all_source_node_complexity_diagnostics.tsv',sep='\t');cand=pd.read_csv(D/'all_candidate_diagnostics.tsv',sep='\t',keep_default_na=False);maps={};maxdiff=0.;checks=0
for a,b in combinations(range(3),2):
 mm={}
 for s in sorted(set(scope)):
  ids=np.flatnonzero(scope==s);cos=np.dot(norm[a][ids],norm[b][ids].T);rows,cols=linear_sum_assignment(-cos)
  for i,j in zip(ids[rows],ids[cols]):mm[int(i)]=int(j)
 maps[a,b]=mm
 for i,j in mm.items():
  row=pairs[(pairs.seed_a==plan['seeds'][a])&(pairs.seed_b==plan['seeds'][b])&(pairs.factor_a==i)].iloc[0];assert row.factor_b==j;value=float(np.dot(norm[a][i],norm[b][j]));maxdiff=max(maxdiff,abs(value-row.cosine));assert np.isclose(value,row.cosine,atol=1e-10)
  ga=set(np.argsort(-models[a]['factors'][i],kind='stable')[:50]);gb=set(np.argsort(-models[b]['factors'][j],kind='stable')[:50]);assert np.isclose(row.top50_jaccard,len(ga&gb)/len(ga|gb));checks+=1
def dummy(v):return pd.get_dummies(pd.Series(v,dtype=str),dtype=float).to_numpy()
def resid(x,d):return x-d@np.linalg.lstsq(d,x,rcond=None)[0]
def rho(x,y):
 tx=1e-12*max(1,float(np.max(np.abs(x))));ty=1e-12*max(1,float(np.max(np.abs(y))))
 return float(spearmanr(x,y).statistic) if np.ptp(x)>tx and np.ptp(y)>ty else np.nan
def eq(x,y):return (np.isnan(x) and np.isnan(y)) or np.isclose(x,y,atol=1e-7)
cache={}
for major in sorted(meta.major.unique()):
 mask=meta.major.eq(major).to_numpy();md=meta.loc[mask];ds=dummy(md.gse_id.to_numpy());dn=dummy((md.gse_id.astype(str)+'::'+md.node_id.astype(str)).to_numpy());values=[m['cell_scores'][mask].astype(float) for m in models]
 cache[major]=dict(scores=values,source_resid=[resid(v,ds) for v in values],node_resid=[resid(v,dn) for v in values])
scorechecks=0
for row in score.itertuples():
 a=plan['seeds'].index(row.seed_a);b=plan['seeds'].index(row.seed_b);cc=cache[row.major]
 assert eq(rho(cc['scores'][a][:,row.factor_a],cc['scores'][b][:,row.factor_b]),row.pooled_score_rho) and eq(rho(cc['source_resid'][a][:,row.factor_a],cc['source_resid'][b][:,row.factor_b]),row.source_residual_score_rho) and eq(rho(cc['node_resid'][a][:,row.factor_a],cc['node_resid'][b][:,row.factor_b]),row.source_node_residual_score_rho);scorechecks+=1
matched=[np.arange(len(scope)),np.array([maps[0,1][i] for i in range(len(scope))]),np.array([maps[0,2][i] for i in range(len(scope))])];weights=np.median(np.stack([norm[i][matched[i]] for i in range(3)]),axis=0);genes=ref['genes'];tops=np.argsort(-weights,axis=1,kind='stable')[:,:50];ivraw=json.loads((R/'all_lineage_state_native_input_audit_20261014_v2/accepted_input_inventory.json').read_text());look={r['major']:r for r in ivraw['files']};support={};dchecks=0
for major in sorted(meta.major.unique()):
 item=look[major];assert sha(item['path'])==item['sha256'];raw=ad.read_h5ad(item['path']);mask=meta.major.eq(major).to_numpy();md=meta.loc[mask];assert np.array_equal(raw.obs_names,md.index);ids=np.flatnonzero((scope=='block_global')|(scope==major));rows=[];cols=[]
 for j,factor in enumerate(ids):
  names=genes[tops[factor]];ii=raw.var_names.get_indexer(names);assert (ii>=0).all();rows.extend(ii.tolist());cols.extend([j]*len(ii))
 incidence=csr_matrix((np.ones(len(rows),np.int32),(rows,cols)),shape=(raw.n_vars,len(ids)));det=((raw.X>0).astype(np.int32)@incidence).toarray();fraction=(raw.X@incidence).toarray()/md.full_RNA_total.to_numpy()[:,None];y=ref['cell_scores'][mask][:,ids].astype(float);ds=dummy(md.gse_id.to_numpy());dn=dummy((md.gse_id.astype(str)+'::'+md.node_id.astype(str)).to_numpy());depth=np.column_stack([np.log1p(md.full_RNA_total),np.log1p(md.full_detected_genes)]);ys=resid(y,ds);yn=resid(y,dn);dr=resid(depth,dn);full=resid(y,np.column_stack([dn,depth]));den=((y-y.mean(axis=0))**2).sum(axis=0);rs=np.divide(den-(ys**2).sum(axis=0),den,out=np.full(den.shape,np.nan),where=den>0);rn=np.divide(den-(yn**2).sum(axis=0),den,out=np.full(den.shape,np.nan),where=den>0);extra=np.divide((yn**2).sum(axis=0)-(full**2).sum(axis=0),den,out=np.full(den.shape,np.nan),where=den>0)
 for j,factor in enumerate(ids):
  row=diag[(diag.major==major)&(diag.anchor_factor==factor)].iloc[0];n=int((det[:,j]>=5).sum());support[int(factor),major]=n;assert row.informative_gene_support_cells==n and np.isclose(row.signature_RNA_fraction_median,np.median(fraction[:,j]));assert eq(rs[j],row.source_r2) and eq(rn[j],row.source_node_r2) and eq(extra[j],row.RNA_detected_increment_after_source_node);assert eq(rho(yn[:,j],dr[:,0]),row.RNA_partial_rank_after_source_node) and eq(rho(yn[:,j],dr[:,1]),row.detected_partial_rank_after_source_node);dchecks+=1
for r in cand.itertuples():
 assert list(genes[tops[r.anchor_factor]])==r.primary50_median_L2_gene_weights.split(';')
 i=r.anchor_factor;js=[int(x[i]) for x in matched];cyc=maps[1,2][js[1]]==js[2];assert cyc==r.cycle_consistent;pp=[p for p in pairs.itertuples() if (p.factor_a==js[plan['seeds'].index(p.seed_a)])];pp=[p for p in pp if p.scope==r.scope];mincos=min(p.cosine for p in pp);minjac=min(p.top50_jaccard for p in pp);assert np.isclose(mincos,r.min_weight_cosine) and np.isclose(minjac,r.min_top50_jaccard);weightok=cyc and mincos>=.9 and minjac>=.5;majorrows=json.loads(r.all_major_score_support_diagnostics);supported=[]
 for x in majorrows:
  values=[]
  for a,b in combinations(range(3),2):
   q=score[(score.seed_a==plan['seeds'][a])&(score.seed_b==plan['seeds'][b])&(score.factor_a==js[a])&(score.major==x['major'])];assert len(q)==1;values.append(q.source_residual_score_rho.iloc[0])
  value=min(values) if np.isfinite(values).all() else np.nan;assert eq(value,x['min_rho']);ok=support[i,x['major']]>=20 and np.isfinite(value) and value>=.7;assert ok==x['supported']
  if ok:supported.append(x['major'])
 assert set(supported)==set(filter(None,r.supported_majors.split(';')));context=bool(weightok and supported);shared=bool(weightok and r.scope=='block_global' and len(supported)>=2);assert r.stable_context_candidate==context and r.stable_shared_candidate==shared and r.stable_shared_or_lineage_candidate==(shared if r.scope=='block_global' else context);assert r.formal_state==False
receipt=dict(status='PASS_ALL_BLOCK_FACTORS_MATCHES_SCORE_SOURCE_NODE_DEPTH_SUPPORT_AND_CANDIDATE_RULES',block=block,observed_UTC=datetime.now(timezone.utc).isoformat(),audit_code_sha256=sha(__file__),producer_inventory_sha256=sha(D/'inventory.json'),factor_pair_checks=checks,score_pair_checks=scorechecks,major_factor_support_and_regression_checks=dchecks,candidate_rules=len(cand),max_cosine_difference=maxdiff,scalar_and_sample_p_values_not_used=True,source_distributions_not_independently_recomputed_here=True,formal_states=False)
A=R/'all_lineage_Spectra_diagnostic_audit_20261014_v1'/block;assert not A.exists();A.mkdir(parents=True);(A/'independent_audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
