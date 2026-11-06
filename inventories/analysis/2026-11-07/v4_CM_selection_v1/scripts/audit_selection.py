"""Independent array-level audit of all source, denominator and factor diagnostics."""
import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261107_v4_CM_selection_v3';base=PROJECT_ROOT/'runs/20261106_v7_selective_CoVarNet_paper_v2';out=r/'diagnostics'
assert json.loads((out/'independent_R_audit.json').read_text())['status'].startswith('PASS_')
f=pd.read_csv(base/'inputs/frequency_PRIVATE.tsv',sep='\t',index_col=0);units=pd.read_csv(base/'inputs/sample_admission_PRIVATE.tsv',sep='\t').set_index('sampleID');src=units.loc[f.columns[1:],'gse_id'].to_numpy();W=pd.read_csv(base/'CM_results/W_PRIVATE.tsv',sep='\t',index_col=0);S=pd.read_csv(base/'CM_results/S_PRIVATE.tsv',sep='\t',index_col=0).to_numpy();H=pd.read_csv(base/'CM_results/H_PRIVATE.tsv',sep='\t',index_col=0).to_numpy();A=W.to_numpy()@S;contribution=A.sum(axis=0)[:,None]*H;act=contribution/contribution.sum(axis=0);assert np.allclose(act.sum(axis=0),1)
savedact=pd.read_csv(out/'scale_invariant_activity_PRIVATE.tsv',sep='\t',index_col=0);assert np.allclose(act,savedact.to_numpy(),atol=1e-12)
scale=np.linspace(.13,12,14);scaled=(A*scale).sum(axis=0)[:,None]*(H/scale[:,None]);assert np.allclose(scaled,contribution)
F=f.iloc[:,1:].to_numpy();pos={x:i for i,x in enumerate(f.index)};den=pd.read_csv(base/'inputs/parent_denominators_PRIVATE.tsv',sep='\t',index_col=0).loc[:,f.columns[1:]];parent=f.majorCluster;edge=pd.read_csv(out/'all_global_edge_source_diagnostics.tsv',sep='\t');quality=pd.read_csv(out/'all_CM_quality_diagnostics.tsv',sep='\t');members=pd.read_csv(base/'CM_results/each_node.tsv',sep='\t');mass=pd.read_csv(out/'all_CM_source_activity.tsv',sep='\t')
def centered(v,groups):
 result=v.copy()
 for group in np.unique(groups):
  ix=groups==group;result[ix]-=np.mean(v[ix])
 return result
def cor(a,b):
 a=a-np.mean(a);b=b-np.mean(b);norm=np.sqrt((a*a).sum()*(b*b).sum());return float((a*b).sum()/norm) if len(a)>=4 and norm>1e-18 else np.nan
def same(a,b):assert (np.isnan(a) and np.isnan(b)) or np.isclose(a,b,atol=1e-10), (a,b)
for z in edge.itertuples():
 a=F[pos[z.node1]];b=F[pos[z.node2]];same(cor(a,b),z.raw_r);same(cor(centered(a,src),centered(b,src)),z.source_centered_r);use=(den.loc[parent[z.node1]].to_numpy()>=20)&(den.loc[parent[z.node2]].to_numpy()>=20);assert use.sum()==z.n_both_parent_ge20;same(cor(a[use],b[use]),z.r_both_parent_ge20);same(cor(centered(a[use],src[use]),centered(b[use],src[use])),z.source_centered_r_parent_ge20);loss=[cor(a[src!=g],b[src!=g]) for g in np.unique(src)];same(float(np.nanmin(loss)),z.LOSO_r_min)
for i,cm in enumerate(W.columns):
 q=quality[quality.cm.eq(cm)].iloc[0];ids=members.loc[members.cm.eq(cm),'subCluster'];ix=[pos[n] for n in ids];same(W.iloc[ix,i].sum()/W.iloc[:,i].sum(),q.network_W_mass_fraction);same(A[ix,i].sum()/A[:,i].sum(),q.network_effective_A_mass_fraction);v=act[i];res=centered(v,src);same(1-(res*res).sum()/((v-v.mean())**2).sum(),q.source_activity_R2);shares=[v[src==g].sum()/v.sum() for g in np.unique(src)];same(max(shares),q.max_source_activity_mass_share)
 for g in np.unique(src):
  row=mass[mass.cm.eq(cm)&mass.source.eq(g)].iloc[0];same(np.mean(v[src==g]),row.activity_mean);same(v[src==g].sum()/v.sum(),row.activity_mass_share)
matches=pd.read_csv(out/'all_CM_refit_matches.tsv',sep='\t');assert len(matches)==252 and matches.cm.nunique()==14
for kind,n in [('seed',5),('LOSO',13)]:assert matches[matches.kind.eq(kind)].groupby('cm').size().eq(n).all()
for job,g in matches.groupby('job'):
 assert g.matched_factor_W.nunique()==14 and g.matched_factor_A.nunique()==14
 folder=r/'fits'/job;w=pd.read_csv(folder/'W_PRIVATE.tsv',sep='\t',index_col=0).reindex(W.index,fill_value=0);a=pd.read_csv(folder/'A_PRIVATE.tsv',sep='\t',index_col=0).reindex(W.index,fill_value=0);a.columns=w.columns
 for z in g.itertuples():
  i=list(W.columns).index(z.cm);v=W.iloc[:,i].to_numpy();u=w[z.matched_factor_W].to_numpy();same(float((v*u).sum()/np.linalg.norm(v)/np.linalg.norm(u)),z.W_cosine);u=a[z.matched_factor_A].to_numpy();v=A[:,i];same(float((v*u).sum()/np.linalg.norm(v)/np.linalg.norm(u)),z.A_cosine_optimal)
receipt={'status':'PASS_independent_all14_source_denominator_activity_and_252_refit_match_checks','UTC':datetime.now(timezone.utc).isoformat(),'n_CM':14,'n_global_edges':len(edge),'n_refit_jobs':18,'n_starts':540,'n_refit_matches':252,'scale_invariant_activity_verified':True,'all_original_CM_definitions_unchanged':True,'no_postselection_or_biological_validation_claim':True,'files':[{'name':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file() and 'PRIVATE' not in p.name]};(out/'FINAL_INDEPENDENT_AUDIT.json').open('x').write(json.dumps(receipt,indent=2));print(json.dumps({k:v for k,v in receipt.items() if k!='files'}))
