"""Match every controlled refit to all original factors; no favorable-only matching."""
import json,hashlib,itertools
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.optimize import linear_sum_assignment

r=PROJECT_ROOT/'runs/20261107_v4_CM_selection_v3';base=PROJECT_ROOT/'runs/20261106_v7_selective_CoVarNet_paper_v2';out=r/'diagnostics';jobs=json.loads((r/'jobs.json').read_text());W=pd.read_csv(base/'CM_results/W_PRIVATE.tsv',sep='\t',index_col=0);A=pd.read_csv(base/'CM_results/A_PRIVATE.tsv',sep='\t',index_col=0);A.columns=W.columns;m=pd.read_csv(base/'CM_results/each_node.tsv',sep='\t');ed=pd.read_csv(base/'CM_results/each_edge.tsv',sep='\t');rows=[];files=[]
def cosine(a,b):
 a=a.to_numpy();b=b.to_numpy();return a.T@b/np.outer(np.linalg.norm(a,axis=0),np.linalg.norm(b,axis=0))
def optimal_value_dp(score):
 # Independent exact subset dynamic program checks scipy's optimal assignment objective.
 n=score.shape[0];v=np.full(1<<n,-np.inf);v[0]=0
 for mask in range(1<<n):
  i=mask.bit_count()
  if i==n:continue
  for j in range(n):
   if not(mask>>j&1):
    t=mask|(1<<j);v[t]=max(v[t],v[mask]+score[i,j])
 return v[-1]
def eset(g):return {tuple(sorted((z.subCluster1,z.subCluster2))) for z in g.itertuples()}
for job in jobs:
 folder=r/'fits'/job['job'];proof=json.loads((folder/'complete.json').read_text());assert proof['nrun']==30 and proof['status']=='PASS_controlled_K14_refit';terminal=json.loads((r/'launch'/(job['job']+'_terminal.json')).read_text());assert terminal['exit_code']==0
 w=pd.read_csv(folder/'W_PRIVATE.tsv',sep='\t',index_col=0);a=pd.read_csv(folder/'A_PRIVATE.tsv',sep='\t',index_col=0);assert a.shape==w.shape;a.columns=w.columns;w=w.reindex(W.index,fill_value=0);a=a.reindex(A.index,fill_value=0);nw=pd.read_csv(folder/'each_node.tsv',sep='\t');ne=pd.read_csv(folder/'each_edge.tsv',sep='\t');cw=cosine(W,w);ca=cosine(A,a);ii,jj=linear_sum_assignment(-cw);ia,ja=linear_sum_assignment(-ca);assert np.array_equal(ii,np.arange(14)) and np.array_equal(ia,np.arange(14));assert abs(cw[ii,jj].sum()-optimal_value_dp(cw))<1e-9 and abs(ca[ia,ja].sum()-optimal_value_dp(ca))<1e-9
 for i,cm in enumerate(W.columns):
  j=jj[i];factor=w.columns[j];orig=set(m.loc[m.cm.eq(cm),'subCluster']);found=set(nw.loc[nw.cm.eq(factor),'subCluster']);origedges=eset(ed[ed.cm.eq(cm)]);foundedges=eset(ne[ne.cm.eq(factor)]);matchedA=w.columns[ja[i]];foundA=set(nw.loc[nw.cm.eq(matchedA),'subCluster']);foundedgesA=eset(ne[ne.cm.eq(matchedA)]);top=set(W[cm].nlargest(5).index);newtop=set(w[factor].nlargest(5).index)
  rows.append({'cm':cm,'job':job['job'],'kind':job['kind'],'held_source':job['held_source'],'matched_factor_W':factor,'W_cosine':cw[i,j],'A_cosine_at_W_match':ca[i,j],'matched_factor_A':matchedA,'A_cosine_optimal':ca[i,ja[i]],'W_A_match_agree':bool(j==ja[i]),'top5_node_Jaccard':len(top&newtop)/len(top|newtop),'network_node_Jaccard':len(orig&found)/len(orig|found) if orig else np.nan,'original_edge_recall':len(origedges&foundedges)/len(origedges) if origedges else np.nan,'original_edge_recall_A_match':len(origedges&foundedgesA)/len(origedges) if origedges else np.nan,'n_matched_network_nodes':len(found),'n_matched_network_edges':len(foundedges),'n_variable_nodes':proof['n_nodes'],'n_samples':proof['n_samples'],'n_starts':proof['nrun'],'rank_reoptimized':False})
 files.append({'job':job['job'],'complete_sha256':hashlib.sha256((folder/'complete.json').read_bytes()).hexdigest(),'model_sha256':hashlib.sha256((folder/'model_PRIVATE.rds').read_bytes()).hexdigest(),'n_nodes':proof['n_nodes'],'removed_constant_nodes':proof['removed_constant_nodes']})
detail=pd.DataFrame(rows);detail.to_csv(out/'all_CM_refit_matches.tsv',sep='\t',index=False);agg=[]
for cm,g in detail.groupby('cm'):
 row={'cm':cm}
 for kind in ['seed','LOSO']:
  h=g[g.kind.eq(kind)];row[kind+'_n']=len(h)
  for name in ['W_cosine','A_cosine_at_W_match','A_cosine_optimal','top5_node_Jaccard','network_node_Jaccard','original_edge_recall','original_edge_recall_A_match']:
   row[kind+'_'+name+'_median']=h[name].median();row[kind+'_'+name+'_min']=h[name].min()
  row[kind+'_exact_network_fraction']=float(h.network_node_Jaccard.eq(1).mean()) if h.network_node_Jaccard.notna().any() else np.nan;row[kind+'_all_edges_recovered_fraction']=float(h.original_edge_recall.eq(1).mean()) if h.original_edge_recall.notna().any() else np.nan
 agg.append(row)
pd.DataFrame(agg).to_csv(out/'all_CM_refit_stability_summary.tsv',sep='\t',index=False);(out/'refit_complete.json').write_text(json.dumps({'status':'PASS_all_18_refits_540_starts_all_14_factors_compared','UTC':datetime.now(timezone.utc).isoformat(),'n_jobs':len(jobs),'n_matches':len(rows),'rank':14,'nrun_per_job':30,'all_assignments_independently_optimal_DP':True,'profiles':'W and A=W*S column-L2 cosine; missing constant nodes zero-filled','primary_match':'W-optimal one-to-one assignment; independently A-optimal assignment also reported','input_binding':files,'limitation':'Controlled fixed-rank sensitivity, not rank re-selection, independent cohort validation or patient-block inference'},indent=2));print(json.dumps({'status':'PASS_all_refits_compared','n_jobs':len(jobs),'n_factors':14,'n_matches':len(rows)}))
