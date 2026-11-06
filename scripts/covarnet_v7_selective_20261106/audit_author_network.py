"""Independently reconstruct Pearson, author P adjustment, specificity and CM edges."""
import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.stats import t as tdist
r=RUN_ROOT;o=r/'CM_results';cfg=json.loads((r/'config/paper_config.json').read_text());proof=json.loads((o/'independent_R_audit.json').read_text());assert proof['status']=='PASS_independent_R_saved_models_CCC_WSH_and_permutations';completion=json.loads((o/'complete.json').read_text())
def adjust(p,method):
 p=np.asarray(p);n=len(p);ix=np.argsort(p,kind='stable');s=p[ix]
 if method=='holm':z=np.minimum(1,np.maximum.accumulate(s*np.arange(n,0,-1)))
 else:z=np.minimum(1,np.minimum.accumulate((s*n/np.arange(1,n+1))[::-1])[::-1])
 a=np.empty(n);a[ix]=z;return a
result={'UTC':datetime.now(timezone.utc).isoformat(),'all_atlas_or_biological_validation_claim':False}
if completion['status']=='NO_ADMISSIBLE_RANK_UNDER_PAPER_CCC_RULE':
 assert proof.get('selected_rank') in [None,[]];result.update(status='PASS_no_admissible_paper_rank_no_CM_forced',n_rank_models=19)
else:
 fq=pd.read_csv(r/'inputs/frequency_PRIVATE.tsv',sep='\t',index_col=0);labels=pd.Index(fq.index);parent=fq.majorCluster.astype(str);X=fq.iloc[:,1:].to_numpy(float);R=np.corrcoef(X);ii,jj=np.triu_indices(len(labels),1);rr=np.clip(R[ii,jj],-1,1);stat=np.abs(rr)*np.sqrt((X.shape[1]-2)/np.maximum(1-rr**2,0));raw=2*tdist.sf(stat,df=X.shape[1]-2);holm=adjust(raw,'holm');fdr=adjust(holm,'BH');rawbh=adjust(raw,'BH');spe=np.array([np.mean(np.r_[R[i,np.arange(len(labels))!=i],R[(np.arange(len(labels))!=i)&(np.arange(len(labels))!=j),j]]<=R[i,j]) for i,j in zip(ii,jj)])
 expected={(labels[i],labels[j]):(rr[z],raw[z],holm[z],fdr[z],rawbh[z],spe[z]) for z,(i,j) in enumerate(zip(ii,jj))};pairs=pd.read_csv(o/'all_pairs_author_and_raw_BH.tsv',sep='\t');assert len(pairs)==len(expected);err=0.
 for z in pairs.to_dict('records'):
  i=labels.get_loc(z['subCluster1']);j=labels.get_loc(z['subCluster2']);a,b=sorted([i,j]);v=expected[(labels[a],labels[b])];err=max(err,abs(v[0]-z['correlation']));assert abs(v[0]-z['correlation'])<1e-10 and abs(v[5]-z['spe'])<1e-12
  for actual,want in [(z['raw_Pearson_p'],v[1]),(z['author_pval'],v[2]),(z['author_pval_fdr'],v[3]),(z['raw_Pearson_BH'],v[4])]:assert np.isclose(actual,want,atol=1e-10,rtol=1e-6),(actual,want)
 cutoff=1-(2*(cfg['top_n']-1)-1)/(2*(len(labels)-1)-1);edges={(labels[i],labels[j]) for z,(i,j) in enumerate(zip(ii,jj)) if rr[z]>cfg['correlation_cutoff'] and fdr[z]<=cfg['FDR'] and spe[z]>=cutoff};save=pd.read_csv(o/'global_edge.tsv',sep='\t');edgekey=lambda a,b:tuple(sorted([a,b],key=lambda n:labels.get_loc(n)));assert {edgekey(z.subCluster1,z.subCluster2) for z in save.itertuples()}==edges
 W=pd.read_csv(o/'W_PRIVATE.tsv',sep='\t',index_col=0);top=pd.read_csv(o/'factor_top15_candidates.tsv',sep='\t');each_node=pd.read_csv(o/'each_node.tsv',sep='\t');each_edge=pd.read_csv(o/'each_edge.tsv',sep='\t');summary=pd.read_csv(o/'CM_network_topology_and_permutation.tsv',sep='\t')
 for cm in W.columns:
  expectedtop=set(W[cm].sort_values(ascending=False,kind='stable').head(cfg['top_n']).index);actual=set(top.loc[top.cm.eq(cm),'subCluster']);assert expectedtop==actual
  want_edges={x for x in edges if set(x)<=expectedtop};want_nodes={node for edge in want_edges for node in edge};assert set(each_node.loc[each_node.cm.eq(cm),'subCluster'])==want_nodes;group=each_edge[each_edge.cm.eq(cm)];assert {edgekey(z.subCluster1,z.subCluster2) for z in group.itertuples()}==want_edges
  q=summary[summary.cm.eq(cm)].iloc[0];assert q.n_network_nodes==len(want_nodes) and q.n_network_edges==len(want_edges);assert q.n_parent_types==parent.loc[list(want_nodes)].nunique();assert q.n_cross_parent_edges==sum(parent[a]!=parent[b] for a,b in want_edges)
 result.update(status='PASS_independent_Pearson_author_Holm_then_BH_specificity_and_all_CM_members_edges',n_rank_models=19,n_input_nodes=len(labels),n_samples=X.shape[1],selected_rank=completion['selected_rank'],n_all_pairs=len(pairs),n_global_edges=len(edges),n_multicellular_networks=int(summary.multicellular_network.sum()),max_absolute_Pearson_error=err,all_factor_top15_and_isolated_node_removal_exact=True,all_permutations_verified_by_independent_R=True)
result['aggregate_files']=[{'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in o.iterdir() if p.is_file() and 'PRIVATE' not in p.name];(o/'FINAL_INDEPENDENT_AUDIT.json').open('x').write(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='aggregate_files'}))
