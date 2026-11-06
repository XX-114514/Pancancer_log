"""Trace all 26 original-CM03 edge opportunities through each LOSO filter."""
import json,hashlib
import numpy as np,pandas as pd
from scipy.stats import t
r=PROJECT_ROOT/'runs/20261107_v4_CM_selection_v3';b=PROJECT_ROOT/'runs/20261106_v7_selective_CoVarNet_paper_v2';F=pd.read_csv(b/'inputs/frequency_PRIVATE.tsv',sep='\t',index_col=0).iloc[:,1:];src=pd.read_csv(b/'inputs/sample_admission_PRIVATE.tsv',sep='\t').set_index('sampleID').loc[F.columns,'gse_id'];match=pd.read_csv(r/'diagnostics/all_CM_refit_matches.tsv',sep='\t');g=match[match.cm.eq('CM03')&match.kind.eq('LOSO')];original=pd.read_csv(b/'CM_results/each_edge.tsv',sep='\t');original=original[original.cm.eq('CM03')];rows=[]
def adj(p,kind):
 n=len(p);ix=np.argsort(p);v=p[ix];v=np.minimum(1,np.maximum.accumulate(v*np.arange(n,0,-1))) if kind=='holm' else np.minimum(1,np.minimum.accumulate((v*n/np.arange(1,n+1))[::-1])[::-1]);out=np.empty(n);out[ix]=v;return out
for z in g.itertuples():
 folder=r/'fits'/z.job;w=pd.read_csv(folder/'W_PRIVATE.tsv',sep='\t',index_col=0);fq=F.loc[w.index,src.ne(z.held_source)];R=np.corrcoef(fq.to_numpy());n=len(w);ii,jj=np.triu_indices(n,1);rho=R[ii,jj];p=2*t.sf(np.abs(rho)*np.sqrt((fq.shape[1]-2)/np.maximum(1-rho*rho,1e-300)),df=fq.shape[1]-2);fdr=adj(adj(p,'holm'),'BH');threshold=1-27/(2*n-3);top=set(w[z.matched_factor_W].sort_values(ascending=False,kind='stable').head(15).index);actual=pd.read_csv(folder/'each_edge.tsv',sep='\t');actual=actual[actual.cm.eq(z.matched_factor_W)];actualkeys={tuple(sorted((a.subCluster1,a.subCluster2))) for a in actual.itertuples()};idx={label:i for i,label in enumerate(w.index)}
 for e in original.itertuples():
  i,j=sorted((idx[e.subCluster1],idx[e.subCluster2]));v=R[i,j];k=np.flatnonzero((ii==i)&(jj==j))[0];bg=np.r_[R[i,np.arange(n)!=i],R[(np.arange(n)!=i)&(np.arange(n)!=j),j]];spe=np.mean(bg<=v);endpoints=e.subCluster1 in top and e.subCluster2 in top;globalpass=v>.2 and fdr[k]<=.05 and spe>=threshold;expected=endpoints and globalpass;observed=tuple(sorted((e.subCluster1,e.subCluster2))) in actualkeys;assert expected==observed
  reasons=[]
  if not endpoints:reasons.append('endpoint_outside_matched_factor_top15')
  if v<=.2:reasons.append('r_le_0.2')
  if fdr[k]>.05:reasons.append('author_FDR_gt_0.05')
  if spe<threshold:reasons.append('specificity_below_threshold')
  rows.append({'held_source':z.held_source,'matched_factor':z.matched_factor_W,'node1':e.subCluster1,'node2':e.subCluster2,'r':v,'author_FDR':fdr[k],'specificity':spe,'specificity_cutoff':threshold,'both_endpoints_in_top15':endpoints,'global_edge_passes':bool(globalpass),'edge_in_matched_network':observed,'failure_reason':';'.join(reasons) or 'retained'})
p=r/'diagnostics/CM03_LOSO_edge_filter_trace.tsv';pd.DataFrame(rows).to_csv(p,sep='\t',index=False);print(json.dumps({'status':'PASS_all26_CM03_edge_filter_stages_independently_reconstructed','files':{'CM03_LOSO_edge_filter_trace.tsv':{'text':p.read_text(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}}}))
