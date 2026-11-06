"""All-factor descriptive diagnostics; preserve original author CM definitions."""
import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.stats import pearsonr

base=PROJECT_ROOT/'runs/20261106_v7_selective_CoVarNet_paper_v2';run=PROJECT_ROOT/'runs/20261107_v4_CM_selection_v3';run.mkdir(exist_ok=False);out=run/'diagnostics';out.mkdir();(run/'scripts').mkdir();(run/'fits').mkdir();(run/'launch').mkdir()
assert json.loads((base/'CM_results/FINAL_INDEPENDENT_AUDIT.json').read_text())['status'].startswith('PASS_')
fq=pd.read_csv(base/'inputs/frequency_PRIVATE.tsv',sep='\t',index_col=0);labels=fq.index;parents=fq.majorCluster;F=fq.iloc[:,1:].astype(float);units=pd.read_csv(base/'inputs/sample_admission_PRIVATE.tsv',sep='\t').set_index('sampleID').loc[F.columns];sources=units.gse_id;D=pd.read_csv(base/'inputs/parent_denominators_PRIVATE.tsv',sep='\t',index_col=0).loc[:,F.columns];counts=pd.read_csv(base/'inputs/counts_all_non_epithelial_PRIVATE.tsv',sep='\t',index_col=0).loc[labels,F.columns];W=pd.read_csv(base/'CM_results/W_PRIVATE.tsv',sep='\t',index_col=0);A=pd.read_csv(base/'CM_results/A_PRIVATE.tsv',sep='\t',index_col=0);H=pd.read_csv(base/'CM_results/H_PRIVATE.tsv',sep='\t',index_col=0);members=pd.read_csv(base/'CM_results/CM_members_after_author_network_filter.tsv',sep='\t');edges=pd.read_csv(base/'CM_results/global_edge.tsv',sep='\t');summary=pd.read_csv(base/'CM_results/CM_network_topology_and_permutation.tsv',sep='\t');X=(F.sub(F.min(axis=1),axis=0)).div(F.max(axis=1)-F.min(axis=1),axis=0)
assert list(A.columns)==['V'+str(k) for k in range(1,15)]
S=pd.read_csv(base/'CM_results/S_PRIVATE.tsv',sep='\t',index_col=0);assert np.allclose(W.to_numpy()@S.to_numpy(),A.to_numpy(),rtol=1e-10,atol=1e-12)
assert list(A.index)==list(W.index) and list(H.index)==list(W.columns)
# S has numeric columns, so R's W %*% S inherits 1..14, in the original factor order.
A.columns=W.columns
def corr(a,b):
 a=np.asarray(a,float);b=np.asarray(b,float)
 return float(np.corrcoef(a,b)[0,1]) if len(a)>=4 and np.std(a)>1e-12 and np.std(b)>1e-12 else np.nan
def center(z,source):return z-pd.Series(z).groupby(np.asarray(source)).transform('mean').to_numpy()
def r2(y,g):
 y=np.asarray(y,float);total=np.sum((y-y.mean())**2)
 return float(1-np.sum(center(y,g)**2)/total) if total>1e-20 else np.nan
edge_rows=[];within_rows=[]
for e in edges.itertuples():
 a,b=e.subCluster1,e.subCluster2;xx=F.loc[a].to_numpy();yy=F.loc[b].to_numpy();keep=(D.loc[parents[a]].to_numpy()>=20)&(D.loc[parents[b]].to_numpy()>=20);both=(D.loc[parents[a]].to_numpy()>0)&(D.loc[parents[b]].to_numpy()>0);loo=[]
 for source in sorted(sources.unique()):
  use=sources.ne(source).to_numpy();loo.append(corr(xx[use],yy[use]));sel=sources.eq(source).to_numpy();within_rows.append({'node1':a,'node2':b,'source':source,'n_observations':int(sel.sum()),'within_source_r':corr(xx[sel],yy[sel]),'both_parent_ge20':int((sel&keep).sum())})
 avail=[x for x in loo if np.isfinite(x)]
 within=[z['within_source_r'] for z in within_rows if z['node1']==a and z['node2']==b and np.isfinite(z['within_source_r'])]
 edge_rows.append({'node1':a,'node2':b,'parent1':parents[a],'parent2':parents[b],'raw_r':corr(xx,yy),'source_centered_r':corr(center(xx,sources),center(yy,sources)),'n_both_parent_observed':int(both.sum()),'r_both_parent_observed':corr(xx[both],yy[both]),'n_both_parent_ge20':int(keep.sum()),'r_both_parent_ge20':corr(xx[keep],yy[keep]),'source_centered_r_parent_ge20':corr(center(xx[keep],sources.to_numpy()[keep]),center(yy[keep],sources.to_numpy()[keep])) if keep.sum()>=4 else np.nan,'LOSO_r_min':min(avail) if avail else np.nan,'LOSO_r_max':max(avail) if avail else np.nan,'LOSO_positive_fraction':float(np.mean(np.array(avail)>0)) if avail else np.nan,'n_estimable_sources':len(within),'n_positive_sources':int(np.sum(np.array(within)>0))})
pd.DataFrame(edge_rows).to_csv(out/'all_global_edge_source_diagnostics.tsv',sep='\t',index=False);pd.DataFrame(within_rows).to_csv(out/'all_global_edge_within_source.tsv',sep='\t',index=False)
# Use sample contributions A-column-L1 times H; invariant to reciprocal factor scaling.
C=H.mul(A.sum(axis=0),axis=0);activity=C.div(C.sum(axis=0),axis=1);activity.to_csv(out/'scale_invariant_activity_PRIVATE.tsv',sep='\t')
node_rows=[]
for node in labels:
 d=D.loc[parents[node]];y=F.loc[node];cc=counts.loc[node];by=cc.groupby(sources).sum();singleton=int((parents==parents[node]).sum())==1
 node_rows.append({'subCluster':node,'majorCluster':parents[node],'n_cells':int(cc.sum()),'n_positive_samples':int(cc.gt(0).sum()),'n_positive_sources':int(by.gt(0).sum()),'largest_source_cell_share':float(by.max()/by.sum()),'largest_source':str(by.idxmax()),'parent_ge20_samples':int(d.ge(20).sum()),'parent_absent_samples':int(d.eq(0).sum()),'source_frequency_R2':r2(y,sources),'single_node_parent':singleton})
pd.DataFrame(node_rows).to_csv(out/'node_coverage_source_diagnostics.tsv',sep='\t',index=False)
module_rows=[];source_rows=[];edge_df=pd.DataFrame(edge_rows)
for q in summary.itertuples():
 cm=q.cm;nodes=members.loc[members.cm.eq(cm),'subCluster'].tolist();ownW=W[cm];myW=ownW.reindex(nodes);myA=A[cm].reindex(nodes);acts=activity.loc[cm];by=acts.groupby(sources).sum();share=by/by.sum();source_mean=acts.groupby(sources).mean();topq=acts>=acts.quantile(.8);tally=sources[topq].value_counts();e=edge_df.loc[edge_df.node1.isin(nodes)&edge_df.node2.isin(nodes)];cross=e[e.parent1!=e.parent2]
 for source in sorted(sources.unique()):
  use=sources.eq(source);source_rows.append({'cm':cm,'source':source,'n_samples':int(use.sum()),'activity_mean':float(acts[use].mean()),'activity_median':float(acts[use].median()),'activity_mass_share':float(share[source]),'n_top_quintile':int((topq&use).sum()),'network_coobserved_samples':int(counts.loc[nodes,use].gt(0).all(axis=0).sum()) if nodes else 0})
 module_rows.append({'cm':cm,'network_W_mass_fraction':float(myW.sum()/ownW.sum()),'network_effective_A_mass_fraction':float(myA.sum()/A[cm].sum()),'minimum_member_relative_W':float(myW.min()/ownW.max()) if nodes else np.nan,'network_nodes':len(nodes),'network_edges':q.n_network_edges,'network_components':q.n_components,'cross_parent_edges':q.n_cross_parent_edges,'network_all_nodes_observed_samples':int(counts.loc[nodes].gt(0).all(axis=0).sum()) if nodes else 0,'network_all_parents_ge20_samples':int(D.loc[parents.loc[nodes].unique()].ge(20).all(axis=0).sum()) if nodes else 0,'source_activity_R2':r2(acts,sources),'max_source_activity_mass_share':float(share.max()),'max_source_activity_mass_source':str(share.idxmax()),'activity_effective_sources':float(1/(share**2).sum()),'top_quintile_largest_source_fraction':float(tally.max()/topq.sum()),'source_centered_edge_r_min':float(e.source_centered_r.min()) if len(e) else np.nan,'source_centered_cross_edge_r_min':float(cross.source_centered_r.min()) if len(cross) else np.nan,'LOSO_cross_edge_r_min':float(cross.LOSO_r_min.min()) if len(cross) else np.nan,'n_cross_edges_positive_source_centered':int(cross.source_centered_r.gt(0).sum()),'n_cross_edges_positive_parent_ge20':int(cross.source_centered_r_parent_ge20.gt(0).sum())})
pd.DataFrame(module_rows).to_csv(out/'all_CM_quality_diagnostics.tsv',sep='\t',index=False);pd.DataFrame(source_rows).to_csv(out/'all_CM_source_activity.tsv',sep='\t',index=False)
# Job plan fixed before inspecting diagnostic effects. Fixed-rank tests, not new optimal rank claims.
jobs=[{'job':'seed_%d'%j,'kind':'seed','held_source':'','seed':2026110700+j} for j in range(1,6)]+[{'job':'LOSO_%02d'%j,'kind':'LOSO','held_source':s,'seed':20261106+14} for j,s in enumerate(sorted(sources.unique()),1)]
(run/'jobs.json').write_text(json.dumps(jobs,indent=2));units[['gse_id']].to_csv(run/'sample_source_PRIVATE.tsv',sep='\t');(out/'diagnostic_complete.json').write_text(json.dumps({'status':'ALL_14_CM_DESCRIPTIVE_DIAGNOSTICS_COMPLETE','UTC':datetime.now(timezone.utc).isoformat(),'n_samples':len(F.columns),'n_sources':sources.nunique(),'n_nodes':len(F),'n_global_edges':len(edges),'n_factors':len(W.columns),'n_followup_jobs':len(jobs),'inference':'Descriptive selection diagnostics, no postselection biological validation or independent-patient inference','input_files':{str(p.relative_to(base)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [base/'inputs/frequency_PRIVATE.tsv',base/'CM_results/W_PRIVATE.tsv',base/'CM_results/H_PRIVATE.tsv',base/'CM_results/A_PRIVATE.tsv']}},indent=2));print(json.dumps({'status':'ALL_CM_DIAGNOSTICS_COMPLETE','jobs':len(jobs),'sources':sources.nunique()}))
