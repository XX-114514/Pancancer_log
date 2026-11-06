"""Identify the omitted source behind every global-edge sensitivity value."""
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261107_v4_CM_selection_v3';b=PROJECT_ROOT/'runs/20261106_v7_selective_CoVarNet_paper_v2';f=pd.read_csv(b/'inputs/frequency_PRIVATE.tsv',sep='\t',index_col=0).iloc[:,1:];u=pd.read_csv(b/'inputs/sample_admission_PRIVATE.tsv',sep='\t').set_index('sampleID');src=u.loc[f.columns,'gse_id'];e=pd.read_csv(b/'CM_results/global_edge.tsv',sep='\t');rows=[]
for z in e.itertuples():
 for s in sorted(src.unique()):
  use=src.ne(s);a=f.loc[z.subCluster1,use];c=f.loc[z.subCluster2,use];v=float(np.corrcoef(a,c)[0,1]) if np.std(a)>1e-12 and np.std(c)>1e-12 else np.nan;rows.append({'node1':z.subCluster1,'node2':z.subCluster2,'held_source':s,'n_retained_samples':int(use.sum()),'r':v})
pd.DataFrame(rows).to_csv(r/'diagnostics/all_global_edge_LOSO_detail.tsv',sep='\t',index=False);print({'all_global_edge_source_leaveouts':len(rows)})
