"""Counterexamples for real statistical failure modes, not implementation mirrors."""
import numpy as np,pandas as pd
from vascular_spatial_controls_v2 import paired_p_greater,fdr,regional_ES,edge_expression_null
assert paired_p_greater([1,2,3,4,5,6])==1/64
assert paired_p_greater([-1,-2,-3,-4,-5,-6])==1
assert paired_p_greater([0,0,0,0])==1
assert np.allclose(fdr([.01,.04,1]),[.03,.06,1])
# A missing cell class must not force its receptor to zero in regional total RNA.
cells=[];counts=[];genes=['CCL19','CCR7','CCL21','CXCL12','CXCR4','ACKR1']
for group in ['AB','A_only','B_only']:
 for j in range(5):
  b=len(cells)//22;types=['Blood_EC']*2+(['Mural']*5 if group!='B_only' else ['Other']*5)+(['B']*5 if group!='A_only' else ['Other']*5)+['Other']*10
  for k,typ in enumerate(types):
   cells.append({'fov':'F','x_um':b*110+20+(k%5),'y_um':20+k//5,'native_total':50,'cell_type':typ});v=np.zeros(6);v[0]=1 if k==2 else 0;v[1]=1 if k==7 else 0;counts.append(v)
r=regional_ES(pd.DataFrame(cells),genes,np.array(counts));z=r[(r.sender=='Mural')&(r.receiver=='B')&(r.ligand=='CCL19')].iloc[0]
assert z.status=='eligible_descriptive_effect' and z.A_only_ES_raw==1 and z.B_only_ES_raw==1 and z.AB_minus_A==0 and z.AB_minus_B==0
# Observed zero does not imply a zero null expectation; depletion must be retained.
c=pd.DataFrame({'fov':['F']*20,'x_um':[10]*5+[110]*5+[10]*5+[110]*5,'y_um':[0]*20,'native_total':[50]*20,'cell_type':['Mural']*10+['B']*10});x=np.zeros((20,6));x[:5,0]=1;x[15:,1]=1
r=edge_expression_null(c,genes,x,nboot=39,seed=7);z=r[(r.sender=='Mural')&(r.receiver=='B')&(r.ligand=='CCL19')].iloc[0];assert z.observed_LR_positive_edges==0 and z.excess_edge_fraction<0 and z.P_MC_add1==1
# Identical coordinates in disconnected FOVs must never create cross-field contacts.
c['fov']=['F1']*10+['F2']*10;r=edge_expression_null(c,genes,x,nboot=9);assert r.n_edges_50um.sum()==0 and r.P_MC_add1.isna().all()
print('PASS: exact patient P, full FDR, non-tautological regional expression, zero-observed negative excess, no crossFOV edges')
