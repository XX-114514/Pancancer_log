"""Counterexamples for compartment exclusion versus complete parent denominators."""
import importlib.util
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts/prepare_TME_immune_only_design_v1.py'
spec=importlib.util.spec_from_file_location('immune_scope',p);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
parents=[{'code':i,'name':n} for i,n in enumerate(['T_NK','B_Plasma','Myeloid','Epithelial','Mesenchymal_stromal','Endothelial','Mast'])]
allowed=['T_NK','B_Plasma','Myeloid','Mast']
# A highly lymphocyte-rich material is not assigned a sorting protocol.
common=np.ones(12,bool);axis=np.array([0,0,0,1,2,2,3,4,5,6,-1,0]);common[11]=False
m=mod.immune_scope_mask(common,axis,parents,allowed)
np.testing.assert_array_equal(np.flatnonzero(m),[0,1,2,3,4,5,9])
# A recorded malignant immune cell is excluded; an immune label alone is not its proof.
known=np.zeros(12,bool);known[0]=True
mm=mod.immune_scope_mask(common,axis,parents,allowed,known)
assert not mm[0] and mm[1] and not mm[11]
# Excluding an unselected discovery child would corrupt the denominator: all T children remain.
node=np.array([0,0,1,2,3,3,4,5,6,7,8,9]);t=(axis==0)&m
assert np.count_nonzero(t)==3 and np.count_nonzero(t&(node==0))/np.count_nonzero(t)==2/3
try: mod.immune_scope_mask(common,axis,parents,['Undefined'])
except AssertionError: pass
else: raise AssertionError('Unknown compartment must not silently enter')
print('PASS: immune domain, malignancy flag, full denominator, unknown parent counterexamples')
