"""Regression: adopted ambient deltas may be fractional; raw counts may not."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_TME_standard_Harmony_v2 import gather_nonnegative_float_delta
from train_TME_scVI_GPU_v1 import gather_counts
class CSR:
 def __init__(self,values):
  self.attrs={'shape':[2,3]};self.fields={'indptr':np.array([0,2,3]),'indices':np.array([0,2,1]),'data':np.asarray(values,float)}
 def __getitem__(self,key):return self.fields[key]
def test_fractional_delta_and_row_order():
 x=gather_nonnegative_float_delta(CSR([.25,1.5,.75]),np.array([1,0,1]))
 np.testing.assert_array_equal(x.toarray(),[[0,.75,0],[.25,0,1.5],[0,.75,0]])
def test_negative_delta_rejected():
 try:gather_nonnegative_float_delta(CSR([.25,-.5,.75]),np.array([0,1]))
 except ValueError:return
 raise AssertionError('Negative delta accepted')
def test_nonfinite_delta_rejected():
 try:gather_nonnegative_float_delta(CSR([.25,float('nan'),.75]),np.array([0,1]))
 except ValueError:return
 raise AssertionError('Nonfinite delta accepted')
def test_fractional_raw_counts_rejected():
 try:gather_counts(CSR([.25,1.5,.75]),np.array([0,1]))
 except ValueError:return
 raise AssertionError('Fractional raw counts accepted')
if __name__=='__main__':
 for test in [test_fractional_delta_and_row_order,test_negative_delta_rejected,test_nonfinite_delta_rejected,test_fractional_raw_counts_rejected]:test()
 print('4 count-domain and fractional-delta regression cases passed')
