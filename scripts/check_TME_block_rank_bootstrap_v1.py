"""Independent expanded-row SciPy checks for unequal block sizes and ties."""
import argparse,datetime,hashlib,json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from TME_block_rank_bootstrap_v1 import weighted_spearman,block_multiplicities

def main(R):
 x=np.array([0,0,.2,.2,.7,1,1.]);y=np.array([.9,.9,.8,.3,.3,0,0]);blocks=np.array([0,1,1,2,2,2,3]);w=block_multiplicities(blocks,64,42);r=weighted_spearman(x,y,w)
 for i in range(64):
  ids=np.repeat(np.arange(len(x)),w[i]);a=spearmanr(x[ids],y[ids]).statistic
  np.testing.assert_allclose(r[i],a,atol=1e-12,equal_nan=True)
 same=weighted_spearman(x,x,w);assert np.isfinite(same).any();np.testing.assert_allclose(same[np.isfinite(same)],np.ones(np.isfinite(same).sum()),atol=1e-12)
 assert np.isnan(weighted_spearman(np.zeros(7),y,w)).all();np.testing.assert_array_equal(w,block_multiplicities(blocks,64,42))
 op=R/'outputs/accepted_block_rank_bootstrap_math_v1.json';op.write_text(json.dumps(dict(schema='TME_block_rank_bootstrap_synthetic_math_v1',all_tests_passed=True,n_synthetic_bootstraps_replayed=64,unequal_block_sizes_and_ties_checked=True,positive_negative_and_constant_checked=True,expanded_rows_match_Scipy=True),indent=2));(R/'provenance/terminal_v1.json').write_text(json.dumps(dict(status='EXACT_INTEGER_BLOCK_BOOTSTRAP_SPEARMAN_MATH_PASS',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=hashlib.sha256(op.read_bytes()).hexdigest()),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  (R/'provenance/failure_v1.json').write_text(json.dumps(dict(traceback=traceback.format_exc()),indent=2));raise
