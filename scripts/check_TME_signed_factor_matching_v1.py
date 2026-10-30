"""Sign/permutation/scale and unmatched-zero counterexamples for factor matching."""
import argparse,datetime,hashlib,json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from audit_TME_state_MOFA_fits_v1 import signed_match

def main(R):
 q=np.linalg.qr(np.random.default_rng(42).normal(size=(12,8)))[0];A=q[:,:3];per=[2,0,1];sign=np.array([-1,1,-1]);B=A[:,per]*sign*3
 matches,ua,ub=signed_match(A,B,.6);assert len(matches)==3 and not ua and not ub
 z=np.random.default_rng(43).normal(size=(50,3));zz=z[:,per]*sign/3
 for i,j,c in matches:
  assert per[j]==i and abs(abs(c)-1)<1e-12;assert spearmanr(z[:,i],zz[:,j]*np.sign(c)).statistic>.999999
 matches,ua,ub=signed_match(np.column_stack([A,np.zeros(12)]),np.column_stack([B,q[:,4]]),.6);assert len(matches)==3 and ua==[3] and ub==[3]
 op=R/'outputs/accepted_signed_factor_matching_math_v1.json';op.write_text(json.dumps(dict(schema='TME_signed_factor_matching_math_v1',all_tests_passed=True,n_tests=4,sign_scale_permutation_and_one_to_one_unmatched_zero_checked=True,not_human_state_validation=True),indent=2));(R/'provenance/terminal_v1.json').write_text(json.dumps(dict(status='SIGNED_FACTOR_MATCHING_SYNTHETIC_PASS',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=hashlib.sha256(op.read_bytes()).hexdigest()),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  (R/'provenance/failure_v1.json').write_text(json.dumps(dict(traceback=traceback.format_exc()),indent=2));raise
