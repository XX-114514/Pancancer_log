"""Counterexamples for collinearity, missingness and the real official MOFA API."""
import argparse,datetime,hashlib,json,sys
from pathlib import Path
import numpy as np
from prepare_TME_state_multiview_inputs_v1 import residual,transform
from run_TME_state_MOFA_v1 import fit_mofa,wait

def main(R):
 cfg=json.loads((R/'config/spec_v1.json').read_text());import time;wait(cfg['vendor_stage'],3600,time.monotonic());V=Path(cfg['vendor_stage']);p=V/'outputs/accepted_official_MOFA_vendor_v1.json'
 assert json.loads((V/'provenance/supervisor_terminal_v1.json').read_text())['actual_child_exit_code']==0
 assert hashlib.sha256(p.read_bytes()).hexdigest()==json.loads((V/'provenance/terminal_v1.json').read_text())['acceptance_SHA256'];v=json.loads(p.read_text());sys.path.insert(0,v['vendor_path'])
 rng=np.random.default_rng(42);n=120;cov=rng.normal(size=(n,3));y=3*cov[:,0]-2*cov[:,1]+rng.normal(size=n)
 r,status=residual(y,cov,10);assert status=='estimable' and np.ptp(r)>0 and r.min()<0<r.max();assert max(abs(cov.T@r))<1e-10
 rr,status=residual(y,np.column_stack([cov[:,0],cov[:,0],cov[:,1]]),10);assert rr is None and status=='collinear_covariates'
 rr,status=residual(y[:4],cov[:4],10);assert rr is None and status=='insufficient_units'
 x=np.array([[0.,np.nan,2.],[1.,np.nan,2.],[2.,0.,2.]])
 zz,pars=transform(x);assert np.isnan(zz[:,1]).all();assert np.isnan(zz[:,2]).all();assert zz[0,0]<0<zz[2,0]
 latent=rng.normal(size=(n,2));values=[latent@rng.normal(size=(2,d))+rng.normal(0,.15,(n,d)) for d in [20,25,18]]
 values[0][::5,:]=np.nan;values[1][1::7,::3]=np.nan;values[2][2::9,:]=np.nan;before=[np.isnan(a).copy() for a in values]
 e,Z,W=fit_mofa(values,['synthetic_A','synthetic_B','synthetic_C'],[['f'+str(i) for i in range(a.shape[1])] for a in values],['ARTIFICIAL_'+str(i) for i in range(n)],dict(factors=3,iterations=300,seed=42),R/'outputs/private_synthetic_MOFA.hdf5')
 assert all(np.array_equal(np.isnan(a),b) for a,b in zip(values,before))
 r2=np.asarray(e.model.calculate_variance_explained(total=True)[0]);assert r2.min()>.65
 op=R/'outputs/accepted_MOFA_math_and_API_synthetic_v1.json';op.write_text(json.dumps(dict(schema='TME_official_MOFA_math_API_synthetic_acceptance_v1',all_tests_passed=True,official_version=v['version'],n_tests=6,source_depth_covariates_not_automatically_dropped=True,missing_mask_preserved=True,negative_residuals_preserved=True,synthetic_three_view_total_R2=r2.tolist(),not_human_science_result=True),indent=2))
 (R/'provenance/terminal_v1.json').write_text(json.dumps(dict(status='OFFICIAL_MOFA_MISSINGNESS_AND_MATH_SYNTHETIC_PASS',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=hashlib.sha256(op.read_bytes()).hexdigest()),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  (R/'provenance/failure_v1.json').write_text(json.dumps(dict(traceback=traceback.format_exc()),indent=2));raise
