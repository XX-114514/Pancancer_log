"""Manual rank and full permutation/bootstrap replay, no model/genes refit.
An independent numerical audit of the fixed heldout test, not new data.
"""
from pathlib import Path
import os,json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
def ranks(a):
 a=np.asarray(a);return np.array([np.sum(a<v)+(np.sum(a==v)+1)/2 for v in a],float)
def rho(x,y):
 a=ranks(x);b=ranks(y);a-=a.mean();b-=b.mean();d=np.sqrt(sum(a*a)*sum(b*b));return sum(a*b)/d if d>0 else np.nan
def main():
 O=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])/'CRC_fine_pair_DIALOGUE_20261014_v1';s=pd.read_csv(O/'private_frozen_projection_v1/all_patient_mean_scores_private.tsv',sep='\t');out=pd.read_csv(O/'heldout_patient_results_v1.tsv',sep='\t');loo=pd.read_csv(O/'heldout_leave_one_patient_v1.tsv',sep='\t');rng=np.random.default_rng(20261014);checks=0
 for row in out.itertuples():
  if row.status!='frozen_holdout_estimated':assert pd.isna(row.permutation_p) and pd.isna(row.BH_two_slot_family);continue
  w=s[(s.partition=='holdout')&(s.program==row.program)].pivot(index='sample',columns='cell_type',values='score');x=w.Tnode9.to_numpy();y=w.Macnode18.to_numpy();assert len(x)==12;obs=rho(x,y);assert abs(obs-row.rho)<1e-12;extreme=0
  for _ in range(20000):extreme+=abs(rho(x,y[rng.permutation(12)]))>=abs(obs)-1e-12;checks+=1
  p=(extreme+1)/20001;assert abs(p-row.permutation_p)<1e-12;boot=[];undefined=0
  for _ in range(2000):
   ix=rng.integers(0,12,12);r=rho(x[ix],y[ix]);checks+=1
   if np.isfinite(r):boot.append(r)
   else:undefined+=1
  np.testing.assert_allclose(np.quantile(boot,[.025,.975]),[row.lower,row.upper],atol=1e-12);assert undefined==row.undefined_bootstrap_draws
  for i in range(12):m=np.arange(12)!=i;actual=loo[(loo.program==row.program)&(loo.left_out_ordinal==i+1)].rho.iloc[0];assert abs(rho(x[m],y[m])-actual)<1e-12;checks+=1
  assert abs(min(1,2*p)-row.BH_two_slot_family)<1e-12
 parity=pd.read_csv(O/'frozen_projection_numeric_audit.tsv',sep='\t');assert parity.max_saved_score_error.max()<1e-6;assert parity.n_discovery_cells.sum()==2225 and parity.n_holdout_cells.sum()==1170
 j=dict(status='PASS_FIXED_HOLDOUT_STATS_AND_DISCOVERY_SCORE_PARITY',observed_UTC=datetime.now(timezone.utc).isoformat(),manual_rank_draws_and_LOO_checked=checks,permutations_replayed=20000,bootstraps_replayed=2000,all_supported_discovery_scores_reconstructed=2225,holdout_scores_frozen_projection=1170,parity_max_error=float(parity.max_saved_score_error.max()),family_slots=2,new_model_or_genes_fit_on_holdout=False,not_independent_project_cohort_validation=True,auditor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());(O/'independent_holdout_audit_v1.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j))
if __name__=='__main__':main()
