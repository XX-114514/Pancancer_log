"""Independent readback of readout sidecars, original means and quantiles."""
from pathlib import Path
import json,sys,time,datetime,hashlib,traceback
import numpy as np
import pandas as pd

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def main(A,R,P):
 deadline=time.monotonic()+3600
 while not (R/'provenance/terminal.json').exists():
  if (R/'provenance/failure.json').exists():raise RuntimeError('Producer failed')
  if time.monotonic()>deadline:raise RuntimeError('Producer timed out')
  time.sleep(15)
 term=load(R/'provenance/terminal.json');assert term['actual_exit_code']==0
 names=load(P/'outputs/score_program_axis_v1.json')['program_ids'];means=[];errors=[];n=0;nv=0
 for rec in load(P/'outputs/source_score_index_v1.json'):
  s=rec['source'];u=pd.read_csv(R/'outputs'/('private_'+s+'_units.tsv.gz'),sep='\t',dtype={'sample':str,'subtype':str});d=np.load(R/'outputs'/('private_'+s+'_readout_diagnostics.npz'))
  original=pd.read_csv(P/'outputs'/('private_'+s+'_unit_subtype_means_v1.tsv.gz'),sep='\t',dtype={'sample':str,'subtype':str});oq=pd.read_csv(P/'outputs'/('private_'+s+'_unit_subtype_q095_v1.tsv.gz'),sep='\t',dtype={'sample':str,'subtype':str})
  assert list(zip(u['sample'],u.subtype))==list(zip(original['sample'],original.subtype))==list(zip(oq['sample'],oq.subtype))
  np.testing.assert_array_equal(u.n_cells,original.n_cells);np.testing.assert_allclose(d['mean'],original[names],atol=1e-7,rtol=1e-6);np.testing.assert_allclose(d['q95'],oq[names],atol=1e-7,rtol=1e-6)
  expected=np.sqrt(np.mean((d['half_mean_A']-d['half_mean_B'])**2,axis=1))/2;np.testing.assert_allclose(d['sampling_se_diagnostic'],expected,equal_nan=True,atol=1e-12)
  for size in [20,40,80]:
   arr=d['downsample_'+str(size)];mask=u.n_cells.to_numpy()>=size;assert np.isfinite(arr[mask]).all() and np.isnan(arr[~mask]).all();assert arr[mask].min(initial=0)>=0 and arr[mask].max(initial=0)<=1.000001
  means.append(d['mean']);errors.append(expected);n+=int(u.n_cells.sum());nv+=d['mean'].size
 assert n==776270
 np.testing.assert_allclose(np.load(R/'outputs/private_S_candidate_means.npy'),np.concatenate(means));np.testing.assert_allclose(np.load(R/'outputs/private_Q_sampling_SE_diagnostic.npy'),np.concatenate(errors),equal_nan=True)
 c=pd.read_csv(R/'outputs/private_C_TNK_scope_counts.tsv.gz',sep='\t',index_col=[0,1]);assert int(c.to_numpy().sum())==n
 pr=pd.read_csv(R/'outputs/private_C_TNK_scope_proportions.tsv.gz',sep='\t',index_col=[0,1]);np.testing.assert_allclose(pr.to_numpy().sum(1),1)
 assert all(not r['module_admitted'] for r in load(R/'outputs/program_applicability_and_admission.json'))
 summary=dict(completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_cells=n,n_mean_values_checked=nv,n_q95_values_crosschecked=nv,producer_terminal_SHA256=sha(R/'provenance/terminal.json'),no_new_module_admission=True,scope='independent_mean_q95_readback_and_sampling_error_mask_composition_checks; statistical_correlations_not_independently_recomputed',actual_exit_code=0)
 save(A/'outputs/acceptance.json',summary);save(A/'provenance/terminal.json',summary)
if __name__=='__main__':
 A,R,P=map(Path,sys.argv[1:])
 try:main(A,R,P)
 except BaseException:save(A/'provenance/failure.json',dict(traceback=traceback.format_exc()));raise
