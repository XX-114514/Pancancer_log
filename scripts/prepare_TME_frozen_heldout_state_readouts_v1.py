"""Fixed programme readouts in source lockboxes after development freeze.

No new features, parent scope or thresholds are selected in the held-out data.
Only released groups are measured; low coverage and absent groups remain NaN.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def bound(b):
 p=Path(b['path']);assert sha(p)==b['SHA256'];return p
def binding(p):return dict(path=str(p),SHA256=sha(p))
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def main(R):
 cfg=load(R/'config/spec_v1.json');A=Path(cfg['reference_auditor']);start=time.monotonic()
 while not (A/'provenance/supervisor_terminal_v1.json').exists():
  assert time.monotonic()-start<cfg['dependency_wait_seconds'];time.sleep(5)
 assert load(A/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 ap=A/'outputs/accepted_development_state_reference_numeric_audit_v1.json';assert sha(ap)==load(A/'provenance/terminal_v1.json')['acceptance_SHA256'];aud=load(ap);model=load(bound(aud['producer_acceptance']));assert model['n_lockbox_score_arrays_opened']==0
 inp=load(bound(model['input_acceptance']));roles=load(bound(model['source_roles']));app=load(bound(model['applicability']));ra=load(bound(model['score_numerical_acceptance']));idx=load(bound(ra['producer_index']));context=np.load(bound(app['context_codes']));applicable=np.load(bound(app['context_program_applicability']));lockbox=np.load(bound(roles['lockbox_mask']));role={r['source']:r['role'] for r in roles['source_roles']};n_programmes=len(idx['program_ids'])
 sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);pc=load(Path(cfg['score_producer'])/'config/spec_v1.json');H=Path(pc['model']);totals=np.load(H/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r');design=Path(cfg['design']);allowed=np.load(design/'outputs/private_all_qualified_tumour_unit_mask_v1.npy').copy()
 closure=Path(cfg['unit_closure']);unitaxis=np.load(bound(load(closure/'outputs/accepted_fit_qualification_unit_bridge_v3.json')['unit_codes']));common=np.load(Path(pc['preparation'])/'outputs/private_common_assessed_cell_mask_v1.npy');commoncounts=np.bincount(unitaxis[common],minlength=len(allowed));allowed&=commoncounts>=cfg['minimum_total_cells']
 fixed={v['context']['code']:{idx['program_ids'].index(n) for n in v['program_ids']} for v in inp['views']};sources=[];total=0;released=0
 for record in idx['sources']:
  if role[record['source']]!='lockbox':continue
  with np.load(bound(record['cell_axis']),allow_pickle=False) as z:rows=z['fit_row'];units=z['qualification_unit_code'];detected=z['detected_native_genes']
  assert lockbox[rows].all();cc=context[rows];y=np.load(bound(record['scores']),mmap_mode='r');measured=np.array([p['native_coverage_pass'] for p in record['coverage_and_interpretation']],bool)
  keys,inv=np.unique(units.astype(np.int64)*app['n_contexts']+cc,return_inverse=True);counts=np.bincount(inv);order=np.argsort(inv,kind='stable');starts=np.r_[0,np.cumsum(counts)];mean=np.full((len(keys),n_programmes),np.nan,np.float32);halfA=np.full_like(mean,np.nan);halfB=np.full_like(mean,np.nan);rna=np.empty(len(keys));ngene=np.empty(len(keys));main=np.zeros(len(keys),bool)
  for g,key in enumerate(keys):
   pos=order[starts[g]:starts[g+1]];c=int(key%app['n_contexts']);u=int(key//app['n_contexts']);rna[g]=np.log1p(totals[rows[pos]]).mean();ngene[g]=detected[pos].mean();main[g]=len(pos)>=cfg['minimum_primary_cells'] and allowed[u]
   if c not in fixed or len(pos)<cfg['minimum_release_cells']:continue
   js=np.array(sorted(fixed[c]),int);js=js[applicable[c,js]&measured[js]]
   if not len(js):continue
   values=np.asarray(y[pos][:,js]);assert np.isfinite(values).all();mean[g,js]=values.mean(0,dtype=np.float64)
   if len(pos)>=cfg['minimum_primary_cells']:
    sample=np.random.default_rng(np.random.SeedSequence([cfg['seed'],u,c])).permutation(len(pos));n=len(sample)//2;halfA[g,js]=values[sample[:n]].mean(0,dtype=np.float64);halfB[g,js]=values[sample[n:]].mean(0,dtype=np.float64);np.testing.assert_allclose((n*halfA[g,js]+(len(pos)-n)*halfB[g,js])/len(pos),mean[g,js],atol=2e-7,rtol=1e-5)
  p=R/'outputs'/('private_'+record['source']+'_frozen_state_lockbox_readouts_v1.npz');np.savez_compressed(p,unit_code=keys//app['n_contexts'],context_code=keys%app['n_contexts'],n_cells=counts,means=mean,split_mean_A=halfA,split_mean_B=halfB,mean_log1p_RNA=rna,mean_detected_native_genes=ngene,primary_material_context_eligible=main)
  sources.append(dict(source=record['source'],n_cells=len(rows),n_context_unit_groups=len(keys),n_primary_context_units=int(main.sum()),readouts=binding(p)));total+=len(rows);released+=int(np.isfinite(mean).any(1).sum());storage_gate(R,cfg)
 assert total==roles['n_lockbox_cells'] and len(sources)==roles['n_lockbox_sources'];op=R/'outputs/accepted_frozen_lockbox_state_readouts_v1.json';save(op,dict(schema='TME_frozen_development_programmes_lockbox_readouts_v1',development_reference_numeric_acceptance=binding(ap),input_acceptance=model['input_acceptance'],source_roles=model['source_roles'],applicability=model['applicability'],score_numerical_acceptance=model['score_numerical_acceptance'],program_ids=idx['program_ids'],n_sources=len(sources),n_cells=total,n_released_context_unit_readouts=released,sources=sources,minimum_release_cells=cfg['minimum_release_cells'],minimum_primary_cells=cfg['minimum_primary_cells'],no_feature_parent_threshold_or_state_model_refit=True,all_readouts_require_native_coverage_and_context_scope=True,clinical_outcomes_used=False,heldout_is_conditional_on_frozen_all_source_atlas=True))
 save(R/'provenance/terminal_v1.json',dict(status='FROZEN_STATE_PROGRAMME_LOCKBOX_READOUTS_PREPARED_INCREMENT_TEST_PENDING',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(op)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
