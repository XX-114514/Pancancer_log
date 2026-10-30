"""Audit all state-membership axes/opportunities and sampled fixed projections.

Pure identity, missing assays and weak-vote branches are checked over all
cells. Native-score/classifier replay is explicitly a stratified sample audit.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np,pandas as pd,joblib
from threadpoolctl import threadpool_limits

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
def wait(root,start,seconds):
 p=Path(root)/'provenance/supervisor_terminal_v1.json'
 while not p.exists():
  assert time.monotonic()-start<seconds;time.sleep(5)
 assert load(p)['actual_child_exit_code']==0 and load(p).get('science_terminal_exists',True)

def main(R):
 cfg=load(R/'config/spec_v1.json');start=time.monotonic()
 for root in [cfg['projection'],cfg['reference_auditor']]:wait(root,start,cfg['dependency_wait_seconds'])
 P=Path(cfg['projection']);pc=load(P/'config/spec_v1.json');threadpool_limits(limits=pc['threads']);ap=P/'outputs/accepted_state_membership_collection_v1.json';assert sha(ap)==load(P/'provenance/terminal_v1.json')['acceptance_SHA256'];collection=load(ap)
 A=Path(cfg['reference_auditor']);aa=A/'outputs/accepted_development_state_reference_numeric_audit_v1.json';assert sha(aa)==load(A/'provenance/terminal_v1.json')['acceptance_SHA256'];audit=load(aa);assert audit['producer_acceptance']==collection['frozen_state_model_acceptance'];model_manifest=load(bound(audit['producer_acceptance']))
 app=load(bound(model_manifest['applicability']));scope=load(bound(app['scope_acceptance']));broad=np.load(bound(scope['codes']));context=np.load(bound(app['context_codes']));numeric=load(bound(model_manifest['score_numerical_acceptance']));idx=load(bound(numeric['producer_index']));names=idx['program_ids'];col={v:i for i,v in enumerate(names)}
 sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);common=np.load(Path(pc['preparation'])/'outputs/private_common_assessed_cell_mask_v1.npy');assert common.sum()==idx['n_cells'];roles=load(bound(model_manifest['source_roles']));role={r['source']:r['role'] for r in roles['source_roles']};design=pd.read_csv(Path(pc['design'])/'outputs/private_analysis_unit_context_and_coverage_v1.tsv',sep='\t');C=Path(pc['unit_closure']);units=np.load(bound(load(C/'outputs/accepted_fit_qualification_unit_bridge_v3.json')['unit_codes']))
 memberships={};lookups={};codes={};opportunities={}
 for s in collection['schemes']:
  m=load(bound(s['acceptance']));key=m['scheme_id'].removeprefix('identity_state_');node=np.load(bound(m['node_codes']));parent=np.load(bound(m['parent_codes']));op=np.load(bound(m['node_measurement_mask']));np.testing.assert_array_equal(node>=0,common);np.testing.assert_array_equal(parent[common],broad[common]);mapping=np.array([n['parent_code'] for n in m['nodes']]);np.testing.assert_array_equal(mapping[node[common]],broad[common]);assert op.shape==(len(m['nodes']),len(design))
  lookups[key]={(n['context_code'],n['state_pattern_integer']):n['code'] for n in m['nodes']};memberships[key]=m;codes[key]=node;opportunities[key]=op
 models={m['context']['code']:joblib.load(bound(m['reference_model'])) for m in model_manifest['models']};primary=collection['primary_scheme'].removeprefix('identity_state_');confidence=np.load(bound(collection['neighbour_confidence']));fallback=np.array([lookups['identity_context_only'][(int(c),None)] for c in context[common]])
 np.testing.assert_array_equal(codes['identity_context_only'][common],fallback)
 toconf=np.array([lookups['primary_confident_projection'][(n['context_code'],n['state_pattern_integer'])] for n in memberships[primary]['nodes']]);fallbackconf=np.array([lookups['primary_confident_projection'][(int(c),None)] for c in context[common]])
 expected=np.where(confidence[common]>=pc['minimum_neighbour_agreement'],toconf[codes[primary][common]],fallbackconf);np.testing.assert_array_equal(expected,codes['primary_confident_projection'][common])
 training=np.load(bound(collection['training_cell_mask']));np.testing.assert_array_equal(training,design.source.map(role).eq('development').to_numpy()[units])
 reports=[];nprobes=0;n_batch_replays=0
 for record in idx['sources']:
  source=record['source'];unitcols=design.source.eq(source).to_numpy();measured=np.array([p['native_coverage_pass'] for p in record['coverage_and_interpretation']],bool)
  with np.load(bound(record['cell_axis']),allow_pickle=False) as z:rows=z['fit_row']
  cc=context[rows];y=np.load(bound(record['scores']),mmap_mode='r');nsource=0
  # Opportunity must be metadata-derived even for zero observed contexts.
  for key,m in memberships.items():
   for node in m['nodes']:
    c=node['context_code'];expected_op=True
    if node['state_pattern_integer'] is not None:expected_op=bool(measured[[col[n] for n in models[c]['program_ids']]].all())
    assert (opportunities[key][node['code'],unitcols]==expected_op).all()
  for c in np.unique(cc):
   c=int(c);pos=np.flatnonzero(cc==c);rr=rows[pos];model=models.get(c);possible=model is not None and measured[[col[n] for n in model['program_ids']]].all()
   if not possible:
    for key in memberships:assert (codes[key][rr]==lookups[key][(c,None)]).all()
    continue
   js=[col[n] for n in model['program_ids']];reference={int(r):i for i,r in enumerate(model['reference_fit_rows'])};probes=set([0,len(pos)//2,len(pos)-1]);inref=np.flatnonzero(np.isin(rr,model['reference_fit_rows']))
   if len(inref):probes.add(int(inref[0]))
   weak=np.flatnonzero(confidence[rr]<pc['minimum_neighbour_agreement'])
   if len(weak):probes.add(int(weak[0]))
   for q in sorted(probes):
    row=int(rr[q]);value=np.asarray(y[pos[q]:pos[q]+1][:,js]);assert np.isfinite(value).all();z=model['PCA'].transform(model['scaler'].transform(value));nsource+=1;nprobes+=1
    for key,classifier in model['classifiers'].items():
     predicted=int(classifier.predict(z)[0]);vote=classifier.predict_proba(z)[0];position={int(v):i for i,v in enumerate(classifier.classes_)}
     if row in reference:predicted=int(model['reference_labels'][key][reference[row]])
     expected=lookups[key][(c,predicted)]
     if codes[key][row]!=expected:
      # Matrix multiplication can round differently for a single probe.
      # Resolve any actual discrete mismatch with the exact producer chunk.
      lo=(q//pc['projection_chunk_size'])*pc['projection_chunk_size'];qq=pos[lo:lo+pc['projection_chunk_size']];allv=np.asarray(y[qq][:,js]);zz=model['PCA'].transform(model['scaler'].transform(allv));prediction=classifier.predict(zz);probability=classifier.predict_proba(zz);predicted=int(prediction[q-lo]);vote=probability[q-lo]
      if row in reference:predicted=int(model['reference_labels'][key][reference[row]])
      expected=lookups[key][(c,predicted)];n_batch_replays+=1
     assert codes[key][row]==expected,'Frozen classifier membership replay mismatch'
     if key==primary:
      agreement=float(vote[position[predicted]])
      if not np.isclose(confidence[row],agreement,rtol=1e-4,atol=1e-5):
       lo=(q//pc['projection_chunk_size'])*pc['projection_chunk_size'];qq=pos[lo:lo+pc['projection_chunk_size']];zz=model['PCA'].transform(model['scaler'].transform(np.asarray(y[qq][:,js])));agreement=float(classifier.predict_proba(zz)[q-lo,position[predicted]]);n_batch_replays+=1
      np.testing.assert_allclose(confidence[row],agreement,rtol=1e-4,atol=1e-5)
  reports.append(dict(source=source,n_native_classifier_probe_cells=nsource,all_cell_missing_opportunity_and_fallback_and_parent_axes_checked=True));storage_gate(R,cfg);save(R/'checkpoints'/('source_'+str(len(reports))+'_audited_v1.json'),dict(n_sources=len(reports),n_probes=nprobes))
 op=R/'outputs/accepted_state_projection_numeric_audit_v1.json';save(op,dict(schema='TME_frozen_state_projection_all_axes_sampled_classifier_audit_v1',producer_acceptance=binding(ap),development_reference_audit=binding(aa),n_sources=len(reports),n_common_cells=int(common.sum()),n_schemes=len(memberships),n_native_classifier_probe_cells=nprobes,n_exact_chunk_replays=n_batch_replays,source_reports=reports,all_cell_parent_and_pure_identity_and_weak_vote_and_missing_assay_branches_checked=True,all_node_unit_measurement_opportunities_reconstructed=True,whole_work_unit_development_mask_checked=True,classifier_replay_is_stratified_sample_not_all_cells=True,heldout_only_fixed_projection_not_retraining=True,not_functional_or_clinical_validation=True))
 save(R/'provenance/terminal_v1.json',dict(status='STATE_PROJECTION_ALL_AXES_AND_STRATIFIED_CLASSIFIER_NUMERIC_AUDIT_COMPLETE',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(op)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
