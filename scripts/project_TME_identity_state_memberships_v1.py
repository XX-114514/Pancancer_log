"""Project frozen development state patterns without refitting held-out sources.

Every common-domain cell keeps one broad parent and one partition membership.
Unmodelled/missing opportunities use explicit identity-retained fallbacks.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np,pandas as pd,joblib

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
 cfg=load(R/'config/spec_v1.json');start=time.monotonic();T=Path(cfg['training']);tp=T/'provenance/supervisor_terminal_v1.json'
 while not tp.exists():
  assert time.monotonic()-start<cfg['dependency_wait_seconds'];time.sleep(5)
 assert load(tp)['actual_child_exit_code']==0
 sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);ap=T/'outputs/accepted_development_identity_state_partitions_v1.json';assert sha(ap)==load(T/'provenance/terminal_v1.json')['acceptance_SHA256'];trained=load(ap);assert trained['n_lockbox_score_arrays_opened']==0
 app=load(bound(trained['applicability']));scope=load(bound(app['scope_acceptance']));parents=np.load(bound(scope['codes']));ctx=np.load(bound(app['context_codes']));N=len(parents)
 audit=load(bound(trained['score_numerical_acceptance']));idx=load(bound(audit['producer_index']));names=idx['program_ids'];col={n:j for j,n in enumerate(names)}
 P=Path(cfg['preparation']);common=np.load(P/'outputs/private_common_assessed_cell_mask_v1.npy');assert sha(P/'outputs/private_common_assessed_cell_mask_v1.npy')==idx['common_cell_mask_SHA256'];assert common.sum()==idx['n_cells']
 C=Path(cfg['unit_closure']);units=np.load(bound(load(C/'outputs/accepted_fit_qualification_unit_bridge_v3.json')['unit_codes']));nu=int(units.max())+1
 design=pd.read_csv(Path(cfg['design'])/'outputs/private_analysis_unit_context_and_coverage_v1.tsv',sep='\t');assert len(design)==nu and design.qualification_unit_code.tolist()==list(range(nu))
 contexts={x['code']:x for x in app['contexts']};models={x['context']['code']:joblib.load(bound(x['reference_model'])) for x in trained['models']}
 keys=sorted(next(iter(models.values()))['classifiers']);assert all(sorted(m['classifiers'])==keys for m in models.values());primary=next(iter(models.values()))['primary_key'];schemes=keys+['identity_context_only','primary_confident_projection']
 node_defs={};lookup={};opportunity={};codes={};confidence=np.full(N,np.nan,np.float32)
 for key in schemes:
  nodes=[];mapping={}
  for c in sorted(contexts):
   parent=scope['compartments'].index(contexts[c]['compartment']);mapping[(c,-1)]=len(nodes);nodes.append(dict(code=len(nodes),name='STATE_N_'+str(len(nodes)),parent_code=parent,context_code=c,state_pattern_integer=None,reference_label=contexts[c]['canonical_identity'] or contexts[c]['compartment'],discovery_eligible=key=='identity_context_only',fallback_preserves_original_identity=True))
   if key=='identity_context_only' or c not in models:continue
   partition=primary if key=='primary_confident_projection' else key
   for label in sorted(models[c]['classifiers'][partition].classes_,key=int):
    mapping[(c,int(label))]=len(nodes);nodes.append(dict(code=len(nodes),name='STATE_N_'+str(len(nodes)),parent_code=parent,context_code=c,state_pattern_integer=int(label),reference_label=contexts[c]['canonical_identity'] or contexts[c]['compartment'],discovery_eligible=True,biological_state_function_unassigned=True))
  node_defs[key]=nodes;lookup[key]=mapping;opportunity[key]=np.ones((len(nodes),nu),bool);codes[key]=np.full(N,-1,np.int32)
 source_reports=[]
 for rec in idx['sources']:
  source=rec['source'];with_rows=np.load(bound(rec['cell_axis']),allow_pickle=False);rows=with_rows['fit_row'];cc=ctx[rows];coverage=np.array([x['native_coverage_pass'] for x in rec['coverage_and_interpretation']],bool)
  score=np.load(bound(rec['scores']),mmap_mode='r');source_report=dict(source=source,n_cells=len(rows),contexts=[])
  # Opportunity is defined from assay metadata even when a context has zero
  # observed cells. An absent predicted class cannot establish measurability.
  for c0,model0 in models.items():
   possible0=bool(coverage[[col[n] for n in model0['program_ids']]].all())
   for key in schemes:
    if key=='identity_context_only':continue
    for node in node_defs[key]:
     if node['context_code']==c0 and node['state_pattern_integer'] is not None:opportunity[key][node['code'],design.source.eq(source).to_numpy()]=possible0
  for c in np.unique(cc):
   pos=np.flatnonzero(cc==c);rr=rows[pos];c=int(c)
   for key in schemes:codes[key][rr]=lookup[key][(c,-1)]
   if c not in models:
    source_report['contexts'].append(dict(context_code=c,n_cells=len(rr),status='no_frozen_development_state_model_identity_retained'));continue
   model=models[c];js=np.array([col[n] for n in model['program_ids']]);possible=bool(coverage[js].all())
   for key in schemes:
    if key=='identity_context_only':continue
    for node in node_defs[key]:
     if node['context_code']==c and node['state_pattern_integer'] is not None:opportunity[key][node['code'],design.source.eq(source).to_numpy()]=possible
   if not possible:
    source_report['contexts'].append(dict(context_code=c,n_cells=len(rr),status='missing_native_programme_opportunity_identity_retained'));continue
   weak=0
   for lo in range(0,len(pos),cfg['projection_chunk_size']):
    pp=pos[lo:lo+cfg['projection_chunk_size']];fitrows=rows[pp];yy=np.asarray(score[pp][:,js]);assert np.isfinite(yy).all();pc=model['PCA'].transform(model['scaler'].transform(yy))
    for key in keys:
     pred=model['classifiers'][key].predict(pc).astype(int)
     refpos=np.searchsorted(np.sort(model['reference_fit_rows']),fitrows);ordered=np.argsort(model['reference_fit_rows']);sortedrows=model['reference_fit_rows'][ordered];inref=(refpos<len(sortedrows))
     inref[inref]=sortedrows[refpos[inref]]==fitrows[inref]
     pred[inref]=model['reference_labels'][key][ordered[refpos[inref]]].astype(int)
     codes[key][fitrows]=np.array([lookup[key][(c,int(k))] for k in pred],np.int32)
     if key==primary:
      probability=model['classifiers'][key].predict_proba(pc);labelpos={int(v):i for i,v in enumerate(model['classifiers'][key].classes_)};vote=probability[np.arange(len(pred)),np.array([labelpos[int(v)] for v in pred])];confidence[fitrows]=vote;good=vote>=cfg['minimum_neighbour_agreement'];weak+=int((~good).sum());codes['primary_confident_projection'][fitrows[good]]=np.array([lookup['primary_confident_projection'][(c,int(k))] for k in pred[good]],np.int32)
   source_report['contexts'].append(dict(context_code=c,n_cells=len(rr),status='fixed_development_model_projected',n_weak_neighbour_votes=weak,confidence_is_neighbour_agreement_not_calibrated_probability=True))
  source_reports.append(source_report);save(R/'checkpoints'/('source_'+str(len(source_reports))+'_complete_v1.json'),dict(n_sources=len(source_reports)));storage_gate(R,cfg)
 outputs=[];eligible=common&(parents>=0);assert np.array_equal(eligible,common)
 for key in schemes:
  np.testing.assert_array_equal(codes[key]>=0,common);nodes=node_defs[key];mapping=np.array([n['parent_code'] for n in nodes]);np.testing.assert_array_equal(mapping[codes[key][common]],parents[common])
  count=np.bincount(codes[key][common].astype(np.int64)*nu+units[common],minlength=len(nodes)*nu).reshape(len(nodes),nu);assert (count[~opportunity[key]]==0).all()
  nc=R/'outputs'/('private_'+key+'_state_node_codes_v1.npy');np.save(nc,codes[key]);pc=R/'outputs'/('private_'+key+'_state_parent_codes_v1.npy');parent_axis=np.where(common,parents,-1).astype(np.int16);np.save(pc,parent_axis);oc=R/'outputs'/('private_'+key+'_state_measurement_v1.npy');np.save(oc,opportunity[key])
  manifest=R/'outputs'/('private_'+key+'_accepted_state_membership_v1.json');save(manifest,dict(schema='accepted_TME_identity_scheme_membership_v1',accepted_purpose=cfg['accepted_purpose'],scheme_id='identity_state_'+key,n_fit=N,final_membership_codes_not_display_names=True,unresolved_fine_identity_has_explicit_fallback_node=True,parents=[dict(code=i,name=n) for i,n in enumerate(scope['compartments'])],nodes=nodes,parent_codes=binding(pc),node_codes=binding(nc),node_measurement_mask=binding(oc),input_axis_SHA256=scope['input_axis_SHA256'],scope_acceptance=app['scope_acceptance'],state_model_acceptance=binding(ap),n_assessed_cells=int(common.sum()),n_nodes=len(nodes),original_continuous_vectors_retained=True,no_clinical_outcome_used=True))
  outputs.append(dict(scheme_id='identity_state_'+key,acceptance=binding(manifest),n_nodes=len(nodes),n_discovery_eligible_nodes=sum(n['discovery_eligible'] for n in nodes)))
 co=R/'outputs/private_primary_projection_neighbour_agreement_v1.npy';np.save(co,confidence)
 roles=load(bound(trained['source_roles']));source_role={r['source']:r['role'] for r in roles['source_roles']};assert set(design.source)<=set(source_role);unit_is_development=design.source.map(source_role).eq('development').to_numpy();train=unit_is_development[units];trainp=R/'outputs/private_development_source_training_mask_v1.npy';np.save(trainp,train)
 op=R/'outputs/accepted_state_membership_collection_v1.json';save(op,dict(schema='accepted_TME_annotation_membership_collection_v1',accepted_purpose=cfg['accepted_purpose'],n_fit=N,schemes=outputs,training_cell_mask=dict(**binding(trainp),frozen_before_module_fitting=True),input_axis_SHA256=scope['input_axis_SHA256'],frozen_state_model_acceptance=binding(ap),source_roles=trained['source_roles'],source_projection_reports=source_reports,primary_scheme='identity_state_'+primary,identity_control='identity_state_identity_context_only',neighbour_confidence=binding(co),heldout_sources_only_projected_no_retraining=True,patterns_computational_not_exclusive_functional_states=True,own_programme_profiles_not_independent_state_validation=True,all_cells_exactly_one_membership_and_complete_broad_parent_denominator=True))
 storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='FROZEN_STATE_PATTERNS_PROJECTED_ALL_SOURCES_WITH_EXPLICIT_OPPORTUNITY',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(op)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
