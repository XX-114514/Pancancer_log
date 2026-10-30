"""Freeze whole-unit state increment inputs on exact development panels.

Reconstruct the original precoverage node pools to apply the same eligibility
predicate to source lockboxes. Receiver targets remain masked, never zero-filled.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np,pandas as pd

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
 for dep in ['projection_auditor','composition_auditor','lockbox_readouts','MOFA_auditor']:wait(cfg[dep],start,cfg['dependency_wait_seconds'])
 sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);PA=Path(cfg['projection_auditor']);pp=PA/'outputs/accepted_state_projection_numeric_audit_v1.json';assert sha(pp)==load(PA/'provenance/terminal_v1.json')['acceptance_SHA256'];pma=load(pp);collection=load(bound(pma['producer_acceptance']));model=load(bound(collection['frozen_state_model_acceptance']));app=load(bound(model['applicability']));scope=load(bound(app['scope_acceptance']));mi=load(bound(model['input_acceptance']));roles=load(bound(model['source_roles']));role={r['source']:r['role'] for r in roles['source_roles']}
 MA=Path(cfg['MOFA_auditor']);mp=MA/'outputs/accepted_state_MOFA_numeric_comparison_v1.json';assert sha(mp)==load(MA/'provenance/terminal_v1.json')['acceptance_SHA256'];mfa=load(mp);mf=load(bound(mfa['producer_acceptance']));assert mf['input_acceptance']==model['input_acceptance']
 LA=Path(cfg['lockbox_readouts']);lp=LA/'outputs/accepted_frozen_lockbox_state_readouts_v1.json';assert sha(lp)==load(LA/'provenance/terminal_v1.json')['acceptance_SHA256'];lock=load(lp);assert lock['input_acceptance']==model['input_acceptance'];candidate=load(bound(mi['candidate_acceptance']));dev=load(bound(candidate['readout_acceptance']));assert dev['source_roles']==model['source_roles']
 CA=Path(cfg['composition_auditor']);cap=CA/'outputs/accepted_comparable_CM_inputs_v1.json';assert sha(cap)==load(CA/'provenance/terminal_v1.json')['acceptance_SHA256'];accepted=load(cap);prep=Path(load(CA/'config/spec_v1.json')['preparation']);spec=load(prep/'config/spec_v1.json');common=np.load(prep/'outputs/private_common_assessed_cell_mask_v1.npy');C=Path(spec['unit_closure']);units=np.load(bound(load(C/'outputs/accepted_fit_qualification_unit_bridge_v3.json')['unit_codes']));D=Path(spec['design']);meta=pd.read_csv(D/'outputs/private_analysis_unit_context_and_coverage_v1.tsv',sep='\t');nu=len(meta);assert meta.qualification_unit_code.tolist()==list(range(nu));parent_axis=np.load(bound(scope['codes']));nparent=len(scope['compartments']);parent_count=np.bincount(parent_axis[common].astype(np.int64)*nu+units[common],minlength=nparent*nu).reshape(nparent,nu);common_count=parent_count.sum(0);unit_role=meta.source.map(role).to_numpy();assert set(unit_role)<={'development','lockbox'}
 schemes={}
 for b in spec['schemes']:
  a=load(bound(b));node=np.load(bound(a['node_codes']));par=np.load(bound(a['parent_codes']));np.testing.assert_array_equal(par[common],parent_axis[common]);count=np.bincount(node[common].astype(np.int64)*nu+units[common],minlength=len(a['nodes'])*nu).reshape(len(a['nodes']),nu);op=np.load(bound(a['node_measurement_mask']));assert (count[~op]==0).all()
  mapping=np.array([n['parent_code'] for n in a['nodes']]);den=parent_count[mapping];frequency=np.divide(count,den,out=np.full(count.shape,np.nan,float),where=den>=spec['minimum_parent_cells'])
  for p in range(nparent):np.testing.assert_array_equal(count[mapping==p].sum(0),parent_count[p])
  schemes[a['scheme_id']]=dict(manifest=a,binding=b,count=count,opportunity=op,frequency=frequency,mapping=mapping)
 primary=collection['primary_scheme'];identity='CT_identity_cluster';assert primary in schemes and identity in schemes
 programme_names=dev['program_ids'];assert lock['program_ids']==programme_names;readout_sources=dev['readout_sources']+lock['sources'];assert len({s['source'] for s in readout_sources})==len(readout_sources)
 mi_view={v['context']['code']:v for v in mi['views']};view_arrays={}
 for v in mf['views']:
  c=int(v['name'].removeprefix('context_'));pids=v['features'];view_arrays[c]=dict(program_ids=pids,means=np.full((nu,len(pids)),np.nan,np.float32),split_A=np.full((nu,len(pids)),np.nan,np.float32),split_B=np.full((nu,len(pids)),np.nan,np.float32),n_cells=np.zeros(nu,np.int32),mean_log1p_RNA=np.full(nu,np.nan),mean_detected_native_genes=np.full(nu,np.nan),parent_code=scope['compartments'].index(mi_view[c]['context']['compartment']))
 for record in readout_sources:
  with np.load(bound(record['readouts']),allow_pickle=False) as z:
   u=z['unit_code'];cc=z['context_code'];n=z['n_cells'];means=z['means'];isdev=role[record['source']]=='development';A=z['split_means_A'][0] if isdev else z['split_mean_A'];B=z['split_means_B'][0] if isdev else z['split_mean_B']
   for c,a in view_arrays.items():
    take=np.flatnonzero(cc==c);ix=u[take];assert (a['n_cells'][ix]==0).all();a['n_cells'][ix]=n[take];a['mean_log1p_RNA'][ix]=z['mean_log1p_RNA'][take];a['mean_detected_native_genes'][ix]=z['mean_detected_native_genes'][take];columns=[programme_names.index(p) for p in a['program_ids']];release=n[take]>=cfg['minimum_receiver_context_cells'];i=ix[release];rr=take[release];a['means'][i]=means[rr][:,columns];a['split_A'][i]=A[rr][:,columns];a['split_B'][i]=B[rr][:,columns]
 # Independent all-parent depth backgrounds are measured on the same common
 # cell domain; these do not consume receiver state programme responses.
 numeric=load(bound(model['score_numerical_acceptance']));score_index=load(bound(numeric['producer_index']));H=Path(cfg['model']);totals=np.load(H/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r');sumrna=np.zeros((nparent,nu));sumgenes=np.zeros_like(sumrna)
 for source in score_index['sources']:
  with np.load(bound(source['cell_axis']),allow_pickle=False) as z:r=z['fit_row'];u=z['qualification_unit_code'];ng=z['detected_native_genes']
  key=parent_axis[r].astype(np.int64)*nu+u;sumrna+=np.bincount(key,weights=np.log1p(totals[r]),minlength=nparent*nu).reshape(nparent,nu);sumgenes+=np.bincount(key,weights=ng,minlength=nparent*nu).reshape(nparent,nu)
 meanrna=np.divide(sumrna,parent_count,out=np.full(sumrna.shape,np.nan),where=parent_count>0);meangenes=np.divide(sumgenes,parent_count,out=np.full(sumgenes.shape,np.nan),where=parent_count>0)
 cp=R/'outputs/private_frozen_state_increment_composition_depth_v1.npz';np.savez_compressed(cp,unit_code=np.arange(nu),parent_counts=parent_count,identity_counts=schemes[identity]['count'],primary_state_counts=schemes[primary]['count'],mean_parent_log1p_RNA=meanrna,mean_parent_detected_native_genes=meangenes,common_counts=common_count,development=unit_role=='development',lockbox=unit_role=='lockbox',candidate_block=meta.patient_candidate_block.to_numpy(dtype=str))
 view_outputs=[]
 for c,a in view_arrays.items():
  p=R/'outputs'/('private_context_'+str(c)+'_receiver_state_targets_v1.npz');np.savez_compressed(p,**{k:v for k,v in a.items() if k not in ['program_ids','parent_code']});view_outputs.append(dict(context=mi_view[c]['context'],parent_code=a['parent_code'],program_ids=a['program_ids'],targets=binding(p)))
 panels=[]
 for panel in spec['panels']:
  base=np.load(D/'outputs'/('private_'+panel['material_panel']+'_unit_mask_v1.npy'))&(common_count>=spec['minimum_total_cells']);parents=[scope['compartments'].index(p) for p in panel['parents']];base&=(parent_count[parents]>=spec['minimum_parent_cells']).all(0);preliminary=base&(unit_role=='development');eligible=base.copy();pools=[]
  for sid,s in schemes.items():
   pool=[]
   for j,node in enumerate(s['manifest']['nodes']):
    if node['parent_code'] not in parents or not node['discovery_eligible']:continue
    observed=preliminary&s['opportunity'][j];f=s['frequency'][j,observed]
    if s['count'][j,observed].sum()>=spec['minimum_node_cells'] and (s['count'][j,observed]>0).sum()>=spec['minimum_node_observed_units'] and observed.sum() and np.ptp(f)>0:pool.append(j)
   if pool:eligible&=s['opportunity'][pool].all(0)
   pools.append(dict(scheme=sid,n_precoverage_nodes=len(pool)))
  train=eligible&(unit_role=='development');held=eligible&(unit_role=='lockbox');primary_input=next(i for i in accepted['inputs'] if i['panel']==panel['id'] and i['scheme']==primary);pa=load(bound(primary_input['acceptance']));pm=pd.read_csv(bound(pa['metadata']),sep='\t');np.testing.assert_array_equal(np.flatnonzero(train),pm.qualification_unit_code.to_numpy())
  fq=pd.read_csv(bound(pa['frequency']),sep='\t',index_col=0);by={n['name']:n['code'] for n in schemes[primary]['manifest']['nodes']};nodecodes=np.array([by[n] for n in fq.index]);np.testing.assert_allclose(schemes[primary]['frequency'][np.ix_(nodecodes,train)],fq.iloc[:,1:].to_numpy(float),rtol=1e-12,atol=1e-14)
  observed=schemes[primary]['frequency'][np.ix_(nodecodes,train)];minimum=observed.min(1);span=observed.max(1)-minimum;assert (span>0).all();p=R/'outputs'/('private_'+panel['id']+'_increment_unit_masks_v1.npz');np.savez_compressed(p,development=train,lockbox=held,requested_parent_codes=parents,primary_fitted_node_codes=nodecodes,primary_minmax_minimum=minimum,primary_minmax_span=span)
  panels.append(dict(panel=panel['id'],n_development_units=int(train.sum()),n_lockbox_units=int(held.sum()),n_development_sources=int(meta.loc[train,'source'].nunique()),n_lockbox_sources=int(meta.loc[held,'source'].nunique()),original_development_panel_exactly_reproduced=True,identical_eligibility_predicate_applied_to_lockbox=True,precoverage_pool_counts=pools,masks=binding(p),primary_CM_input_acceptance=primary_input['acceptance']))
 op=R/'outputs/accepted_state_increment_prediction_inputs_v1.json';save(op,dict(schema='TME_state_increment_prediction_inputs_v1',projection_audit=binding(pp),composition_audit=binding(cap),continuous_state_audit=binding(mp),lockbox_readout_acceptance=binding(lp),n_work_units=nu,source_roles=model['source_roles'],metadata=binding(D/'outputs/private_analysis_unit_context_and_coverage_v1.tsv'),composition_and_depth=binding(cp),identity_membership=schemes[identity]['binding'],primary_state_membership=schemes[primary]['binding'],parents=scope['compartments'],receiver_contexts=view_outputs,panels=panels,minimum_receiver_context_cells=cfg['minimum_receiver_context_cells'],missing_targets_retained_as_NaN=True,all_identity_children_including_low_abundance_fallbacks_retained=True,no_receiver_state_values_enter_sender_composition=True,heldout_targets_not_used_to_select_features_or_parameters=True,clinical_outcomes_used=False,predictive_comparison_not_yet_fitted=True))
 storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='FROZEN_STATE_MODULE_INCREMENT_PREDICTION_INPUTS_PREPARED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(op)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
