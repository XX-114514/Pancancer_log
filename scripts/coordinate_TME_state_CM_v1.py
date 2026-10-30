"""Finite state-partition CoVarNet pipeline after frozen development training.

The coordinator creates immutable, individually supervised scientific stages.
It never reruns an existing root or changes rank/arms from observed effects.
"""
import argparse,datetime,hashlib,json,os,subprocess,sys,time
from pathlib import Path

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def binding(p):return dict(path=str(p),SHA256=sha(p))
def bound(b):
 p=Path(b['path']);assert sha(p)==b['SHA256'];return p
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def wait(root,deadline):
 p=Path(root)/'provenance/supervisor_terminal_v1.json'
 while not p.exists():
  assert time.monotonic()<deadline;time.sleep(5)
 t=load(p);assert t['actual_child_exit_code']==0 and t.get('science_terminal_exists',True) and not (Path(root)/'provenance/failure_v1.json').exists(),'A dependency failed; no automatic retry'
def stage(R,name,entry,files,cfg,wall,deadline):
 root=Path(R.parent)/name;assert not root.exists(),'Observe existing root; never overwrite or restart'
 root.mkdir(mode=0o700)
 for n in ['scripts','config','outputs','checkpoints','provenance']:(root/n).mkdir()
 hashes={}
 for n in sorted(set(files+['supervise_atlas_stage_v1.py'])):
  p=root/'scripts'/n;p.write_bytes((R/'scripts'/n).read_bytes());hashes['scripts/'+n]=sha(p)
 save(root/'config/spec_v1.json',cfg);save(root/'config/supervisor_request_v1.json',dict(entrypoint=entry,arguments=['--output',str(root)],wall_seconds=wall,script_SHA256=hashes))
 with (root/'provenance/supervisor_stdout.txt').open('x') as o,(root/'provenance/supervisor_stderr.txt').open('x') as e:p=subprocess.Popen([sys.executable,'-B',str(root/'scripts/supervise_atlas_stage_v1.py'),'--output',str(root)],stdout=o,stderr=e,start_new_session=True,env=os.environ.copy())
 ticks=(Path('/proc')/str(p.pid)/'stat').read_text().rsplit(')',1)[1].split()[19];save(root/'provenance/coordinator_launch_v1.json',dict(UTC=now(),PID=p.pid,startticks=ticks,coordinator=str(R)))
 wait(root,deadline);return root

def main(R):
 cfg=load(R/'config/spec_v1.json');deadline=time.monotonic()+cfg['pipeline_wall_seconds'];wait(cfg['projection'],deadline)
 sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);P=Path(cfg['projection']);pp=P/'outputs/accepted_state_membership_collection_v1.json';assert sha(pp)==load(P/'provenance/terminal_v1.json')['acceptance_SHA256'];m=load(pp)
 assert m['heldout_sources_only_projected_no_retraining'] and m['all_cells_exactly_one_membership_and_complete_broad_parent_denominator']
 original=load(Path(cfg['identity_memberships'])/'outputs/accepted_annotation_membership_collection_v1.json');anchor=next(x for x in original['schemes'] if x['scheme_id']=='CT_identity_cluster')
 schemes=[x['acceptance'] for x in m['schemes']]+[anchor['acceptance']];pc=dict(cfg['panel_base_config']);pc.update(schemes=schemes,training_cell_mask=m['training_cell_mask'])
 prep=stage(R,cfg['prefix']+'_inputs_v1','prepare_TME_CM_comparable_panels_v1.py',['prepare_TME_CM_comparable_panels_v1.py','TME_CM_comparison_contract_v1.py'],pc,7200,deadline)
 ac=dict(cfg['audit_base_config'],preparation=str(prep),binding_mode='frozen_state_projection',state_projection=str(P))
 audit=stage(R,cfg['prefix']+'_inputs_audit_v1','audit_TME_CM_comparable_inputs_v2.py',['audit_TME_CM_comparable_inputs_v2.py'],ac,7200,deadline)
 ap=audit/'outputs/accepted_comparable_CM_inputs_v1.json';assert sha(ap)==load(audit/'provenance/terminal_v1.json')['acceptance_SHA256'];inputs=load(ap)['inputs']
 if not inputs:
  op=R/'outputs/accepted_state_CM_pipeline_result_v1.json';save(op,dict(status='NO_COMPLETE_COMPARABLE_STATE_CM_PANEL',preparation=binding(prep/'outputs/comparable_CM_panel_preparation_v1.json'),auditor=binding(ap),no_missing_to_zero_or_unapproved_scope_fallback=True))
  save(R/'provenance/terminal_v1.json',dict(status='STATE_CM_PANELS_NONESTIMABLE_REASONS_RETAINED',completed_UTC=now(),acceptance_SHA256=sha(op)));return
 # The first primary fit must independently accept before finite expansion.
 primary=next(x for x in inputs if x['panel']==cfg['first_panel'] and x['scheme']==m['primary_scheme']);fc=json.loads(json.dumps(cfg['base_fit_config']));fc['input_acceptance']=primary['acceptance'];fc['R_parameters']['maximum_rank']=primary['maximum_rank']
 first_name=cfg['prefix']+'_first_primary_v1';fit_files=['execute_TME_CM_scheme_v1.py','run_TME_CoVarNet_scheme_v1.R','frozen_author_CoVarNet_v1.R']
 first=stage(R,first_name,'execute_TME_CM_scheme_v1.py',fit_files,fc,14400,deadline)
 requests=[];runs_by_panel={};skipped=[]
 for i,item in enumerate(inputs):
  name=first_name if item is primary else cfg['prefix']+'_own_'+str(i)+'_v1';runs_by_panel.setdefault(item['panel'],[]).append(str(R.parent/name))
  if item is not primary:requests.append(dict(stage_name=name,panel=item['panel'],scheme=item['scheme'],rank_mode='own_survey_rank'))
 common=[]
 for i,item in enumerate(inputs):
  if item['panel']!=cfg['first_panel']:continue
  if item['maximum_rank']<cfg['common_rank']:
   skipped.append(dict(panel=item['panel'],scheme=item['scheme'],reason='frozen_common_rank_not_identifiable'));continue
  name=cfg['prefix']+'_commonK'+str(cfg['common_rank'])+'_'+str(i)+'_v1';common.append(str(R.parent/name));requests.append(dict(stage_name=name,panel=item['panel'],scheme=item['scheme'],rank_mode='fixed_common_rank'))
 qc=dict(cfg['queue_base_config']);qc.update(initial_fits=[str(first)],input_auditor=str(audit),requests=requests,runs_parent=str(R.parent),base_fit_config=cfg['base_fit_config'],fixed_common_rank=cfg['common_rank'])
 queue=stage(R,cfg['prefix']+'_finite_fit_queue_v1','run_TME_CM_frozen_fit_queue_v1.py',['run_TME_CM_frozen_fit_queue_v1.py']+fit_files,qc,86400,deadline)
 comparisons=[]
 if len(common)==len(runs_by_panel.get(cfg['first_panel'],[])):runs_by_panel['prespecified_common_rank']=common
 for i,(panel,fits) in enumerate(sorted(runs_by_panel.items())):
  cc=dict(cfg['comparison_base_config'],preparation=str(prep),fits=fits,comparison_role='development_state_patterns_vs_identity_'+panel)
  comp=stage(R,cfg['prefix']+'_comparison_'+str(i)+'_v1','compare_TME_CM_membership_modules_v1.py',['compare_TME_CM_membership_modules_v1.py','TME_CM_comparison_contract_v1.py'],cc,7200,deadline)
  comparisons.append(binding(comp/'outputs/accepted_membership_module_comparison_v1.json'))
 op=R/'outputs/accepted_state_CM_pipeline_result_v1.json';save(op,dict(status='STATE_PARTITION_COVARNET_FITS_AND_COMPARISONS_COMPLETE',projection=binding(pp),input_auditor=binding(ap),first_fit=binding(first/'outputs/numeric_readback_acceptance_v1.json'),queue_acceptance=binding(queue/'outputs/completed_frozen_fit_queue_v1.json'),comparisons=comparisons,skipped_fixed_rank=skipped,source_holdout_and_incremental_validation_not_replaced=True,clinical_outcomes_used=False))
 storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='STATE_PATTERN_CM_EXECUTION_COMPLETE_INCREMENT_VALIDATION_PENDING',completed_UTC=now(),acceptance_SHA256=sha(op)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
