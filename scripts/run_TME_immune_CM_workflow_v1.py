"""Finite immune-only CM workflow; frozen TME labels and author fitting unchanged."""
import argparse,datetime,hashlib,json,os,subprocess,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
from TME_CM_factor_label_contract_v1 import official_factor_labels


def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bound(b):
 p=Path(b['path']);assert sha(p)==b['SHA256'];return p

def save(p,x):
 with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False)
def binding(p):return dict(path=str(p),SHA256=sha(p))
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()


def main(R):
 cfg=load(R/'config/spec_v1.json');P=R.parent;deadline=time.monotonic()+cfg['wall_seconds']
 sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg)
 def wait(root):
  p=root/'provenance/supervisor_terminal_v1.json'
  while not p.exists():
   assert time.monotonic()<deadline;time.sleep(5)
  assert load(p)['actual_child_exit_code']==0 and load(p)['science_terminal_exists'],str(root)
 def launch(suffix,entry,c,files,source=None,wall=23000):
  root=P/(cfg['prefix']+'_'+suffix+'_v1');assert not root.exists(),'Inspect existing exact stage; never overwrite/restart'
  root.mkdir(mode=0o700)
  for name in ['scripts','config','outputs','checkpoints','provenance']:(root/name).mkdir()
  source=source or R
  for name in files:
   (root/'scripts'/name).write_bytes((source/'scripts'/name).read_bytes())
  if 'supervise_atlas_stage_v1.py' not in files:(root/'scripts/supervise_atlas_stage_v1.py').write_bytes((R/'scripts/supervise_atlas_stage_v1.py').read_bytes())
  save(root/'config/spec_v1.json',c)
  hashes={'scripts/'+p.name:sha(p) for p in (root/'scripts').iterdir() if p.is_file()}
  save(root/'config/supervisor_request_v1.json',dict(entrypoint=entry,arguments=['--output',str(root)],wall_seconds=wall,script_SHA256=hashes))
  env=os.environ.copy();env.update(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2',PYTHONDONTWRITEBYTECODE='1')
  with (root/'provenance/supervisor_stdout.txt').open('x') as o,(root/'provenance/supervisor_stderr.txt').open('x') as e:
   proc=subprocess.Popen([sys.executable,'-B',str(root/'scripts/supervise_atlas_stage_v1.py'),'--output',str(root)],env=env,stdout=o,stderr=e,start_new_session=True)
  ticks=(Path('/proc')/str(proc.pid)/'stat').read_text().rsplit(')',1)[1].split()[19]
  save(root/'provenance/workflow_launch_v1.json',dict(PID=proc.pid,startticks=ticks,UTC=now(),parent=str(R)))
  (R/'checkpoints/current_stage_v1.json').write_text(json.dumps(dict(UTC=now(),root=str(root),entry=entry)))
  return root
 base=dict(joined=cfg['joined'],threads=2,quota_group=cfg['quota_group'],output_cap_bytes=128*1024**2,workspace_reserve_bytes=cfg['workspace_reserve_bytes'])
 original=Path(cfg['original_preparation']);pc=load(original/'config/spec_v1.json')
 dc=dict(base,original_preparation=str(original),original_design=pc['design'],scope_policy=cfg['scope_policy'],scope_policy_SHA256=cfg['scope_policy_SHA256'])
 design=launch('design','prepare_TME_immune_only_design_v1.py',dc,['prepare_TME_immune_only_design_v1.py']);wait(design)
 pp=dict(pc,design=str(design),panels=cfg['scope_policy']['panels'],output_cap_bytes=256*1024**2)
 prep=launch('preparation','prepare_TME_CM_comparable_panels_v1.py',pp,['prepare_TME_CM_comparable_panels_v1.py','TME_CM_comparison_contract_v1.py']);wait(prep)
 ac=dict(base,preparation=str(prep),SingleR_auditor=pc['SingleR_auditor'])
 auditor=launch('input_audit','audit_TME_CM_comparable_inputs_v2.py',ac,['audit_TME_CM_comparable_inputs_v2.py']);wait(auditor)
 ap=auditor/'outputs/accepted_comparable_CM_inputs_v1.json';assert sha(ap)==load(auditor/'provenance/terminal_v1.json')['acceptance_SHA256'];inputs=load(ap)['inputs']
 assert inputs,'No immune panel meets frozen coverage; report without fitting'
 main_panel='Immune_T_NK_B_Plasma_Myeloid';ct='CT_identity_cluster'
 primary=next((x for x in inputs if x['panel']==main_panel and x['scheme']==ct),None)
 assert primary is not None,'Main immune panel nonestimable: preserve scope and inspect coverage'
 fitcfg=load(Path(cfg['baseline_fit'])/'config/spec_v1.json');fitcfg['input_acceptance']=primary['acceptance'];fitcfg['R_parameters']['maximum_rank']=primary['maximum_rank'];fitcfg['R_parameters'].pop('rank_override',None)
 fitfiles=['execute_TME_CM_scheme_v1.py','run_TME_CoVarNet_scheme_v1.R','frozen_author_CoVarNet_v1.R']
 pilot=launch('main_CT_identity_own','execute_TME_CM_scheme_v1.py',fitcfg,fitfiles);wait(pilot)
 assert (pilot/'outputs/numeric_readback_acceptance_v1.json').exists()
 # Requests were frozen before outcomes; coverage alone determines whether a request can run.
 requests=[];not_estimable=[];mainfits=[pilot];commonfits=[]
 for req in cfg['fit_requests']:
  if req['suffix']=='main_CT_identity_own':continue
  item=next((x for x in inputs if x['panel']==req['panel'] and x['scheme']==req['scheme']),None)
  if item is None or (req['rank_mode']=='fixed_common_rank' and item['maximum_rank']<4):
   not_estimable.append(dict(request=req,reason='predeclared_input_coverage_or_rank_bound'));continue
  target=P/(cfg['prefix']+'_'+req['suffix']+'_v1')
  requests.append(dict(stage_name=target.name,panel=req['panel'],scheme=req['scheme'],rank_mode=req['rank_mode']))
  if req['panel']==main_panel:(commonfits if req['rank_mode']=='fixed_common_rank' else mainfits).append(target)
 queuecfg=dict(base,base_fit_config=fitcfg,initial_fits=[str(pilot)],input_auditor=str(auditor),requests=requests,runs_parent=str(P),queue_wall_seconds=cfg['wall_seconds']-1800,per_fit_supervisor_wall_seconds=23000,maximum_concurrent_fits=2,fixed_common_rank=4,frozen_before_module_result_review=True,clinical_outcomes_used=False)
 queue=launch('finite_queue','run_TME_CM_frozen_fit_queue_v1.py',queuecfg,['run_TME_CM_frozen_fit_queue_v1.py',*fitfiles],wall=cfg['wall_seconds']-900);wait(queue)
 roots=[pilot]+[P/x['stage_name'] for x in requests]
 comparisons=[]
 for tag,group in [('own_rank',mainfits),('commonK4',commonfits)]:
  if len(group)!=5:
   not_estimable.append(dict(comparison=tag,reason='not_all5schemes_comparable'));continue
  cc=dict(base,accepted_purpose=pc['accepted_purpose'],preparation=str(prep),fits=[str(x) for x in group],minimum_cosine_similarity=.6,comparison_role='immune_only_frozen_scope_'+tag,dependency_wait_seconds=3600,output_cap_bytes=512*1024**2)
  comp=launch('comparison_'+tag,'compare_TME_CM_membership_modules_v1.py',cc,['compare_TME_CM_membership_modules_v1.py','TME_CM_comparison_contract_v1.py'],wall=4000);wait(comp);comparisons.append(str(comp))
 catalogue=[];fit_summaries=[]
 for root in roots:
  wait(root);fc=load(root/'config/spec_v1.json');a=load(bound(fc['input_acceptance']));m=load(bound(a['annotation_acceptance']));nodes={v['name']:v for v in m['nodes']};parents=m['parents'];coverage=a['coverage']
  fq=pd.read_csv(bound(a['frequency']),sep='\t',index_col=0);meta=pd.read_csv(bound(a['metadata']),sep='\t')
  assert set(fq.majorCluster)<=set(cfg['scope_policy']['allowed_parents'])
  fit_summaries.append(dict(root=str(root),panel=a['requested_parent_panel']['id'],scheme=m['scheme_id'],rank=load(root/'outputs/CoVarNet_execution_report_v1.json')['selected_rank'],n_work_units=len(meta),n_sources=int(meta.source.nunique()),n_nodes=len(fq)))
  for seed in fc['R_parameters']['final_seeds']:
   d=root/'outputs'/('seed_'+str(seed));W=pd.read_csv(d/'private_W_v1.tsv',sep='\t',index_col=0);A=pd.read_csv(d/'private_A_v1.tsv',sep='\t',index_col=0);H=pd.read_csv(d/'private_H_v1.tsv',sep='\t',index_col=0);labels=official_factor_labels(W,A,H);assert H.columns.tolist()==meta.sampleID.tolist()
   try:edges=pd.read_csv(d/'author_compatibility_each_edge_v1.tsv',sep='\t')
   except pd.errors.EmptyDataError:edges=pd.DataFrame(columns=['cm','subCluster1','subCluster2'])
   assert set(edges.cm.astype(str))<=set(labels)
   for label in labels:
    e=edges.loc[edges.cm.astype(str)==label];pairs=sorted({tuple(sorted([r['subCluster1'],r['subCluster2']])) for r in e.to_dict('records')});members=sorted({v for p in pairs for v in p});descriptions=[]
    for name in members:
     v=nodes[name];parent=parents[v['parent_code']]['name'];assert parent in cfg['scope_policy']['allowed_parents'];descriptions.append(dict(node=name,parent=parent,reference_label=v.get('reference_label'),numeric_cluster=v.get('numeric_cluster'),parent_marker=v.get('parent_cluster_marker'),marker_borrowed_from_parent=v.get('parent_marker_borrowed',False),own_node_DEG_not_claimed=v.get('own_node_DEG_marker_not_yet_assigned',True)))
    catalogue.append(dict(version='2026-11-01_qualified_TME_immune_only',root=root.name,scheme=m['scheme_id'],panel=a['requested_parent_panel']['id'],seed=seed,module=label,edges=pairs,members=descriptions,n_work_units=len(meta),n_sources=int(meta.source.nunique()),empty_author_network=not pairs,shared_edges_not_independent_evidence=True))
 # Source/context evidence reuses the accepted block-bootstrap implementation.
 sr=P/'20261030_TME_CM_source_context_recurrence_v1';sc=load(sr/'config/spec_v1.json');sc['schemes']=[dict(id=x['panel']+'__'+x['scheme'],root=x['root']) for x in fit_summaries if 'commonK4' not in x['root']]
 sc.update(dependency_wait_seconds=cfg['wall_seconds']-900,output_cap_bytes=256*1024**2)
 files=[p.name for p in (sr/'scripts').iterdir() if p.is_file() and p.suffix=='.py']
 recurrence=launch('source_context','measure_TME_CM_source_context_recurrence_v1.py',sc,files,source=sr,wall=cfg['wall_seconds']-900);wait(recurrence)
 cp=R/'outputs/private_immune_CM_factor_member_catalogue_v1.json';save(cp,catalogue)
 accepted=R/'outputs/accepted_immune_CM_workflow_v1.json';save(accepted,dict(schema='TME_immune_only_CM_complete_workflow_v1',completed_UTC=now(),scope_policy=cfg['scope_policy'],immune_design=binding(design/'outputs/accepted_immune_only_design_v1.json'),input_audit=binding(ap),fits=fit_summaries,n_formal_fits=len(roots),n_formal_seed_outputs=2*len(roots),n_factor_records=len(catalogue),all_network_members_immune_only=True,not_estimable=not_estimable,comparisons=comparisons,source_recurrence=str(recurrence),catalogue=binding(cp),same_author_algorithm_and_previous_rank_variant_preserved=True,no_new_clinical_LR_or_state_function_claim=True))
 storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='IMMUNE_ONLY_COVARNET_WORKFLOW_NUMERICALLY_ACCEPTED',completed_UTC=now(),acceptance_SHA256=sha(accepted)))


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',{'traceback':traceback.format_exc()})
  raise
