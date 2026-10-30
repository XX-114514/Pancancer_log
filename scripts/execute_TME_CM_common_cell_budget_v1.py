"""Supervise official nsNMF fits on audited identical-cell budget inputs."""
import argparse,datetime,hashlib,json,os,subprocess,sys,time
from pathlib import Path

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def main(R):
 cfg=load(R/'config/spec_v1.json');start=time.monotonic();P=Path(cfg['inputs'])
 while not (P/'provenance/supervisor_terminal_v1.json').exists():
  assert time.monotonic()-start<cfg['dependency_wait_seconds'];time.sleep(5)
 assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 ip=P/'outputs/accepted_common_cell_budget_inputs_v1.json';assert sha(ip)==load(P/'provenance/terminal_v1.json')['acceptance_SHA256'];idx=load(ip)
 assert idx['identical_draws_across_annotation_schemes'] and idx['all_fallbacks_retained_in_complete_parent_denominator'];digests={str(ip):sha(ip)}
 for rep in idx['repeats']:
  for b in [rep['same_cell_draw']]+[s['frequency'] for s in rep['schemes']]:assert sha(b['path'])==b['SHA256'];digests[b['path']]=b['SHA256']
 sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);request=dict(cfg,input_index=str(ip),math_source=str(R/'scripts/run_TME_CM_block_source_robustness_v1.R'),output=str(R/'outputs'));rp=R/'config/frozen_R_request_v1.json';save(rp,request)
 cmd=[cfg['Rscript'],'--vanilla',str(R/'scripts/run_TME_CM_common_cell_budget_v1.R'),str(rp)]
 with (R/'provenance/private_R_stdout_v1.txt').open('x') as o,(R/'provenance/private_R_stderr_v1.txt').open('x') as e:
  c=subprocess.Popen(cmd,stdout=o,stderr=e,start_new_session=True);proc=Path('/proc')/str(c.pid)/'stat';ticks=proc.read_text().rsplit(')',1)[1].split()[19];save(R/'provenance/R_child_launch_v1.json',dict(PID=c.pid,startticks=ticks))
  try:code=c.wait(timeout=cfg['R_wall_seconds'])
  except subprocess.TimeoutExpired:
   import signal
   assert proc.exists() and proc.read_text().rsplit(')',1)[1].split()[19]==ticks;os.killpg(c.pid,signal.SIGTERM);code=c.wait(timeout=30)
  save(R/'provenance/R_child_terminal_v1.json',dict(actual_exit_code=code))
 assert code==0
 assert all(sha(p)==h for p,h in digests.items());report=R/'outputs/CM_same_cell_budget_report_v1.json';assert load(report)['status']=='SAME_CELL_DRAWS_ALL_SCHEMES_COVARNET_ROBUSTNESS_COMPLETE';storage_gate(R,cfg)
 save(R/'provenance/terminal_v1.json',dict(status='COMMON_CELL_BUDGET_COVARNET_REFITS_AND_MATCHING_COMPLETE',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(report),all_input_files_unchanged=True))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
