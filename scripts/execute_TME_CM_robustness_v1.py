"""Bounded official nsNMF resampling with exact frozen input bindings."""
import argparse,datetime,hashlib,json,os,subprocess,sys
from pathlib import Path

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def main(R):
 cfg=load(R/'config/spec_v1.json');sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);inputs={}
 test=R/'outputs/robustness_math_synthetic_acceptance_v1.json'
 check=subprocess.run([cfg['Rscript'],'--vanilla',str(R/'scripts/check_TME_CM_robustness_math_v1.R'),str(R/'scripts/run_TME_CM_block_source_robustness_v1.R'),str(test)],capture_output=True,text=True,timeout=180)
 save(R/'provenance/math_synthetic_terminal_v1.json',dict(actual_exit_code=check.returncode,stderr=check.stderr[-4000:]))
 assert check.returncode==0 and load(test)['all_tests_passed'],'Robustness math/API test failed before human fits'
 for item in cfg['schemes']:
  root=Path(item['root']);assert load(root/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
  a=load(root/'config/spec_v1.json')['input_acceptance'];assert sha(a['path'])==a['SHA256'];m=load(a['path'])
  for name in ['frequency','metadata']:
   b=m[name];assert sha(b['path'])==b['SHA256'];inputs[b['path']]=b['SHA256']
 req=dict(cfg,output=str(R/'outputs'));p=R/'config/frozen_R_request_v1.json';save(p,req)
 cmd=[cfg['Rscript'],'--vanilla',str(R/'scripts/run_TME_CM_block_source_robustness_v1.R'),str(p)]
 with (R/'provenance/private_R_stdout_v1.txt').open('x') as o,(R/'provenance/private_R_stderr_v1.txt').open('x') as e:
  child=subprocess.Popen(cmd,stdout=o,stderr=e,start_new_session=True);pr=Path('/proc')/str(child.pid)/'stat';ticks=pr.read_text().rsplit(')',1)[1].split()[19];save(R/'provenance/R_child_launch_v1.json',dict(PID=child.pid,startticks=ticks))
  try:code=child.wait(timeout=cfg['R_wall_seconds'])
  except subprocess.TimeoutExpired:
   import signal
   assert pr.exists() and pr.read_text().rsplit(')',1)[1].split()[19]==ticks;os.killpg(child.pid,signal.SIGTERM);code=child.wait(timeout=30)
  save(R/'provenance/R_child_terminal_v1.json',dict(actual_exit_code=code))
 assert code==0
 for p,h in inputs.items():assert sha(p)==h
 report=R/'outputs/CM_block_source_robustness_report_v1.json';assert load(report)['status']=='FROZEN_MATRIX_CM_BLOCK_AND_SOURCE_REFITS_COMPLETE'
 storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='CM_FROZEN_INPUT_BLOCK_AND_SOURCE_ROBUSTNESS_COMPLETE',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(report),input_files_unchanged=True))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
