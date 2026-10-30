"""Task-owned direct-GPU gate around the unchanged official scVI recipe."""
import argparse,fcntl,json,os,subprocess,time
from pathlib import Path
import train_TME_scVI_GPU_v3 as standard
def main(R):
 cfg=standard.load(R/'config/spec_v1.json')
 for key in ['producer','auditor']:
  P=Path(cfg[key]);end=P/'provenance/supervisor_terminal_v1.json'
  while not end.exists():time.sleep(5)
  assert standard.load(end)['actual_child_exit_code']==0
 for name,digest in cfg['standard_adapter_code_SHA256'].items():assert standard.sha(R/'scripts'/name)==digest
 result=subprocess.run(['nvidia-smi','--query-compute-apps=pid','--format=csv,noheader,nounits'],capture_output=True,text=True,check=True,timeout=15);other=[int(s.strip()) for s in result.stdout.splitlines() if s.strip() and int(s.strip())!=os.getpid()];assert not other,'GPU has another process; do not interrupt it or silently share'
 standard.main(R)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:
  cfg=standard.load(R/'config/spec_v1.json')
  with open(cfg['GPU_lock_path'],'a+') as lock:
   fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);main(R)
 except BaseException:
  import traceback
  standard.save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
