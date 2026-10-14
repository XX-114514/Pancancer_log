from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,subprocess
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);C=Path(os.environ['PANCANCER_STATE_INPUT_CODE']);M=R/'state_singlecell_inputs_20261014_v1.execution.json';L=R/'state_singlecell_inputs_20261014_v1.stdout.log'
assert not M.exists() and not L.exists()
def now():return datetime.now(timezone.utc).isoformat()
def write(x):
 t=Path(str(M)+'.atomic_tmp');t.write_text(json.dumps(x,indent=2)+'\n');os.replace(t,M)
x=dict(status='RUNNING_INPUT_ASSEMBLY',started_UTC=now(),worker_pid=os.getpid(),worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),producer_sha256=hashlib.sha256((C/'assemble_state_singlecell_inputs_v1.py').read_bytes()).hexdigest(),state_model_trained=False,state_labels_assigned=False)
write(x);code=1
try:
 env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
 with L.open('xb') as f:
  p=subprocess.Popen([os.sys.executable,str(C/'assemble_state_singlecell_inputs_v1.py')],env=env,stdout=f,stderr=subprocess.STDOUT);x['child_pid']=p.pid;write(x);code=p.wait()
except Exception as e:x['error_type']=type(e).__name__
finally:
 x.update(status='TERMINAL_SUCCESS_INPUTS_PENDING_READBACK_AUDIT' if code==0 else 'TERMINAL_FAILED_INPUT_ASSEMBLY',finished_UTC=now(),exit_code=code);write(x)
print(json.dumps(x),flush=True)
raise SystemExit(code)
