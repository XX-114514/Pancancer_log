from pathlib import Path
from datetime import datetime,timezone
import os,json,subprocess,hashlib
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);C=Path(os.environ['PANCANCER_STATE_CODE']);manifest=R/'cell_state_programme_readiness_20261013_v1.execution.json';log=R/'cell_state_programme_readiness_20261013_v1.stdout.log'
assert not manifest.exists() and not log.exists()
def now():return datetime.now(timezone.utc).isoformat()
def write(j):
 t=Path(str(manifest)+'.atomic_tmp');t.write_text(json.dumps(j,indent=2)+'\n');os.replace(t,manifest)
j=dict(status='RUNNING',started_UTC=now(),worker_pid=os.getpid(),role='state_programme_reference_readiness_not_annotation',worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),child_script_sha256=hashlib.sha256((C/'state_programme_readiness_v1.py').read_bytes()).hexdigest(),scientific_release=False,state_labels_assigned=False)
write(j);code=1
try:
 env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
 with log.open('xb') as output:
  child=subprocess.Popen([os.sys.executable,str(C/'state_programme_readiness_v1.py')],env=env,stdout=output,stderr=subprocess.STDOUT);j['child_pid']=child.pid;write(j);code=child.wait()
except Exception as e:j['error_type']=type(e).__name__
finally:
 j.update(status='TERMINAL_SUCCESS_READINESS_ONLY' if code==0 else 'TERMINAL_FAILED_READINESS',finished_UTC=now(),exit_code=code);write(j)
print(json.dumps(j),flush=True)
raise SystemExit(code)
