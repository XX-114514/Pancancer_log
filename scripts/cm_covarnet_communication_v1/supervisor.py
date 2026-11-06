"""Resumable bounded producer; lock prevents duplicate job submission."""
from pathlib import Path
import os,sys,subprocess,json,time,datetime,fcntl,hashlib
r=Path(__file__).resolve().parents[1];lock=(r/'logs/supervisor.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert json.loads((r/'audit/pilot_acceptance.json').read_text())['status']=='technical_pass'
cpdb=os.environ['CPDB_PYTHON'];scanpy=os.environ['SCANPY_PYTHON'];env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
(r/'logs/supervisor_pid.json').write_text(json.dumps({'pid':os.getpid(),'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
def stage(name,command):
 receipt=r/'logs'/(name+'_terminal.json')
 if receipt.exists() and json.loads(receipt.read_text())['exit_code']==0:return
 start=datetime.datetime.now(datetime.timezone.utc).isoformat()
 with (r/'logs'/(name+'.log')).open('a') as f:
  try:code=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=86400).returncode
  except subprocess.TimeoutExpired:code=124
 receipt.write_text(json.dumps({'start_UTC':start,'end_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'exit_code':code,'command':command},indent=2))
 print(name,code,flush=True)
stage('global_formal',[cpdb,str(r/'scripts/run_cpdb.py'),'global'])
stage('summary_global',[scanpy,str(r/'scripts/summarize_cpdb.py')])
deadline=time.time()+3600
while not (r/'audit/full_input_complete.json').exists() and time.time()<deadline:time.sleep(10)
if (r/'audit/full_input_complete.json').exists():
 stage('cytokine',[scanpy,str(r/'scripts/cytokine_response.py')])
 stage('comparisons',[cpdb,str(r/'scripts/run_comparisons.py')])
 stage('summary_final',[scanpy,str(r/'scripts/summarize_cpdb.py')])
else:print('full_input_not_completed_within_wait_window',flush=True)
stage('figures',[scanpy,str(r/'scripts/render_figures.py')])
