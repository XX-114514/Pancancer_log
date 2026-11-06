"""Continue the existing comparison producer without killing a multi-day family at 24h.
Resource/lifecycle correction only: same inputs, seeds, four workers and scientific settings.
"""
from pathlib import Path
import os,sys,time,json,subprocess,fcntl,hashlib,datetime
r=Path(__file__).resolve().parents[1]
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
lock=(r/'logs/supervisor.lock').open('a');(r/'logs/guardian_waiting.json').write_text(json.dumps({'pid':os.getpid(),'UTC':now()}));fcntl.flock(lock,fcntl.LOCK_EX)
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1');cpdb=env['CPDB_PYTHON'];scanpy=env['SCANPY_PYTHON'];script=str(r/'scripts/run_comparisons.py')
(r/'logs/guardian_active.json').write_text(json.dumps({'pid':os.getpid(),'UTC':now(),'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}))
def matching_producers():
 found=[]
 for p in Path('/proc').glob('[0-9]*/cmdline'):
  try:
   args=p.read_bytes().split(b'\0')
   if script.encode() in args:found.append(int(p.parent.name))
  except (OSError,ProcessLookupError):pass
 return found
while matching_producers():time.sleep(30)
# Verify existing producer receipts before using skip/resume behavior.
for p in (r/'results').glob('CM*/complete.json'):
 for name,sha in json.loads(p.read_text())['sha256'].items():assert hashlib.sha256((p.parent/name).read_bytes()).hexdigest()==sha
started=now()
with (r/'logs/comparisons_guardian.log').open('a') as f:code=subprocess.run([cpdb,script],env=env,stdout=f,stderr=subprocess.STDOUT).returncode
(r/'logs/comparisons_guardian_terminal.json').write_text(json.dumps({'start_UTC':started,'end_UTC':now(),'exit_code':code},indent=2))
if code:sys.exit(code)
for name in ['audit_comparisons.py','summarize_cpdb.py','render_figures_reviewed.py']:
 with (r/'logs'/('guardian_'+name+'.log')).open('a') as f:
  result=subprocess.run([scanpy,str(r/'scripts'/name)],env=env,stdout=f,stderr=subprocess.STDOUT)
 if result.returncode:sys.exit(result.returncode)
(r/'logs/guardian_terminal.json').write_text(json.dumps({'end_UTC':now(),'exit_code':0,'B_all18_complete_and_readback':True,'patient_inference':'not_established'},indent=2))
