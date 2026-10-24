
from pathlib import Path
import json,sys,subprocess,time,datetime,hashlib
root=Path(sys.argv[sys.argv.index('--output')+1]);cfg=json.loads((root/'config/audit_followthrough_spec_v1.json').read_text());pending=list(cfg['stages']);launched=[];began=time.monotonic()
def save(path,value):
 with path.open('x') as f:json.dump(value,f,indent=2);f.write('\n')
def live(stage):
 child=json.loads((stage/'provenance/child_launch_v1.json').read_text());stat=Path('/proc')/str(child['PID'])/'stat'
 return stat.exists() and stat.read_text().rsplit(')',1)[1].split()[19]==child['startticks']
while pending:
 for spec in list(pending):
  stage=Path(spec['producer']);terminal=stage/'provenance/supervisor_terminal_v1.json'
  if terminal.exists():
   result=json.loads(terminal.read_text())
   if result['actual_child_exit_code'] or result['wall_timeout']:raise RuntimeError('Full audit predecessor failed; science history preserved')
   launcher=root/spec['launcher']
   if hashlib.sha256(launcher.read_bytes()).hexdigest()!=spec['launcher_sha256']:raise ValueError('Frozen audit launcher changed')
   if Path(spec['audit']).exists():raise ValueError('Audit stage already exists; refusing duplicate launch')
   proc=subprocess.run([sys.executable,'-B',str(launcher)],capture_output=True,text=True,timeout=90)
   (root/'provenance'/(spec['name']+'_launcher_stdout.txt')).write_text(proc.stdout)
   (root/'provenance'/(spec['name']+'_launcher_stderr.txt')).write_text(proc.stderr)
   if proc.returncode:raise RuntimeError('Prepared full audit launcher failed')
   receipt=json.loads(proc.stdout);launched.append({'name':spec['name'],'receipt':receipt,'producer_terminal_sha256':hashlib.sha256(terminal.read_bytes()).hexdigest()})
   save(root/'checkpoints'/(spec['name']+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'name':spec['name'],'status':'FULL_AUDITOR_LAUNCHED_RESULT_PENDING'})
   pending.remove(spec)
  elif not live(stage):
   time.sleep(2)
   if not terminal.exists():raise RuntimeError('Exact full audit predecessor handle unavailable')
 if pending:
  if time.monotonic()-began>cfg['maximum_wait_seconds']:raise TimeoutError('Bounded full audit followthrough wait exhausted')
  time.sleep(15)
save(root/'outputs/audit_launch_receipts_v1.json',launched)
save(root/'provenance/terminal_v1.json',{'status':'FULL_AUDIT_CHILDREN_LAUNCHED_NOT_ACCEPTANCE','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_launched':len(launched),'audit_results_accepted':False,'new_final_identity':False})
