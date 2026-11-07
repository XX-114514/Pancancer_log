"""Finite recovery of audit typing and unexecuted spatial models; reuse all prepared cells."""
from pathlib import Path
from datetime import datetime,timezone
import sys,json,subprocess,time
r=Path(sys.argv[1]);rb=sys.argv[2];project=r.parents[1]
def now():return datetime.now(timezone.utc).isoformat()
def run(tag,command):
 z={'stage':tag,'UTC_start':now()}
 with (r/'launch'/(tag+'.log')).open('xb') as f:p=subprocess.run(command,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT)
 z.update(UTC_end=now(),exit_code=p.returncode);(r/'launch'/(tag+'_terminal.json')).open('x').write(json.dumps(z));return p.returncode

def py(name,args=()):return [sys.executable,'-B','-c',"import runpy,sys;from pathlib import Path;sys.path.insert(0,"+repr(str(r/'scripts'))+");sys.argv=["+repr(name)+"]+"+repr(list(args))+";runpy.run_path("+repr(str(r/'scripts'/name))+",init_globals={'PROJECT_ROOT':Path("+repr(str(project))+")})"]
for _ in range(2880):
 if (r/'launch/bulk_models_audit_v1_terminal.json').exists():break
 time.sleep(15)
results={'bulk_measurement_audit_v2':run('bulk_measurement_audit_v2',py('audit_V4_bulk_measurement_v2.py'))}
failed=r/'launch/imaging_I001_broad_context_terminal.json'
if failed.exists() and json.loads(failed.read_text())['exit_code']!=0 and not (r/'imaging/I001/context_complete.json').exists():results['I001_int64_context_recovery']=run('imaging_I001_int64_context_recovery_v2',py('measure_imaging_broad_CM_context.py',['0']))
results['spatial_recovery']=[]
for job in json.loads((r/'config/imaging_jobs_PRIVATE.json').read_text()):
 unit=job['unit_code'];ready=False
 for _ in range(2880):
  term=r/'launch'/('imaging_'+unit+'_CellChat_terminal.json');prep=r/'launch'/('imaging_'+unit+'_prepare_terminal.json')
  if term.exists():ready=True;break
  if prep.exists() and json.loads(prep.read_text())['exit_code']!=0:break
  time.sleep(15)
 if ready and (r/'imaging'/unit/'CellChat_v1').exists() and not (r/'imaging'/unit/'CellChat_v2').exists():
  z={'unit_code':unit,'fit':run('imaging_'+unit+'_geometry_interface_recovery_v2',[rb,'--vanilla',str(r/'scripts/run_spatial_CellChat_V4_v2.R'),str(r),unit])}
  if z['fit']==0:z['audit']=run('imaging_'+unit+'_geometry_interface_audit_v2',[rb,'--vanilla',str(r/'scripts/audit_spatial_CellChat_V4_v2.R'),str(r),unit])
  results['spatial_recovery'].append(z)
 (r/'launch/interface_recovery_progress_PRIVATE.json').write_text(json.dumps({'UTC':now(),'results':results},indent=2))
(r/'launch/interface_recovery_all_terminal_PRIVATE.json').open('x').write(json.dumps({'UTC':now(),'results':results},indent=2))
