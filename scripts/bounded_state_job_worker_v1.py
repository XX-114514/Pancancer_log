"""Single authorized versioned job; status is execution only, not scientific release."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,subprocess
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);job=os.environ['PANCANCER_STATE_JOB_ID'];script=Path(os.environ['PANCANCER_STATE_PRODUCER']);assert script.is_file()
manifest=R/(job+'.execution.json');log=R/(job+'.stdout.log');assert not manifest.exists() and not log.exists()
record=dict(status='RUNNING',job_id=job,started_UTC=datetime.now(timezone.utc).isoformat(),worker_pid=os.getpid(),producer=script.name)
manifest.write_text(json.dumps(record,indent=2)+'\n')
with log.open('xb') as f:
 p=subprocess.Popen([os.sys.executable,str(script)],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT);record['child_pid']=p.pid;manifest.write_text(json.dumps(record,indent=2)+'\n');code=p.wait()
record.update(status='TERMINAL_SUCCESS_PENDING_SCIENTIFIC_ACCEPTANCE' if code==0 else 'TERMINAL_FAILED',exit_code=code,finished_UTC=datetime.now(timezone.utc).isoformat());manifest.write_text(json.dumps(record,indent=2)+'\n');raise SystemExit(code)
