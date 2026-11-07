"""Finite dependency continuation; one full-cohort extension, no publication."""
import sys,time,subprocess,json
from pathlib import Path
from datetime import datetime,timezone
r=Path(sys.argv[1]);script=Path(sys.argv[2]);out=Path(sys.argv[3]);deadline=time.monotonic()+72*3600
while time.monotonic()<deadline:
 if len(list((r/'sections').glob('S*/independent_audit.json')))==53 or (r/'launch/all_scores_terminal.json').exists():
  code=subprocess.call([sys.executable,str(script),str(r),str(out)])
  out.parent.joinpath('keymodule_full_cohort_terminal.json').open('x').write(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'exit_code':code,'scope':'all accepted sections after original full scoring termination; failed sections retained in denominator53'}));sys.exit(code)
 time.sleep(120)
raise TimeoutError('72h finite dependency deadline; original jobs unchanged')
