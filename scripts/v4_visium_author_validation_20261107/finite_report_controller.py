"""Finite, private-output report rendering after the real cohort audit succeeds."""
import os, sys, time, json, subprocess
from pathlib import Path
from datetime import datetime, timezone

r = Path(os.environ['RUN_ROOT'])
start = time.monotonic()
while not (r / 'audit/FINAL_COHORT_INDEPENDENT_AUDIT.json').exists():
    if time.monotonic() - start > 432000:
        raise TimeoutError('Five-day bounded final audit wait expired')
    terminal = r / 'launch/audit_Visium_cohort.py_terminal.json'
    if terminal.exists() and json.loads(terminal.read_text())['exit_code']:
        raise RuntimeError('Cohort audit failed; report must not render')
    time.sleep(30)
begin = datetime.now(timezone.utc).isoformat()
code = "import os,runpy;from pathlib import Path;runpy.run_path(str(Path(os.environ['RUN_ROOT'])/'scripts/render_Visium_report.py'),init_globals={'PROJECT_ROOT':Path(os.environ['PROJECT_ROOT'])})"
with (r / 'launch/render_Visium_report.log').open('xb') as log:
    p = subprocess.run([sys.executable, '-B', '-c', code], stdout=log, stderr=subprocess.STDOUT, timeout=7200)
(r / 'launch/render_Visium_report_terminal.json').open('x').write(json.dumps({
    'UTC_start': begin, 'UTC_end': datetime.now(timezone.utc).isoformat(),
    'exit_code': p.returncode, 'automatic_publication': False,
}, indent=2))
sys.exit(p.returncode)
