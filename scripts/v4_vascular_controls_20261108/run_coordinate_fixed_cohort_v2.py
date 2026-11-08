"""Finite server execution after verified uploads; all large outputs stay there."""
from pathlib import Path
import json, os, sys, time, traceback, subprocess
from datetime import datetime, timezone
from run_native_vascular_controls_v6 import run
from vascular_spatial_controls_v2 import patient_summary

root = Path(sys.argv[1]).resolve()
manifest = json.loads((root / 'config/cohort_PRIVATE.json').read_text())
cases = manifest['cases']
assert len(cases) == 6
assert len({c['unit'] for c in cases}) == 6
assert len({c['donor'] for c in cases}) == 6
out = root / 'outputs_direct_v1'
out.mkdir(exist_ok=False)
events = out / 'events_PRIVATE.jsonl'
deadline = time.monotonic() + 18 * 3600
done = {}


def event(**data):
    with events.open('a') as handle:
        handle.write(json.dumps({'UTC': datetime.now(timezone.utc).isoformat(), **data}) + '\n')


event(status='COORDINATE_FIXED_SERVER_STAGE_STARTED_WAITING_VERIFIED_INPUTS', requested=6)
while len(done) < len(cases) and time.monotonic() < deadline:
    for case in cases:
        unit = case['unit']
        if unit in done or not (root / 'inputs_direct_v1' / unit / 'ready_PRIVATE.json').exists():
            continue
        try:
            event(unit=unit, status='START_NATIVE_REPLAY_WITH_EXPLICIT_COORDINATE_UNITS')
            result = run(out, unit, root / 'inputs_direct_v1' / unit, root / 'inputs_direct_v1/mapping',
                         expected_member=case['expected_member'],
                         panel_gene_names=set(manifest['panel_genes']))
            event(unit=unit, status='COMPUTED_NOT_BIOLOGICALLY_OR_VISUALLY_ACCEPTED',
                  coordinate_contract=result['coordinate_contract'])
            done[unit] = 'complete'
            subprocess.run([sys.executable, str(root / 'scripts/render_vascular_identity_review_v3.py'),
                            str(out / 'spatial_controls_v1' / unit)], check=True,
                           env={**os.environ, 'MPLCONFIGDIR': str(root / 'config/mpl_cache')})
        except Exception as error:
            done[unit] = 'failed'
            (out / (unit + '_failure_PRIVATE.txt')).write_text(traceback.format_exc())
            event(unit=unit, status='FAILED_RETAINED', exception_type=type(error).__name__)
    if len(done) < len(cases):
        time.sleep(30)
accepted = [c for c in cases if done.get(c['unit']) == 'complete']
if accepted:
    patient_summary([out / 'spatial_controls_v1' / c['unit'] for c in accepted],
                    [c['donor'] for c in accepted], out / 'cohort_summary_PRIVATE')
(out / 'terminal_PRIVATE.json').write_text(json.dumps({
    'UTC': datetime.now(timezone.utc).isoformat(), 'requested': 6,
    'computed': len(accepted), 'failed_or_pending': 6 - len(accepted),
    'biological_and_visual_review_required': True,
    'status': 'FINITE_REPLAY_ENDED_CHECK_INDIVIDUAL_RECEIPTS',
}, indent=2) + '\n')
