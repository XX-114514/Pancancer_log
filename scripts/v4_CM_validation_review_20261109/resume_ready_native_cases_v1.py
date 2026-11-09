"""Continue the original six-case analysis after its input-wait deadline.

Reuse completed cases and original scientific code. No new discovery or cutoff.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import os
import subprocess
import sys
import traceback

root = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(root / 'scripts'))
from run_native_vascular_controls_v6 import run
from vascular_spatial_controls_v2 import patient_summary

config = json.loads((root / 'config/cohort_PRIVATE.json').read_text())
cases = config['cases']
assert len(cases) == len({c['donor'] for c in cases}) == 6
old = root / 'outputs_direct_v1'
assert (old / 'terminal_PRIVATE.json').exists(), 'Do not overlap an active prior consumer'
assert all((root / 'inputs_direct_v1' / c['unit'] / 'ready_PRIVATE.json').exists() for c in cases)
out = root / 'outputs_resume_v1'
out.mkdir(exist_ok=False)
events = out / 'events_PRIVATE.jsonl'
directories, donors, failures = [], [], []


def event(**values):
    with events.open('a') as handle:
        handle.write(json.dumps({'UTC': datetime.now(timezone.utc).isoformat(), **values}) + '\n')


for case in cases:
    unit = case['unit']
    prior = old / 'spatial_controls_v1' / unit
    try:
        if (prior / 'complete.json').exists():
            receipt = json.loads((prior / 'complete.json').read_text())
            assert receipt['native_source_file_sha256'] == case['receipt']['members'][0]['sha256']
            assert receipt['coordinate_contract']['status'].startswith('PASS_')
            assert all((prior / n).exists() for n in ['neighborhood.tsv', 'regional_ES.tsv', 'conditional_LR_edges.tsv'])
            directory = prior
            event(stage=unit, status='REUSE_PREVIOUSLY_COMPUTED_CASE_NO_REFIT')
        else:
            event(stage=unit, status='START_ORIGINAL_PARAMETERS_AFTER_VERIFIED_INPUT_READY')
            receipt = run(out, unit, root / 'inputs_direct_v1' / unit,
                          root / 'inputs_direct_v1/mapping',
                          expected_member=case['expected_member'],
                          panel_gene_names=set(config['panel_genes']))
            directory = out / 'spatial_controls_v1' / unit
            subprocess.run([sys.executable, str(root / 'scripts/render_vascular_identity_review_v3.py'),
                            str(directory)], check=True,
                           env={**os.environ, 'MPLCONFIGDIR': str(root / 'config/mpl_cache')})
            event(stage=unit, status='COMPUTED_REQUIRES_NUMERIC_IDENTITY_AND_VISUAL_REVIEW')
        directories.append(directory)
        donors.append(case['donor'])
    except Exception as error:
        failures.append(unit)
        (out / (unit + '_failure_PRIVATE.txt')).write_text(traceback.format_exc())
        event(stage=unit, status='FAILED_RETAINED', exception_type=type(error).__name__)
if directories:
    patient_summary(directories, donors, out / 'cohort_summary_PRIVATE')
terminal = {'UTC': datetime.now(timezone.utc).isoformat(), 'requested_cases': 6,
            'computed_or_reused_cases': len(directories), 'failed_cases': failures,
            'original_parameters_and_labels_unchanged': True,
            'biological_numeric_and_visual_acceptance_pending': True}
(out / 'terminal_PRIVATE.json').write_text(json.dumps(terminal, indent=2) + '\n')
print(json.dumps(terminal), flush=True)
