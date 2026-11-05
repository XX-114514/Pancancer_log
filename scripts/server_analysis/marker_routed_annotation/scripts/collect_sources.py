"""Validate source outputs and export TME-only reference candidates."""
import argparse
import gzip
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    args = p.parse_args()
    plan = json.loads(args.plan.read_text())
    run = args.plan.parent
    reports, failures, outputs = [], [], []
    columns = ['source']
    for source in plan['sources']:
        s = source['source']; out = Path(source['output'])
        try:
            summary = json.loads((out / 'summary.json').read_text())
            if summary['pilot'] or summary['source'] != s:
                raise ValueError('Pilot or wrong source')
            if summary['source_cells'] != source.get('expected_cells', summary['source_cells']):
                raise ValueError('Source cell count differs from the input inventory')
            if summary['output_cells'] != summary['source_cells']:
                raise ValueError('Full-source prediction output is truncated')
            check = subprocess.run([
                sys.executable, str(run / 'code/check_output.py'), str(out),
                '--malignancy', str(Path(plan['malignancy_dir']) / (s + '.tsv.gz'))],
                capture_output=True, text=True)
            if check.returncode:
                raise ValueError('Output readback failed: ' + check.stderr[-500:])
            # A changed malignancy sidecar invalidates cached reference candidates.
            import hashlib
            current_sha = hashlib.sha256((Path(plan['malignancy_dir']) / (s + '.tsv.gz')).read_bytes()).hexdigest()
            if current_sha != summary['malignancy_sha256']:
                raise ValueError('Malignancy sidecar changed')
            path = out / 'candidate_labels.tsv.gz'
            cols = pd.read_csv(path, sep='\t', nrows=0).columns.tolist()
            columns.extend(c for c in cols if c not in columns)
            outputs.append((s, path))
            counts = summary['route_counts']
            reports.append({'source': s, 'source_cells': summary['source_cells'],
                'cells_tme': summary['output_cells'] - counts.get('excluded_malignant', 0),
                **counts})
        except Exception as e:
            failures.append({'source': s, 'error': str(e)[:700]})
    pd.DataFrame(reports).to_csv(run / 'source_annotation_summary.tsv', sep='\t', index=False)
    # Publish the final filename only after every expected source passes.
    complete = len(outputs) == len(plan['sources']) and not failures
    name = 'tme_candidate_labels.tsv.gz' if complete else 'tme_candidate_labels.INCOMPLETE.tsv.gz'
    total = 0
    with (run / name).open('xb') as raw:
        with gzip.open(raw, 'wt') as handle:
            header = True
            for s, path in outputs:
                for block in pd.read_csv(path, sep='\t', chunksize=50000):
                    block = block.loc[block.malignancy_call != 'malignant'].copy()
                    block.insert(0, 'source', s)
                    total += len(block)
                    block.reindex(columns=columns).to_csv(handle, sep='\t', index=False, header=header)
                    header = False
    summary = {'expected_sources': len(plan['sources']), 'validated_sources': len(outputs),
               'complete': complete, 'tme_cells': total, 'failures': failures,
               'candidate_labels': name, 'final_identity_accepted': False}
    (run / 'collection_summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)
    return 0 if complete else 2


if __name__ == '__main__':
    raise SystemExit(main())
