"""Synthetic readback: complete export versus explicit incomplete export."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import pandas as pd


class CollectionTests(unittest.TestCase):
    def run_case(self, missing=False, changed_sidecar=False):
        with tempfile.TemporaryDirectory(prefix='tme-collect-test-') as tmp:
            root = Path(tmp); (root / 'code').mkdir(); (root / 'malig').mkdir()
            for name in ['check_output.py', 'marker_policy.py']:
                shutil.copy2(Path(__file__).with_name(name), root / 'code' / name)
            sources = []
            for source in ['S1', 'S2']:
                out = root / source; out.mkdir()
                sources.append({'source': source, 'output': str(out), 'expected_cells': 3})
                if missing and source == 'S2':
                    continue
                m = pd.DataFrame({'cell_id': [source + '-i', source + '-f', source + '-m'],
                                  'malignancy_call': ['non_malignant', 'uncertain', 'malignant']})
                mp = root / 'malig' / (source + '.tsv.gz'); m.to_csv(mp, sep='\t', index=False)
                d = m.copy()
                d['marker_parent'] = ['immune', 'fibroblast', 'epithelial']
                d['prediction_route'] = ['immune_celltypist', 'singler_hpca', 'excluded_malignant']
                d['reference_status'] = ['candidate_immune_parent_supported',
                                          'candidate_hpca_parent_supported', 'excluded_malignant']
                for model in ['Immune_All_High', 'Immune_All_Low']:
                    d[model + '_candidate'] = ['T cells', None, None]
                    d[model + '_max_probability'] = [0.9, None, None]
                    d[model + '_margin'] = [0.8, None, None]
                d['hpca_main_raw'] = [None, 'Fibroblasts', None]
                d.to_csv(out / 'candidate_labels.tsv.gz', sep='\t', index=False)
                summary = {'source': source, 'pilot': False, 'source_cells': 3,
                           'output_cells': 3, 'route_counts': d.prediction_route.value_counts().to_dict(),
                           'malignancy_sha256': hashlib.sha256(mp.read_bytes()).hexdigest()}
                if changed_sidecar and source == 'S2':
                    summary['malignancy_sha256'] = 'invalidated'
                (out / 'summary.json').write_text(json.dumps(summary))
            plan = {'malignancy_dir': str(root / 'malig'), 'sources': sources}
            (root / 'plan.json').write_text(json.dumps(plan))
            result = subprocess.run([sys.executable, '-I', str(Path(__file__).with_name('collect_sources.py')),
                                     '--plan', str(root / 'plan.json')], capture_output=True, text=True)
            incomplete = missing or changed_sidecar
            self.assertEqual(result.returncode, 2 if incomplete else 0, result.stderr)
            summary = json.loads((root / 'collection_summary.json').read_text())
            self.assertEqual(summary['complete'], not incomplete)
            self.assertEqual(summary['tme_cells'], 2 if incomplete else 4)
            self.assertEqual((root / 'tme_candidate_labels.tsv.gz').exists(), not incomplete)
            table = pd.read_csv(root / summary['candidate_labels'], sep='\t')
            self.assertFalse((table.malignancy_call == 'malignant').any())
            self.assertEqual(len(table), summary['tme_cells'])

    def test_complete_export_contains_only_tme(self):
        self.run_case()

    def test_missing_source_is_incomplete_and_nonzero(self):
        self.run_case(missing=True)

    def test_changed_malignancy_sidecar_invalidates_source(self):
        self.run_case(changed_sidecar=True)


if __name__ == '__main__':
    unittest.main()
