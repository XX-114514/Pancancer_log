import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from common_measured_feature_contract_v1 import restrict_candidates


class MeasurementContractTests(unittest.TestCase):
    def binding(self, directory, genes):
        path = Path(directory) / 'metadata.json'
        path.write_text(json.dumps({'common_nonheldout_feature_candidates': genes,
            'no_gene_absence_filled_as_zero': True, 'counts_read': False, 'feature_selection_fitted': False}))
        return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}

    def test_missing_axis_is_removed_not_zero_created(self):
        with tempfile.TemporaryDirectory() as directory:
            binding = self.binding(directory, ['A', 'C', 'EXTERNAL_ONLY'])
            self.assertEqual(restrict_candidates(['C', 'B', 'A'], ['STATE'], binding, 2), ['C', 'A'])

    def test_measurement_contract_mutation_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            binding = self.binding(directory, ['A'])
            Path(binding['path']).write_text('{}')
            with self.assertRaises(ValueError):
                restrict_candidates(['A'], [], binding, 1)

    def test_holdout_leak_and_insufficient_axis_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            binding = self.binding(directory, ['A', 'STATE'])
            with self.assertRaises(ValueError):
                restrict_candidates(['A'], ['STATE'], binding, 1)
            binding = self.binding(directory, ['A'])
            with self.assertRaises(ValueError):
                restrict_candidates(['A', 'B'], [], binding, 2)


if __name__ == '__main__':
    unittest.main()
