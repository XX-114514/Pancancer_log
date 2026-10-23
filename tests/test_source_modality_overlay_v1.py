from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_source_modality_overlay_v1 import modality_decision


class SourceModality(unittest.TestCase):
    def test_contradictory_evidence_not_resolved_by_majority(self):
        self.assertEqual(modality_decision({'cell': 20}, '', ['nucleus']), 'conflicting_modality_evidence_requires_review')
        self.assertEqual(modality_decision({'cell': 20}, 'snRNA_explicit_existing_note', []), 'conflicting_modality_evidence_requires_review')

    def test_mixed_or_partially_unknown_values_do_not_become_homogeneous(self):
        self.assertEqual(modality_decision({'cell': 100, 'nucleus': 1}, '', []), 'mixed_cell_and_nucleus_requires_per_record_branch')
        self.assertEqual(modality_decision({'nucleus': 100, 'unknown': 1}, '', ['nucleus']), 'snRNA_source_evidence_native_suspension_incomplete')

    def test_platform_name_is_not_suspension_evidence(self):
        self.assertEqual(modality_decision({}, '10x_3prime_v3', []), 'modality_unresolved_not_inferred_from_platform_name')
        self.assertEqual(modality_decision({'nucleus': 100}, '', ['nucleus']), 'snRNA_native_and_public_dataset_concordant')


if __name__ == '__main__':
    unittest.main()
