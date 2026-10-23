"""Reference applicability guards; expression absence is not identity absence."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import numpy as np
from full_archive_reference_evidence_v1 import reference_scopes, pack_support
from native_annotation_review_v1 import LINEAGE_PANELS, TYPE_PANELS


class ScopeTests(unittest.TestCase):
    def calls(self, n=4):
        return {key: np.zeros(n, bool) for key in {**LINEAGE_PANELS, **TYPE_PANELS}}

    def test_nonimmune_cells_not_forced_into_immune_reference(self):
        calls = self.calls(); calls['Epithelial'][0] = True; calls['T'][1] = True
        immune, organ = reference_scopes(calls, ['Epithelial', 'T_NK', 'unresolved_context', 'unresolved_context'], np.ones(4, bool), True)
        np.testing.assert_array_equal(immune, [False, True, False, False])
        np.testing.assert_array_equal(organ, [True, False, False, False])

    def test_zero_RNA_and_unavailable_organ_remain_out_of_scope(self):
        calls = self.calls(); calls['Epithelial'][0] = True; calls['T'][1] = True
        immune, organ = reference_scopes(calls, ['Epithelial', 'T_NK', 'T_NK', 'unresolved_context'], [True, False, True, True], False)
        np.testing.assert_array_equal(immune, [False, False, True, False]); self.assertFalse(organ.any())

    def test_mixed_lineage_evidence_retained_in_bits(self):
        bits = pack_support({'a': [True, True, False], 'b': [False, True, False]}, ['a', 'b'])
        np.testing.assert_array_equal(bits, [1, 3, 0])


if __name__ == '__main__':
    unittest.main()
