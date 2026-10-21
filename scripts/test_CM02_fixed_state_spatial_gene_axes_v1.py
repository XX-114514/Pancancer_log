"""Required biological branches, missing/ambiguous targets, and axis semantics."""
import unittest
from qualify_CM02_fixed_state_spatial_gene_axes_v1 import qualify, DEFINITIONS


class Tests(unittest.TestCase):
    def test_complete_symbols_not_ensembl(self):
        names = sorted(set(sum(DEFINITIONS.values(), [])))
        self.assertTrue(all(qualify(names)['complete'].values()))
        self.assertFalse(any(qualify(['ENSG00000000001'])['complete'].values()))

    def test_OR_branch_absence_not_silently_changed(self):
        names = sorted(set(sum(DEFINITIONS.values(), [])) - {'UBE2C'})
        result = qualify(names)
        self.assertFalse(result['complete']['E_cycling'])
        self.assertTrue(result['complete']['R_LAG3'])

    def test_ambiguity_not_missing_or_zero(self):
        names = sorted(set(sum(DEFINITIONS.values(), []))) + ['lag3']
        result = qualify(names)
        self.assertFalse(result['complete']['R_LAG3'])
        self.assertEqual(result['ambiguous'], ['LAG3'])
        self.assertEqual(result['missing'], [])


if __name__ == '__main__':
    unittest.main()
