"""Independent dense checks for coalesced CSR and native expression evidence."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
try:
    import numpy as np
    import scipy.sparse as sp
    import zarr
    from anndata._io.specs import write_elem
    from native_panel_aggregate_v1 import gather_csr_rows, MarkerEvidence
    AVAILABLE = True
except ImportError:
    AVAILABLE = False


@unittest.skipUnless(AVAILABLE, 'Run in actual scientific environment')
class NativeEvidenceTests(unittest.TestCase):
    def test_coalesced_gather_restores_duplicates_gaps_and_order(self):
        dense = np.random.default_rng(3).integers(0, 4, (30, 7))
        dense[dense < 2] = 0
        z = zarr.group(store=zarr.MemoryStore())
        write_elem(z, 'X', sp.csr_matrix(dense))
        rows = [29, 0, 4, 4, 19, 2, 18, 1]
        for span in [1, 3, 8192]:
            np.testing.assert_array_equal(gather_csr_rows(z['X'], rows, max_span=span).toarray(), dense[rows])
        self.assertEqual(gather_csr_rows(z['X'], []).shape, (0, 7))
        for bad in [[-1], [30]]:
            with self.assertRaises(IndexError):
                gather_csr_rows(z['X'], bad)
        with self.assertRaises(ValueError):
            gather_csr_rows(z['X'], [0.5])

    def test_native_denominator_missingness_and_zero_RNA(self):
        counts = np.array([[2, 0, 0], [4, 3, 0], [0, 0, 0], [1, 0, 0]])
        totals = np.array([10, 20, 0, 4])
        codes = np.array([0, 0, 1, 1])
        e = MarkerEvidence(['T', 'unresolved'], ['A', 'B', 'MISSING'], [True, True, False])
        e.update(sp.csr_matrix(counts[:2]), totals[:2], codes[:2])
        e.update(sp.csr_matrix(counts[2:]), totals[2:], codes[2:])
        r = e.records()
        self.assertEqual(r[0]['marker_evidence']['A']['mean_native_CP10K'], 2000)
        self.assertAlmostEqual(r[0]['marker_evidence']['B']['mean_native_log1p_CP10K'], np.log1p(1500) / 2)
        self.assertEqual(r[1]['n_archive_cells'], 2)
        self.assertEqual(r[1]['n_native_RNA_evaluable'], 1)
        self.assertEqual(r[1]['marker_evidence']['A']['mean_native_CP10K'], 2500)
        self.assertIsNone(r[0]['marker_evidence']['MISSING']['n_detected'])
        self.assertEqual(r[1]['marker_evidence']['B']['n_detected'], 0)

    def test_invalid_count_mask_and_denominator_fail(self):
        for counts, totals in [([[1, 1]], [5]), ([[4, 0]], [3]), ([[-1, 0]], [4])]:
            e = MarkerEvidence(['x'], ['A', 'MISSING'], [True, False])
            with self.assertRaises(ValueError):
                e.update(sp.csr_matrix(counts), np.array(totals), np.array([0]))
        e = MarkerEvidence(['x'], ['A'], [True])
        e.update(sp.csr_matrix([[0]]), np.array([0]), np.array([0]))
        self.assertIsNone(e.records()[0]['marker_evidence']['A']['detection_fraction'])


if __name__ == '__main__':
    unittest.main()
