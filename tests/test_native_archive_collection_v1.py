"""Archive/context/roster/QC boundaries and immutable source namespaces."""
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
try:
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from freeze_native_archive_manifest_v1 import create_manifest
    from native_archive_collection_v1 import NativeArchiveCollection
    from test_full_context_atlas_v1 import FullContextTests
    AVAILABLE = True
except ImportError:
    AVAILABLE = False


@unittest.skipUnless(AVAILABLE, 'Run counterexamples in the scientific environment')
class NativeArchiveTests(unittest.TestCase):
    def setUp(self):
        base = FullContextTests()
        base.setUp()
        self.root, self.context = base.root, base.manifest
        old = ad.AnnData(sp.csr_matrix([[2, 3, 0], [0, 1, 0], [7, 0, 0], [0, 0, 0], [8, 2, 0]], dtype=np.int32),
            obs=pd.DataFrame({'gse_id': ['OLD1'] * 4 + ['EXCLUDED'],
                              'sampleID': ['raw0', 'raw1', 'raw0', 'raw2', 'raw3'],
                              'sample_uid': ['unit0', 'unit1', 'unit0', 'unit2', 'unit3'],
                              'majorCluster': ['T', 'T', 'Myeloid', 'unresolved', 'T'],
                              'ev_author_label': ['author_T', 'author_T', 'author_M', 'unknown', 'author_T'],
                              'expression_v1_exclusion_reason': ['', '', '', 'identity_unresolved', 'source_excluded']},
                             index=['o0', 'o1', 'o2', 'o3', 'o4']),
            var=pd.DataFrame(index=['G1', 'G2', 'G3']))
        old.write_zarr(self.root / 'all_old.zarr')
        (self.root / 'all_native').mkdir()
        old[:4, :2].copy().write_h5ad(self.root / 'all_native/OLD1_symbol_counts.h5ad')
        self.context['legacy'].update(canonical_path=str(self.root / 'all_old.zarr'),
                                      native_source_root=str(self.root / 'all_native'))
        new = ad.read_zarr(self.root / 'new.zarr')
        new.obs['annotation_exclude_from_primary'] = [False, False, True]
        new.write_zarr(self.root / 'flagged_new.zarr')
        primary = json.loads((self.root / 'primary.json').read_text())
        for p in primary['partitions']:p['path'] = str(self.root / 'flagged_new.zarr')
        (self.root / 'flagged_primary.json').write_text(json.dumps(primary))
        self.context['primary_manifest'] = str(self.root / 'flagged_primary.json')
        self.frozen_path = self.root / 'analysis_frozen.json'
        self.frozen_path.write_text(json.dumps(self.context))
        self.path = self.root / 'archive.json'
        self.spec = create_manifest(self.frozen_path)
        self.save()

    def save(self):
        self.path.write_text(json.dumps(self.spec))

    def test_context_two_roster_three_archive_four_source_exclusion_retained(self):
        r = NativeArchiveCollection(self.path)
        np.testing.assert_array_equal(r.context_rows('OLD1'), [2, 0])
        np.testing.assert_array_equal(r.rows('OLD1'), [0, 1, 2, 3])
        m = r.membership('OLD1')
        np.testing.assert_array_equal(m['in_frozen_context'], [True, False, True, False])
        np.testing.assert_array_equal(m['in_legacy_roster'], [True, True, True, False])
        np.testing.assert_array_equal(m['below_legacy_expression_unit_gate'], [False, True, False, False])
        self.assertNotIn('EXCLUDED', r.partitions)
        self.assertEqual(self.spec['native_store_source_scope'][0]['outside_selected_source_scope'], {'EXCLUDED': 1})
        np.testing.assert_array_equal(r.frozen_labels('OLD1')['fine'], ['F1', 'F0', 'F0', 'unresolved_outside_frozen_roster'])
        self.assertEqual(r.author_labels('OLD1')['values'][0], 'author_T')

    def test_zero_RNA_missing_features_duplicate_read_order_and_units(self):
        r = NativeArchiveCollection(self.path)
        p = r.read_panel('OLD1', [3, 0, 3, 1], ['G1', 'G3'])
        np.testing.assert_array_equal(p['counts'].toarray(), [[0, 0], [2, 0], [0, 0], [0, 0]])
        np.testing.assert_array_equal(p['measured'], [True, False])
        np.testing.assert_array_equal(p['native_total_counts'], [0, 5, 0, 1])
        np.testing.assert_array_equal(r.frozen_expression_units('OLD1'), ['unit0', 'unit1', 'unit0', 'unit2'])
        self.assertNotEqual(r.field('OLD1', 'sampleID')[1], r.frozen_expression_units('OLD1')[1])

    def test_primary_QC_exclusion_is_retained_in_archive_not_model_approval(self):
        r = NativeArchiveCollection(self.path)
        self.assertEqual(len(r.rows('NEW1')), 2)
        np.testing.assert_array_equal(r.field('NEW1', 'annotation_exclude_from_primary'), [False, True])
        self.assertFalse(r.membership('NEW1')['legacy_roster_applicable'].any())
        self.assertEqual(r.membership('NEW1')['new_model_eligibility'], 'unassessed_not_inferred_from_membership')

    def test_native_id_unit_or_source_scope_mutation_rejected(self):
        for key in ['native_rows_sha256', 'native_ids_sha256', 'native_units_sha256']:
            original = self.spec['partitions'][0][key]
            self.spec['partitions'][0][key] = '0' * 64;self.save()
            with self.assertRaises(ValueError):
                r = NativeArchiveCollection(self.path)
                r.rows('OLD1');r.native_ids('OLD1');r.frozen_expression_units('OLD1')
            self.spec['partitions'][0][key] = original
        self.spec['partitions'].pop();self.save()
        with self.assertRaises(ValueError):NativeArchiveCollection(self.path)

    def test_unit_namespace_mismatch_rejected_even_if_contract_rebound(self):
        # A valid digest cannot substitute for agreement with the old exact-cell roster.
        cfg = self.context['legacy']
        with np.load(cfg['roster_path'], allow_pickle=False) as f:
            arrays = {k: f[k].copy() for k in f.files}
        arrays['sample_categories'] = np.array(['unit0', 'WRONG'])
        np.savez(self.root / 'wrong_units.npz', **arrays)
        (self.root / 'wrong_samples.tsv').write_text('sampleID\tgse_id\teligible_min100_cells\nunit0\tOLD1\tTrue\nWRONG\tOLD1\tFalse\n')
        cfg.update(roster_path=str(self.root / 'wrong_units.npz'), sample_metadata_path=str(self.root / 'wrong_samples.tsv'))
        self.frozen_path.write_text(json.dumps(self.context))
        self.spec = create_manifest(self.frozen_path);self.save()
        with self.assertRaises(ValueError):NativeArchiveCollection(self.path).frozen_expression_units('OLD1')


if __name__ == '__main__':
    unittest.main()
