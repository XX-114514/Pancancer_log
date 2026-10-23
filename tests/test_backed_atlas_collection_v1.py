"""Small numerical fixtures for source selection and missing-measurement masks."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

try:
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import zarr
    from anndata._io.specs import write_elem
    DEPS = True
except ImportError:
    DEPS = False


@unittest.skipUnless(DEPS, 'Requires source scientific runtime')
class BackedCollectionTests(unittest.TestCase):
    def setUp(self):
        spec = importlib.util.spec_from_file_location('collection_fixture_reader', Path(__file__).parents[1] / 'scripts/backed_atlas_collection_v1.py')
        self.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.module)
        self.tmp = tempfile.TemporaryDirectory(prefix='atlas-collection-fixture-')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        path = self.root / 'source.zarr'
        z = zarr.open_group(str(path), mode='w')
        self.dense = np.asarray([[1, 2, 8, 0], [3, 0, 5, 1], [0, 4, 1, 2], [7, 1, 0, 1]], dtype=np.int32)
        write_elem(z, 'layers', {'counts': sp.csr_matrix(self.dense), 'ambient_removed_counts': sp.csr_matrix(self.dense * .1)})
        obs = pd.DataFrame({'dataset': ['S', 'S', 'S', 'OTHER'], 'major': ['T', 'Fibroblast', 'Myeloid', 'T']}, index=['r0', 'r1', 'r2', 'r3'])
        var = pd.DataFrame({'canonical_gene_identity': ['ENSG00000000001', 'ENSG00000000002', 'ENSG00000000003', 'AMBIGUOUS::S::D'], 'gene_symbol': ['A', 'B', 'C', 'D'], 'feature_identity_status': ['stable_Ensembl_ID'] * 3 + ['ambiguous_source_feature'], 'source_gene_aliases': ['A|ONE', 'B|TWO', 'C', 'D']}, index=['f0', 'f1', 'f2', 'f3'])
        # Match released stores where the canonical ID is the DataFrame index.
        var = var.set_index('canonical_gene_identity')
        write_elem(z, 'obs', obs); write_elem(z, 'var', var)
        write_elem(z, 'uns', {'source_order': np.asarray(['S'])})
        write_elem(z, 'varm', {'measured_in_source': np.asarray([[1], [1], [0], [1]], bool)})
        self.manifest = {'schema': 'backed_atlas_collection_v1', 'partitions': [{'dataset': 'S', 'path': str(path), 'counts_key': 'layers/counts', 'ambient_removed_key': 'layers/ambient_removed_counts', 'n_context_rows': 3, 'n_immune_rows': 2, 'selector': {'kind': 'field_equals', 'field': 'dataset', 'value': 'S'}, 'immune_selector': {'field': 'major', 'values': ['T', 'Myeloid']}}]}
        self.file = self.root / 'manifest.json'; self.file.write_text(json.dumps(self.manifest))
        self.collection = self.module.BackedAtlasCollection(self.file)

    def test_order_duplicates_missing_and_ambiguous_features(self):
        result = self.collection.read_panel('S', [1, 0, 1], ['ONE', 'B', 'C', 'D', 'missing'])
        np.testing.assert_array_equal(result['counts'].toarray(), [[0, 4, 0, 0, 0], [1, 2, 0, 0, 0], [0, 4, 0, 0, 0]])
        np.testing.assert_array_equal(result['measured'], [True, True, False, False, False])
        np.testing.assert_array_equal(result['source_row_positions'], [2, 0, 2])

    def test_applied_correction_and_context_view(self):
        result = self.collection.read_panel('S', [1, 0], ['A', 'B'], view='context', layer='analysis_counts')
        np.testing.assert_allclose(result['counts'].toarray(), self.dense[[1, 0], :2] * .9)

    def test_unknown_view_out_of_range_and_alias_collision_rejected(self):
        with self.assertRaises(ValueError): self.collection.rows('S', 'unreviewed')
        with self.assertRaises(IndexError): self.collection.read_panel('S', [2], ['A'])
        with self.assertRaises(ValueError): self.collection.read_panel('S', [0], ['A', 'ONE'])

    def test_duplicate_source_and_changed_partition_rejected(self):
        self.manifest['partitions'].append(dict(self.manifest['partitions'][0]))
        self.file.write_text(json.dumps(self.manifest))
        with self.assertRaises(ValueError): self.module.BackedAtlasCollection(self.file)
        self.manifest['partitions'].pop(); self.manifest['partitions'][0]['n_context_rows'] = 99
        self.file.write_text(json.dumps(self.manifest))
        with self.assertRaises(ValueError): self.module.BackedAtlasCollection(self.file).rows('S')


if __name__ == '__main__':
    unittest.main()
