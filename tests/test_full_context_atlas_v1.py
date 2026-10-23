"""Meaningful cell-order, source-measurement and denominator counterexamples."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
try:
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from full_context_atlas_v1 import FullContextAtlas
    AVAILABLE = True
except ImportError:
    AVAILABLE = False


@unittest.skipUnless(AVAILABLE, 'Scientific dependencies unavailable; run in actual scientific environment')
class FullContextTests(unittest.TestCase):
    def setUp(self):
        # Task-owned fixtures outside the source workspace; retain failed inputs.
        self.root = Path(tempfile.mkdtemp(prefix='full-context-counterexample-'))
        a = ad.AnnData(sp.csr_matrix([[2, 3, 0], [0, 1, 0], [7, 0, 0]], dtype=np.int32),
                       obs=pd.DataFrame({'gse_id': ['OLD1'] * 3,
                                         'majorCluster': ['T cell', 'T cell', 'Myeloid']},
                                        index=['o0', 'o1', 'o2']),
                       var=pd.DataFrame(index=['G1', 'G2', 'G3']))
        a.write_zarr(self.root / 'old.zarr')
        native = a[:, ['G1', 'G2']].copy()
        (self.root / 'native').mkdir()
        native.write_h5ad(self.root / 'native/OLD1_symbol_counts.h5ad')
        np.savez(self.root / 'roster.npz', global_position=np.array([2, 0, 1]),
                 sample_code=np.array([0, 0, 1]), label_code=np.array([0, 1, 0]),
                 sample_categories=np.array(['unit0', 'unit1']), label_categories=np.array(['F0', 'F1']))
        (self.root / 'samples.tsv').write_text('sampleID\tgse_id\teligible_min100_cells\nunit0\tOLD1\tTrue\nunit1\tOLD1\tFalse\n')
        b = ad.AnnData(sp.csr_matrix([[11, 2, 0], [4, 1, 0], [5, 6, 0]], dtype=np.int32),
                       obs=pd.DataFrame({'dataset': ['NEW1', 'REF1', 'NEW1'], 'major': ['T', 'T', 'B'],
                                         'node_id': ['N0', 'N1', 'N2'],
                                         'author': ['Tline', 'Epithelial', 'Bline']}, index=['n0', 'n1', 'n2']),
                       var=pd.DataFrame({'gene_symbol': ['G1', 'G2', 'G3'],
                                         'canonical_gene_identity': ['E1', 'E2', 'E3'],
                                         'feature_identity_status': ['stable_Ensembl_ID'] * 3,
                                         'source_gene_aliases': ['G1|ONE', 'G2', 'G3']}))
        b.layers['counts'] = b.X.copy()
        b.uns['source_order'] = ['NEW1', 'REF1']
        b.varm['measured_in_source'] = np.array([[True, True], [True, True], [False, False]])
        b.write_zarr(self.root / 'new.zarr')
        parts = []
        for ds, n, annotation in [('NEW1', 2, {'major_field': 'major', 'subcluster_field': 'node_id'}),
                                   ('REF1', 1, {'field': 'author', 'mapping': {'Tline': 'T'}})]:
            parts.append({'dataset': ds, 'path': str(self.root / 'new.zarr'),
                          'selector': {'kind': 'field_equals', 'field': 'dataset', 'value': ds},
                          'n_context_rows': n, 'counts_key': 'layers/counts', 'annotation': annotation})
        (self.root / 'primary.json').write_text(json.dumps({'schema': 'backed_atlas_collection_v1', 'partitions': parts}))
        self.manifest = {'schema': 'full_context_atlas_collection_v1', 'bindings': [],
                         'primary_manifest': str(self.root / 'primary.json'),
                         'legacy': {'canonical_path': str(self.root / 'old.zarr'), 'counts_key': 'X',
                                    'roster_path': str(self.root / 'roster.npz'),
                                    'sample_metadata_path': str(self.root / 'samples.tsv'),
                                    'native_source_root': str(self.root / 'native')},
                         'partitions': [{'dataset': 'OLD1', 'provider': 'legacy_frozen_native_adapter', 'n_context_rows': 2}] +
                                       [{'dataset': p['dataset'], 'provider': 'backed_existing_primary_collection',
                                         'n_context_rows': p['n_context_rows']} for p in parts]}
        self.path = self.root / 'full.json'
        self.save()

    def save(self):
        self.path.write_text(json.dumps(self.manifest))

    def test_frozen_roster_order_and_excluded_sample(self):
        r = FullContextAtlas(self.path)
        np.testing.assert_array_equal(r.rows('OLD1'), [2, 0])
        np.testing.assert_array_equal(r.frozen_labels('OLD1')['fine'], ['F0', 'F1'])
        np.testing.assert_array_equal(r.rows('NEW1'), [0, 2])

    def test_duplicates_order_missingness_and_full_RNA_denominator(self):
        r = FullContextAtlas(self.path)
        p = r.read_panel('OLD1', [0, 1, 0], ['G1', 'G3', 'ABSENT'])
        np.testing.assert_array_equal(p['counts'].toarray(), [[7, 0, 0], [2, 0, 0], [7, 0, 0]])
        np.testing.assert_array_equal(p['measured'], [True, False, False])
        np.testing.assert_array_equal(p['native_total_counts'], [7, 5, 7])
        q = r.read_panel('NEW1', [1, 0, 1], ['G2', 'G3'])
        np.testing.assert_array_equal(q['counts'].toarray(), [[6, 0], [2, 0], [6, 0]])
        np.testing.assert_array_equal(q['native_total_counts'], [11, 13, 11])

    def test_unmapped_author_context_is_unresolved_not_immune(self):
        r = FullContextAtlas(self.path)
        labels = r.frozen_labels('REF1')
        self.assertEqual(labels['major'].tolist(), ['unresolved'])
        self.assertEqual(labels['fine'].tolist(), ['Epithelial'])

    def test_alias_collision_and_invalid_row_requests_fail(self):
        r = FullContextAtlas(self.path)
        with self.assertRaises(ValueError):
            r.read_panel('NEW1', [0], ['G1', 'ONE'])
        for positions in [[-1], [2]]:
            with self.assertRaises(IndexError):
                r.read_panel('OLD1', positions, ['G1'])
        for positions in [[0.5], list(range(4097))]:
            with self.assertRaises(ValueError):
                r.read_panel('OLD1', positions, ['G1'])
        self.assertEqual(r.read_panel('OLD1', [], [])['counts'].shape, (0, 0))

    def test_changed_bound_input_and_source_count_fail(self):
        p = self.root / 'samples.tsv'
        self.manifest['bindings'] = [{'path': str(p), 'role': 'frozen scope',
                                      'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}]
        self.save()
        p.write_text(p.read_text() + '\n')
        with self.assertRaises(ValueError):
            FullContextAtlas(self.path)
        self.manifest['bindings'] = []
        self.manifest['partitions'][0]['n_context_rows'] = 3
        self.save()
        with self.assertRaises(ValueError):
            FullContextAtlas(self.path).rows('OLD1')

    def test_native_source_cell_order_must_match_core(self):
        path = self.root / 'native/OLD1_symbol_counts.h5ad'
        native = ad.read_h5ad(path)
        native = native[[1, 0, 2]].copy()
        other = self.root / 'reordered'
        other.mkdir()
        native.write_h5ad(other / path.name)
        self.manifest['legacy']['native_source_root'] = str(other)
        self.save()
        with self.assertRaises(ValueError):
            FullContextAtlas(self.path).feature_lookup('OLD1')


if __name__ == '__main__':
    unittest.main()
