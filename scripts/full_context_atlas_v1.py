"""Read frozen legacy and newer atlas partitions without copying native counts.

Native rows and per-source measured features are checked separately. Returned
panel zeros must always be interpreted together with the measurement mask.
"""
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import zarr
from anndata._io.specs import read_elem

from backed_atlas_collection_v1 import BackedAtlasCollection, read_csr_rows


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _strings(values):
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in values], dtype=object)


def _h5_index(path, axis):
    with h5py.File(path, 'r') as root:
        name = root[axis].attrs.get('_index', '_index')
        if isinstance(name, bytes):
            name = name.decode()
        return _strings(root[axis][name][:])


class FullContextAtlas:
    """A provenance-bound collection, not a joint embedding or annotation."""

    def __init__(self, manifest):
        self.manifest_path = Path(manifest)
        self.manifest = json.loads(self.manifest_path.read_text())
        if self.manifest.get('schema') != 'full_context_atlas_collection_v1':
            raise ValueError('Unsupported full-context collection')
        self.partitions = {p['dataset']: p for p in self.manifest['partitions']}
        if len(self.partitions) != len(self.manifest['partitions']):
            raise ValueError('Repeated source key')
        for entry in self.manifest['bindings']:
            if sha256_file(entry['path']) != entry['sha256']:
                raise ValueError('Changed bound input: ' + entry['role'])
        self.primary = BackedAtlasCollection(self.manifest['primary_manifest'])
        self.legacy_config = self.manifest['legacy']
        self.legacy = zarr.open_group(self.legacy_config['canonical_path'], mode='r')
        cfg = self.legacy_config
        with np.load(cfg['roster_path'], allow_pickle=False) as roster:
            self._positions = roster['global_position'].copy()
            self._sample_codes = roster['sample_code'].copy()
            self._label_codes = roster['label_code'].copy()
            self._sample_categories = roster['sample_categories'].copy()
            self._label_categories = roster['label_categories'].copy()
        if not (len(self._positions) == len(self._sample_codes) == len(self._label_codes)):
            raise ValueError('Frozen roster axes differ')
        shape = tuple(self.legacy[cfg['counts_key']].attrs['shape'])
        if (len(self._positions) != len(np.unique(self._positions)) or
                (len(self._positions) and (self._positions.min() < 0 or self._positions.max() >= shape[0]))):
            raise ValueError('Invalid frozen native row positions')
        for codes, categories in [(self._sample_codes, self._sample_categories),
                                  (self._label_codes, self._label_categories)]:
            if len(codes) and (codes.min() < 0 or codes.max() >= len(categories)):
                raise ValueError('Frozen category code outside axis')
        meta = pd.read_csv(cfg['sample_metadata_path'], sep='\t', dtype=str,
                           keep_default_na=False, usecols=['sampleID', 'gse_id', 'eligible_min100_cells'])
        if not meta.sampleID.is_unique:
            raise ValueError('Repeated frozen sample key')
        kept = meta.loc[meta.eligible_min100_cells.eq('True')].set_index('sampleID')
        self._sample_source = kept.gse_id.reindex(self._sample_categories).fillna('').to_numpy()
        self._source_per_roster = self._sample_source[self._sample_codes]
        self._rows, self._fields, self._features = {}, {}, {}
        var = self.legacy['var']
        self._legacy_var = _strings(read_elem(var[var.attrs.get('_index', '_index')]))
        if len(np.unique(self._legacy_var)) != len(self._legacy_var):
            raise ValueError('Legacy feature index is not unique')
        self._legacy_gene_position = {gene: i for i, gene in enumerate(self._legacy_var)}

    def is_legacy(self, dataset):
        provider = self.partitions[dataset]['provider']
        if provider not in {'legacy_frozen_native_adapter', 'backed_existing_primary_collection'}:
            raise ValueError('Unsupported partition provider')
        return provider == 'legacy_frozen_native_adapter'

    def store(self, dataset):
        return self.legacy if self.is_legacy(dataset) else self.primary.store(dataset)

    def rows(self, dataset):
        if dataset not in self._rows:
            if self.is_legacy(dataset):
                selected = self._source_per_roster == dataset
                rows = self._positions[selected]
                observed = np.asarray(read_elem(self.legacy['obs/gse_id']), dtype=object)[rows]
                if not np.all(observed == dataset):
                    raise ValueError('Legacy roster/source mismatch')
            else:
                rows = self.primary.rows(dataset, view='context')
            if len(rows) != self.partitions[dataset]['n_context_rows']:
                raise ValueError('Changed frozen partition count')
            self._rows[dataset] = np.asarray(rows, dtype=np.int64)
        return self._rows[dataset]

    def field_names(self, dataset):
        return list(self.store(dataset)['obs'].attrs['column-order'])

    def frozen_expression_units(self, dataset):
        """Use the frozen roster's keys, preserving the raw obs keys separately."""
        if self.is_legacy(dataset):
            selected = self._source_per_roster == dataset
            units = self._sample_categories[self._sample_codes[selected]]
        else:
            unit_field = self.partitions[dataset].get('unit_field', 'library_id')
            units = self.field(dataset, unit_field)
        if len(units) != len(self.rows(dataset)):
            raise ValueError('Frozen expression-unit axis mismatch')
        expected = self.partitions[dataset].get('n_context_units')
        if expected is not None and len(set(units)) != expected:
            raise ValueError('Changed frozen expression-unit count')
        return units

    def field(self, dataset, name):
        """Return a selected metadata axis; caller keeps individual rows private."""
        z = self.store(dataset)
        if name not in z['obs']:
            raise KeyError(name)
        key = (self.partitions[dataset]['provider'], str(z.store.path), name)
        if key not in self._fields:
            self._fields[key] = np.asarray(read_elem(z['obs/' + name]), dtype=object)
        return self._fields[key][self.rows(dataset)]

    def frozen_labels(self, dataset):
        """Carry original evidence; no label is newly promoted by this reader."""
        if self.is_legacy(dataset):
            selected = self._source_per_roster == dataset
            return {'major': self.field(dataset, 'majorCluster'),
                    'fine': self._label_categories[self._label_codes[selected]],
                    'status': 'frozen_legacy_identity_not_revalidated'}
        annotation = self.primary.partitions[dataset]['annotation']
        if 'mapping' in annotation:
            author = self.field(dataset, annotation['field'])
            major = np.asarray([annotation['mapping'].get(str(v), 'unresolved') for v in author], dtype=object)
            return {'major': major, 'fine': author,
                    'status': 'author_mapping_immune_parents_only_other_context_unresolved'}
        return {'major': self.field(dataset, annotation['major_field']),
                'fine': self.field(dataset, annotation['subcluster_field']),
                'status': 'frozen_immune_v1_identity_not_revalidated'}

    def feature_lookup(self, dataset):
        if dataset in self._features:
            return self._features[dataset]
        if not self.is_legacy(dataset):
            lookup = self.primary.feature_lookup(dataset)
        else:
            path = Path(self.legacy_config['native_source_root']) / (dataset + '_symbol_counts.h5ad')
            genes = _h5_index(path, 'var')
            if len(np.unique(genes)) != len(genes) or any(g not in self._legacy_gene_position for g in genes):
                raise ValueError('Invalid native legacy feature membership')
            gse = np.asarray(read_elem(self.legacy['obs/gse_id']), dtype=object)
            source_rows = np.flatnonzero(gse == dataset)
            obs = self.legacy['obs']
            ids = _strings(read_elem(obs[obs.attrs.get('_index', '_index')]))[source_rows]
            if not np.array_equal(_h5_index(path, 'obs'), ids):
                raise ValueError('Native/core source cell order mismatch')
            lookup = {g: {self._legacy_gene_position[g]} for g in genes}
        self._features[dataset] = lookup
        return lookup

    def read_panel(self, dataset, positions, genes):
        original_positions = np.asarray(positions)
        if original_positions.size and (not np.issubdtype(original_positions.dtype, np.integer)):
            raise ValueError('Selected row positions must be integer indices')
        positions = original_positions.astype(np.int64)
        if positions.ndim != 1 or len(positions) > 4096:
            raise ValueError('Read blocks of at most 4096 selected rows')
        if len(genes) != len(set(genes)):
            raise ValueError('Repeated requested feature')
        rows = self.rows(dataset)
        if len(positions) and (positions.min() < 0 or positions.max() >= len(rows)):
            raise IndexError('Position outside selected source')
        lookup = self.feature_lookup(dataset)
        cols = [next(iter(lookup[g])) if len(lookup.get(g, set())) == 1 else None for g in genes]
        resolved = [c for c in cols if c is not None]
        if len(set(resolved)) != len(resolved):
            raise ValueError('Requested aliases map to one physical feature')
        z = self.store(dataset)
        key = self.legacy_config['counts_key'] if self.is_legacy(dataset) else self.primary.partitions[dataset]['counts_key']
        full = read_csr_rows(z[key], rows[positions])
        native_total = np.asarray(full.sum(axis=1, dtype=np.float64)).ravel()
        columns = [full[:, c] if c is not None else sp.csr_matrix((len(positions), 1), dtype=full.dtype) for c in cols]
        panel = sp.hstack(columns, format='csr') if columns else sp.csr_matrix((len(positions), 0))
        return {'counts': panel, 'measured': np.asarray([c is not None for c in cols]),
                'genes': list(genes), 'native_total_counts': native_total,
                'source_row_positions': rows[positions],
                'count_origin': key, 'native_total_includes_all_stored_features': True}
