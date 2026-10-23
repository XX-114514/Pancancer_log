"""Read small RNA panels from an explicit collection of existing sparse stores.

This is a collection index, not a self-contained AnnData/Zarr replacement and
not a fitted integration model. Missing features have an explicit measured mask.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import zarr
from anndata._io.specs import read_elem


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csr_rows(group, rows):
    """Preserve requested order, including duplicates, without loading all X."""
    rows = np.asarray(rows, dtype=np.int64)
    shape = tuple(map(int, group.attrs['shape']))
    if not len(rows):
        return sp.csr_matrix((0, shape[1]), dtype=group['data'].dtype)
    if rows.min() < 0 or rows.max() >= shape[0]:
        raise IndexError('Source row outside CSR shape')
    splits = np.flatnonzero(np.diff(rows) != 1) + 1
    chunks = []
    for run in np.split(rows, splits):
        ptr = np.asarray(group['indptr'][int(run[0]):int(run[-1]) + 2], dtype=np.int64)
        lo, hi = int(ptr[0]), int(ptr[-1])
        chunks.append(sp.csr_matrix((np.asarray(group['data'][lo:hi]),
                                    np.asarray(group['indices'][lo:hi]), ptr - lo),
                                   shape=(len(run), shape[1])))
    return sp.vstack(chunks, format='csr')


class BackedAtlasCollection:
    def __init__(self, manifest):
        self.manifest_path = Path(manifest)
        self.manifest = json.loads(self.manifest_path.read_text())
        if self.manifest.get('schema') != 'backed_atlas_collection_v1':
            raise ValueError('Unsupported collection schema')
        self.partitions = {p['dataset']: p for p in self.manifest['partitions']}
        if len(self.partitions) != len(self.manifest['partitions']):
            raise ValueError('Repeated dataset partition')
        self._stores, self._arrays, self._rows, self._features = {}, {}, {}, {}

    def store(self, dataset):
        path = self.partitions[dataset]['path']
        if path not in self._stores:
            self._stores[path] = zarr.open_group(path, mode='r')
        return self._stores[path]

    def _array(self, dataset, key):
        cache_key = (self.partitions[dataset]['path'], key)
        if cache_key not in self._arrays:
            self._arrays[cache_key] = np.asarray(read_elem(self.store(dataset)[key]))
        return self._arrays[cache_key]

    def rows(self, dataset, view='context'):
        key = (dataset, view)
        if key in self._rows:
            return self._rows[key]
        p = self.partitions[dataset]
        selector = p['selector']
        if selector['kind'] == 'npy':
            if file_sha256(selector['path']) != selector['sha256']:
                raise ValueError('Changed row selector')
            rows = np.load(selector['path'], allow_pickle=False)
        elif selector['kind'] == 'field_equals':
            values = self._array(dataset, 'obs/' + selector['field']).astype(str)
            rows = np.flatnonzero(values == selector['value'])
        else:
            raise ValueError('Unsupported row selector')
        rows = np.asarray(rows, dtype=np.int64)
        if len(rows) != p['n_context_rows'] or len(np.unique(rows)) != len(rows):
            raise ValueError('Changed partition count or duplicate rows')
        field = p.get('dataset_field', 'dataset')
        if not np.all(self._array(dataset, 'obs/' + field)[rows].astype(str) == dataset):
            raise ValueError('Source identity mismatch')
        if view == 'immune':
            rule = p['immune_selector']
            values = self._array(dataset, 'obs/' + rule['field'])[rows].astype(str)
            rows = rows[np.isin(values, rule['values'])]
            if len(rows) != p['n_immune_rows']:
                raise ValueError('Changed immune-view count')
        elif view != 'context':
            raise ValueError('View must be context or immune')
        self._rows[key] = rows
        return rows

    def feature_lookup(self, dataset):
        if dataset in self._features:
            return self._features[dataset]
        z = self.store(dataset)
        var = read_elem(z['var'])
        source_order = np.asarray(read_elem(z['uns/source_order'])).astype(str).tolist()
        if source_order.count(dataset) != 1:
            raise ValueError('Missing or duplicated source measurement profile')
        measured = np.asarray(read_elem(z['varm/measured_in_source']), dtype=bool)
        if measured.shape != (len(var), len(source_order)):
            raise ValueError('Feature mask shape mismatch')
        available = measured[:, source_order.index(dataset)]
        lookup = {}
        for col, (index, row) in enumerate(var.iterrows()):
            if 'canonical_gene_identity' in row:
                canonical = str(row['canonical_gene_identity'])
            elif var.index.name == 'canonical_gene_identity':
                canonical = str(index)
            else:
                raise ValueError('No explicit canonical feature identity axis')
            if not available[col] or canonical.startswith('AMBIGUOUS::'):
                continue
            labels = {canonical, str(row['gene_symbol'])}
            if str(row['feature_identity_status']) == 'stable_Ensembl_ID':
                labels.update(str(row['source_gene_aliases']).split('|'))
            for name in labels - {'', 'nan', 'None'}:
                lookup.setdefault(name, set()).add(col)
        self._features[dataset] = lookup
        return lookup

    def read_panel(self, dataset, positions, genes, view='immune', layer='counts'):
        positions = np.asarray(positions, dtype=np.int64)
        if positions.ndim != 1 or len(positions) > 4096:
            raise ValueError('Use one-dimensional blocks of at most 4096 rows')
        if len(genes) != len(set(genes)):
            raise ValueError('Repeated requested feature')
        rows = self.rows(dataset, view)
        if len(positions) and (positions.min() < 0 or positions.max() >= len(rows)):
            raise IndexError('Requested position outside selected view')
        lookup = self.feature_lookup(dataset)
        indices = [next(iter(lookup[g])) if len(lookup.get(g, set())) == 1 else None for g in genes]
        resolved = [i for i in indices if i is not None]
        if len(resolved) != len(set(resolved)):
            raise ValueError('Requested aliases resolve to the same physical feature')
        measured = np.asarray([i is not None for i in indices], dtype=bool)
        z = self.store(dataset)
        native = read_csr_rows(z[self.partitions[dataset]['counts_key']], rows[positions])
        if layer == 'analysis_counts':
            key = self.partitions[dataset].get('ambient_removed_key')
            if not key or key not in z:
                raise ValueError('No explicitly stored applied correction layer')
            native = native.astype(np.float64) - read_csr_rows(z[key], rows[positions])
        elif layer != 'counts':
            raise ValueError('Layer must be counts or analysis_counts')
        columns = [native[:, c] if c is not None else sp.csr_matrix((len(positions), 1), dtype=native.dtype) for c in indices]
        panel = sp.hstack(columns, format='csr') if columns else sp.csr_matrix((len(positions), 0))
        return {'counts': panel, 'measured': measured, 'genes': list(genes),
                'source_row_positions': rows[positions], 'view': view, 'layer': layer,
                'missing_policy': 'Zero placeholders require measured=False; never biological zero'}
