"""Native marker evidence with bounded CSR reads and explicit missingness.

This module measures evidence for annotation review. It assigns no cell types
and does not infer malignant status or functional states from a single marker.
"""
import numpy as np
import scipy.sparse as sp

from backed_atlas_collection_v1 import read_csr_rows


LINEAGE_PANELS = {
    'T': ['CD3D', 'CD3E', 'TRAC', 'CD247'],
    'NK_related': ['NKG7', 'GNLY', 'KLRD1', 'FCGR3A'],
    'B': ['CD79A', 'MS4A1', 'CD37', 'CD74'],
    'Plasma': ['MZB1', 'JCHAIN', 'SDC1', 'DERL3'],
    'Myeloid': ['LST1', 'LYZ', 'FCER1G', 'CTSS', 'AIF1'],
    'DC_related': ['FCER1A', 'CD1C', 'CLEC9A', 'CLEC10A'],
    'Neutrophil_related': ['CSF3R', 'FCGR3B', 'S100A8', 'S100A9'],
    'Mast': ['TPSAB1', 'TPSB2', 'KIT', 'CPA3'],
    'Endothelial': ['PECAM1', 'VWF', 'KDR', 'EMCN'],
    'Lymphatic_related': ['PROX1', 'FLT4', 'PDPN', 'CCL21'],
    'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM'],
    'Mural': ['RGS5', 'CSPG4', 'PDGFRB', 'MCAM'],
    'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19'],
    'Glial_Oligodendrocyte_related': ['OLIG1', 'OLIG2', 'MBP', 'PLP1', 'SOX10'],
    'Astrocyte_related': ['GFAP', 'AQP4', 'SLC1A3', 'ALDOC'],
    'Neuron_related': ['SNAP25', 'SYT1', 'RBFOX3', 'STMN2'],
    'Erythroid_related': ['HBB', 'HBA1', 'HBA2', 'AHSP'],
    'Progenitor_related': ['CD34', 'SPINK2', 'GATA2', 'MEIS1'],
}
STATE_PANELS = {
    'cytotoxic_related': ['GZMB', 'PRF1', 'GNLY', 'NKG7'],
    'interferon_related': ['IFIT1', 'ISG15', 'MX1', 'IFIT3'],
    'cycling_related': ['MKI67', 'TOP2A', 'STMN1', 'UBE2C'],
    'regulatory_related': ['FOXP3', 'IL2RA', 'CTLA4', 'TIGIT'],
    'naive_memory_related': ['TCF7', 'LEF1', 'CCR7', 'IL7R'],
    'dysfunction_related': ['PDCD1', 'LAG3', 'HAVCR2', 'TOX'],
}
IDENTITY_CANDIDATE_PANELS = {
    'gamma_delta_T_related': ['TRDC', 'TRGC1', 'TRGC2', 'CD3D'],
}
GENES = sorted({gene for panel in [LINEAGE_PANELS, STATE_PANELS, IDENTITY_CANDIDATE_PANELS]
                for genes in panel.values() for gene in genes})


def gather_csr_rows(group, rows, max_span=8192):
    """Coalesce nearby rows, then restore arbitrary order and duplicates."""
    requested = np.asarray(rows)
    if requested.ndim != 1 or (requested.size and not np.issubdtype(requested.dtype, np.integer)):
        raise ValueError('One-dimensional integer native row positions required')
    if max_span < 1 or max_span > 16384 or len(requested) > 4096:
        raise ValueError('CSR read budget exceeded')
    requested = requested.astype(np.int64)
    shape = tuple(map(int, group.attrs['shape']))
    if not len(requested):
        return sp.csr_matrix((0, shape[1]), dtype=group['data'].dtype)
    if requested.min() < 0 or requested.max() >= shape[0]:
        raise IndexError('Native row outside source matrix')
    unique, inverse = np.unique(requested, return_inverse=True)
    pieces = []
    low = 0
    while low < len(unique):
        high = np.searchsorted(unique, unique[low] + max_span, side='left')
        native = read_csr_rows(group, np.arange(unique[low], unique[high - 1] + 1))
        pieces.append(native[unique[low:high] - unique[low]])
        low = high
    return sp.vstack(pieces, format='csr')[inverse]


def native_marker_block(reader, dataset, positions, genes=GENES):
    rows = reader.rows(dataset)[positions]
    z = reader.store(dataset)
    key = reader.legacy_config['counts_key'] if reader.is_legacy(dataset) else reader.primary.partitions[dataset]['counts_key']
    raw = gather_csr_rows(z[key], rows)
    if not np.isfinite(raw.data).all() or (raw.data < 0).any() or not np.equal(raw.data, np.floor(raw.data)).all():
        raise ValueError('Native input does not meet finite nonnegative count assumption')
    lookup = reader.feature_lookup(dataset)
    cols = [next(iter(lookup[g])) if len(lookup.get(g, set())) == 1 else None for g in genes]
    physical = [c for c in cols if c is not None]
    if len(physical) != len(set(physical)):
        raise ValueError('Marker aliases collide')
    arrays = [raw[:, c] if c is not None else sp.csr_matrix((raw.shape[0], 1), dtype=raw.dtype) for c in cols]
    return {'counts': sp.hstack(arrays, format='csr'),
            'measured': np.asarray([c is not None for c in cols]),
            'native_total_counts': np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel(),
            'full_native_nnz_read': raw.nnz, 'native_rows_read': len(rows)}


class MarkerEvidence:
    def __init__(self, group_names, genes, measured):
        self.group_names = list(group_names)
        self.genes = list(genes)
        self.measured = np.asarray(measured, dtype=bool)
        n, p = len(group_names), len(genes)
        if len(self.measured) != p:
            raise ValueError('Marker measurement axis mismatch')
        self.n_cells = np.zeros(n, dtype=np.int64)
        self.n_evaluable = np.zeros(n, dtype=np.int64)
        self.detected = np.zeros((n, p), dtype=np.int64)
        self.cp10k_sum = np.zeros((n, p), dtype=np.float64)
        self.logcp10k_sum = np.zeros((n, p), dtype=np.float64)

    def update(self, counts, native_totals, codes):
        counts = counts.tocsr()
        totals = np.asarray(native_totals, dtype=np.float64)
        codes = np.asarray(codes)
        if counts.shape != (len(totals), len(self.genes)) or len(codes) != len(totals):
            raise ValueError('Marker/group/denominator axes differ')
        if len(codes) and (not np.issubdtype(codes.dtype, np.integer) or codes.min() < 0 or codes.max() >= len(self.group_names)):
            raise ValueError('Invalid evidence group codes')
        if not np.isfinite(totals).all() or (totals < 0).any():
            raise ValueError('Invalid native totals')
        if not np.isfinite(counts.data).all() or (counts.data < 0).any():
            raise ValueError('Invalid marker values')
        if counts[:, ~self.measured].nnz:
            raise ValueError('Unmeasured marker has nonzero values')
        if np.any(np.asarray(counts.sum(axis=1)).ravel() > totals + 1e-6):
            raise ValueError('Marker sum exceeds complete native RNA denominator')
        valid = totals > 0
        for code in np.unique(codes):
            member = codes == code
            self.n_cells[code] += int(member.sum())
            include = member & valid
            self.n_evaluable[code] += int(include.sum())
            if not include.any():
                continue
            x = counts[include].astype(np.float64)
            self.detected[code] += np.asarray((x > 0).sum(axis=0)).ravel()
            normalized = x.multiply((10000.0 / totals[include])[:, None]).tocsr()
            self.cp10k_sum[code] += np.asarray(normalized.sum(axis=0)).ravel()
            normalized.data = np.log1p(normalized.data)
            self.logcp10k_sum[code] += np.asarray(normalized.sum(axis=0)).ravel()

    def records(self):
        records = []
        for i, group in enumerate(self.group_names):
            n = int(self.n_evaluable[i])
            values = {}
            for j, gene in enumerate(self.genes):
                measured = bool(self.measured[j])
                values[gene] = {'measured': measured,
                                'n_detected': int(self.detected[i, j]) if measured and n else None,
                                'detection_fraction': float(self.detected[i, j] / n) if measured and n else None,
                                'mean_native_CP10K': float(self.cp10k_sum[i, j] / n) if measured and n else None,
                                'mean_native_log1p_CP10K': float(self.logcp10k_sum[i, j] / n) if measured and n else None}
            records.append({'frozen_group': str(group), 'n_archive_cells': int(self.n_cells[i]),
                            'n_native_RNA_evaluable': n, 'marker_evidence': values,
                            'identity_promoted': False, 'state_or_malignancy_promoted': False})
        return records
