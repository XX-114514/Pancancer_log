"""Native marker review and immutable, scope-qualified reference proposals.

Reference probabilities are one-vs-rest scores, not calibrated identity truth.
Missing model features are omitted with explicit coverage. Native lineage/type
rules are annotation evidence; RNA-state and malignancy calls remain separate.
"""
import hashlib
import re

import numpy as np
from scipy.special import expit

from native_panel_aggregate_v1 import LINEAGE_PANELS, STATE_PANELS


TYPE_PANELS = {
    'Treg_related': ['FOXP3', 'IL2RA', 'CTLA4'],
    'CD8_T_related': ['CD8A', 'CD8B', 'CD3D', 'CD3E'],
    'CD4_T_related': ['CD4', 'IL7R', 'CD3D', 'CD3E'],
    'gamma_delta_T_related': ['TRDC', 'TRGC1', 'TRGC2', 'CD3D'],
    'Macrophage_related': ['C1QA', 'C1QB', 'APOE', 'CSF1R'],
    'Monocyte_related': ['FCN1', 'S100A8', 'S100A9', 'LYZ'],
    'cDC1_related': ['CLEC9A', 'XCR1', 'BATF3'],
    'cDC2_related': ['CD1C', 'FCER1A', 'CLEC10A'],
    'pDC_related': ['CLEC4C', 'IL3RA', 'TCF4', 'GZMB'],
}
REVIEW_GENES = sorted({g for panels in [LINEAGE_PANELS, STATE_PANELS, TYPE_PANELS]
                       for panel in panels.values() for g in panel})


def model_fingerprint(model):
    digest = hashlib.sha256()
    digest.update('\0'.join(model.features).encode())
    digest.update('\0'.join(model.cell_types).encode())
    for v in [model.classifier.coef_, model.classifier.intercept_, model.scaler.mean_, model.scaler.scale_]:
        x = np.asarray(v);digest.update(str(x.dtype).encode());digest.update(str(x.shape).encode());digest.update(x.tobytes())
    return digest.hexdigest()


def predict_reference(log_expression, genes, model, minimum_feature_fraction=.8):
    """Reproduce installed CellTypist scaling and scores without mutating it.

    Input already uses complete native RNA totals; never renormalize this model
    subset. Only measured, uniquely resolved physical columns are supplied.
    """
    if len(genes) != len(set(genes)) or log_expression.shape[1] != len(genes):
        raise ValueError('Reference input gene axis invalid')
    values = log_expression.data if hasattr(log_expression, 'tocsr') else np.asarray(log_expression)
    if not np.isfinite(values).all() or (values < 0).any() or (values > np.log1p(10000) + 1e-5).any():
        raise ValueError('Expected native log1p CP10K')
    positions = {g: j for j, g in enumerate(model.features)}
    query_columns = [i for i, g in enumerate(genes) if g in positions]
    model_columns = np.array([positions[genes[i]] for i in query_columns], dtype=np.int64)
    n = log_expression.shape[0]
    coverage = len(model_columns) / len(model.features)
    result = {'feature_fraction': coverage, 'n_model_features_measured': len(model_columns),
              'best_code': np.full(n, -1, np.int16), 'unique_code': np.full(n, -1, np.int16),
              'n_positive_labels': np.full(n, -1, np.int16), 'maximum_probability': np.full(n, np.nan),
              'status': 'insufficient_measured_model_features'}
    if coverage < minimum_feature_fraction:return result
    x = log_expression[:, query_columns]
    if hasattr(x, 'toarray'):x = x.toarray()
    means = model.scaler.mean_[model_columns] if model.scaler.with_mean else 0
    x = (x - means) / model.scaler.scale_[model_columns]
    x[x > 10] = 10
    score = x @ model.classifier.coef_[:, model_columns].T + model.classifier.intercept_
    if score.shape[1] == 1:score = np.column_stack([-score[:, 0], score[:, 0]])
    prob = expit(score)
    positive = (prob > .5).sum(axis=1)
    best = score.argmax(axis=1)
    result.update(status='reference_proposal_not_identity_truth', best_code=best.astype(np.int16),
                  unique_code=np.where(positive == 1, best, -1).astype(np.int16),
                  n_positive_labels=positive.astype(np.int16), maximum_probability=prob.max(axis=1),
                  decision_scores=score, probability_scores=prob)
    return result


def native_measured_genes(reader, dataset):
    """All uniquely measured native gene names, one per physical column."""
    from anndata._io.specs import read_elem
    z = reader.store(dataset);var = read_elem(z['var']);lookup = reader.feature_lookup(dataset)
    names, columns = [], []
    for col, index in enumerate(var.index):
        name = str(var.iloc[col]['gene_symbol']) if 'gene_symbol' in var.columns else str(index)
        if lookup.get(name, set()) == {col}:
            names.append(name);columns.append(col)
    if len(names) != len(set(names)):raise ValueError('Native gene-name ambiguity escaped lookup')
    return names, np.asarray(columns, dtype=np.int64)


def native_group_evidence(raw, genes, totals, group_codes, n_groups):
    """Complete measured-native-gene descriptive statistics, no DE p-values."""
    import scipy.sparse as sp
    totals = np.asarray(totals, dtype=float);group_codes = np.asarray(group_codes)
    if raw.shape != (len(totals), len(genes)) or len(group_codes) != len(totals):raise ValueError('Native group axis differs')
    if len(set(genes)) != len(genes) or n_groups < 1:raise ValueError('Invalid gene or group axis')
    if len(group_codes) and (not np.issubdtype(group_codes.dtype, np.integer) or group_codes.min() < 0 or group_codes.max() >= n_groups):
        raise ValueError('Invalid group code')
    if not np.isfinite(totals).all() or (totals < 0).any() or not np.isfinite(raw.data).all() or (raw.data < 0).any():
        raise ValueError('Invalid native counts or totals')
    if np.any(np.asarray(raw.sum(axis=1)).ravel() > totals + 1e-6):raise ValueError('Gene sum exceeds native denominator')
    good = totals > 0
    group = sp.csr_matrix((np.ones(good.sum()), (group_codes[good], np.flatnonzero(good))), shape=(n_groups, len(totals)))
    scale = np.zeros(len(totals), float);scale[good] = 10000 / totals[good]
    cp = raw.astype(float).multiply(scale[:, None]).tocsr()
    log = cp.copy();log.data = np.log1p(log.data)
    return {'n': np.bincount(group_codes[good], minlength=n_groups),
            'cp_sum': (group @ cp).toarray(), 'log_sum': (group @ log).toarray(),
            'detected': (group @ (raw > 0)).toarray()}


def marker_calls(detected, genes, measured):
    """Per-cell multi-gene evidence; no lineage diagnosis from a single gene."""
    lookup = {g: j for j, g in enumerate(genes)}
    out = {}
    for name, panel in {**LINEAGE_PANELS, **TYPE_PANELS}.items():
        cols = [lookup[g] for g in panel if g in lookup and measured[lookup[g]]]
        out[name] = (detected[:, cols].sum(axis=1) >= 2) if len(cols) >= 2 else np.zeros(len(detected), bool)
    def positive(g):return detected[:, lookup[g]] if g in lookup and measured[lookup[g]] else np.zeros(len(detected), bool)
    # CD37/CD74 alone also occur outside B cells. A B-lineage anchor is needed.
    out['B'] &= positive('CD79A') | positive('MS4A1')
    out['Treg_related'] = positive('FOXP3') & (positive('IL2RA') | positive('CTLA4'))
    out['gamma_delta_T_related'] = positive('TRDC') & (positive('TRGC1') | positive('TRGC2')) & positive('CD3D')
    out['CD8_T_related'] = positive('CD8A') & positive('CD8B') & (positive('CD3D') | positive('CD3E'))
    return out


def classify_group(panel_fraction, detection, reference_unique_counts, n_cells, minimum=30):
    """Conservative local group evidence calls, retained separately from truth."""
    families = {'T': 'T', 'NK': 'NK_related', 'B': 'B', 'Plasma': 'Plasma', 'Myeloid': 'Myeloid',
                'Mast': 'Mast', 'Endothelial': 'Endothelial', 'Fibroblast': 'Fibroblast',
                'Mural': 'Mural', 'Epithelial': 'Epithelial', 'Glial': 'Glial_Oligodendrocyte_related',
                'Astrocyte': 'Astrocyte_related', 'Neuron': 'Neuron_related', 'Erythroid': 'Erythroid_related'}
    supported = {name for name, panel in families.items() if panel_fraction.get(panel, 0) >= .2}
    if panel_fraction.get('DC_related', 0) >= .2:supported.add('Myeloid')
    if panel_fraction.get('Neutrophil_related', 0) >= .2:supported.add('Myeloid')
    if panel_fraction.get('Lymphatic_related', 0) >= .2:supported.add('Endothelial')
    if 'T' in supported:supported.discard('NK')
    if 'Plasma' in supported:supported.discard('B')
    out = {'local_major': 'unresolved', 'local_celltype': 'unresolved',
           'annotation_evidence_status': 'insufficient_native_lineage_evidence',
           'native_supported_families': sorted(supported), 'confidence': 'unresolved',
           'malignant_status': 'not_inferred_from_lineage', 'RNA_state_separate': True}
    if n_cells < minimum:out['annotation_evidence_status'] = 'small_group_review';return out
    if len(supported) > 1:out['annotation_evidence_status'] = 'cross_lineage_marker_conflict_review';return out
    if not supported:return out
    major = next(iter(supported))
    out.update(local_major=major, local_celltype=major, confidence='native_multi_marker_supported',
               annotation_evidence_status='local_native_lineage_supported_fine_review')
    def reference_fraction(predicate):
        return sum(n for label, n in reference_unique_counts.items() if predicate(label.lower())) / n_cells
    if major == 'T':
        out['gamma_delta_RNA_candidate_fraction_for_review'] = panel_fraction.get('gamma_delta_T_related', 0)
        if panel_fraction.get('gamma_delta_T_related', 0) >= .05 and reference_fraction(lambda s: 'gamma-delta' in s) >= .5:
            out.update(local_celltype='T_gamma_delta_RNA_candidate', confidence='RNA_candidate_requires_orthogonal_review')
        elif panel_fraction.get('Treg_related', 0) >= .05 and reference_fraction(lambda s: 'regulatory t' in s or 'treg' in s) >= .5:
            out.update(local_celltype='Treg', confidence='native_and_reference_supported')
        elif panel_fraction.get('CD8_T_related', 0) >= .05 and reference_fraction(lambda s: ('cytotoxic t' in s and 'cd4' not in s) or s.startswith('cd8')) >= .5:
            out.update(local_celltype='CD8_T', confidence='native_and_reference_supported')
        elif detection.get('CD4', 0) >= .1 and reference_fraction(lambda s: 'helper t' in s or 'cd4' in s) >= .5:
            out.update(local_celltype='CD4_T', confidence='native_and_reference_supported')
    elif major == 'Myeloid':
        if panel_fraction.get('Neutrophil_related', 0) >= .2 and reference_fraction(lambda s: 'neutrophil' in s) >= .5:
            out.update(local_celltype='Neutrophil', confidence='native_and_reference_supported')
        elif panel_fraction.get('pDC_related', 0) >= .2 and reference_fraction(lambda s: s=='pdc') >= .5:
            out.update(local_celltype='pDC', confidence='native_and_reference_supported')
        elif panel_fraction.get('cDC1_related', 0) >= .05 and reference_fraction(lambda s: s=='dc1') >= .5:
            out.update(local_celltype='cDC1', confidence='native_and_reference_supported')
        elif panel_fraction.get('cDC2_related', 0) >= .2 and reference_fraction(lambda s: s in {'dc','dc2','dc3','transitional dc'}) >= .5:
            out.update(local_celltype='cDC2_related', confidence='native_and_reference_supported')
        elif panel_fraction.get('Macrophage_related', 0) >= .2 and reference_fraction(lambda s: 'macrophage' in s or s=='mono-mac') >= .5:
            out.update(local_celltype='Macrophage', confidence='native_and_reference_supported')
        elif panel_fraction.get('Monocyte_related', 0) >= .2 and reference_fraction(lambda s: 'monocyte' in s) >= .5:
            out.update(local_celltype='Monocyte', confidence='native_and_reference_supported')
    return out


def descriptive_markers(genes, group_cp, group_det, other_cp, other_det, minimum_detection=.1, limit=20):
    """Expression-specificity fallback; never label these as inferential DEGs."""
    genes = np.asarray(genes)
    effect = np.log2((group_cp + .1) / (other_cp + .1))
    keep = ((group_det >= minimum_detection) & (group_det > other_det) & (effect > 0) &
            np.asarray([not re.match(r'^(MT-|RPL|RPS|ENSG|IG[HKL][VDJC])', g) for g in genes]))
    rank = effect * np.maximum(group_det - other_det, 0)
    order = sorted(np.flatnonzero(keep), key=lambda j: (-rank[j], str(genes[j])))[:limit]
    return [{'gene': str(genes[j]), 'native_mean_CP10K': float(group_cp[j]),
             'native_detection_fraction': float(group_det[j]), 'comparator_detection_fraction': float(other_det[j]),
             'descriptive_log2_mean_ratio': float(effect[j]), 'marker_evidence_type': 'own_group_expression_specificity_not_inferential_DEG'}
            for j in order]
