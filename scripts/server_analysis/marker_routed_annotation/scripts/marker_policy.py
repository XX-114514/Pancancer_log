"""Exploratory marker routing, independent of reference-model predictions."""
import numpy as np
import pandas as pd

PANELS = {
    'T': ['CD3D', 'CD3E', 'TRAC', 'CD2'],
    'NK': ['NKG7', 'GNLY', 'KLRD1'],
    'B': ['MS4A1', 'CD79A', 'CD19'],
    'plasma': ['JCHAIN', 'MZB1', 'IGHG1'],
    'myeloid': ['LST1', 'TYROBP', 'FCER1G', 'CSF1R'],
    'mast': ['TPSAB1', 'TPSB2', 'CPA3', 'KIT'],
    'pDC': ['LILRA4', 'IL3RA', 'GZMB'],
    'fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM'],
    'endothelial': ['PECAM1', 'VWF', 'CDH5', 'EMCN'],
    'mural': ['RGS5', 'PDGFRB', 'CSPG4', 'NOTCH3', 'MYH11'],
    'epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19'],
    'hepatocyte': ['ALB', 'APOA1', 'TTR', 'APOA2'],
    'melanocytic': ['PMEL', 'MLANA', 'TYR', 'DCT'],
    'glial': ['GFAP', 'OLIG1', 'OLIG2', 'SOX10'],
}
IMMUNE_PANELS = {'T', 'NK', 'B', 'plasma', 'myeloid', 'mast', 'pDC'}
MIN_SCORE, MIN_GAP, MIN_MARKERS = 1.0, 0.25, 2
IMMUNE_HIGH = {'B cells', 'B-cell lineage', 'DC', 'DC precursor',
              'Double-negative thymocytes', 'Double-positive thymocytes', 'ETP',
              'Granulocytes', 'HSC/MPP', 'ILC', 'ILC precursor', 'MNP',
              'Macrophages', 'Mast cells', 'Mono-mac', 'Monocyte precursor',
              'Monocytes', 'Myelocytes', 'Plasma cells', 'Promyelocytes',
              'T cells', 'pDC', 'pDC precursor'}
EXPLICIT_NONIMMUNE = {'Epithelial cells', 'Endothelial cells', 'Fibroblasts',
                     'Erythrocytes', 'Erythroid', 'Megakaryocytes/platelets',
                     'Megakaryocyte precursor', 'Early MK', 'Platelets'}


def score_clusters(marker_expression, groups):
    means = marker_expression.groupby(groups, sort=True).mean()
    panel = pd.DataFrame(index=means.index)
    evidence = pd.DataFrame(index=means.index)
    for name, genes in PANELS.items():
        present = [g for g in genes if g in means]
        panel[name] = means[present].mean(1) if len(present) >= MIN_MARKERS else 0.
        evidence[name] = (means[present] >= .5).sum(1)
    # Competition occurs between broad compartments, not between T and NK.
    broad = panel.drop(columns=list(IMMUNE_PANELS)).copy()
    broad['immune'] = panel[sorted(IMMUNE_PANELS)].max(1)
    winner = broad.idxmax(1)
    immune_winner = panel[sorted(IMMUNE_PANELS)].idxmax(1)
    best = np.sort(broad.to_numpy(), axis=1)
    result = pd.DataFrame({'marker_parent': winner, 'marker_score': best[:, -1],
                           'marker_gap': best[:, -1] - best[:, -2]})
    support = [evidence.loc[c, immune_winner[c] if w == 'immune' else w]
               for c, w in winner.items()]
    result['supporting_markers'] = support
    result['marker_panel'] = [immune_winner[c] if w == 'immune' else w
                              for c, w in winner.items()]
    ok = ((result.marker_score >= MIN_SCORE) & (result.marker_gap >= MIN_GAP)
          & (result.supporting_markers >= MIN_MARKERS))
    result['marker_gate_status'] = np.where(ok, 'provisional_marker_supported',
                                          'weak_or_competing_marker_evidence')
    result.loc[~ok, 'marker_parent'] = 'unresolved'
    return result.join(panel.add_prefix('panel_mean_'))


def routes(parent, malignancy):
    parent, malignancy = np.asarray(parent), np.asarray(malignancy)
    if not np.isin(malignancy, ['malignant', 'uncertain', 'non_malignant']).all():
        raise ValueError('Unknown or missing malignancy calls')
    result = np.full(len(parent), 'unresolved_no_compatible_reference', dtype=object)
    result[parent == 'immune'] = 'immune_celltypist'
    result[np.isin(parent, ['endothelial', 'mural'])] = 'vascular_celltypist'
    result[np.isin(parent, ['fibroblast', 'epithelial', 'hepatocyte', 'glial'])] = 'singler_hpca'
    result[malignancy == 'malignant'] = 'excluded_malignant'
    return result


def vascular_family(label):
    label = label.lower()
    if label.endswith('_ec') or '_ec_' in label or label.endswith('_lec'):
        return 'endothelial'
    if label.endswith('_pc') or label.endswith('_smc') or label == 'smc_pc_intermediate':
        return 'mural'
    return 'unsupported'


def counterexamples():
    parents = ['immune', 'fibroblast', 'endothelial', 'mural', 'unresolved',
               'immune', 'melanocytic', 'epithelial']
    calls = ['non_malignant'] * 5 + ['malignant', 'uncertain', 'uncertain']
    assert routes(parents, calls).tolist() == [
        'immune_celltypist', 'singler_hpca', 'vascular_celltypist',
        'vascular_celltypist', 'unresolved_no_compatible_reference',
        'excluded_malignant', 'unresolved_no_compatible_reference', 'singler_hpca']
    assert vascular_family('ven_ec_1') == 'endothelial'
    assert vascular_family('smc_pc_intermediate') == 'mural'
    assert vascular_family('Fibroblasts') == 'unsupported'
    # A missing panel, weak signal and competing compartments must not be forced.
    genes = sorted({g for gs in PANELS.values() for g in gs})
    x = pd.DataFrame(0., index=['immune', 'fibro', 'weak', 'mixed'], columns=genes)
    x.loc['immune', PANELS['T']] = 3.
    x.loc['fibro', PANELS['fibroblast']] = 3.
    x.loc['weak', 'CD3D'] = .4
    x.loc['mixed', PANELS['T'] + PANELS['fibroblast']] = 3.
    g = score_clusters(x, x.index)
    assert g.loc['immune', 'marker_parent'] == 'immune'
    assert g.loc['fibro', 'marker_parent'] == 'fibroblast'
    assert g.loc['weak', 'marker_parent'] == 'unresolved'
    assert g.loc['mixed', 'marker_parent'] == 'unresolved'
    return {'routing_cases': len(parents), 'marker_cases': len(x), 'passed': True}
