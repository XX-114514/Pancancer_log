"""Conservative source-text channels, not an annotation model."""
import re
import numpy as np
from NSCLC_external_lineage_role_rules_v1 import author_role, PARENTS

CHANNELS = ['author_exact_VCAN_mac', 'author_gamma_with_mac',
            'author_gamma_with_myeloid_parent', 'author_marker_with_myeloid_parent',
            'frozen_marker_grouping']
T_PARENTS = PARENTS['author_CD8'] | PARENTS['author_Treg']
M_PARENTS = PARENTS['author_Macrophage']


def gamma_text(text):
    raw = str(text).strip().casefold()
    s = raw.replace('γ', 'gamma ').replace('δ', 'delta ').replace('α', 'alpha ').replace('β', 'beta ')
    s = re.sub(r'[_\-]+', ' ', s)
    s = re.sub(r'\s+', ' ', s)
    explicit = bool(re.search(r'gamma\s*delta|\b(?:gdt?|tgd|tcrgdt?)\b', s))
    mixed = explicit and bool(re.search(r'\b(?:nk|nkt)\b|alpha\s*beta', s))
    marker = bool(re.search(r'\btrd[cv]\d*\b|\btrgc?\d*\b|\bv\s*delta\s*\d+\b', s))
    return 2 if mixed else 1 if explicit else 3 if marker else 0


def vcan_text(text):
    s = re.sub(r'[_\-]+', ' ', str(text).strip().casefold())
    return bool(re.search(r'\bvcan\b', s))


def fixture():
    cases = {'T_gdT_TRDV1': 1, 'gamma-delta T cells': 1, 'Tγδ': 1,
             'γδ T cells': 1, 'Tgd': 1, 'TCRgd': 1, 'gdT/NK': 2,
             'gamma delta / alpha beta T': 2, 'T_TRDC_001': 3,
             'TRDV1_high': 3, 'Vδ1_T': 3, 'NKT': 0, 'CD8_GZMK': 0,
             'Treg_FOXP3': 0, 'unknown': 0, '': 0}
    for label, expected in cases.items():
        assert gamma_text(label) == expected, label
    assert vcan_text('Mφ_VCAN') and vcan_text('Macrophage_VCAN_2')
    assert not vcan_text('Macro_VCANlike')
    assert author_role('VCAN monocyte/macrophage') != 2
    return dict(status='PASS_EXPLICIT_MARKER_MIXED_AND_UNRELATED_LABEL_FIXTURE', examples=len(cases))


def masks(chunk, catalog):
    author = chunk['author']
    classes = np.array([gamma_text(v) for v in author.categories], np.int8)
    vcan = np.array([vcan_text(v) for v in author.categories], bool)
    macro = np.array([author_role(v) == 2 for v in author.categories], bool)
    c = np.zeros(len(author), np.int8)
    mv = np.zeros(len(author), bool)
    mm = np.zeros(len(author), bool)
    valid = author.codes >= 0
    c[valid] = classes[author.codes[valid]]
    mv[valid] = vcan[author.codes[valid]]
    mm[valid] = macro[author.codes[valid]]
    major = np.asarray(chunk['major'])[chunk['node'].codes]
    T = np.isin(major, list(T_PARENTS)); M = np.isin(major, list(M_PARENTS))
    explicit = (c == 1) & T
    marker = np.isin(c, [1, 3]) & T
    mac = mm & M
    names = np.asarray(chunk['node'].categories).astype(str)
    labels = catalog.set_index('node_id').loc[names, 'label'].astype(str)
    fg = np.array([gamma_text(v) in [1, 3] for v in labels], bool)[chunk['node'].codes] & T
    fm = np.array([vcan_text(v) for v in labels], bool)[chunk['node'].codes] & M
    pairs = [(mac & mv, explicit), (mac, explicit), (M, explicit), (M, marker), (fm, fg)]
    diagnostics = dict(missing_author_label=~valid, mixed_gamma_label=c == 2,
                       explicit_gamma_parent_discordant=(c == 1) & ~T,
                       marker_gamma_parent_discordant=(c == 3) & ~T)
    return pairs, diagnostics, classes, vcan, macro
