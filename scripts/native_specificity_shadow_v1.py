"""Diagnostic family-support sensitivity; never releases an identity or QC call.

Only two newly proposed restrictions are compared. Existing B, erythroid and
neutrophil anchors remain exactly as supplied by the accepted evidence layer.
Fewer conflicting families is not evidence of better annotation accuracy.
"""

FAMILY_COMPONENTS = {
    'T_NK': ('T', 'NK_related'),
    'B_Plasma': ('B', 'Plasma'),
    'Myeloid': ('Myeloid', 'DC_related', 'Neutrophil_related'),
    'Mast': ('Mast',),
    'Endothelial': ('Endothelial', 'Lymphatic_related'),
    'Fibroblast': ('Fibroblast',),
    'Mural': ('Mural',),
    'Epithelial': ('Epithelial',),
    'Glial': ('Glial_Oligodendrocyte_related',),
    'Astrocyte': ('Astrocyte_related',),
    'Neuron': ('Neuron_related',),
    'Erythroid': ('Erythroid_related',),
}
ANCHORS = {'NK_related': ('GNLY', 'KLRD1'),
           'Myeloid': ('LST1', 'LYZ', 'AIF1')}
SCENARIOS = ('baseline', 'NK_anchor', 'myeloid_anchor', 'both_anchors')


def shadow_support(calls, detected_genes, measured_genes, scenario):
    """Restrict existing evidence only when every diagnostic anchor is measured.

    Original panel calls are frozen inputs, not independently validated labels.
    A T signal preserves T/NK-family evidence when weak NK evidence is removed;
    a DC/neutrophil signal likewise preserves myeloid-family evidence. Missing
    measurements never count as absence. No expression threshold is optimized.
    """
    if scenario not in SCENARIOS:
        raise ValueError('Unreviewed diagnostic scenario')
    required={p for ps in FAMILY_COMPONENTS.values() for p in ps}
    if not required.issubset(calls) or any(type(calls[p]) is not bool for p in required):
        raise ValueError('Complete Boolean accepted panel calls required')
    detected, measured=set(detected_genes),set(measured_genes)
    if not detected.issubset(measured):
        raise ValueError('Detection outside measured gene scope')
    adjusted=dict(calls);review={}
    for panel,anchors in ANCHORS.items():
        complete=set(anchors).issubset(measured)
        present=bool(set(anchors)&detected)
        weak=bool(calls[panel] and complete and not present)
        enabled=(scenario=='both_anchors' or scenario=='NK_anchor' and panel=='NK_related'
                 or scenario=='myeloid_anchor' and panel=='Myeloid')
        if enabled and weak:adjusted[panel]=False
        review[panel]={'all_anchors_measured':complete,
                       'any_anchor_detected':present,
                       'original_support_without_measured_anchor_detection':weak,
                       'restriction_applied':bool(enabled and weak)}
    family=lambda x:{f for f,ps in FAMILY_COMPONENTS.items() if any(x[p] for p in ps)}
    original,shadow=family(calls),family(adjusted)
    if not shadow.issubset(original):raise AssertionError('Diagnostic created new family support')
    return {'original_families':original,'shadow_families':shadow,
            'anchor_review':review,'identity_assigned':False,
            'doublet_or_ambient_diagnosed':False,'cells_removed':False}
