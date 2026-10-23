"""Version 2: retain broad family when dropout reference opposes group major.

Evidence levels describe provenance, not calibrated probabilities or independent
truth. Coherent minority-cell evidence is retained separately from group labels;
small discovery branches do not erase direct cell-type evidence.
"""
from native_lineage_routing_v1 import FAMILIES, local_family

FINE_TYPES = ['Treg', 'CD4_T', 'CD8_T', 'T_gamma_delta_RNA_candidate',
              'Neutrophil', 'pDC', 'cDC1', 'cDC2_related', 'Macrophage', 'Monocyte']
FINE_MAJOR = {t: ('T' if t in FINE_TYPES[:4] else 'Myeloid') for t in FINE_TYPES}
PARENT_TYPE = {'Treg': 'CD4_T'}
FINE_PANEL = {'Treg': 'Treg_related', 'CD8_T': 'CD8_T_related',
              'T_gamma_delta_RNA_candidate': 'gamma_delta_T_related',
              'Neutrophil': 'Neutrophil_related', 'pDC': 'pDC_related',
              'cDC1': 'cDC1_related', 'cDC2_related': 'cDC2_related',
              'Macrophage': 'Macrophage_related', 'Monocyte': 'Monocyte_related'}


def low_reference_identity(label):
    """Frozen Low reference vocabulary -> coarse major and eligible fine cues.

Unspecified/precursor/cycling labels do not acquire invented mature fine types.
Coarse family scope is already checked by the complete reference dictionary.
"""
    s = str(label).lower();major = None;fine = set()
    if not s:return major, fine
    if 'nk cells' in s or s == 'transitional nk':major = 'NK'
    elif s.startswith('ilc'):major = 'T_NK'
    elif 't cells' in s or 'thymocytes' in s or s in {'cd8a/a', 'cd8a/b(entry)', 't(agonist)', 'treg(diff)'}:major = 'T'
    elif s in {'plasma cells', 'plasmablasts'}:major = 'Plasma'
    elif 'b cells' in s or 'pre-b' in s or 'pro-b' in s:major = 'B'
    elif any(x in s for x in ['macrophage', 'monocyte', 'neutrophil']) or s.startswith(('dc', 'pdc')) or s in {'mono-mac','mnp','hofbauer cells','kupffer cells','granulocytes','myelocytes','promyelocytes','transitional dc','migratory dcs'}:major = 'Myeloid'
    elif s == 'mast cells':major = 'Mast'
    elif 'erythroid' in s or s == 'erythrocytes':major = 'Erythroid'
    if 'gamma-delta' in s:fine.add('T_gamma_delta_RNA_candidate')
    elif s in {'regulatory t cells', 'treg(diff)'}:fine.update(['Treg', 'CD4_T'])
    elif 'helper t' in s or 'cd4' in s:fine.add('CD4_T')
    elif ('cytotoxic t' in s and 'cd4' not in s) or s.startswith('cd8'):fine.add('CD8_T')
    if s == 'neutrophils':fine.add('Neutrophil')
    elif s == 'pdc':fine.add('pDC')
    elif s == 'dc1':fine.add('cDC1')
    elif s in {'dc','dc2','dc3','transitional dc'}:fine.add('cDC2_related')
    elif 'macrophage' in s or s in {'hofbauer cells','kupffer cells'}:fine.add('Macrophage')
    elif 'monocyte' in s and 'precursor' not in s:fine.add('Monocyte')
    return major, fine


def prune_ancestors(types):
    result = set(types)
    for child, parent in PARENT_TYPE.items():
        if child in result:result.discard(parent)
    return result


def compatible_major(a, b):
    if not a or not b:return True
    if a == b:return True
    return a in {'T_NK', 'B_Plasma'} and local_family(b) == a or b in {'T_NK', 'B_Plasma'} and local_family(a) == b


def decide_cell(native_families, reference_families, native_major, reference_major,
                group_major, group_type, own_fine, reference_fine, positive=True,
                fine_fit=True):
    """Return a guarded candidate label; group memberships remain untouched."""
    native_families, reference_families = set(native_families), set(reference_families)
    if not native_families.issubset(FAMILIES) or not reference_families.issubset(FAMILIES):
        raise ValueError('Unreviewed coarse-family vocabulary')
    own_fine, reference_fine = set(own_fine), set(reference_fine)
    if not own_fine.issubset(FINE_TYPES) or not reference_fine.issubset(FINE_TYPES):
        raise ValueError('Unreviewed fine-type vocabulary')
    result = {'consensus_major': 'unresolved', 'consensus_celltype': 'unresolved',
              'coarse_evidence': 'insufficient_cell_evidence', 'fine_evidence': 'coarse_unresolved',
              'group_major_disagreement': False, 'group_type_disagreement': False,
              'fine_discovery_fit_available': bool(fine_fit), 'direct_fine_candidate': 'unresolved'}
    if not positive:result['coarse_evidence'] = 'zero_native_RNA';return result
    if len(native_families) > 1:result['coarse_evidence'] = 'native_cross_lineage_conflict';return result
    if len(reference_families) > 1:result['coarse_evidence'] = 'reference_cross_lineage_conflict';return result
    if native_families and reference_families and native_families != reference_families:
        result['coarse_evidence'] = 'native_reference_lineage_conflict';return result
    family = next(iter(native_families or reference_families), None)
    if family is None:return result
    if native_major and local_family(native_major) != family:raise ValueError('Native major/family inconsistency')
    if reference_major and local_family(reference_major) != family:raise ValueError('Reference major/family inconsistency')
    group_family = local_family(group_major)
    if native_families:
        major = native_major or family
        evidence = 'native_reference_concordant' if reference_families else 'native_supported_reference_unresolved'
    elif group_family == family:
        major = reference_major or group_major or family
        evidence = 'reference_group_inferred_native_dropout'
        if not compatible_major(reference_major, group_major):
            major = family;evidence = 'within_family_reference_group_conflict_native_dropout'
    else:
        result['coarse_evidence'] = 'reference_only_without_group_corroboration';return result
    if not compatible_major(native_major, reference_major):
        major = family;evidence = 'within_family_native_reference_conflict'
    result.update(consensus_major=major, consensus_celltype=major, coarse_evidence=evidence,
                  fine_evidence='coarse_identity_retained',
                  group_major_disagreement=bool(group_family and (group_family != family or not compatible_major(group_major,major))))
    if evidence in {'within_family_native_reference_conflict','within_family_reference_group_conflict_native_dropout'}:return result
    eligible_fine = {t for t in FINE_TYPES if compatible_major(FINE_MAJOR[t],major)}
    direct = prune_ancestors(own_fine & reference_fine & eligible_fine)
    if len(direct) == 1 and native_families:
        fine = next(iter(direct));result.update(consensus_celltype=fine, direct_fine_candidate=fine,
            fine_evidence='own_native_and_unique_reference_supported',
            group_type_disagreement=bool(group_type not in {'unresolved',fine,FINE_MAJOR[fine]}))
        if fine == 'T_gamma_delta_RNA_candidate':result['fine_evidence'] = 'own_gamma_RNA_and_reference_candidate_requires_orthogonal_review'
        return result
    if len(direct) > 1:result['fine_evidence'] = 'multiple_direct_fine_candidates';return result
    if result['group_major_disagreement'] or group_major == 'unresolved':
        result['fine_evidence'] = 'group_cannot_override_cell_major';return result
    if group_type not in eligible_fine:return result
    if not fine_fit:result['fine_evidence'] = 'coarse_only_branch_without_direct_fine_evidence';return result
    ref_specific = prune_ancestors(reference_fine & eligible_fine)
    compatible_reference = group_type in reference_fine
    # A broad CD4 reference is compatible with a marker-supported Treg proposal.
    broad_reference = group_type in PARENT_TYPE and ref_specific == {PARENT_TYPE[group_type]}
    if ref_specific and not compatible_reference and not broad_reference:
        result['fine_evidence'] = 'group_reference_fine_disagreement';return result
    if group_type in own_fine and native_families:
        result.update(consensus_celltype=group_type,fine_evidence='own_native_and_group_supported_reference_unresolved')
        if group_type == 'T_gamma_delta_RNA_candidate':result['fine_evidence'] = 'own_gamma_RNA_and_group_candidate_requires_orthogonal_review'
    elif compatible_reference and group_type != 'T_gamma_delta_RNA_candidate':
        result.update(consensus_celltype=group_type,fine_evidence='reference_group_inferred_without_complete_native_support')
    elif group_type == 'T_gamma_delta_RNA_candidate':result['fine_evidence'] = 'gamma_candidate_requires_own_RNA_evidence'
    else:result['fine_evidence'] = 'group_fine_insufficient_cell_evidence'
    return result
