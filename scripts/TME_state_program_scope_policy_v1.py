"""Exact resource-parent scope with legacy-default intersection sensitivity.

This is an applicability rule for diagnostic expression-program scores, not
an identity classifier or proof of functional cell state.
"""
IMMUNE_COMPARTMENTS={'T_NK','B_Plasma','Myeloid','Mast'}
CD4={'CD4_T_reference_related'}
CD8={'CD8_T_reference_related'}
TREG={'Treg_reference_related'}
DC={'DC_reference_related','cDC1_reference_related','cDC2_reference_related','DC3_reference_related'}

def default_scope(default,context):
    canonical=context.get('canonical_identity');parent=context.get('reference_parent')
    if default=='global':return True
    if default=='B_memory':return canonical=='Memory_B_reference_related' and parent=='B'
    if default=='CD4_T':return canonical in CD4 and parent=='T'
    if default=='CD8_T':return canonical in CD8 and parent=='T'
    if default=='Treg':return canonical in TREG and parent=='T'
    if default=='DC':return canonical in DC and parent=='Myeloid'
    if default=='pDC':return canonical=='pDC_reference_related' and parent=='Myeloid'
    if default=='MDC':return parent=='Myeloid'
    if default=='mast':return parent=='Mast'
    raise ValueError('Unreviewed default resource identity scope: '+str(default))

def resource_parent_scope(parent,context):
    biological=context.get('reference_parent');canonical=context.get('canonical_identity')
    if parent=='all-cells':return True
    if parent=='leukocyte':return context['compartment'] in IMMUNE_COMPARTMENTS
    if parent=='B':return biological=='B'
    if parent=='CD4-T':return biological=='T' and canonical in CD4
    if parent=='CD8-T':return biological=='T' and canonical in CD8
    if parent=='Treg':return biological=='T' and canonical in TREG
    if parent=='TNK':return biological in {'T','NK'}
    if parent=='abT':return biological=='T' and canonical in CD4|CD8|TREG
    if parent=='Mac':return biological=='Myeloid' and canonical=='Macrophage_reference_related'
    if parent=='DC':return biological=='Myeloid' and canonical in DC|{'pDC_reference_related'}
    if parent=='p-DC':return biological=='Myeloid' and canonical=='pDC_reference_related'
    if parent=='mast':return biological=='Mast'
    # The abbreviated source hierarchy node M is not silently expanded into
    # either monocytes or all myeloid cells before its definition is verified.
    if parent=='M':return False
    raise ValueError('Unreviewed Cytopus resource parent: '+str(parent))

def applicable(default,metadata,context,mode='legacy_default_intersection'):
    if not metadata.get('exact_Cytopus_gene_match',False):return False,'resource_exact_match_pending'
    parents=metadata['reference_parents']
    if not parents:return False,'resource_parent_unresolved'
    if parents==['M']:return False,'abbreviated_resource_parent_definition_pending'
    original=default_scope(default,context)
    biological=any(resource_parent_scope(p,context) for p in parents)
    if mode=='legacy_default_intersection':
        ok=original and biological
        return ok,('default_and_resource_parent_intersection' if ok else 'outside_reviewed_reference_scope')
    if mode=='resource_hierarchy_primary':
        return biological,('resource_parent_within_legacy_default' if biological and original else
            'explicit_resource_parent_beyond_legacy_default' if biological else 'outside_reviewed_reference_scope')
    raise ValueError('Unreviewed programme applicability mode: '+str(mode))
