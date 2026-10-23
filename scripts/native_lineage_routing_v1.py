"""Provisional coarse routing with explicit native/reference disagreement.

References and local clusters are proposals, not independent ground truth.
Immune-reference nonhematopoietic predictions cannot label tumor/stromal cells.
"""
import numpy as np

FAMILIES=['T_NK','B_Plasma','Myeloid','Mast','Endothelial','Fibroblast','Mural',
          'Epithelial','Glial','Astrocyte','Neuron','Erythroid']
UNRESOLVED='unresolved_context'
STATUS=['zero_native_RNA','native_cross_lineage_review','reference_cross_lineage_review',
        'native_reference_disagreement','native_reference_concordant',
        'native_local_group_concordant_reference_unresolved','native_only_provisional',
        'reference_local_group_concordant_native_dropout','insufficient_coarse_evidence']


def reference_family(model,label):
    """Explicit frozen-model vocabulary; unexpected classes require review."""
    s=label.lower()
    if model.startswith('immune'):
        ignored={'cmp','elp','etp','early mk','early lymphoid/t lymphoid','gmp','hsc/mpp','memp',
                 'megakaryocyte precursor','megakaryocyte-erythroid-mast cell progenitor',
                 'megakaryocytes/platelets','cycling cells','endothelial cells','epithelial cells','fibroblasts'}
        if s in ignored:return None
        if ('t cells' in s or 'thymocytes' in s or 'nk cells' in s or s.startswith('ilc') or
            s in {'cd8a/a','cd8a/b(entry)','t(agonist)','treg(diff)','transitional nk'}):return 'T_NK'
        if ('b cells' in s or 'pre-b' in s or 'pro-b' in s or s in {'b-cell lineage','plasma cells','plasmablasts'}):return 'B_Plasma'
        if s=='mast cells':return 'Mast'
        if ('macrophage' in s or 'monocyte' in s or s.startswith('dc') or s.startswith('pdc') or
            s.endswith(' dcs') or s=='transitional dc' or s in {'hofbauer cells','kupffer cells','mnp','mono-mac',
            'granulocytes','myelocytes','promyelocytes','neutrophils','neutrophil-myeloid progenitor'}):return 'Myeloid'
        if 'erythroid' in s or s=='erythrocytes':return 'Erythroid'
    elif model=='lung':
        epithelial={'AT0','AT1','AT2','AT2 proliferating','Basal resting','Club (nasal)','Club (non-nasal)',
          'Deuterosomal','Goblet (bronchial)','Goblet (nasal)','Goblet (subsegmental)','Hillock-like','Ionocyte',
          'Multiciliated (nasal)','Multiciliated (non-nasal)','Neuroendocrine','SMG duct','SMG mucous',
          'SMG serous (bronchial)','SMG serous (nasal)','Suprabasal','Tuft','pre-TB secretory'}
        if label in epithelial:return 'Epithelial'
        if label in {'Hematopoietic stem cells','Mesothelium'}:return None
        if 'fibroblast' in s:return 'Fibroblast'
        if s.startswith('ec ') or s.startswith('lymphatic ec'):return 'Endothelial'
        if s.startswith('sm ') or s.startswith('smooth muscle') or s=='pericytes':return 'Mural'
        if 't cells' in s or s=='nk cells':return 'T_NK'
        if s in {'b cells','plasma cells'}:return 'B_Plasma'
        if s=='mast cells':return 'Mast'
        if ('mph' in s or 'macrophage' in s or 'monocyte' in s or s.startswith('dc') or s.endswith(' dcs')):return 'Myeloid'
    elif model=='breast':
        if s.startswith(('cd4-','cd8-')) or s in {'gd','nk','nk-ilcs','nkt','t_prol'}:return 'T_NK'
        if s.startswith(('b_naive','bmem_','plasma_')):return 'B_Plasma'
        if s.startswith(('macro-','mono-')) or s in {'neutrophil','cdc1','cdc2','mdc','mye-prol','pdc'}:return 'Myeloid'
        if s=='mast':return 'Mast'
        if s.startswith('fibro-'):return 'Fibroblast'
        if s.startswith(('lymph-','vas-')):return 'Endothelial'
        if s in {'pericytes','vsmc'}:return 'Mural'
        if s.startswith(('lummhr-','lumsec-')) or s=='basal':return 'Epithelial'
    raise ValueError('Unreviewed frozen reference vocabulary: '+model+' / '+label)


def local_family(major):
    if major in {'T','NK'}:return 'T_NK'
    if major in {'B','Plasma'}:return 'B_Plasma'
    return major if major in FAMILIES else None


def native_family_support(calls):
    return {'T_NK':calls['T']|calls['NK_related'], 'B_Plasma':calls['B']|calls['Plasma'],
            'Myeloid':calls['Myeloid']|calls['DC_related']|calls['Neutrophil_related'],
            'Mast':calls['Mast'],'Endothelial':calls['Endothelial']|calls['Lymphatic_related'],
            'Fibroblast':calls['Fibroblast'],'Mural':calls['Mural'],'Epithelial':calls['Epithelial'],
            'Glial':calls['Glial_Oligodendrocyte_related'],'Astrocyte':calls['Astrocyte_related'],
            'Neuron':calls['Neuron_related'],'Erythroid':calls['Erythroid_related']}


def route_cells(native,reference,local_codes,positive):
    """No majority override of a coherent minority or native contradiction.

Confidence names describe evidence provenance, not calibrated probability.
Native-only calls are low-confidence routing candidates, never final identities.
"""
    native,reference=np.asarray(native,bool),np.asarray(reference,bool)
    local_codes,positive=np.asarray(local_codes),np.asarray(positive,bool)
    if native.shape!=reference.shape or native.shape!=(len(positive),len(FAMILIES)) or len(local_codes)!=len(positive):
        raise ValueError('Coarse routing axes differ')
    n=len(positive);route=np.full(n,-1,np.int16);reason=np.zeros(n,np.uint8)
    native_n=native.sum(axis=1);reference_n=reference.sum(axis=1)
    native_code=np.where(native_n==1,native.argmax(axis=1),-1)
    reference_code=np.where(reference_n==1,reference.argmax(axis=1),-1)
    for i in range(n):
        if not positive[i]:status='zero_native_RNA'
        elif native_n[i]>1:status='native_cross_lineage_review'
        elif reference_n[i]>1:status='reference_cross_lineage_review'
        elif native_n[i]==1:
            if reference_n[i]==1 and native_code[i]!=reference_code[i]:status='native_reference_disagreement'
            else:
                route[i]=native_code[i]
                if reference_n[i]==1:status='native_reference_concordant'
                elif local_codes[i]==native_code[i]:status='native_local_group_concordant_reference_unresolved'
                else:status='native_only_provisional'
        elif reference_n[i]==1 and local_codes[i]==reference_code[i]:
            route[i]=reference_code[i];status='reference_local_group_concordant_native_dropout'
        else:status='insufficient_coarse_evidence'
        reason[i]=STATUS.index(status)
    return route,reason
