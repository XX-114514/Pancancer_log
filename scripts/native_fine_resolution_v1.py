"""Positive multi-gene patterns for a resolution audit, not an identity release.

Each pattern is an AND of OR clauses. Unknown measurement stays unknown unless
a measured clause establishes support or absence of the required pattern. These
are diagnostic patterns; lack of this stronger pattern does not disprove an
existing label, and a positive pattern does not diagnose an identity or doublet.
"""
import numpy as np

FINE_ANCHOR_PATTERNS = {
    'Treg': [('FOXP3',), ('IL2RA', 'CTLA4'), ('CD3D', 'CD3E')],
    'CD4_T': [('CD4',), ('CD3D', 'CD3E')],
    'CD8_T': [('CD8A',), ('CD8B',), ('CD3D', 'CD3E')],
    'T_gamma_delta_RNA_candidate': [('TRDC',), ('TRGC1', 'TRGC2'), ('CD3D',)],
    'Neutrophil': [('CSF3R',), ('FCGR3B',), ('S100A8', 'S100A9')],
    'pDC': [('CLEC4C',), ('IL3RA',), ('TCF4', 'GZMB')],
    'cDC1': [('CLEC9A',), ('XCR1',)],
    'cDC2_related': [('CD1C',), ('FCER1A',)],
    'Macrophage': [('C1QA',), ('C1QB',), ('APOE', 'CSF1R')],
    'Monocyte': [('FCN1',), ('LYZ',), ('S100A8', 'S100A9')],
    'Lymphatic_endothelial_related': [('PROX1', 'FLT4'), ('PDPN', 'CCL21')],
}
ANCHOR_FAMILIES = {t: ('T_NK' if t in list(FINE_ANCHOR_PATTERNS)[:4] else
                     'Endothelial' if t == 'Lymphatic_endothelial_related' else
                     'Myeloid') for t in FINE_ANCHOR_PATTERNS}


def evaluate_patterns(packed, genes, measured, patterns=FINE_ANCHOR_PATTERNS):
    """Return positive/unknown bits; zero means the required pattern is absent.

    Zero is not a biological-negative identity. Genes missing from the axis or
    flagged unmeasured cannot supply evidence of a required gene's absence.
    """
    x=np.asarray(packed);m=np.asarray(measured)
    if (x.ndim!=2 or x.dtype!=np.uint8 or m.dtype!=bool or m.shape!=(len(genes),)
            or x.shape[1]!=(len(genes)+7)//8 or len(set(genes))!=len(genes)):
        raise ValueError('Packed detection, unique gene axis and measurement mask differ')
    if len(patterns)>16:raise ValueError('Diagnostic pattern bit capacity exceeded')
    lookup={g:i for i,g in enumerate(genes)};n=len(x)
    pos=np.zeros(n,np.uint16);unknown=np.zeros(n,np.uint16);coverage={}
    for j,(name,clauses) in enumerate(patterns.items()):
        if not clauses or any(not c or len(set(c))!=len(c) for c in clauses):
            raise ValueError('Every pattern needs nonempty unique-gene OR clauses')
        states=[];description=[]
        for clause in clauses:
            present=[g for g in clause if g in lookup and m[lookup[g]]]
            detected=np.zeros(n,bool)
            for g in present:
                k=lookup[g];detected|=((x[:,k//8]>>(k%8))&1).astype(bool)
            fully_measured=len(present)==len(clause)
            state=np.where(detected,1,0 if fully_measured else 2).astype(np.uint8)
            states.append(state);description.append({'requested':list(clause),'measured':present,'fully_measured':fully_measured})
        truth=np.stack(states)
        positive=(truth==1).all(axis=0)
        unresolved=(truth!=0).all(axis=0)&(truth==2).any(axis=0)
        assert not (positive&unresolved).any()
        pos|=positive.astype(np.uint16)<<j;unknown|=unresolved.astype(np.uint16)<<j
        coverage[name]={'AND_of_OR_clauses':description,'all_requested_genes_measured':all(c['fully_measured'] for c in description)}
    return {'positive_anchor_bitmask':pos,'unassessable_anchor_bitmask':unknown,'coverage':coverage,'pattern_order':list(patterns)}
