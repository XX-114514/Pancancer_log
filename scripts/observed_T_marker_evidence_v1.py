"""Frozen molecule-positive evidence rules, without assay-specific QC/geometry.

Eligibility is provided by the caller. These are RNA observation patterns,
not formal cell identities or biological functions.
"""
import numpy as np

GENES=['CD3D','CD3E','TRAC','CD4','CD8A','FOXP3','IL2RA','CTLA4','TNFRSF9','LAG3','PDCD1','TOX','HAVCR2','MKI67','TOP2A','UBE2C','TCF7','IL7R','MS4A1','CD79A','EPCAM','SOX2']
METRICS=['T_support','regulatory_only','CD8_only','dual_ambiguous','neither_parent',
         'R_TNFRSF9','R_LAG3','R_cycling','E_dysfunction','E_cycling','E_memory',
         'B_cycling','epithelial_cycling']


def vector_evidence(counts,eligible,parent):
    assert parent in ['dual_regulatory_evidence','costate_regulatory_sensitivity']
    x={g:counts[:,i] for i,g in enumerate(GENES)}
    t=(x['CD3D']>0)|(x['CD3E']>0)|(x['TRAC']>0)
    r=t&(x['FOXP3']>0)&((x['IL2RA']>0) if parent=='dual_regulatory_evidence' else ((x['IL2RA']>0)|(x['CTLA4']>0)))&eligible
    e=t&(x['CD8A']>0)&eligible;overlap=r&e;r&=~overlap;e&=~overlap
    cycle=(x['MKI67']>0)&((x['TOP2A']>0)|(x['UBE2C']>0))
    return np.column_stack([t&eligible,r,e,overlap,eligible&~r&~e&~overlap,
        r&(x['TNFRSF9']>0),r&(x['LAG3']>0),r&cycle,
        e&(x['PDCD1']>0)&((x['TOX']>0)|(x['HAVCR2']>0)),e&cycle,
        e&((x['TCF7']>0)|(x['IL7R']>0)),
        (x['MS4A1']>0)&(x['CD79A']>0)&eligible&~t&cycle,
        (x['EPCAM']>0)&eligible&~t&cycle])


def scalar_evidence(counts,eligible,parent):
    """Independent row-by-row Boolean truth table on all measured observations."""
    rows=[]
    for values,okay in zip(counts,eligible):
        x=dict(zip(GENES,values))
        t=any(x[g]>=1 for g in ['CD3D','CD3E','TRAC'])
        r=bool(okay and t and x['FOXP3']>=1 and
            (x['IL2RA']>=1 if parent=='dual_regulatory_evidence' else
             any(x[g]>=1 for g in ['IL2RA','CTLA4'])))
        e=bool(okay and t and x['CD8A']>=1)
        both=r and e
        if both:r=False;e=False
        cycle=x['MKI67']>=1 and any(x[g]>=1 for g in ['TOP2A','UBE2C'])
        rows.append([okay and t,r,e,both,okay and not(r or e or both),
            r and x['TNFRSF9']>=1,r and x['LAG3']>=1,r and cycle,
            e and x['PDCD1']>=1 and any(x[g]>=1 for g in ['TOX','HAVCR2']),
            e and cycle,e and any(x[g]>=1 for g in ['TCF7','IL7R']),
            okay and not t and x['MS4A1']>=1 and x['CD79A']>=1 and cycle,
            okay and not t and x['EPCAM']>=1 and cycle])
    return np.asarray(rows,bool)
