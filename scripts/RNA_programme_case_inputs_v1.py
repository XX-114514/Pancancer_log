"""RAM-only reference/score aggregation with all-current-T denominators.

Author reference absence never becomes biological absence. Overlap removal,
source-key authentication and fixed parent labels precede this pure helper.
Nothing in this module writes arrays, cell rows, or patient identifiers.
"""
import numpy as np
import pandas as pd

CONTEXTS=['current_T_author_reference','CM02_network_member_union_author_reference']
PARENTS=['author_regulatory_reference','author_conventional_CD4_reference','author_CD8_reference']
REGULATORY=['4-1BB+ Treg','HLA-DR+ Treg','LAG3+ Treg','Prolif Treg']
CONVENTIONAL=['CXCL13+ CD4T','IFN-stim CD4T','Tn/Tmem']
CD8=['CXCL13+Tex','Tdys','GZMK+T','Tmem','Tprolif','Tstr']


def parent(label):
    if label in REGULATORY:return PARENTS[0]
    if label in CONVENTIONAL:return PARENTS[1]
    if label in CD8:return PARENTS[2]
    raise ValueError('UNDECLARED_AUTHOR_REFERENCE_LABEL')


def aggregate(tab,cur,unique,nodes,scores,total,detected):
    registry=tab.set_index('library_id');frame=unique.copy()
    assert registry.index.is_unique and frame.position.is_unique
    assert frame.major.eq('T').all() and frame.position.is_monotonic_increasing
    assert len(frame)==len(scores)==len(total)==len(detected) and np.asarray(scores).shape==(len(frame),2)
    assert np.isfinite(scores).all() and (np.asarray(total)>=0).all() and (np.asarray(detected)>=0).all()
    assert frame.library_id.isin(registry.index).all() and cur.library_id.isin(registry.index).all()
    assert frame.library_id.map(registry.timepoint_registry).eq(frame.timepoint).all()
    assert not tab.duplicated(['patient_registry','timepoint_registry']).any()
    frame['parent']=frame.label.map(parent)
    frame['s.genes']=np.asarray(scores)[:,0];frame['g2m.genes']=np.asarray(scores)[:,1]
    frame['mean_log1p_UMI']=np.log1p(total);frame['mean_log1p_detected']=np.log1p(detected)
    output={};coverage=[]
    for context in CONTEXTS:
        current=cur.loc[cur.major.eq('T')];selected=frame
        if context==CONTEXTS[1]:
            current=current.loc[current.node_id.isin(nodes)];selected=frame.loc[frame.node_id.isin(nodes)]
        denominators=current.groupby('library_id').size()
        values=selected.groupby(['library_id','parent']).agg(
            n=('position','size'),S=('s.genes','mean'),G2M=('g2m.genes','mean'),
            mean_log1p_UMI=('mean_log1p_UMI','mean'),mean_log1p_detected=('mean_log1p_detected','mean'))
        for (library,role),r in values.iterrows():
            meta=registry.loc[library];denom=int(denominators.get(library,0));assert denom>=int(r.n)>0
            key=(context,meta.timepoint_registry,meta.patient_registry,role)
            assert key not in output
            output[key]=dict(n=int(r.n),T_total=denom,fraction=float(r.n/denom),
                S=float(r.S),G2M=float(r.G2M),mean_log1p_UMI=float(r.mean_log1p_UMI),
                mean_log1p_detected=float(r.mean_log1p_detected),arm=meta.neo_type,HPV=meta['HPV/p16'])
        for phase in ['pre','post','post2']:
            phase_tab=tab.loc[tab.timepoint_registry.eq(phase)]
            for left in PARENTS[:2]:
                reasons={};eligible=0
                for case in phase_tab.patient_registry:
                    l=output.get((context,phase,case,left));r=output.get((context,phase,case,PARENTS[2]))
                    reason='PARENT_REFERENCE_BELOW20_OR_UNCAPTURED'
                    if l is not None and r is not None and min(l['n'],r['n'])>=20:
                        assert l['arm']==r['arm'] and l['HPV']==r['HPV']
                        reason='UNKNOWN_ARM_OR_HPV'
                        if l['arm'] in ['D','D+T'] and l['HPV'] in ['positive','negative']:
                            reason='ELIGIBLE_REFERENCE_AND_BACKGROUND';eligible+=1
                    reasons[reason]=reasons.get(reason,0)+1
                assert sum(reasons.values())==len(phase_tab)
                coverage.append(dict(context=context,phase=phase,left_parent=left,right_parent=PARENTS[2],
                    reviewed_cases=len(phase_tab),eligible_cases=eligible,reasons=reasons,
                    missing_reference_not_biological_absence=True))
    return output,coverage


def pair_arrays(units,slot):
    context=slot['context'];phase=slot['phase'];left=slot['left_parent'];right=slot['right_parent']
    cases=sorted({case for co,ph,case,role in units if co==context and ph==phase})
    cases=[case for case in cases if (context,phase,case,left) in units and (context,phase,case,right) in units]
    pairs=[];accepted=[]
    for case in cases:
        l=units[context,phase,case,left];r=units[context,phase,case,right]
        assert l['arm']==r['arm'] and l['HPV']==r['HPV']
        if min(l['n'],r['n'])>=20 and l['arm'] in ['D','D+T'] and l['HPV'] in ['positive','negative']:
            pairs.append((l,r));accepted.append(case)
    field={'s.genes':'S','g2m.genes':'G2M'}
    x=np.array([l[field[slot['left_programme']]] for l,r in pairs])
    y=np.array([r[field[slot['right_programme']]] for l,r in pairs])
    controls=np.array([[l['fraction'],r['fraction'],np.log1p(l['n']),np.log1p(r['n']),
        l['mean_log1p_UMI'],r['mean_log1p_UMI'],l['mean_log1p_detected'],r['mean_log1p_detected']] for l,r in pairs],float).reshape(-1,8)
    categories=np.array([[l['arm'],l['HPV']] for l,r in pairs],str).reshape(-1,2)
    return x,y,controls,categories,accepted
