"""Finite descriptive follow-up; individual values stay in caller RAM.

LOO never deletes observations from the closed primary analyses. Half-signature
positions follow the already fixed external list, with no data-driven genes.
Fine-label summaries suppress expression distributions below five qualified
cases and exclude post2 from expression descriptions. No new P or CI is made.
"""
from collections import Counter
import numpy as np
import pandas as pd
import RNA_programme_partial_rank_v1 as stats
import audit_RNA_programme_partial_rank_v1 as audit
import RNA_programme_case_inputs_v1 as inputs


def split_programmes(programmes):
    result={k:list(v) for k,v in programmes.items()}
    genes=result['g2m.genes'];assert len(set(genes))==len(genes)
    assert not {'MKI67','TOP2A','UBE2C'}.intersection(genes)
    result['g2m.odd_positions']=genes[0::2];result['g2m.even_positions']=genes[1::2]
    assert set(result['g2m.odd_positions']).isdisjoint(result['g2m.even_positions'])
    assert set(result['g2m.odd_positions'])|set(result['g2m.even_positions'])==set(genes)
    return result


def describe_fit(x,y,c,cat,levels):
    if len(x)<20:return dict(status='NOT_ESTIMABLE',reason='INSUFFICIENT_REVIEWED_CASES')
    try:fit=stats.fit(x,y,c,cat,levels)
    except ValueError as error:
        try:audit.reference_fit(x,y,c,cat,levels)
        except ValueError as independently:
            assert str(error)==str(independently)
            return dict(status='NOT_ESTIMABLE',reason=str(error))
        raise AssertionError('independent_gate_disagreement')
    reference=audit.reference_fit(x,y,c,cat,levels)
    np.testing.assert_allclose(fit['rho'],reference[0],atol=1e-10,rtol=1e-10)
    return dict(status='DESCRIPTIVE_ESTIMABLE',rho=fit['rho'],unadjusted_rho=fit['unadjusted'],
                condition=fit['condition'],residual_df=fit['residual_df'])


def leave_one_case(x,y,c,cat):
    levels=stats.categorical_spec(cat);effects=[];failure=Counter()
    for omit in range(len(x)):
        keep=np.arange(len(x))!=omit
        result=describe_fit(x[keep],y[keep],c[keep],cat[keep],levels)
        if result['status']=='DESCRIPTIVE_ESTIMABLE':effects.append(result['rho'])
        else:failure[result['reason']]+=1
    return dict(total_case_omissions=len(x),estimable_omissions=len(effects),
        positive=sum(v>0 for v in effects),negative=sum(v<0 for v in effects),zero=sum(v==0 for v in effects),
        rho_min=float(min(effects)) if effects else None,rho_max=float(max(effects)) if effects else None,
        nonestimable_reasons=dict(failure),all_omissions_independently_replayed=True,
        no_P_CI_or_case_deletion_from_primary=True)


def fine_profiles(tab,cur,frame,nodes,scores):
    """All fixed author labels; equal-case summaries after >=20 fine cells."""
    assert scores.shape==(len(frame),4) and frame.position.is_unique and frame.major.eq('T').all()
    registry=tab.set_index('library_id');f=frame.copy()
    assert not tab.duplicated(['patient_registry','timepoint_registry']).any()
    assert f.library_id.map(registry.timepoint_registry).eq(f.timepoint).all()
    assert set(f.label)<=set(inputs.REGULATORY+inputs.CONVENTIONAL+inputs.CD8)
    f['case']=f.library_id.map(registry.patient_registry)
    for k,name in enumerate(['S','G2M','G2M_odd','G2M_even']):f[name]=scores[:,k]
    output=[]
    for context in inputs.CONTEXTS:
        selected=f if context==inputs.CONTEXTS[0] else f.loc[f.node_id.isin(nodes)]
        current=cur.loc[cur.major.eq('T')]
        if context==inputs.CONTEXTS[1]:current=current.loc[current.node_id.isin(nodes)]
        denominator=current.groupby('library_id').size()
        for phase in ['pre','post']:
            phase_tab=tab.loc[tab.timepoint_registry.eq(phase)]
            for label in inputs.REGULATORY+inputs.CONVENTIONAL+inputs.CD8:
                chosen=selected.loc[selected.timepoint.eq(phase)&selected.label.eq(label)]
                rows=chosen.groupby(['library_id','case']).agg(n=('position','size'),S=('S','mean'),G2M=('G2M','mean'),G2M_odd=('G2M_odd','mean'),G2M_even=('G2M_even','mean'))
                assert rows.index.get_level_values('case').is_unique
                qualified=rows.loc[rows.n>=20]
                summary={}
                if len(qualified)>=5:
                    for programme in ['S','G2M','G2M_odd','G2M_even']:
                        summary[programme]=dict(case_mean_min_q25_median_q75_max=np.quantile(qualified[programme],[0,.25,.5,.75,1]).tolist())
                fractions=[float(r.n/denominator.loc[library]) for (library,case),r in qualified.iterrows()]
                output.append(dict(context=context,phase=phase,author_fine_label=label,reference_parent=inputs.parent(label),
                    reviewed_phase_cases=len(phase_tab),cases_with_unique_reference_cells=len(rows),cases_atleast20_fine_cells=len(qualified),
                    unique_reference_cell_observations=len(chosen),expression_summaries=summary,
                    fine_fraction_all_currentT_case_quantiles=np.quantile(fractions,[0,.25,.5,.75,1]).tolist() if len(qualified)>=5 else None,
                    distribution_status='DESCRIPTIVE_EQUAL_CASE_MEANS' if len(qualified)>=5 else 'BELOW_FIVE_QUALIFIED_CASES_SCORE_SUPPRESSED',
                    missing_reference_not_lineage_absence=True))
    assert len(output)==52
    return output
