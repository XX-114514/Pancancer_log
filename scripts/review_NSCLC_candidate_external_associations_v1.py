"""Complete accepted external family summaries and one full primary forest.

Private aggregate tables/figure only. No new tests, scoring or endpoint changes.
"""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,8*1024**2)
    src=root/'NSCLC_candidate_external_association_results_20261018_v1'
    proof=src/'private_independent_external_candidate_result_audit.json'
    assert sha(proof)=='dc9e287dea77c2044f03f4b7087c3cc1462436f2d913088315e33736dd0a818b'
    audit=json.loads(proof.read_text());assert audit['producer_terminal_sha256']==sha(src/'private_terminal.json')
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    path=root/'NSCLC_candidate_external_association_scope_20261018_v1.json'
    assert sha(path)=='d4141c94fe11d2fbccd764fbc2bff503d8a068e7e32132cb6f5c790bb158d752'
    scope=json.loads(path.read_text());metadata=pd.DataFrame(scope['candidates'])
    effects=pd.read_csv(src/'private_all_estimable_external_effects_and_sensitivities.tsv.gz',sep='\t')
    eligibility=pd.read_csv(src/'private_complete_planned_eligibility.tsv.gz',sep='\t')
    assert len(effects)==audit['all_estimable_models'] and len(eligibility)==audit['all_planned_positions']
    out=root/'NSCLC_candidate_external_association_review_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),audit_sha256=sha(proof),budget=budget,
        scope='All planned statuses and all primary/reviewed-family effects, no favorable-source selection.',
        same_cases='Match exact active_case_rows_json between families; same underlying units do not count as independent replication.',
        plot='Every estimable all-unit primary row in original test_index order; hollow markers identify unavailable CIs, not zero effects. All source/cancer/material/phase labels retained.',
        no_new_P_or_primary_changes=True)
    (out/'plan_before_complete_external_result_review.json').write_text(json.dumps(plan,indent=2)+'\n')
    data=effects.merge(metadata,on='record',how='left',validate='many_to_one')
    assert data.ligand.notna().all()
    data['CI_available']=data[['CI_low','CI_high']].notna().all(axis=1)
    data['positive_BH_MC_and_CI']=data.partial_rank_rho.gt(0)&data.BH_q.lt(.05)&data.MC_upper_BH_q.lt(.05)&data.CI_low.gt(0)
    data['same_case_set_as_other_family']=False
    for _,g in data.groupby(['record','context_id'],sort=False):
        if len(g)==2:
            same=set(json.loads(g.active_case_rows_json.iloc[0]))==set(json.loads(g.active_case_rows_json.iloc[1]))
            data.loc[g.index,'same_case_set_as_other_family']=same
    rows=[]
    for c in metadata.itertuples():
        for family in scope['families']:
            g=data[data.record.eq(c.record)&data.family.eq(family)]
            rows.append(dict(record=c.record,ligand=c.ligand,sender=c.sender_view,receiver=c.receiver_view,family=family,
                planned_contexts=int((eligibility.record.eq(c.record)&eligibility.family.eq(family)).sum()),estimated_contexts=len(g),
                distinct_estimated_GSE=g.dataset.nunique(),positive_effects=int(g.partial_rank_rho.gt(0).sum()),negative_effects=int(g.partial_rank_rho.lt(0).sum()),
                CI_unavailable=int((~g.CI_available).sum()),positive_BH_MC_and_CI=int(g.positive_BH_MC_and_CI.sum()),
                BY_below0p05=int(g.BY_q.lt(.05).sum()),all_LOO_valid_same_sign=int((g.LOO_all_valid&g.LOO_same_sign).sum()),
                same_case_set_as_other_family_rows=int(g.same_case_set_as_other_family.sum())))
    data.to_csv(out/'private_all_accepted_effects_with_endpoint_and_case_dependence.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(rows).to_csv(out/'private_all_candidate_family_summary.tsv',sep='\t',index=False)
    failures=eligibility.assign(reason=eligibility.status.str.split(';')).explode('reason').groupby(['record','family','reason'],dropna=False).size().rename('records').reset_index()
    failures.to_csv(out/'private_all_planned_failure_reasons.tsv',sep='\t',index=False)
    g=data[data.family.eq('all_units')].sort_values('test_index');y=np.arange(len(g))
    fig,ax=plt.subplots(figsize=(22,max(5,.53*len(g)+3)))
    labels=[];mapping=[]
    for pos,row in enumerate(g.itertuples()):
        sender=row.sender_view.replace('author_','');receiver=row.receiver_view.replace('author_','')
        label=f'{sender} -> {receiver} | {row.ligand} | {row.dataset} | {row.cancer_label}\n{row.material_class}; {row.phase_key}; units={row.active_units}, reviewed={row.active_reviewed_units}'
        labels.append(label)
        if row.CI_available:
            ax.plot([row.CI_low,row.CI_high],[pos,pos],color='#466c91',lw=1.5)
            ax.plot(row.partial_rank_rho,pos,'o',color='#245b87',ms=6)
        else:ax.plot(row.partial_rank_rho,pos,'o',markerfacecolor='white',markeredgecolor='#245b87',ms=6)
        ax.text(1.03,pos,f'BH={row.BH_q:.4g}; BY={row.BY_q:.4g}',fontsize=8,va='center')
        mapping.append(dict(model=row.model,point=row.partial_rank_rho,CI_low=None if not row.CI_available else row.CI_low,CI_high=None if not row.CI_available else row.CI_high))
    ax.axvline(0,color='#666',lw=.9);ax.set_xlim(-1,1);ax.set_ylim(len(g)-.4,-.6)
    ax.set_yticks(y,labels,fontsize=8);ax.set_xlabel('Prior-oriented partial rank correlation and fixed case-bootstrap 95% CI')
    ax.set_title('All estimable primary external models: full frozen all-unit family\nHollow point = unavailable CI; reviewed-only repeated units are not independent replications',fontsize=12,pad=18)
    fig.subplots_adjust(left=.47,right=.82,bottom=.11,top=.89)
    fig.text(.02,.025,'Frozen subtype-pool assay extension; not direct binding, active ligand, clinical prediction or pan-cancer validation. All planned missing models remain in the full eligibility table.',fontsize=8)
    path=out/'private_all_primary_external_models_forest.png';fig.savefig(path,dpi=150);plt.close(fig)
    (out/'private_exact_primary_forest_mapping.json').write_text(json.dumps(mapping,indent=2,allow_nan=False)+'\n')
    files=[p for p in out.iterdir() if p.is_file()];assert sum(p.stat().st_size for p in files)<8*1024**2
    terminal=dict(status='COMPLETE_ALL_EXTERNAL_RESULT_SUMMARIES_AND_FULL_PRIMARY_FOREST_PIXEL_REVIEW_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),accepted_result_audit_sha256=sha(proof),
        primary_rows_plotted=len(g),all_estimated_rows_retained=len(data),all_planned_rows_accounted=len(eligibility),
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},new_tests_or_independent_validation=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],completed_UTC=terminal['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
