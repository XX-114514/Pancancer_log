"""Post-result all major backgrounds, including sender/receiver-disjoint controls.

Recompute baseline on exactly the same eligible case subset for each background.
No new P/CI, sender replacement or candidate rescue; all absent/failed designs kept.
"""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_external_lineage_inputs_v2 import load
from NSCLC_native_MFA_interpretation_inputs_v1 import sha,VIEWS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from test_NSCLC_candidate_external_associations_v1 import inputs
from prepare_NSCLC_candidate_external_associations_v1 import QCOLS
import NSCLC_candidate_external_association_math_v1 as math
from audit_NSCLC_external_programme_coordination_v3 import calculate
from measure_NSCLC_supported_external_source_identity_v1 import OUT


def checked(x,y,q,b):
    p=math.prepare(x,y,q,b);r=math.base.effect(p)
    oracle=calculate(x,y,q,b);assert np.array_equal(p['keep'],oracle['keep'])
    assert abs(r-oracle['rho'])<1e-10
    return r,len(p['a']),p['design_rank'],p['residual_df']


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,8*1024**2)
    source=root/OUT;proof=source/'private_independent_full_selected_cell_identity_audit.json'
    assert sha(proof)=='b6069f69ece39ee092df7e82b9927fe5c45f5c71b37d245fbb0b1a358732675e'
    terminal=json.loads((source/'private_terminal.json').read_text())
    name='private_all_library_node_source_lineage_counts.tsv.gz'
    assert sha(source/name)==terminal['files'][name]['sha256']
    plan=json.loads((source/'plan_before_source_positive_cell_expression.json').read_text())
    table=pd.read_csv(source/name,sep='\t');data=load()
    _,_,contract,_,cases,_=inputs(root)
    models=pd.read_csv(root/'NSCLC_candidate_external_association_results_20261018_v1/private_all_estimable_external_effects_and_sensitivities.tsv.gz',sep='\t').set_index('model')
    out=root/'NSCLC_supported_external_other_lineage_background_20261018_v2'
    assert not out.exists();out.mkdir(mode=0o700)
    save=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
              identity_audit_sha256=sha(proof),source_table_sha256=sha(source/name),
              scope='All supported original main leads, ALL captured frozen major labels in BOTH full-major and sender/receiver-pool-disjoint scope, one background at a time. Post-result descriptive extension after v1 total-T part-whole reversal; no tests or main changes.',
              eligibility='Same original paired library set, background present in EVERY paired library and at least20 background cells across the case; absent is unmeasurable, never zero ligand.',
              features='Equal original paired-library mean of background per-cell logCP10K; log1p total background case cells; equal-library mean background fraction among all canonical annotated cells. All three appended to original12 quality/fraction controls plus full blocks.',
              comparison='Recompute baseline on exact same background-eligible cases and same non-singleton blocks. At least20 units/df8 and original full-span gates; retain failures. Baseline and adjusted rho additionally checked by independent least-squares projection.',
              caveat='Major labels containing sender or receiver are explicitly marked part-whole/overlapping, not independent third-lineage controls. Background attenuation can represent real shared ecology; it does not prove ambient contamination or causality.')
    (out/'plan_before_all_other_lineage_background_effects.json').write_text(json.dumps(save,indent=2)+'\n')
    candidate_lookup=pd.DataFrame(contract['candidates']).set_index('record')
    result=[];missing=[]
    for lead in plan['selected_records']:
        original=models.loc[lead['model']];allcases=cases.iloc[json.loads(lead['active_case_rows_json'])]
        g=table[table.model.eq(lead['model'])]
        for background_scope in ['full_major','sender_receiver_pool_disjoint']:
            candidate=candidate_lookup.loc[lead['record']]
            excluded=[VIEWS.index(candidate.sender_view),VIEWS.index(candidate.receiver_view)]
            for major in sorted(g.major.unique()):
                part=g[g.major.eq(major)]
                if background_scope=='sender_receiver_pool_disjoint':part=part[~part.frozen_role.isin(excluded)]
                grouped=part.groupby('library_index')[['n_cells','source_sum_logCP10K']].sum()
                rows=[];extra=[]
                for case in allcases.itertuples():
                    li=np.array(json.loads(case.library_indices_json),int)
                    bg=grouped.reindex(li);present=bg.n_cells.notna().all() and bg.n_cells.gt(0).all()
                    total=float(bg.n_cells.sum()) if present else np.nan
                    ok=bool(present and total>=20)
                    missing.append(dict(model=lead['model'],background_scope=background_scope,major=major,case_row=case.Index,all_original_libraries_present=bool(present),background_case_cells=total,eligible=ok))
                    if not ok:continue
                    rows.append(case.Index)
                    extra.append([float((bg.source_sum_logCP10K/bg.n_cells).mean()),float(np.log1p(total)),float((bg.n_cells.to_numpy()/data['lib'].loc[li,'annotated_scope_cells'].to_numpy(float)).mean())])
                use=cases.loc[rows];effect=adjusted=np.nan;n=rank=df=n2=rank2=df2=0
                try:
                    x,y=use.ligand.to_numpy(float),use.response.to_numpy(float);q=use[QCOLS].to_numpy(float);b=use.block_key.to_numpy(str)
                    effect,n,rank,df=checked(x,y,q,b)
                    adjusted,n2,rank2,df2=checked(x,y,np.column_stack([q,np.array(extra,float)]),b)
                    assert n==n2;status='ESTIMABLE'
                except ValueError as e:status=str(e)
                result.append(dict(model=lead['model'],record=lead['record'],dataset=lead['dataset'],background_scope=background_scope,major=major,
                    overlaps_original_four_pools=bool(part.frozen_role.ge(0).any()),original_units=len(allcases),background_eligible_cases=len(rows),baseline_active_units=n,
                    baseline_same_subset_rho=effect,background_adjusted_rho=adjusted,baseline_rank=rank,baseline_df=df,adjusted_rank=rank2,adjusted_df=df2,status=status,
                    original_full_case_rho=original.partial_rank_rho,same_sign_as_baseline=bool(np.isfinite(adjusted) and np.sign(effect)==np.sign(adjusted))))
    combined=pd.DataFrame(result)
    previous=root/'NSCLC_supported_external_other_lineage_background_20261018_v1'
    pt=json.loads((previous/'private_terminal.json').read_text())
    pn='private_all_major_background_effect_diagnostics.tsv'
    assert sha(previous/pn)==pt['files'][pn]['sha256']
    prior=pd.read_csv(previous/pn,sep='\t')
    full=combined[combined.background_scope.eq('full_major')].reset_index(drop=True)
    assert full.major.tolist()==prior.major.tolist() and full.status.tolist()==prior.status.tolist()
    for c in prior.select_dtypes(include='number').columns:assert np.allclose(full[c],prior[c],atol=1e-10,rtol=1e-10,equal_nan=True)
    combined.to_csv(out/'private_all_major_background_effect_diagnostics.tsv',sep='\t',index=False)
    pd.DataFrame(missing).to_csv(out/'private_all_case_background_coverage.tsv.gz',sep='\t',index=False,compression='gzip')
    files=[p for p in out.iterdir() if p.is_file()]
    t=dict(status='COMPLETE_ALL_MAJOR_SAME_CASE_BACKGROUND_NO_NEW_P',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),models=len(result),all_accepted_fits_OLS_crosschecked=True,
           files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},biological_validation_or_primary_change=False)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n')
    print(json.dumps(dict(status=t['status'],completed_UTC=t['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))

if __name__=='__main__':main()
