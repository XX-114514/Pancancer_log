"""Verify complete aggregate scope and real stream coverage after producer exit.

Numerical acceptance does not assert independent patients or causal biology.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import numpy as np


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--previous-task',type=Path,required=True);a=parser.parse_args();d=a.task/'statistics_v1'
    read=lambda n:json.loads((d/n).read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    terminal=read('private_terminal.json');plan=json.loads((a.task/'private_frozen_statistics_v1.json').read_text())
    for name,digest in plan['code_SHA256'].items():assert sha(Path(__file__).parent/name)==digest
    assert terminal['status'] in ['COMPLETE_SUBTYPE_OTHER_CD4_FAMILY_WORKING_NULL_ACCEPTED','COMPLETE_FAMILY_NULL_UNACCEPTED_OR_NO_MAIN_FITS']
    components=read('private_all112_components.json');primary=read('private_primary30.json');fine=read('private_within_subtype_clinical9.json')
    points=read('private_all702_point_diagnostics.json');audits=read('private_all_real_random_stream_R_audits.json');nulls=read('private_all_fixed_working_null_diagnostics.json')
    coverage=read('private_all_fixed_cohort_coverage.json')
    assert [len(v) for v in [components,primary,fine,points]]==[112,30,9,702]
    assert len({r['label'] for r in components})==112
    assert {(r['R_state'],r['E_state']) for r in primary}=={tuple(v) for v in plan['all30_pairs']}
    assert {(r['subtype'],r['marker']) for r in fine}=={(s,g) for s in plan['author_lineage_names']['Treg'] for g in plan['R_markers']}
    assert len({(r['setting'],r['component']) for r in points})==702
    fitted=[r for r in components if r['status']=='ESTIMABLE'];lookup={r['label']:r for r in components}
    audit={r['label']:r for r in audits};assert len(audit)==len(audits)
    expected={'NSCLC_state_vs_subtype_v1|'+r['label'] for r in fitted}
    assert expected<={*audit};assert {r['label'] for r in nulls}=={r['label'] for r in fitted}
    for r in fitted:
        rec=audit['NSCLC_state_vs_subtype_v1|'+r['label']]
        assert rec['status']=='PASS_ALL_SUBMITTED_REAL_RANDOM_STREAMS_R_QR'
        assert rec['full_wild_statistics_checked']==19999 and rec['full_bootstrap_statistics_checked']==2000
        assert sum(r['bootstrap_status_counts'].values())==2000
        assert r['bootstrap_status_counts'].get('ESTIMABLE',0)==r['bootstrap_valid']
        assert (r['CI95'] is not None)==(r['bootstrap_valid']>=1900)
        if r['CI95'] is not None:assert len(r['CI95'])==2 and -1<=r['CI95'][0]<=r['CI95'][1]<=1
        loo=r['leave_one_author_case_out'];assert loo['estimable']+sum(loo['missing_status_counts'].values())==r['n']
        assert loo['same_direction']<=loo['estimable']
    for r in primary:
        label=r['R_state']+'|'+r['E_state'];c=lookup['coordination|'+label];x=lookup['conditional_R|'+label]
        assert r['regulatory_addition_IUT_P']==max(c.get('raw_P',1.),x.get('raw_P',1.))
    for r in fine:
        if r['status']!='ESTIMABLE':assert r['within9_BH']>=0 and r.get('rho') is None
    for r in nulls:
        assert {s['scenario'] for s in r['scenarios']}=={'normal','heteroskedastic_normal'}
        assert all(s['trials']==256 and s['passed_finite_alarm']==(s['rejection_rate']<=.10) for s in r['scenarios'])
    accepted=bool(nulls) and all(s['passed_finite_alarm'] for r in nulls for s in r['scenarios'])
    assert accepted==(terminal['status']=='COMPLETE_SUBTYPE_OTHER_CD4_FAMILY_WORKING_NULL_ACCEPTED')
    # All30 same-case baseline correlations must match the previously closed
    # matched-reference R relation, rather than being attributed to new data.
    old_acceptance=a.previous_task/'private_complete_source_family_acceptance_v1.json'
    measurement=json.loads((a.task/'private_frozen_measurement_v1.json').read_text())
    assert sha(old_acceptance)==measurement['old_family_acceptance_SHA256']
    old_manifest=json.loads(old_acceptance.read_text())
    reference=a.previous_task/'state_family_v1/private_all30_matched_other_CD4_references.json'
    assert sha(reference)==old_manifest['source_files_SHA256'][reference.name]
    previous=json.loads(reference.read_text())
    base={r['component']:r for r in points if r['setting']=='base'};errors=[]
    for r in previous:
        x=base['coordination|'+r['R_state']+'|'+r['E_state']]
        assert x['n']==r['n'] and x['status']==r['R_status']=='ESTIMABLE'
        np.testing.assert_allclose(x['rho'],r['R_rho'],atol=1e-12,rtol=0);errors.append(abs(x['rho']-r['R_rho']))
    within={r['label']:r for r in components if r['kind']=='within_subtype_clinical'}
    assert len(within)==9
    mainn=coverage['main']['retained_for_current_question']
    assert all(r['n']==mainn for r in components if r['kind']!='within_subtype_clinical')
    for r in within.values():assert r['n']==coverage[r['subtype']]['retained_for_current_question']
    positive=next((r['label'] for r in fitted if r['rho']>0),None)
    negative=next((r['label'] for r in fitted if r['rho']<0),None)
    missing=next((r['label'] for r in components if r['status']!='ESTIMABLE'),None)
    result=dict(status='ACCEPT_COMPLETE_NUMERICAL_SOURCE_FAMILY',observed_UTC=datetime.now(timezone.utc).isoformat(),
        scientific_working_null_accepted=accepted,all112_and702_and_fixed_families_complete=True,
        all_real_fitted_random_streams_R_accepted=True,actual_fitted_components=len(fitted),actual_R_audits=len(audits),
        max_R_numeric_error=max(r['max_numeric_error'] for r in audits) if audits else None,
        all30_same_case_baseline_correlations_match_previous_accepted_reference=True,
        max_previous_baseline_error=max(errors),main_author_material_units=mainn,reviewed_independent_patients=0,
        positive_example=positive,negative_example=negative,unestimable_example=missing,
        source_files_SHA256={p.name:sha(p) for p in d.glob('*.json') if p.name!='private_checkpoint.json'},
        private_statistics_contract_SHA256=sha(a.task/'private_frozen_statistics_v1.json'),
        no_direct_LR_space_function_prediction_or_independent_validation_claim=True)
    target=a.task/'private_full_family_postrun_acceptance_v1.json';assert not target.exists()
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','scientific_working_null_accepted','actual_fitted_components','actual_R_audits','max_R_numeric_error','max_previous_baseline_error']}))


if __name__=='__main__':main()
