"""Complete source-output and fixed-family audit; no individual data read.

Raw numeric/individual invariants are dual-implemented by the producer.
This audit checks the complete exported family and clarifies scope coverage.
"""
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import math

PHASES=['pre','post']
PARENTS=['dual_regulatory_evidence','costate_regulatory_sensitivity']
CONTEXTS=['full_current_source','CM02_network_member_union']
DEPTHS=['all_native_filtered','low_native_capture','mid_native_capture','high_native_capture']
R=['R_TNFRSF9','R_LAG3','R_cycling']
E=['E_dysfunction','E_cycling','E_memory']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def correlation(v):
    assert v is None or math.isfinite(v) and -1<=v<=1


def main():
    parser=argparse.ArgumentParser();parser.add_argument('source',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();root=args.source;out=args.output
    assert not out.exists();out.mkdir(mode=0o700)
    terminal=json.loads((root/'private_terminal.json').read_text())
    assert terminal['status']=='COMPLETE_ALL288_OBSERVED_T_PAIR_CAPTURE_DESCRIPTIONS'
    for name,proof in terminal['files'].items():assert sha(root/name)==proof['SHA256']
    numeric=json.loads((root/'private_numeric_proofs.json').read_text())
    assert all(numeric.values())
    rows=json.loads((root/'private_all288_observed_T_pair_capture_descriptions.json').read_text())
    profiles=json.loads((root/'private_all32_source_scope_and_quality_profiles.json').read_text())
    fields=['phase','parent','context','capture_scope']
    key=lambda r:tuple(r[k] for k in fields)
    expected=set(product(PHASES,PARENTS,CONTEXTS,DEPTHS))
    lookup={key(r):r for r in profiles}
    assert len(profiles)==32 and len(lookup)==32 and set(lookup)==expected
    full={(phase,parent):lookup[phase,parent,CONTEXTS[0],DEPTHS[0]]['source_phase_cases'] for phase,parent in product(PHASES,PARENTS)}
    assert all(full[p,PARENTS[0]]==full[p,PARENTS[1]] for p in PHASES)
    qualified=[]
    for r in profiles:
        authority=full[r['phase'],r['parent']]
        observed=r['source_phase_cases']
        assert 0<=r['eligible_reviewed_cases']<=observed<=authority
        assert sum(r['exclusions'].values())+r['eligible_reviewed_cases']==observed
        assert sum(x['reviewed_cases'] for x in r['eligible_by_arm'])==r['eligible_reviewed_cases']
        for values in r['background_capture_median_profiles'].values():
            if r['eligible_reviewed_cases']<5:assert values is None
            else:assert values==sorted(values) and all(math.isfinite(v) and v>=0 for v in values)
        # Raw source_phase_cases counts observed cases in this scope; preserve
        # that file and add the full authority denominator explicitly.
        qualified.append(dict(r,authority_source_phase_cases=authority,
            context_phase_cases_with_observations=observed,
            source_cases_without_context_observations=authority-observed))
    expected_rows={p+(r,e) for p in expected for r,e in product(R,E)}
    keys=[key(r)+(r['regulatory_state'],r['effector_state']) for r in rows]
    assert len(rows)==288 and len(set(keys))==288 and set(keys)==expected_rows
    for r in rows:
        n=r['reviewed_cases'];assert n==lookup[key(r)]['eligible_reviewed_cases']
        values=[r[k] for k in ['capture_adjusted_partial_Spearman','abundance_treatment_HPV_partial_Spearman','unadjusted_Spearman']]
        values+=list(r['sensitivity_partial_Spearman'].values())
        for value in values:correlation(value)
        if n<20:assert all(v is None for v in values)
        for prefix,point in [('capture',r['capture_adjusted_partial_Spearman']),('base',r['abundance_treatment_HPV_partial_Spearman'])]:
            count=r[prefix+'_valid_bootstraps'];interval=r[prefix+'_pointwise_CI']
            assert 0<=count<=2000
            if count<1800:assert interval is None
            if interval is not None:
                assert point is not None and len(interval)==2 and interval==sorted(interval)
                for value in interval:correlation(value)
            loo=r[prefix+'_valid_leave_one'];assert 0<=loo<=n
            if n-1<20:assert loo==0
            bounds=r[prefix+'_leave_one_range']
            if not loo:assert bounds is None
            else:
                assert bounds==sorted(bounds)
                for value in bounds:correlation(value)
            sign=r[prefix+'_leave_one_same_direction']
            if point is None:assert sign is None
            else:assert 0<=sign<=loo
        assert r['no_P_FDR_or_confirmatory_promotion']
    target=out/'private_all32_qualified_scope_denominators.json'
    target.write_text(json.dumps(qualified,indent=2)+'\n')
    receipt=dict(status='PASS_ALL288_FIXED_OPPORTUNITIES_SOURCE_SCOPE_GATES_AND_UNCERTAINTY_BOUNDS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),auditor_SHA256=sha(Path(__file__)),
        producer_terminal_SHA256=sha(root/'private_terminal.json'),all_source_opportunities=288,
        qualified_scope_denominator_SHA256=sha(target),observed_scope_vs_full_authority_denominator_explicit=True,
        all_patient_count_and_covariate_gates_retained=True,all_small_source_quality_distributions_suppressed=True,
        no_independent_patient_array_or_resampling_recomputation_claim=True,
        producer_all_raw_counts_and_every_statistic_rank_QR_SVD_dual_implementation_verified_by_bound_terminal=True)
    target=out/'private_source_audit_receipt.json';target.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(status=receipt['status'],receipt_SHA256=sha(target))))


if __name__=='__main__':main()
