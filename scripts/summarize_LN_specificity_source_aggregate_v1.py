"""Interpret audited source aggregates only; no patient/cell/graph arrays."""
from pathlib import Path
import hashlib
import json
import math
import os

ENDPOINTS = ('E0_PD1_TIM3', 'E1_CD8_MKI67', 'E2_CCR7_CD127')
DIRECT = ('noncycling', 'cycling', 'combined_MH', 'cycling_minus_noncycling')
NULL = 'local_composition_degree_and_cycle_conditional_null'


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def independent_adjustment(p, by=False):
    q = [v if finite(v) else 1.0 for v in p]
    order = sorted(range(len(q)), key=q.__getitem__)
    scale = sum(1 / k for k in range(1, len(q)+1)) if by else 1.0
    adjusted = [None] * len(q)
    ceiling = 1.0
    for rank in range(len(q), 0, -1):
        index = order[rank-1]
        ceiling = min(ceiling, q[index] * len(q) / rank * scale)
        adjusted[index] = ceiling if finite(p[index]) else None
    return adjusted


def validate_family(production):
    rows = production['primary_rows']
    expected = {(e, est) for e in ENDPOINTS for est in (*DIRECT, NULL)}
    assert len(rows) == 15
    keys = [(r['endpoint'], r['estimand']) for r in rows]
    assert len(set(keys)) == 15 and set(keys) == expected
    assert all(r.get('grid') == 4 for r in rows if r['estimand'] == NULL)
    p = [r['primary_p'] for r in rows]
    for field, by in [('BH', False), ('BY', True)]:
        reference = independent_adjustment(p, by)
        for a, row in zip(reference, rows):
            b = row[field]
            assert finite(a) == finite(b)
            if finite(a):
                assert abs(a-b) < 1e-11 and 0 <= b <= 1
    for row in rows:
        assert row['cases'] >= 0 and row['eligible_images'] >= 0
        if finite(row['primary_p']):
            assert row['cases'] >= 5 and 0 < row['primary_p'] <= 1
            assert finite(row['median_effect']) and finite(row['bootstrap_lower']) and finite(row['bootstrap_upper'])
        if finite(row['bootstrap_lower']) and finite(row['bootstrap_upper']):
            assert row['bootstrap_lower'] <= row['bootstrap_upper']
    return {(r['endpoint'], r['estimand']): r for r in rows}


def positive_support(row, correction='BH'):
    return (row['status'] == 'ESTIMABLE' and finite(row['median_effect']) and row['median_effect'] > 0
            and finite(row['bootstrap_lower']) and row['bootstrap_lower'] > 0
            and finite(row[correction]) and row[correction] < .05)


def negative_support(row, correction='BH'):
    return (row['status'] == 'ESTIMABLE' and finite(row['median_effect']) and row['median_effect'] < 0
            and finite(row['bootstrap_upper']) and row['bootstrap_upper'] < 0
            and finite(row[correction]) and row[correction] < .05)


def supported_direction(row, correction='BH'):
    if positive_support(row, correction):
        return 'positive'
    if negative_support(row, correction):
        return 'negative'
    return 'not_established' if estimable(row) else 'not_estimable'


def estimable(row):
    return row['status'] == 'ESTIMABLE' and finite(row['primary_p'])


def interpret_endpoint(records, endpoint):
    main = records[endpoint, 'combined_MH']
    local = records[endpoint, NULL]
    contrast = records[endpoint, 'cycling_minus_noncycling']
    main_yes, local_yes = positive_support(main), positive_support(local)
    main_negative, local_negative = negative_support(main), negative_support(local)
    if not estimable(main):
        label = 'WITHIN_MKI67_COMPARISON_NOT_ESTIMABLE'
    elif not estimable(local):
        label = 'LOCAL_NULL_NOT_ESTIMABLE_WITHIN_MKI67_SUPPORT_REPORTED_SEPARATELY'
    elif main_yes and local_yes:
        label = 'SOURCE_SPECIFIC_STATE_PAIR_RETAINED_UNDER_LOCAL_COMPOSITION_NULL'
    elif main_negative and local_negative:
        label = 'SOURCE_SPECIFIC_NEGATIVE_STATE_ASSOCIATION_RETAINED_UNDER_LOCAL_COMPOSITION_NULL'
    elif (main_yes and local_negative) or (main_negative and local_yes):
        label = 'SUPPORTED_DIRECTIONS_DIFFER_ACROSS_ESTIMANDS_AND_POTENTIALLY_DIFFERENT_ELIGIBLE_SETS'
    elif main_yes:
        label = 'WITHIN_MKI67_SUPPORT_LOCAL_CONDITIONAL_TEST_NOT_SUPPORTED'
    elif local_yes:
        label = 'LOCAL_CONDITIONAL_SUPPORT_WITHOUT_MAIN_COMBINED_SUPPORT'
    elif main_negative:
        label = 'WITHIN_MKI67_NEGATIVE_ASSOCIATION_WITHOUT_LOCAL_CONDITIONAL_SUPPORT'
    elif local_negative:
        label = 'LOCAL_CONDITIONAL_NEGATIVE_ASSOCIATION_WITHOUT_MAIN_COMBINED_SUPPORT'
    else:
        label = 'NO_SUPPORTED_DIRECTION_UNDER_CURRENT_SECONDARY_FAMILY'
    # Differences between two P values never establish a difference between strata.
    contrast_supported = (estimable(contrast) and finite(contrast['BH']) and contrast['BH'] < .05
                          and finite(contrast['bootstrap_lower']) and finite(contrast['bootstrap_upper'])
                          and (contrast['bootstrap_lower'] > 0 or contrast['bootstrap_upper'] < 0))
    return dict(endpoint=endpoint, evidence_label=label,
        within_MKI67_combined_supported_direction_BH=supported_direction(main),
        local_composition_supported_direction_BH=supported_direction(local),
        within_MKI67_combined_support_BH=main_yes, local_composition_support_BH=local_yes,
        within_MKI67_combined_support_BY=positive_support(main, 'BY'),
        local_composition_support_BY=positive_support(local, 'BY'),
        within_MKI67_combined_negative_support_BH=main_negative,
        local_composition_negative_support_BH=local_negative,
        within_MKI67_combined_negative_support_BY=negative_support(main, 'BY'),
        local_composition_negative_support_BY=negative_support(local, 'BY'),
        MKI67_remaining_stratum_support=positive_support(records[endpoint, 'noncycling']),
        MKI67_high_stratum_support=positive_support(records[endpoint, 'cycling']),
        direct_between_stratum_difference_supported=contrast_supported,
        direct_between_stratum_difference_direction=('high_stronger' if contrast['median_effect'] > 0 else 'remaining_stronger') if contrast_supported else 'not_established',
        patient_leave_one_out_positive_main=(finite(main.get('leave_one_case_median_min')) and main['leave_one_case_median_min'] > 0),
        patient_leave_one_out_positive_local=(finite(local.get('leave_one_case_median_min')) and local['leave_one_case_median_min'] > 0),
        patient_leave_one_out_negative_main=(finite(main.get('leave_one_case_median_max')) and main['leave_one_case_median_max'] < 0),
        patient_leave_one_out_negative_local=(finite(local.get('leave_one_case_median_max')) and local['leave_one_case_median_max'] < 0),
        note='Direct and local-null rows may use different eligible images. No attenuation fraction, cross-grid effect difference, biological repulsion or directional reversal is inferred.')


def main():
    task = Path(os.environ['PANCANCER_LN_STRUCTURE_TASK']).resolve()
    repository = Path(__file__).resolve().parents[1]
    assert task != repository and repository not in task.parents
    complete = task/'private_LN_specificity_complete_source_v1.json'
    exited = task/'private_LN_specificity_actual_exit_v1.json'
    assert complete.exists() and exited.exists(), 'ACTUAL_AUDITED_SOURCE_EXIT_NOT_READY'
    source = json.loads(complete.read_text())
    exit_receipt = json.loads(exited.read_text())
    frozen = json.loads((task/'private_frozen_LN_specificity_execution_v1.json').read_text())
    family = frozen['analysis_config']['family_contract_SHA256']
    assert exit_receipt['actual_exit_code'] == 0 and source['source_R_exit_code'] == 0
    assert exit_receipt['source_result_SHA256'] == hashlib.sha256(complete.read_bytes()).hexdigest()
    assert source['source_payload_RAM_released'] and source['no_dynamic_followup_nodes']
    assert source['source_object_SHA256'] == frozen['expected_file_SHA256']
    assert source['R_code_SHA256'] == frozen['R_code_SHA256']
    assert source['frozen_secondary_family_SHA256'] == family
    summary = source['summary'];production = summary['production'];audit = summary['audit']
    assert summary['status'] == 'COMPLETE_FROZEN_SECONDARY_SPECIFICITY_AND_FULL_NUMERIC_AUDIT'
    assert audit['status'] == 'PASS_FULL_SECONDARY_SPECIFICITY_NUMERIC_AUDIT'
    assert production['family_contract_SHA256'] == audit['family_contract_SHA256'] == family
    records = validate_family(production)
    report = dict(status='AUDITED_SOURCE_SECONDARY_EVIDENCE_SUMMARY',family_contract_SHA256=family,
        actual_source_result_SHA256=hashlib.sha256(complete.read_bytes()).hexdigest(),
        endpoint_interpretations=[interpret_endpoint(records, e) for e in ENDPOINTS],
        primary_rows=production['primary_rows'], scale_point_checks=production['sensitivities'],
        edge_margin_quality=production['edge_margin_quality'], coverage_reasons=production['coverage_reasons'],
        local_exchangeability_coverage=production['mobility'],
        interpretation=['Data-informed secondary analysis of the same source; not independent replication.',
                        'MKI67-high and remaining are relative marker groups, not proof of cycling or quiescence.',
                        'Conditional randomization fixes local counts, lineage, degree and receiver geometry; fine-scale sender autocorrelation is not fully preserved.',
                        'A conditional-null effect is a centered log enrichment factor, not the uncentered MH odds ratio.',
                        'Eligible sets can differ. No attenuation percentages or cross-scale hypothesis tests are added.',
                        'BH covers the15 planned secondary tests, not the full exploration history.',
                        'Supported negative associations remain visible; spatial depletion is not proof of active repulsion. Opposite directions from different estimands or eligible sets do not establish reversal.',
                        'No endpoint is labelled direct communication, immune suppression, resistance or clinical prediction.'])
    target = task/'private_LN_specificity_source_evidence_summary_v1.json'
    assert not target.exists()
    target.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],family_contract_SHA256=family,
                         endpoint_evidence_labels=[r['evidence_label'] for r in report['endpoint_interpretations']])))


if __name__ == '__main__':
    main()
