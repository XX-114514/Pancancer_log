"""Recheck full result review tables, original rules, counts and plot scope."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from review_NSCLC_native_interface_associations_v1 import bound_wait


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'NSCLC_native_interface_association_review_20261018_v1'
    terminal = out / 'private_terminal.json'
    bound_wait(root, 'review_NSCLC_native_interface_associations_v1', terminal)
    resource(root, 1024**2)
    manifest = json.loads(terminal.read_text())
    for name, item in manifest['files'].items():
        assert sha(out / name) == item['sha256']
    src = root / 'NSCLC_native_interface_association_results_20261018_v1'
    proof = src / 'private_independent_native_interface_result_audit.json'
    assert sha(proof) == manifest['full_result_audit_sha256']
    assert json.loads(proof.read_text())['producer_terminal_sha256'] == sha(src / 'private_terminal.json')
    original = pd.read_csv(src / 'private_all_primary_and_sensitivity_effects.tsv.gz', sep='\t', keep_default_na=False, na_values=['NA'])
    full = pd.read_csv(out / 'private_all1338_primary_sensitivity_and_coverage_review.tsv.gz', sep='\t', keep_default_na=False, na_values=['NA'])
    assert len(full) == len(original) == 1338
    pd.testing.assert_frame_equal(full[original.columns], original, check_dtype=False, atol=1e-14, rtol=1e-12)
    candidate = (original.partial_rank_rho > 0) & (original.BH_q < .05) & (original.MC_upper_BH_q < .05) & (original.CI_low > 0) & original.LOO_same_sign & (original.discovery_rho > 0)
    assert np.array_equal(candidate, original.internal_candidate_under_prespecified_rules)
    tests = [('PRIOR_ORIENTED_RHO_POSITIVE', original.partial_rank_rho > 0),
             ('PRIMARY_BH_BELOW_0P05', original.BH_q < .05),
             ('MC_UPPER_BH_BELOW_0P05', original.MC_upper_BH_q < .05),
             ('CI_STRICTLY_POSITIVE', original.CI_low > 0),
             ('ALL_LEAVE_ONE_SAME_SIGN', original.LOO_same_sign),
             ('DISCOVERY_RHO_POSITIVE', original.discovery_rho > 0)]
    failures = []
    for i in range(len(original)):
        failures.append(';'.join('NOT_' + name for name, test in tests if not test.iloc[i]))
    assert full.original_rule_failures.tolist() == failures
    failure_counts = pd.read_csv(out / 'private_original_rule_failure_counts.tsv', sep='\t')
    assert failure_counts.records.tolist() == [int((~test).sum()) for _, test in tests]
    assert failure_counts.rule_failure.tolist() == ['NOT_' + name for name, _ in tests]
    for name, test in tests:
        assert np.array_equal(full['passes_' + name], test)
    for label in ['original_receiver_programme', 'sender_programme_adjusted', 'native_immune_abundance_adjusted', 'complete_author_subtype_adjusted']:
        value = original[label + '_rho']
        for suffix, expected in [('_signed_change_from_primary', value - original.partial_rank_rho), ('_absolute_change_from_primary', abs(value) - abs(original.partial_rank_rho))]:
            assert np.allclose(full[label + suffix], expected, atol=1e-14, rtol=1e-12, equal_nan=True)
        assert np.array_equal(full[label + '_same_sign_as_primary'], value.notna() & (np.sign(value) == np.sign(original.partial_rank_rho)))
    triplets = pd.read_csv(out / 'private_all_literal_triplet_memberships.tsv.gz', sep='\t')
    assert triplets.literal_triplet_id.is_unique
    keys = ['sender_view', 'ligand', 'receiver_view']
    assert len(triplets) == len(original[keys].drop_duplicates())
    copied = full.set_index('record')
    all_members = []
    for row in triplets.itertuples():
        group = original[original.sender_view.eq(row.sender_view) & original.ligand.eq(row.ligand) & original.receiver_view.eq(row.receiver_view)]
        assert json.loads(row.member_records) == group.record.tolist()
        all_members.extend(group.record.tolist())
        assert copied.loc[group.record, 'literal_triplet_id'].eq(row.literal_triplet_id).all()
        assert row.all_dependent_records == len(group) and row.factors_tested == group.factor.nunique()
        selected = group[group.internal_candidate_under_prespecified_rules]
        assert row.candidate_records == len(selected) and row.candidate_factors == selected.factor.nunique()
        assert not row.independent_axes_or_direct_LR_validation
    assert sorted(all_members) == sorted(original.record.tolist())
    expected = original[candidate].sort_values('test_index')
    saved = pd.read_csv(out / 'private_all_original_internal_candidates.tsv.gz', sep='\t', keep_default_na=False, na_values=['NA'])
    pd.testing.assert_frame_equal(saved[original.columns].reset_index(drop=True), expected.reset_index(drop=True), check_dtype=False)
    ext = root / 'NSCLC_native_interface_external_opportunities_20261018_v1'
    ep = ext / 'private_independent_external_opportunity_audit.json'
    assert sha(ep) == manifest['external_opportunity_audit_sha256']
    e = json.loads((ext / 'private_terminal.json').read_text())
    name = 'private_all173940_interface_context_opportunities.tsv.gz'
    assert sha(ext / name) == e['files'][name]['sha256']
    coverage = pd.read_csv(ext / name, sep='\t', usecols=['record', 'dataset', 'all_units_coverage_opportunity', 'reviewed_units_coverage_opportunity'])
    for label in ['all_units', 'reviewed_units']:
        selected = coverage[coverage[label + '_coverage_opportunity']]
        sources = {record: set(g.dataset) for record, g in selected.groupby('record')}
        assert full[label + '_external_distinct_GSE_coverage_only'].tolist() == [len(sources.get(record, set())) for record in full.record]
    scope = pd.read_csv(root / 'NSCLC_continuous_prior_native_association_inputs_20261016_v1/private_all4800_pre_effect_eligibility.tsv.gz', sep='\t')
    summary = pd.read_csv(out / 'private_all120_factor_direction_summaries.tsv', sep='\t')
    assert len(summary) == 120 and not summary.duplicated(['factor', 'sender_view', 'receiver_view']).any()
    for row in summary.itertuples():
        choose = lambda d: d.factor.eq(row.factor) & d.sender_view.eq(row.sender_view) & d.receiver_view.eq(row.receiver_view)
        g = original[choose(original)]
        counts = [len(scope[choose(scope)]), len(g), int((g.partial_rank_rho > 0).sum()), int((g.partial_rank_rho < 0).sum()), int((g.BH_q < .05).sum()), int((g.BY_q < .05).sum()), int((g.MC_upper_BH_q < .05).sum()), int((g.CI_low > 0).sum()), int((g.CI_high < 0).sum()), int(g.CI_low.isna().sum()), int(g.internal_candidate_under_prespecified_rules.sum())]
        assert list(row)[4:] == counts
    maps = json.loads((out / 'private_exact_plot_record_mapping.json').read_text())
    assert maps[0]['records'] == maps[1]['records'] == original.record.tolist()
    assert [record for page in maps[2:] for record in page['records']] == expected.record.tolist()
    for page in maps:
        assert sha(out / page['name']) == page['sha256']
    report = dict(status='PASS_ALL1338_REVIEW_DERIVATIONS120_SUMMARIES_AND_EXACT_PLOT_SCOPE_NOT_PIXEL_QA',
        finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        producer_terminal_sha256=sha(terminal), primary_rows_checked=1338,
        new_models_P_thresholds_or_direct_communication_validation=False,
        actual_visual_review_completed=False)
    proof = out / 'private_independent_full_review_audit.json'
    assert not proof.exists()
    proof.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(status=report['status'], finished_UTC=report['finished_UTC'], audit_sha256=sha(proof))))


if __name__ == '__main__':
    main()
