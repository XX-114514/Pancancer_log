"""Complete fixed-family interpretation after the bound full numerical audit.

Derive descriptive review tables/figures only; no new model, P, candidate
threshold, endpoint or source selection. All scientific payload remains private.
"""
import itertools
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha, VIEWS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

SENSITIVITIES = ['original_receiver_programme', 'sender_programme_adjusted',
                 'native_immune_abundance_adjusted', 'complete_author_subtype_adjusted']


def bound_wait(root, task, proof):
    receipt = root / 'HNSCC_receiver_followup_job_control_20261016_v2' / task / 'launch_receipt.json'
    launch = json.loads(receipt.read_text())
    deadline = time.monotonic() + 43200
    while not proof.exists():
        proc = Path('/proc') / str(launch['pid']) / 'stat'
        assert proc.exists(), 'Bound dependency missing without proof; inspect original log'
        fields = proc.read_text().rsplit(')', 1)[1].split()
        assert fields[0] != 'Z' and fields[19] == str(launch['proc_start_ticks'])
        assert time.monotonic() < deadline, 'Observation expired; no termination inferred'
        time.sleep(5)
    assert json.loads(proof.read_text())['source_sha256'] == launch['source_sha256']
    return launch


def save_plot(fig, path, figures, records):
    fig.savefig(path, dpi=140)
    plt.close(fig)
    figures.append(dict(name=path.name, records=records, sha256=sha(path)))


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 32 * 1024**2)
    src = root / 'NSCLC_native_interface_association_results_20261018_v1'
    proof = src / 'private_independent_native_interface_result_audit.json'
    task = 'audit_NSCLC_native_interface_association_results_v1'
    receipt = root / 'HNSCC_receiver_followup_job_control_20261016_v2' / task / 'launch_receipt.json'
    out = root / 'NSCLC_native_interface_association_review_20261018_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        budget=budget, bound_full_audit_launch_sha256=sha(receipt),
        scope='All4800 planned positions, all1338 pre-effect eligible tests, all10 factors and12 directed pool pairs. Read effects only after bound full numerical audit.',
        classification='Exact original internal_candidate flag and each of its original rule failures. No new inferential category, ranking threshold or hypothesis.',
        sensitivity='Full receiver, add sender own programme, add both native immune pool fractions, add all author subtype fractions: existing same-case effects/status only, signed and absolute changes. No extra P, cutoff or artifact/causal diagnosis.',
        coverage='Join independently accepted effect-blind external opportunities for every record; count distinct GSE, not phase/arm records. Coverage is not replication or complete model qualification.',
        figures='All10x12 planned/eligible/candidate counts; all primary-vs-sensitivity points; every original internal candidate in fixed test_index order, paginated24 per forest. Full table retains all noncandidates and negative associations.',
        semantics='Original native internal author cases, not untouched external independent patients. New native F01..F10 are not2026-10 CM01..CM09. No direct LR, receptor-specific response, module dominance or clinical validation.',
        new_models_P_or_thresholds=False, effects_loaded_before_full_audit=False)
    (out / 'plan_before_full_native_result_review.json').write_text(json.dumps(plan, indent=2) + '\n')
    bound_wait(root, task, proof)
    audit = json.loads(proof.read_text())
    assert audit['status'] == 'PASS_COMPLETE_FROZEN_NATIVE_INTERFACE_FAMILY_ALL_REAL_EFFECTS_STREAMS_SENSITIVITIES_AND_CORRECTIONS'
    assert audit['producer_terminal_sha256'] == sha(src / 'private_terminal.json')
    terminal = json.loads((src / 'private_terminal.json').read_text())
    for name, item in terminal['files'].items():
        assert sha(src / name) == item['sha256']
    pre = root / 'NSCLC_continuous_prior_native_association_inputs_20261016_v1'
    scope = pd.read_csv(pre / 'private_all4800_pre_effect_eligibility.tsv.gz', sep='\t', keep_default_na=False)
    result = pd.read_csv(src / 'private_all_primary_and_sensitivity_effects.tsv.gz', sep='\t')
    active = scope[scope.test_index.ge(0)]
    assert len(result) == len(active) == audit['primary_tests_checked'] == 1338
    assert result.record.tolist() == active.record.tolist()
    checks = {
        'PRIOR_ORIENTED_RHO_POSITIVE': result.partial_rank_rho.gt(0),
        'PRIMARY_BH_BELOW_0P05': result.BH_q.lt(.05),
        'MC_UPPER_BH_BELOW_0P05': result.MC_upper_BH_q.lt(.05),
        'CI_STRICTLY_POSITIVE': result.CI_low.gt(0),
        'ALL_LEAVE_ONE_SAME_SIGN': result.LOO_same_sign,
        'DISCOVERY_RHO_POSITIVE': result.discovery_rho.gt(0)}
    assert all(x.dtype == bool for x in checks.values())
    passes = np.column_stack(list(checks.values()))
    assert np.array_equal(passes.all(axis=1), result.internal_candidate_under_prespecified_rules)
    result['original_rule_failures'] = [';'.join('NOT_' + k for k, ok in zip(checks, row) if not ok) for row in passes]
    for name, values in checks.items():
        result['passes_' + name] = values
    for label in SENSITIVITIES:
        values = result[label + '_rho']
        result[label + '_signed_change_from_primary'] = values - result.partial_rank_rho
        result[label + '_absolute_change_from_primary'] = values.abs() - result.partial_rank_rho.abs()
        result[label + '_same_sign_as_primary'] = values.notna() & (np.sign(values) == np.sign(result.partial_rank_rho))
    keys = ['sender_view', 'ligand', 'receiver_view']
    result['literal_triplet_id'] = ['NATIVE_RNA_TRIPLET_' + hashlib.sha256(
        json.dumps(list(values), separators=(',', ':')).encode()).hexdigest()[:16]
        for values in result[keys].itertuples(index=False, name=None)]
    assert result[keys].drop_duplicates().shape[0] == result.literal_triplet_id.nunique()
    triplets = []
    for key, g in result.groupby('literal_triplet_id', sort=True):
        assert len(g[keys].drop_duplicates()) == 1
        selected = g[g.internal_candidate_under_prespecified_rules]
        triplets.append(dict(literal_triplet_id=key, **g.iloc[0][keys].to_dict(),
            all_dependent_records=len(g), factors_tested=g.factor.nunique(),
            candidate_records=len(selected), candidate_factors=selected.factor.nunique(),
            member_records=json.dumps(g.record.tolist(), separators=(',', ':')),
            independent_axes_or_direct_LR_validation=False))
    pd.DataFrame(triplets).to_csv(out / 'private_all_literal_triplet_memberships.tsv.gz',
                                sep='\t', index=False, compression='gzip')
    ext = root / 'NSCLC_native_interface_external_opportunities_20261018_v1'
    ep = ext / 'private_independent_external_opportunity_audit.json'
    assert sha(ep) == 'ee906fff50c57729431a7e018a47403c5b7f05ad24343ceb1d229611ccdb5bdf'
    assert json.loads(ep.read_text())['producer_terminal_sha256'] == sha(ext / 'private_terminal.json')
    em = json.loads((ext / 'private_terminal.json').read_text())
    name = 'private_all173940_interface_context_opportunities.tsv.gz'
    assert sha(ext / name) == em['files'][name]['sha256']
    contexts = pd.read_csv(ext / name, sep='\t', usecols=['record', 'dataset', 'all_units_coverage_opportunity', 'reviewed_units_coverage_opportunity'])
    for label in ['all_units', 'reviewed_units']:
        n = contexts[contexts[label + '_coverage_opportunity']].groupby('record').dataset.nunique()
        result[label + '_external_distinct_GSE_coverage_only'] = result.record.map(n).fillna(0).astype(int)
    result.to_csv(out / 'private_all1338_primary_sensitivity_and_coverage_review.tsv.gz', sep='\t', index=False, compression='gzip', na_rep='NA')
    result[result.internal_candidate_under_prespecified_rules].to_csv(out / 'private_all_original_internal_candidates.tsv.gz', sep='\t', index=False, compression='gzip', na_rep='NA')
    pd.DataFrame([dict(rule_failure='NOT_' + k, records=int((~v).sum())) for k, v in checks.items()]).to_csv(out / 'private_original_rule_failure_counts.tsv', sep='\t', index=False)
    factors = ['NSCLC_MFA_20261016_F%02d' % k for k in range(1, 11)]
    directions = list(itertools.permutations(VIEWS, 2))
    summaries = []
    for factor, (sender, receiver) in itertools.product(factors, directions):
        g = result[result.factor.eq(factor) & result.sender_view.eq(sender) & result.receiver_view.eq(receiver)]
        planned = scope[scope.factor.eq(factor) & scope.sender_view.eq(sender) & scope.receiver_view.eq(receiver)]
        summaries.append(dict(factor=factor, sender_view=sender, receiver_view=receiver,
            planned_positions=len(planned), eligible_primary_tests=len(g),
            prior_oriented_positive=int(g.partial_rank_rho.gt(0).sum()),
            prior_oriented_negative=int(g.partial_rank_rho.lt(0).sum()),
            BH_supported=int(g.BH_q.lt(.05).sum()), BY_supported=int(g.BY_q.lt(.05).sum()),
            MC_upper_BH_supported=int(g.MC_upper_BH_q.lt(.05).sum()),
            positive_CI=int(g.CI_low.gt(0).sum()), negative_CI=int(g.CI_high.lt(0).sum()),
            CI_unavailable=int(g.CI_low.isna().sum()),
            original_internal_candidates=int(g.internal_candidate_under_prespecified_rules.sum())))
    summary = pd.DataFrame(summaries)
    assert summary.planned_positions.sum() == 4800 and summary.eligible_primary_tests.sum() == 1338
    summary.to_csv(out / 'private_all120_factor_direction_summaries.tsv', sep='\t', index=False)
    figures = []
    names = {v: s for v, s in zip(VIEWS, ['CD8', 'Treg', 'Mac', 'B/Plasma'])}
    labels = [names[a] + '→' + names[b] for a, b in directions]
    fig, axes = plt.subplots(1, 2, figsize=(18, 7), layout='constrained')
    for ax, field, title in zip(axes, ['eligible_primary_tests', 'original_internal_candidates'], ['Pre-effect eligible records', 'Original internal candidate rules met']):
        matrix = summary[field].to_numpy().reshape(10, 12)
        ax.imshow(matrix, aspect='auto', cmap='Blues', vmin=0, vmax=max(1, matrix.max()))
        for (i, j), value in np.ndenumerate(matrix):
            ax.text(j, i, str(value), ha='center', va='center', fontsize=9)
        ax.set_xticks(range(12), labels, rotation=60, ha='right')
        ax.set_yticks(range(10), ['F%02d' % k for k in range(1, 11)])
        ax.set_title(title)
    fig.suptitle('Native internal ligand–programme associations; counts are dependent records, not independent axes\nZero eligible records means untestable coverage; no direct communication or clinical validation')
    save_plot(fig, out / 'private_all120_record_count_overview.png', figures, result.record.tolist())
    fig, axes = plt.subplots(2, 2, figsize=(12, 12), layout='constrained')
    candidate = result.internal_candidate_under_prespecified_rules.to_numpy(bool)
    for ax, field in zip(axes.flat, SENSITIVITIES):
        valid = np.isfinite(result[field + '_rho'])
        ax.scatter(result.loc[valid, 'partial_rank_rho'], result.loc[valid, field + '_rho'], c=np.where(candidate[valid], '#d55e00', '#667788'), s=8, alpha=.5)
        ax.plot([-1, 1], [-1, 1], color='black', linewidth=.6)
        ax.axhline(0, color='gray', linewidth=.4); ax.axvline(0, color='gray', linewidth=.4)
        ax.set(xlim=(-1, 1), ylim=(-1, 1), xlabel='Primary partial-rank rho', ylabel='Existing descriptive comparison rho',
               title=field.replace('_', ' ') + '\nfinite %d / %d' % (valid.sum(), len(result)))
    fig.suptitle('All1338 records; orange=original internal candidate; no new P or rescue by sensitivity')
    save_plot(fig, out / 'private_all1338_existing_sensitivity_comparisons.png', figures, result.record.tolist())
    candidates = result[result.internal_candidate_under_prespecified_rules].sort_values('test_index')
    for start in range(0, len(candidates), 24):
        g = candidates.iloc[start:start + 24]
        fig, ax = plt.subplots(figsize=(15, max(4, .34 * len(g) + 2)), layout='constrained')
        labels = []
        for j, row in enumerate(g.itertuples()):
            labels.append('%s %s→%s %s (sign %+d)' % (row.factor.rsplit('_', 1)[1], names[row.sender_view], names[row.receiver_view], row.ligand, row.direction))
            ax.plot([row.CI_low, row.CI_high], [j, j], color='#2166ac')
            ax.scatter(row.partial_rank_rho, j, color='#2166ac', s=20)
        ax.set_yticks(range(len(g)), labels, fontsize=9)
        ax.set(xlim=(-1, 1), xlabel='Primary prior-oriented partial-rank rho and95% case-bootstrap CI',
               title='Every original internal candidate, frozen test order; page%d\nRNA association only; receptor-specific causality and external replication unproven' % (start // 24 + 1))
        ax.axvline(0, color='gray', linewidth=.7)
        ax.invert_yaxis()
        save_plot(fig, out / ('private_all_internal_candidates_page_%03d.png' % (start // 24 + 1)), figures, g.record.tolist())
    (out / 'private_exact_plot_record_mapping.json').write_text(json.dumps(figures, indent=2) + '\n')
    report = dict(status='COMPLETE_FULL_NATIVE_RESULT_REVIEW_DERIVATIONS_AUDIT_AND_VISUAL_REVIEW_PENDING',
        finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        full_result_audit_sha256=sha(proof), original_primary_tests=1338, planned_positions=4800,
        external_opportunity_audit_sha256=sha(ep),
        files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
        new_models_P_or_thresholds=False, direct_communication_clinical_or_external_validation=False)
    assert sum(v['bytes'] for v in report['files'].values()) < budget['output_cap_bytes']
    (out / 'private_terminal.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(status=report['status'], finished_UTC=report['finished_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
