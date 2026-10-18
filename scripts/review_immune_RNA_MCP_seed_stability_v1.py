"""Evaluate every frozen programme match; never select a model by effect/P."""
import itertools
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import rankdata, spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def centred_rank_correlation(left, right):
    """Within-source standardized ranks, followed by one pooled correlation."""
    paired = left.merge(right, on=['sample', 'dataset'], suffixes=('_A', '_B'), validate='one_to_one')
    if len(paired) != 38:
        return np.nan, {}, len(paired)
    vectors = [[], []]; source_rho = {}
    for source, g in paired.groupby('dataset', sort=True):
        if len(g) < 10:
            return np.nan, {}, len(paired)
        x = g.score_A.to_numpy(float); y = g.score_B.to_numpy(float)
        if not np.isfinite(np.r_[x, y]).all() or min(np.ptp(x), np.ptp(y)) == 0:
            return np.nan, {}, len(paired)
        ranks = [rankdata(v) for v in [x, y]]
        for i, v in enumerate(ranks):
            vectors[i].extend((v - v.mean()) / v.std(ddof=1))
        source_rho[source] = float(spearmanr(x, y).statistic)
    return float(np.corrcoef(vectors)[0, 1]), source_rho, len(paired)


def signed_weight_cosine(left, right):
    a = left.set_index('gene').weight; b = right.set_index('gene').weight
    if not len(a) or not len(b):
        return np.nan
    common = a.index.union(b.index)
    x = a.reindex(common, fill_value=0).to_numpy(float)
    y = b.reindex(common, fill_value=0).to_numpy(float)
    denominator = np.linalg.norm(x) * np.linalg.norm(y)
    return float(x @ y / denominator) if denominator > 0 else np.nan


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); budget = resource(root, 8 * 1024**2)
    contract_path = root / 'immune_RNA_MCP_postfit_contract_20261018_v1/plan_before_learned_outputs_or_projection.json'
    assert sha(contract_path) == '1172d887cffb7c3b0c28a21583154fbe3718a1f13d79edb85500343ead5e2b09'
    contract = json.loads(contract_path.read_text())
    audit_path = root / 'immune_RNA_MCP_frozen_projection_audit_20261018_v1/private_terminal.json'
    audit = json.loads(audit_path.read_text())
    assert audit['status'] == 'ALL_PLANNED_SLOTS_AND_ALL_SUPPORTED_PROJECTION_VALUES_ACCEPTED'
    producer = root / 'immune_RNA_MCP_frozen_projection_20261018_v1'
    assert sha(producer / 'private_terminal.json') == audit['producer_terminal_sha256']
    p = producer / 'private_R_projection'
    sp = p / 'private_all_case_mean_scores.tsv'
    scores = pd.read_csv(sp, sep='\t', dtype={'sample': str}) if sp.exists() else pd.DataFrame(
        columns=['fit_id','view','program','partition','sample','dataset','score'])
    scores = scores[scores.partition.eq('discovery')].copy()
    wp = p / 'private_frozen_gene_weights.tsv'
    weights = pd.read_csv(wp, sep='\t') if wp.exists() else pd.DataFrame(
        columns=['fit_id','view','program','gene','weight'])
    slots = pd.read_csv(p / 'private_all_planned_projection_slots.tsv', sep='\t')
    assert len(slots) == 33
    output = root / 'immune_RNA_MCP_seed_stability_20261018_v1'
    assert not output.exists(); output.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        contract_sha256=sha(contract_path), projection_acceptance_sha256=sha(audit_path), budget=budget,
        scope='All anchor2 by comparator2/2/3/2 programmes, all3 views and both source descriptions; no heldout scores, clinical outcomes or scientific P selection.',
        engineering_rule='Both other K2 seeds must align to anchor with each-view rank correlation>=0.70 and signed-weight cosine>=0.50; not a calibrated biological/clinical significance threshold.')
    (output / 'plan_before_seed_matching.json').write_text(json.dumps(plan, indent=2) + '\n')
    views = list(contract['view_to_official_name']); anchor = contract['anchor']
    rows = []; assignments = []; classifications = []
    for item in contract['planned_fits']:
        fit_id = 'K%d_seed%d_%s' % (item['k'], item['seed'], item['source_mode'])
        if fit_id == anchor:
            continue
        for a, b, view in itertools.product(['MCP1', 'MCP2'],
                ['MCP%d' % i for i in range(1, item['k'] + 1)], views):
            sa = scores[scores.fit_id.eq(anchor) & scores.view.eq(view) & scores.program.eq(a)]
            sb = scores[scores.fit_id.eq(fit_id) & scores.view.eq(view) & scores.program.eq(b)]
            wa = weights[weights.fit_id.eq(anchor) & weights.view.eq(view) & weights.program.eq(a)]
            wb = weights[weights.fit_id.eq(fit_id) & weights.view.eq(view) & weights.program.eq(b)]
            rho, source, n = centred_rank_correlation(sa, sb)
            cosine = signed_weight_cosine(wa, wb)
            rows.append(dict(comparator=fit_id, comparator_role=item['role'], anchor_program=a,
                comparator_program=b, view=view, n_discovery_cases=n,
                within_source_centred_rank_correlation=rho, signed_weight_cosine=cosine,
                HNSCC_rank_correlation=source.get('GSE286827', np.nan),
                ccRCC_rank_correlation=source.get('GSE294109', np.nan),
                status='AVAILABLE' if np.isfinite(rho) and np.isfinite(cosine) else 'UNAVAILABLE_SUPPORT_OR_VARIATION'))
        block = pd.DataFrame([r for r in rows if r['comparator'] == fit_id])
        matrix = np.full((2, item['k']), np.nan)
        for i, j in itertools.product(range(2), range(item['k'])):
            g = block[block.anchor_program.eq('MCP%d' % (i + 1)) & block.comparator_program.eq('MCP%d' % (j + 1))]
            if len(g) == 3 and g.status.eq('AVAILABLE').all():
                matrix[i, j] = g.within_source_centred_rank_correlation.mean()
        left, right = linear_sum_assignment(np.where(np.isfinite(matrix), -np.abs(matrix), 1e6))
        for i, j in zip(left, right):
            value = matrix[i, j]; valid = np.isfinite(value)
            sign = 1 if not valid or value >= 0 else -1
            g = block[block.anchor_program.eq('MCP%d' % (i + 1)) & block.comparator_program.eq('MCP%d' % (j + 1))]
            for r in g.to_dict('records'):
                r.update(global_orientation=sign, aligned_correlation=sign * r['within_source_centred_rank_correlation'],
                    aligned_weight_cosine=sign * r['signed_weight_cosine'], matching_available=bool(valid))
                assignments.append(r)
    full = pd.DataFrame(rows); matches = pd.DataFrame(assignments)
    assert len(full) == 54 and len(matches) == 24
    for program in ['MCP1', 'MCP2']:
        g = matches[matches.anchor_program.eq(program) & matches.comparator_role.eq('seed_stability')]
        assert len(g) == 6
        available = g.matching_available.all() and np.isfinite(g.aligned_correlation).all() and np.isfinite(g.aligned_weight_cosine).all()
        passes = available and (g.aligned_correlation >= .70).all() and (g.aligned_weight_cosine >= .50).all()
        classifications.append(dict(program=program, planned_seed_view_checks=6,
            all_seed_views_available=bool(available), prespecified_engineering_stability_pass=bool(passes),
            status='ENGINEERING_STABLE' if passes else 'UNSTABLE_OR_UNAVAILABLE',
            no_clinical_or_heldout_selection=True))
    full.to_csv(output / 'private_all54_seed_rank_source_pair_comparisons.tsv', sep='\t', index=False)
    matches.to_csv(output / 'private_all24_joint_matches_by_view.tsv', sep='\t', index=False)
    pd.DataFrame(classifications).to_csv(output / 'private_all2_anchor_engineering_classification.tsv', sep='\t', index=False)
    files = {p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file()}
    (output / 'private_terminal.json').write_text(json.dumps(dict(
        status='COMPLETE_ALL_FROZEN_MATCHES_INDEPENDENT_NUMERIC_REVIEW_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), files=files,
        planned_pair_view_comparisons=54, matched_pair_view_rows=24, anchor_programmes=2,
        heldout_or_clinical_tests_used_for_selection=False), indent=2) + '\n')
    print('ALL_FROZEN_STABILITY_MATCHES_DONE_NUMERIC_REVIEW_PENDING')


if __name__ == '__main__':
    main()
