"""Independently rebuild niche/CM activity summaries from saved spot outputs."""
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd

r = PROJECT_ROOT / 'runs/20261107_v4_Visium_author_validation_v2'
jobs = json.loads((r / 'inputs/section_jobs_PRIVATE.json').read_text())
summary = pd.read_csv(r / 'niches/all_CM_by_source_niche_alignment_summary.tsv', sep='\t')
data = pd.read_csv(r / 'niches/all_section_CM_niche_alignment_PRIVATE.tsv', sep='\t')
assert len(summary) == 14 * 2 * len({j['dataset'] for j in jobs})
assert not data.duplicated(['section_code', 'cm', 'posterior']).any()
rechecked = 0
for key in sorted(data.section_code.unique()):
    proof = json.loads((r / 'niches' / key / 'independent_audit.json').read_text())
    assert proof['status'].startswith('PASS_')
    group = pd.read_csv(r / 'niches' / key / 'spot_niches_PRIVATE.tsv.gz', sep='\t', index_col=0).niche
    for posterior in ['means', 'q05']:
        act = pd.read_csv(r / 'sections' / key / (posterior + '_CM_activity_PRIVATE.tsv.gz'), sep='\t', index_col=0)
        assert act.index.equals(group.index)
        for cm in act.columns:
            row = data[data.section_code.eq(key) & data.posterior.eq(posterior) & data.cm.eq(cm)].iloc[0]
            if row.n_members:
                y = act[cm]
                grand = y.mean()
                ss_total = np.square(y - grand).sum()
                group_means = y.groupby(group).mean()
                ss_between = sum(int((group == name).sum()) * (value - grand)**2 for name, value in group_means.items())
                expected = ss_between / ss_total if ss_total > 0 else np.nan
                assert (pd.isna(expected) and pd.isna(row.activity_eta_squared)) or np.isclose(expected, row.activity_eta_squared, atol=1e-10)
            else:
                assert pd.isna(row.activity_eta_squared)
            rechecked += 1
for row in summary.itertuples():
    selected = data[data.dataset.eq(row.dataset) & data.cm.eq(row.cm) & data.posterior.eq(row.posterior)]
    assert len(selected) == row.n_sections_with_both_workflows
    assert len([j for j in jobs if j['dataset'] == row.dataset]) == row.n_requested_sections
    actual = selected.activity_eta_squared.median()
    assert (pd.isna(actual) and pd.isna(row.median_activity_eta_squared)) or np.isclose(actual, row.median_activity_eta_squared, atol=1e-12)
proof = {'status': 'PASS_independent_all_requested_source_14_CM_niche_activity_variance_and_summary_readback',
    'UTC': datetime.now(timezone.utc).isoformat(), 'n_requested_sections': len(jobs),
    'n_sections_both_workflows': data.section_code.nunique(), 'n_section_CM_activity_rows_recomputed': rechecked,
    'n_source_CM_posterior_rows': len(summary),
    'no_independent_subtype_truth_or_contact_or_causal_or_significance_claim': True}
(r / 'niches/FINAL_NICHE_COHORT_INDEPENDENT_AUDIT.json').open('x').write(json.dumps(proof, indent=2))
print(json.dumps(proof))
