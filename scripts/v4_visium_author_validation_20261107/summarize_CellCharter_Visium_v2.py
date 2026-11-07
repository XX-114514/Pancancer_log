"""All-request niche alignment diagnostics; spatial domains are section-specific."""
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd

r = PROJECT_ROOT / 'runs/20261107_v4_Visium_author_validation_v2'
base = PROJECT_ROOT / 'runs/20261106_v7_selective_CoVarNet_paper_v2'
jobs = json.loads((r / 'inputs/section_jobs_PRIVATE.json').read_text())
niche_terminal = json.loads((r / 'launch/all_niches_accepted_v2_terminal.json').read_text())
main_terminal = json.loads((r / 'launch/all_scores_terminal.json').read_text())
ns = {z['section_code']: z['status'] for z in niche_terminal['results']}
branches = {z['section_code']: z['niche_branch'] for z in niche_terminal['results']}
ms = {z['section_code']: z['status'] for z in main_terminal['results']}
members = pd.read_csv(base / 'CM_results/CM_members_after_author_network_filter.tsv', sep='\t')
rows = []
for job in jobs:
    key = job['section_code']
    if ns[key] != 'PASS_CELLCHARTER_FIT_AND_INDEPENDENT_AUDIT' or ms[key] != 'PASS_MODEL_AUTHOR_METRICS_AND_INDEPENDENT_AUDIT':
        continue
    p = r / 'sections' / key
    spots = pd.read_csv(r / branches[key] / key / 'spot_niches_PRIVATE.tsv.gz', sep='\t', index_col=0)
    coords = pd.read_csv(p / 'coordinates_PRIVATE.tsv.gz', sep='\t', index_col=0)
    assert spots.index.isin(coords.index).all()
    coords = coords.loc[spots.index]
    assert np.allclose(spots[['x', 'y']], coords[['x', 'y']])
    group = spots.niche
    for posterior in ['means', 'q05']:
        activity = pd.read_csv(p / (posterior + '_CM_activity_PRIVATE.tsv.gz'), sep='\t', index_col=0).loc[spots.index]
        abundance = pd.read_csv(p / (posterior + '_abundance_PRIVATE.tsv.gz'), sep='\t', index_col=0).loc[spots.index]
        assert activity.index.equals(spots.index) and abundance.index.equals(spots.index)
        means = abundance.groupby(group).mean()
        for cm in ['CM%02d' % i for i in range(1, 15)]:
            nodes = members.loc[members.cm.eq(cm), 'subCluster'].tolist()
            eta = enrichment = np.nan
            same_peak = False
            if nodes:
                values = activity[cm]
                total = float(np.square(values - values.mean()).sum())
                within = values - values.groupby(group).transform('mean')
                eta = 1 - float(np.square(within).sum()) / total if total > 0 else np.nan
                peak = values.groupby(group).mean().idxmax()
                denominator = abundance[nodes].mean()
                ratios = means.loc[peak, nodes] / denominator.replace(0, np.nan)
                enrichment = float(ratios.min()) if np.isfinite(ratios).all() else np.nan
                same_peak = bool(means[nodes].idxmax().nunique() == 1)
            rows.append({'section_code': key, 'dataset': job['dataset'], 'cm': cm, 'posterior': posterior,
                'n_members': len(nodes), 'activity_eta_squared': eta,
                'minimum_member_enrichment_at_activity_peak_niche': enrichment,
                'all_members_share_highest_mean_niche': same_peak,
                'n_niches': group.nunique()})
out = r / 'niches_final_v2'
out.mkdir(exist_ok=False)
data = pd.DataFrame(rows)
data.to_csv(out / 'all_section_CM_niche_alignment_PRIVATE.tsv', sep='\t', index=False)
assert len(data), 'No real section passed both independent workflows'
summary = []
for dataset in sorted({j['dataset'] for j in jobs}):
    requested = [j for j in jobs if j['dataset'] == dataset]
    accepted = [j for j in requested if ns[j['section_code']] == 'PASS_CELLCHARTER_FIT_AND_INDEPENDENT_AUDIT']
    for cm in ['CM%02d' % i for i in range(1, 15)]:
        for posterior in ['means', 'q05']:
            g = data[data.dataset.eq(dataset) & data.cm.eq(cm) & data.posterior.eq(posterior)]
            summary.append({'dataset': dataset, 'cm': cm, 'posterior': posterior,
                'n_requested_sections': len(requested), 'n_independently_audited_niche_sections': len(accepted),
                'n_sections_with_both_workflows': len(g),
                'n_evaluable_activity_eta_squared': int(g.activity_eta_squared.notna().sum()),
                'median_activity_eta_squared': g.activity_eta_squared.median(),
                'median_minimum_member_enrichment_at_activity_peak_niche': g.minimum_member_enrichment_at_activity_peak_niche.median(),
                'n_sections_all_members_share_highest_mean_niche': int(g.all_members_share_highest_mean_niche.sum()),
                'interpretation': 'Descriptive conditional projections onto expression-position niches; no permutation P, contact, causality or independent subtype accuracy claim'})
pd.DataFrame(summary).to_csv(out / 'all_CM_by_source_niche_alignment_summary.tsv', sep='\t', index=False)
receipt = {'status': 'ALL14_CM_SOURCE_NICHE_ALIGNMENT_SUMMARIZED_PENDING_COHORT_READBACK',
    'UTC': datetime.now(timezone.utc).isoformat(), 'n_requested_sections': len(jobs),
    'n_sections_both_workflows': data.section_code.nunique(), 'n_rows': len(summary),
    'n_original_CM': 14, 'posterior': ['means', 'q05'], 'spot_or_patient_independence_claim': False,
    'niches_fitted_without_CM_labels_or_abundance': True, 'public_private_output_requires_review': True}
(out / 'summary_complete.json').open('x').write(json.dumps(receipt, indent=2))
print(json.dumps(receipt))

qc=[]
for dataset in sorted({j['dataset'] for j in jobs}):
 selected=[j for j in jobs if j['dataset']==dataset];accepted=[j for j in selected if ns[j['section_code']]=='PASS_CELLCHARTER_FIT_AND_INDEPENDENT_AUDIT'];before=after=0
 for job in accepted:
  rec=json.loads((r/branches[job['section_code']]/job['section_code']/'complete.json').read_text());after+=rec['n_spots'];before+=rec.get('n_spots_before_niche_information_filter',rec['n_spots'])
 qc.append({'dataset':dataset,'n_requested_sections':len(selected),'n_audited_niche_sections':len(accepted),'n_failed_or_unassessable_sections':len(selected)-len(accepted),'n_spots_before_niche_information_filter_accepted_sections':before,'n_spots_niche_eligible_accepted_sections':after,'n_spots_zero_HVG_unassessable_for_niches':before-after})
pd.DataFrame(qc).to_csv(out/'niche_information_QC_by_source_summary.tsv',sep='\t',index=False)
