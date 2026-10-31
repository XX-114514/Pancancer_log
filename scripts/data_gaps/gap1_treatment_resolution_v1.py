# Gap 1: resolve treatment status for qualified units whose recorded treatment_phase is unverified.
# Evidence order: recorded phase (kept) > per-sample GEO field > sample-name label > series/paper
# statement for the whole source. Mixed sources without per-sample evidence stay unresolved.
# Usage: gap1_treatment_resolution_v1.py MATERIAL.tsv QUALIFIED_UNITS.tsv UNIT_GEO.tsv EVIDENCE.tsv PRIVATE_OUT.tsv AGG_OUT.tsv
import sys, re
import pandas as pd

mat, qual, geo, ev, priv, agg = sys.argv[1:7]
c = pd.read_csv(mat, sep='\t')
q = pd.read_csv(qual, sep='\t')[['qualification_unit_key_sha256', 'qualification_unit_id']]
c = c.merge(q, left_on='qualification_unit_key_SHA256', right_on='qualification_unit_key_sha256', how='left')
g = pd.read_csv(geo, sep='\t')
c = c.merge(g[['qualification_unit_id', 'treatment']], on='qualification_unit_id', how='left')
e = pd.read_csv(ev, sep='\t').set_index('source')
NAIVE = {'0', 'no', 'untreated', 'no treatment', 'treatment naïve', 'treatment naive'}


def resolve(r):
    if '待核实' not in str(r.treatment_phase):
        return 'recorded', str(r.treatment_phase)
    t = str(r.treatment).strip().lower() if r.treatment == r.treatment else ''
    if t:
        return ('gsm_field', 'untreated') if t in NAIVE else ('gsm_field', 'prior_or_other:' + str(r.treatment))
    if r.source == 'GSE182109':
        lab = re.search(r'(ndGBM|rGBM|LGG)', str(r.qualification_unit_id))
        if lab:
            return 'sample_label', {'ndGBM': 'untreated_newly_diagnosed', 'rGBM': 'recurrent_after_standard_therapy',
                                    'LGG': 'LGG_status_unstated'}[lab.group(1)]
    if r.source in e.index and e.loc[r.source, 'value'] in ('naive', 'naive_chemo_RT'):
        return e.loc[r.source, 'rule'], e.loc[r.source, 'value']
    return 'unresolved', ''


c[['gap1_rule', 'gap1_status']] = c.apply(lambda r: pd.Series(resolve(r)), axis=1)
c[['qualification_unit_key_SHA256', 'source', 'gap1_rule', 'gap1_status']].to_csv(priv, sep='\t', index=False)
a = c[c.gap1_rule != 'recorded'].groupby(['source', 'gap1_rule', 'gap1_status']).size().rename('units')
a.to_csv(agg, sep='\t')
print(c.gap1_rule.value_counts().to_string())
print('unresolved sources:', c[c.gap1_rule == 'unresolved'].groupby('source').size().to_dict())
