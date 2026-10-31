# Gap 2: derive a patient key for qualified units that have no recorded patient basis, using
# per-source naming rules checked against GEO sample titles and the publications (2026-11-01).
# Units with a recorded patient basis keep their existing patient_candidate_block.
# Usage: gap2_patient_keys_v1.py MATERIAL.tsv QUALIFIED_UNITS.tsv UNIT_GEO.tsv GSM_TEXT_DIR PRIVATE_OUT.tsv AGG_OUT.tsv
import sys, re
import pandas as pd

mat, qual, geo, gdir, priv, agg = sys.argv[1:7]
c = pd.read_csv(mat, sep='\t')
q = pd.read_csv(qual, sep='\t')[['qualification_unit_key_sha256', 'qualification_unit_id']]
c = c.merge(q, left_on='qualification_unit_key_SHA256', right_on='qualification_unit_key_sha256', how='left')
g = pd.read_csv(geo, sep='\t')[['qualification_unit_id', 'gsm']]
c = c.merge(g, on='qualification_unit_id', how='left')
titles = {}
for s in c.source.unique():
    try:
        cur = None
        for line in open(f'{gdir}/{s}_gsm.txt', encoding='utf-8', errors='replace'):
            if line.startswith('^SAMPLE'):
                cur = line.split('=')[1].strip()
            elif line.startswith('!Sample_title'):
                titles[cur] = line.split('=', 1)[1].strip()
    except FileNotFoundError:
        pass
c['title'] = c.gsm.map(titles).fillna('')


def stage274229(t):
    return 'mCRPC' if 'mCRPC' in t else 'mHSPC' if 'mHSPC' in t else 'localized' if 'Localized' in t else None


# source -> (field, regex or function, evidence)
RULES = {
    'GSE127465': ('id', r'human_p(\d+)t', 'sample name pXtY = patient X tumour sample Y; GEO lists 7 patients'),
    'GSE138709': ('id', r'ICC_(\d+)_', 'ICC_<patient>_Tumor[n]'),
    'GSE155698': ('title', r'PDAC_TISSUE_(\d+)', 'PDAC_TISSUE_11A/11B share patient 11'),
    'GSE178481': ('title', r'RCC-([A-Z]+\d+)-', 'RCC-<patient>-PTumor[n]: multi-region'),
    'GSE294109': ('title', r'\b(sc\d+)', 'sc<patient>, region suffix: multi-region'),
    'GSE274229': ('fn', lambda r: (lambda st, m: f'{st}_{m.group(1)}' if st and m else None)(
        stage274229(r.title), re.search(r'_S\d+_(\d+)_GEX', str(r.qualification_unit_id))),
        'stage + second file-name number; 6 mCRPC + 14 mHSPC + 5 localized = 25 patients as in PMID 39436703'),
}
UNKNOWN = {'GSE163558': 'primary and metastatic samples; pairing to patients not in GEO',
           'GSE228501': '4 samples 19355X1-4; patient count not stated in GEO'}


def key(r):
    if str(r.patient_basis_available) == 'True':
        return 'recorded', str(r.patient_candidate_block)
    if r.source in RULES:
        field, rule, _ = RULES[r.source]
        if field == 'fn':
            k = rule(r)
        else:
            m = re.search(rule, r.qualification_unit_id if field == 'id' else r.title)
            k = m.group(1) if m else None
        return ('source_rule', f'{r.source}:{k}') if k else ('rule_no_match', '')
    if r.source in UNKNOWN:
        return 'unknown', ''
    return 'assumed_one_unit_per_patient', f'{r.source}:{r.qualification_unit_id}'


c[['gap2_rule', 'patient_key']] = c.apply(lambda r: pd.Series(key(r)), axis=1)
c[['qualification_unit_key_SHA256', 'source', 'gap2_rule', 'patient_key']].to_csv(priv, sep='\t', index=False)
x = c[c.gap2_rule != 'recorded']
a = x.groupby(['source', 'gap2_rule']).agg(units=('source', 'size'), patients=('patient_key', lambda s: s[s != ''].nunique()))
a['evidence'] = [RULES[s][2] if s in RULES else UNKNOWN.get(s, '') for s, _ in a.index]
a.to_csv(agg, sep='\t')
print(a.to_string())
print(x.gap2_rule.value_counts().to_dict())
