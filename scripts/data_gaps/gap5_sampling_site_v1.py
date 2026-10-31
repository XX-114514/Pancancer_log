# Gap 5: classify the sampling site of qualified units whose material_class is
# tumor_stage_unspecified, from GEO source name / tissue fields. Sampling site (organ of origin vs
# metastatic site) is kept apart from disease stage (e.g. prostate biopsy in mHSPC, advanced TNBC).
# Usage: gap5_sampling_site_v1.py MATERIAL.tsv QUALIFIED_UNITS.tsv UNIT_GEO.tsv GSM_TEXT_DIR PRIVATE_OUT.tsv AGG_OUT.tsv
import sys, re
import pandas as pd

mat, qual, geo, gdir, priv, agg = sys.argv[1:7]
c = pd.read_csv(mat, sep='\t')
q = pd.read_csv(qual, sep='\t')[['qualification_unit_key_sha256', 'qualification_unit_id']]
c = c.merge(q, left_on='qualification_unit_key_SHA256', right_on='qualification_unit_key_sha256')
c = c.merge(pd.read_csv(geo, sep='\t')[['qualification_unit_id', 'gsm']], on='qualification_unit_id', how='left')
site = {}
for s in c.source.unique():
    try:
        cur = None
        for line in open(f'{gdir}/{s}_gsm.txt', encoding='utf-8', errors='replace'):
            if line.startswith('^SAMPLE'):
                cur = line.split('=')[1].strip(); site[cur] = []
            elif line.startswith('!Sample_source_name_ch1'):
                site[cur].append(line.split('=', 1)[1].strip())
            elif line.startswith('!Sample_characteristics_ch1') and re.match(r'\s*(tissue|sample type)\s*:', line.split('=', 1)[1], re.I):
                site[cur].append(line.split(':', 1)[1].strip())
    except FileNotFoundError:
        pass
METASTATIC = r'peritoneum|peritonium|metasta|lymph node'
# organ-of-origin words per source cancer (the sampled organ is where the tumour arises)
PRIMARY = (r'lung|stomach|primary gastric|pancrea|kidney|nasopharyn|prostate|thyroid|ovarian|colorectal|'
           r'cervical|testis|esophag|primary hepatocellular|breast|gastric cancer|primary tumou?r|primary tissue|human tumor tissue')
# sources where the sampled site is primary by biology or design, without a per-sample field
SOURCE_PRIMARY = {'GSE182109': 'diffuse glioma: CNS primary (gliomas do not metastasise outside the CNS)',
                  'GSE138665': 'GEO title: primary uveal melanomas',
                  'GSE138709': 'GEO summary: freshly resected ICC; unit names ICC_<n>_Tumor',
                  'GSE222315': 'GEO characteristics: bladder cancer tissues; unit names p<n>_BCa'}


def classify(r):
    if r.material_class != 'tumor_stage_unspecified':
        return 'recorded', r.material_class
    txt = ' | '.join(site.get(r.gsm, [])) if r.gsm == r.gsm else ''
    if re.search(METASTATIC, txt, re.I):
        return 'geo_site', 'metastatic_site'
    if re.search(PRIMARY, txt, re.I):
        return 'geo_site', 'primary_organ_site'
    if r.source in SOURCE_PRIMARY:
        return 'source_rule', 'primary_organ_site'
    return 'unresolved', ''


c[['gap5_rule', 'sampling_site']] = c.apply(lambda r: pd.Series(classify(r)), axis=1)
c[['qualification_unit_key_SHA256', 'source', 'gap5_rule', 'sampling_site']].to_csv(priv, sep='\t', index=False)
a = c[c.gap5_rule != 'recorded'].groupby(['source', 'gap5_rule', 'sampling_site']).size().rename('units')
a.to_csv(agg, sep='\t')
print(c[c.gap5_rule != 'recorded'].groupby(['gap5_rule', 'sampling_site']).size().to_string())
print('unresolved:', c[c.gap5_rule == 'unresolved'].groupby('source').size().to_dict())
