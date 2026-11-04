# Source / cohort / sample-tissue registry for the 72 sources of integration v2 (2026-11-04).
# Inputs: GEO family SOFT files (fetched 2026-11-04, one per GSE), integration v2 obs (units actually used),
# PubMed esummary JSON (pmids from GEO, else from the project's candidate tables).
# Per source: GEO title, submission date, PubMed (journal, year, first/last author, title), contact
#   institute / city / country, contributors, platform(s), library chemistry keywords found in the sample
#   protocols, GEO sample count vs units used, units bridged to a GSM, tissue/source descriptions of the used
#   units, and flags (units whose GEO description says adjacent normal / blood / PBMC / lymph node / organoid /
#   cell line; same contact group as another source -> check patient overlap).
# Unit -> GSM: GSM accession in the unit id; else a unique GSM whose title, bracketed name or supplementary file
#   name matches a unit-id token; else a unique GSM whose title is a substring of the unit id.
# Outputs: OUT_PRIVATE/unit_registry_v1.tsv (per unit; server only), OUT_REPO/source_registry_v1.tsv,
#   OUT_REPO/unit_flags_v1.tsv (GSM accession + flag only). Sample descriptions in the source table are
#   aggregated, with patient numbers in GEO titles masked ('patient #')
# Usage: build_source_registry_v1.py SOFT_DIR V2_ZARR PUBMED_JSON PMID_FALLBACK_TSV OUT_PRIVATE OUT_REPO
import sys, os, re, gzip, json, glob, collections
import pandas as pd, zarr
from anndata.io import read_elem

soft, v2, pubj, pfall, opriv, orepo = sys.argv[1:7]
os.makedirs(opriv, exist_ok=True); os.makedirs(orepo, exist_ok=True)
CHEM = [('10x 5prime', r"5['’ ]?\s*(?:prime|end|v[12]|gene expression|reagent|library|\+ ?vdj|single cell 5)|single cell 5|chromium (?:next gem )?single cell 5"),
        ('10x 3prime', r"3['’ ]?\s*(?:prime|end|v[23]|gene expression|reagent|library)|single cell 3|chromium (?:next gem )?single cell 3"),
        ('10x (version unstated)', r'chromium|10x genomics|10 ?x'), ('Singleron GEXSCOPE', r'gexscope|singleron'),
        ('BD Rhapsody', r'rhapsody'), ('Smart-seq2', r'smart-?seq'), ('Drop-seq', r'drop-?seq'), ('DNBelab C4', r'dnbelab|c4 '),
        ('SeekOne', r'seekone|seekgene'), ('Seq-Well', r'seq-?well'), ('inDrop', r'indrop')]
FLAGS = [('adjacent_or_normal', r'adjacent|para-?(?:tumou?r|cancer|carcinoma)|non-?tumou?r|normal (?:tissue|mucosa|lung|liver)|healthy|\bnat\b|distant normal'),
         ('blood_or_pbmc', r'\bblood\b|pbmc|peripheral'), ('lymph_node', r'lymph node(?! stage)|\bln\b'), ('organoid_or_line', r'organoid|cell line|xenograft|pdx'),
         ('precancer', r'adenoma|dysplasia|precancer|pre-?malignan|\bais\b|\bmia\b|in situ')]
TISSUE_KEYS = ('tissue', 'organ', 'site', 'sample type', 'tissue type', 'location', 'region', 'disease', 'cell type', 'source')


comp = lambda t: re.sub(r'[^a-z0-9]', '', t.lower())
heads = lambda x: [x['title'], x['title'].split(':')[0].strip(), x['title'].split(' [')[0].strip()] + re.findall(r'\[([^\]]+)\]', x['title'])


def parse(gse):
    S = {'samples': [], 'platforms': {}}
    cur = None; plat = None; sec = None
    with gzip.open(f'{soft}/{gse}_family.soft.gz', 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('^'):
                sec = line[1:].split('=')[0].strip().upper(); acc = line.split('=', 1)[1].strip() if '=' in line else ''
                cur = {'gsm': acc, 'title': '', 'source_name': '', 'ch': [], 'supp': [], 'proto': []} if sec == 'SAMPLE' else None
                if cur: S['samples'].append(cur)
                plat = acc if sec == 'PLATFORM' else None
                continue
            if '=' not in line:
                continue
            k, v = line.split('=', 1); k = k.strip(); v = v.strip()
            if sec == 'SERIES':
                S.setdefault(k, []).append(v)
            elif sec == 'PLATFORM' and k == '!Platform_title':
                S['platforms'][plat] = v
            elif cur is not None:
                if k == '!Sample_title': cur['title'] = v
                elif k.startswith('!Sample_source_name'): cur['source_name'] = v
                elif k.startswith('!Sample_characteristics'): cur['ch'].append(v)
                elif k.startswith('!Sample_supplementary_file'): cur['supp'].append(os.path.basename(v))
                elif k.startswith(('!Sample_extract_protocol', '!Sample_data_processing', '!Sample_description', '!Sample_library')): cur['proto'].append(v)
                elif k == '!Sample_platform_id': cur['platform'] = v
    return S


g = zarr.open(v2, 'r')
obs = read_elem(g['obs'])[['source', 'unit_id', 'patient_id', 'sampling_site', 'material_detail', 'analysis_main', 'oncotree_code']]
U = obs.groupby('unit_id', observed=True).agg(source=('source', 'first'), cells=('source', 'size'), patient_id=('patient_id', 'first'),
                                               sampling_site=('sampling_site', 'first'), material_detail=('material_detail', 'first'),
                                               oncotree_code=('oncotree_code', 'first')).reset_index()
pub = json.load(open(pubj))['result'] if os.path.exists(pubj) else {}
fb = pd.read_csv(pfall, sep='\t', dtype=str).dropna().drop_duplicates('accession').set_index('accession')['pmid'].to_dict() if os.path.exists(pfall) else {}
urows, srows = [], []
for src, grp in U.groupby('source'):
    gse = src.split('_')[0]
    S = parse(gse)
    G = S['samples']
    for _, r in grp.iterrows():
        u = r['unit_id']; tail = u.split('|', 1)[-1]
        toks = (set(re.split(r'[|:\s]+', tail)) | set(re.split(r'[|_:\s]+', tail))) - {gse, src, ''}
        m = [x for x in G if x['gsm'] in u]; how = 'gsm_accession'
        core = re.sub(rf'^(?:{gse}_+)+', '', tail.split('|')[-1])  # unit tail without the GSE prefix
        if not m and core:
            wb = re.compile(rf'(?<![A-Za-z0-9]){re.escape(core)}(?![A-Za-z0-9])')
            for how, test in (('title_exact_norm', lambda x: comp(core) in {comp(h) for h in heads(x)}),
                              ('supp_file_norm', lambda x: any(re.sub(r'^GSM\d+_', '', fn).startswith((core + '_', core + '.')) for fn in x['supp'])),
                              ('title_word', lambda x: bool(wb.search(x['title']))),
                              ('characteristic_value', lambda x: comp(core) in {comp(c.split(':', 1)[-1]) for c in x['ch']})):
                m = [x for x in G if test(x)]
                if len(m) == 1:
                    break
            if len(m) != 1:
                m = []
        if not m:
            for how, test in (('title', lambda x: x['title'] in toks),
                              ('title_bracket', lambda x: any(b in toks for b in re.findall(r'\[([^\]]+)\]', x['title']))),
                              ('supp_file', lambda x: any(any(fn.split('_', 1)[-1].startswith(t + '_') or fn.split('_', 1)[-1].startswith(t + '.') for t in toks if len(t) >= 2) for fn in x['supp'])),
                              ('title_substring', lambda x: x['title'] and len(x['title']) >= 3 and x['title'] in tail)):
                m = [x for x in G if test(x)]
                if m:
                    break
        row = {'source': src, 'unit_id': u, 'cells_v2': r['cells'], 'patient_id': r['patient_id'], 'sampling_site': r['sampling_site'],
               'material_detail': r['material_detail'], 'match': how if len(m) == 1 else ('ambiguous_%d' % len(m) if m else 'none')}
        if len(m) == 1:
            x = m[0]
            tis = [c for c in x['ch'] if c.split(':', 1)[0].strip().lower().startswith(TISSUE_KEYS)]
            txt = ' | '.join([x['title'], x['source_name']] + x['ch'])
            row.update({'gsm': x['gsm'], 'gsm_title': x['title'], 'source_name': x['source_name'], 'tissue_characteristics': '; '.join(tis),
                        'geo_sample_description': ' | '.join(t for t in [re.sub(r'(?i)\b(patient|pt|case)[ _-]*[A-Za-z]*\d+\w*', r'\1 #', x['title']), x['source_name']] + tis if t),
                        'all_characteristics': '; '.join(x['ch']),
                        'geo_flags': ','.join(f for f, rx in FLAGS if re.search(rx, txt, re.I))})
        urows.append(row)
    proto = ' '.join(p for x in G for p in x['proto']).lower()
    chem = [n for n, rx in CHEM if re.search(rx, proto, re.I)]
    if '10x 5prime' in chem or '10x 3prime' in chem:
        chem = [c for c in chem if c != '10x (version unstated)']
    pmids = S.get('!Series_pubmed_id', []) or ([fb[gse]] if gse in fb else [])
    pm = pub.get(pmids[0], {}) if pmids else {}
    au = [a['name'] for a in pm.get('authors', [])]
    srows.append({'source': src, 'gse': gse, 'geo_title': (S.get('!Series_title') or [''])[0], 'geo_submission': (S.get('!Series_submission_date') or [''])[0],
                  'pmid': ';'.join(pmids), 'pmid_from': 'GEO' if S.get('!Series_pubmed_id') else ('project_table' if pmids else 'none'),
                  'journal': pm.get('source', ''), 'pub_year': pm.get('pubdate', '')[:4], 'first_author': au[0] if au else '', 'last_author': au[-1] if au else '',
                  'paper_title': pm.get('title', ''), 'contact_institute': (S.get('!Series_contact_institute') or [''])[0],
                  'contact_city': (S.get('!Series_contact_city') or [''])[0], 'contact_country': (S.get('!Series_contact_country') or [''])[0],
                  'geo_contributors_first_last': ' ... '.join(dict.fromkeys([(S.get('!Series_contributor') or [''])[0], (S.get('!Series_contributor') or [''])[-1]])),
                  'platforms': '; '.join(f'{k}: {v}' for k, v in S['platforms'].items()), 'chemistry_from_protocols': '; '.join(chem) or 'not stated',
                  'geo_samples_total': len(G), 'units_used': len(grp)})
P = pd.DataFrame(urows)
P.to_csv(f'{opriv}/unit_registry_v1.tsv', sep='\t', index=False)
R = pd.DataFrame(srows).set_index('source')
R['units_bridged_to_gsm'] = P.groupby('source')['gsm'].count()
OK = ['gsm_accession', 'title_exact_norm', 'supp_file_norm', 'title_word', 'characteristic_value', 'title', 'title_bracket', 'supp_file', 'title_substring']
R['units_unmatched'] = P.groupby('source')['match'].agg(lambda s: int((~s.isin(OK)).sum()))
R['geo_sample_descriptions_of_used_units'] = P.groupby('source')['geo_sample_description'].agg(lambda s: '; '.join(f'{k} ({v})' for k, v in s.fillna('unmatched').value_counts().items()))
dist = lambda s: '; '.join(f'{k} ({v})' for k, v in s.fillna('').replace('', 'unmatched').value_counts().items())
R['geo_source_name_of_used_units'] = P.groupby('source')['source_name'].agg(dist)
R['geo_tissue_of_used_units'] = P.groupby('source')['tissue_characteristics'].agg(dist)
R['our_sampling_site'] = P.groupby('source')['sampling_site'].agg(dist)
R['units_flagged'] = P.groupby('source')['geo_flags'].agg(lambda s: int((s.fillna('') != '').sum()))
R['flag_types'] = P.groupby('source')['geo_flags'].agg(lambda s: '; '.join(f'{k} ({v})' for k, v in collections.Counter(f for x in s.dropna() for f in x.split(',') if f).items()))
grpkey = R['last_author'].where(R['last_author'] != '', R['contact_institute'])
R['same_group_as'] = [', '.join(sorted(set(R.index[(grpkey == k) & (R.index != s)]))) if k else '' for s, k in grpkey.items()]
R.to_csv(f'{orepo}/source_registry_v1.tsv', sep='\t')
P.loc[P['geo_flags'].fillna('') != '', ['source', 'gsm', 'match', 'sampling_site', 'geo_flags']].to_csv(f'{orepo}/unit_flags_v1.tsv', sep='\t', index=False)
print(R[['gse', 'pmid', 'journal', 'pub_year', 'first_author', 'contact_country', 'chemistry_from_protocols', 'geo_samples_total', 'units_used', 'units_bridged_to_gsm', 'units_flagged', 'same_group_as']].to_string())
print(P['match'].value_counts().to_dict())
