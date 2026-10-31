# Second-pass screen of GEO candidates from search_tissue_sources_v1.py against the frozen spec
# (methods/tissue-source-search-spec-v1-20261101.md). Read-only; the user makes the final call.
# Pass 1 (series level): human-only, 3-300 samples (plate-based series list one GSM per cell), not a
#   known source, not obviously non-tissue or functional/bulk, has a publication, offers count-like
#   supplementary files, title names tumour tissue.
# Pass 2 (sample level, top TOP per target by sample count): GSM records are read for species,
#   per-sample patient / tissue / treatment fields, and sorting / blood / non-scRNA wording.
# Usage: screen_tissue_candidates_v1.py CANDIDATES.tsv OUT_DIR [TOP]
import sys, os, re, time, json, urllib.request, collections
import pandas as pd

cand, out = sys.argv[1:3]
TOP = int(sys.argv[3]) if len(sys.argv) > 3 else 15
d = pd.read_csv(cand, sep='\t', dtype=str).fillna('')
g = d[d['repository'] == 'GEO'].drop_duplicates(['target', 'accession']).copy()
g['n'] = pd.to_numeric(g['n_samples'], errors='coerce')
text = g['title']
TUMOUR = r'tumou?r|cancer|carcinoma|glioma|glioblastoma|sarcoma|melanoma|seminoma|neoplasm|malignan'
g['p1_human_known_ok'] = ~g['flags'].str.contains('KNOWN_SOURCE|nontissue|multi_species')
g['p1_n3'] = g['n'] >= 3
g['p1_pub'] = g['pubmed'] != ''
g['p1_counts'] = g['supp'].str.contains(r'H5|MTX|TSV|CSV|TXT|RDS|H5AD|LOOM|TAR', case=False)
g['p1_tumour'] = text.str.contains(TUMOUR, case=False)
# plate-based series list one GSM per cell, and functional / bulk series are not tissue atlases
g['p1_size_3_300'] = g['n'].between(3, 300)
g['p1_not_functional'] = ~text.str.contains(r'knock-?down|knock-?out|cell lines?|organoid|xenograft|PDX|bulk|\[RNA-?seq\]|miRNA|methylation|in vitro|CRISPR|mouse|murine', case=False)
g['pass1'] = g[[c for c in g if c.startswith('p1_')]].all(axis=1)


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode('utf-8', 'replace').replace('\r', '')
        except Exception:
            time.sleep(2 ** i)
    return ''


KEY = {'patient': r'patient|donor|individual|subject|case', 'tissue': r'tissue|site|source|location|organ',
       'treatment': r'treat|therap|naive|chemo|prior|neoadjuvant|timepoint|time point'}
SORT = r'CD45\+|CD45-pos|sorted|FACS|CD3\+|CD8\+|CD4\+|CD14\+|MACS|enrich'
BLOOD = r'PBMC|peripheral blood|blood'
NONSC = r'ATAC|bulk|TCR|BCR|spatial|Visium|CITE|ChIP|methyl|WES|WGS|exome'
rows = []
for target, grp in g[g['pass1']].groupby('target'):
    for _, r in grp.sort_values('n', ascending=False).head(TOP).iterrows():
        t = get(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={r.accession}&targ=gsm&form=text&view=brief')
        time.sleep(1.0)
        samples, cur = [], None
        for line in t.splitlines():
            if line.startswith('^SAMPLE'):
                cur = {'org': '', 'title': '', 'src': '', 'keys': set(), 'chars': []}; samples.append(cur)
            elif cur is None:
                continue
            elif line.startswith('!Sample_organism_ch1'):
                cur['org'] = line.split('=', 1)[1].strip()
            elif line.startswith('!Sample_title'):
                cur['title'] = line.split('=', 1)[1].strip()
            elif line.startswith('!Sample_source_name_ch1'):
                cur['src'] = line.split('=', 1)[1].strip()
            elif line.startswith('!Sample_characteristics_ch1') and ':' in line:
                k, v = line.split('=', 1)[1].split(':', 1)
                cur['keys'].add(k.strip().lower()); cur['chars'].append(v.strip())
        n = len(samples) or 1
        hum = sum(s['org'] == 'Homo sapiens' for s in samples)
        blob = [' '.join([s['title'], s['src']] + s['chars']) for s in samples]
        frac = lambda rx: round(sum(bool(re.search(rx, b, re.I)) for b in blob) / n, 2)
        has = {k: round(sum(any(re.search(rx, kk) for kk in s['keys']) for s in samples) / n, 2) for k, rx in KEY.items()}
        rows.append({**r[['target', 'accession', 'title', 'n_samples', 'pubmed', 'date', 'supp', 'flags']].to_dict(),
                     'gsm_read': len(samples), 'human_frac': round(hum / n, 2),
                     'patient_field_frac': has['patient'], 'tissue_field_frac': has['tissue'],
                     'treatment_field_frac': has['treatment'], 'sorted_frac': frac(SORT),
                     'blood_frac': frac(BLOOD), 'non_scRNA_frac': frac(NONSC)})
s = pd.DataFrame(rows)


def verdict(r):
    issues = []
    if r.human_frac < 1: issues.append('non-human samples')
    if r.blood_frac >= 0.5: issues.append('mostly blood')
    if r.sorted_frac >= 0.5: issues.append('mostly sorted')
    if r.non_scRNA_frac >= 0.5: issues.append('mostly non-scRNA')
    if r.patient_field_frac == 0: issues.append('no patient field')
    if r.treatment_field_frac == 0: issues.append('no treatment field')
    return 'review' if not issues else 'review_with_issues: ' + ', '.join(issues)


if len(s):
    s['auto_note'] = s.apply(verdict, axis=1)
os.makedirs(out, exist_ok=True)
g.to_csv(os.path.join(out, 'pass1_series_screen.tsv'), sep='\t', index=False)
s.to_csv(os.path.join(out, 'pass2_shortlist.tsv'), sep='\t', index=False)
summary = g.groupby('target').agg(geo_hits=('accession', 'size'), pass1=('pass1', 'sum')).join(
    s.groupby('target').agg(shortlist=('accession', 'size'),
                            clean=('auto_note', lambda x: int((x == 'review').sum()))) if len(s) else None)
summary.to_csv(os.path.join(out, 'screen_summary.tsv'), sep='\t')
print(summary.to_string())
