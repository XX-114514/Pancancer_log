# Bridge qualified work units to GEO GSM records and pull per-sample characteristics that close
# data gaps: tumour grade/histology (gap 4 fine subtype) and treatment status (gap 1).
# Matching: unit id contains the GSM accession; else exactly one GSM whose title equals a unit-id token,
# or whose bracketed sample name does, or (last) whose title is a substring of the unit id.
# Writes a private per-unit table (server only) and a per-source aggregate (safe to commit).
# Usage: unit_geo_characteristics_v1.py QUALIFIED_UNITS.tsv GSM_TEXT_DIR PRIVATE_OUT.tsv AGG_OUT.tsv
import sys, re
import pandas as pd

units, gdir, priv, agg = sys.argv[1:5]
FIELDS = {  # characteristic key (lower case, prefix) -> output column
    'tumor grade': 'grade', 'histology-reduced': 'histology', 'histologic type': 'histology',
    'histolgical type': 'histology', 'tumor subtype': 'histology', 'tumor type': 'histology',
    'cancer type': 'histology', 'treatment': 'treatment', 'pre-surgery treatment': 'treatment',
    'treatment received prior to surgery': 'treatment', 'type of treatment': 'treatment',
    'time': 'timepoint', 'timepoint': 'timepoint', 'disease state': 'disease_state',
}


def gsms(g):
    out, cur = [], None
    for line in open(f'{gdir}/{g}_gsm.txt', encoding='utf-8', errors='replace'):
        if line.startswith('^SAMPLE'):
            cur = {'gsm': line.split('=')[1].strip(), 'title': '', 'ch': {}}
            out.append(cur)
        elif cur and line.startswith('!Sample_title'):
            cur['title'] = line.split('=', 1)[1].strip()
        elif cur and line.startswith('!Sample_characteristics_ch1'):
            kv = line.split('=', 1)[1].strip()
            if ':' in kv:
                k, v = kv.split(':', 1)
                k = k.strip().lower()
                for key, col in FIELDS.items():
                    if k.startswith(key):
                        cur['ch'].setdefault(col, v.strip())
                        break
    return out


q = pd.read_csv(units, sep='\t')
rows = []
for src, grp in q.groupby('source'):
    if not src.startswith('GSE'):
        continue
    G = gsms(src)
    for u in grp.qualification_unit_id:
        tail = u.split('|', 1)[-1]
        m = [x for x in G if x['gsm'] in u]
        how = 'gsm_accession'
        if not m:
            toks = set(re.split(r'[|:\s]+', tail)) | set(re.split(r'[|_:\s]+', tail))
            toks -= {src, ''}
            # exact title, then the bracketed sample name in the title, then a unique substring
            for how, test in (('title', lambda x: x['title'] in toks),
                              ('title_bracket', lambda x: any(b in toks for b in re.findall(r'\[([^\]]+)\]', x['title']))),
                              ('title_substring', lambda x: x['title'] and x['title'] in tail)):
                m = [x for x in G if test(x)]
                if m:
                    break
        r = {'source': src, 'qualification_unit_id': u, 'match': how if len(m) == 1 else ('ambiguous' if m else 'none')}
        if len(m) == 1:
            r['gsm'] = m[0]['gsm']
            r.update(m[0]['ch'])
        rows.append(r)
p = pd.DataFrame(rows)
p.to_csv(priv, sep='\t', index=False)
cols = [c for c in ('grade', 'histology', 'treatment', 'timepoint') if c in p]
a = p.groupby('source').agg(units=('match', 'size'), matched=('match', lambda s: int(s.str.startswith(('gsm', 'title')).sum())))
for c in cols:
    a[c] = p.groupby('source')[c].agg(lambda s: '; '.join(f'{k}: {v}' for k, v in s.dropna().value_counts().items()))
a.to_csv(agg, sep='\t')
print(a.to_string())
