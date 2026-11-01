# Detailed check of shortlisted GEO candidates against the frozen spec
# (methods/tissue-source-search-spec-v1-20261101.md). Read-only on public records; downloads at most
# one small count file per series into OUT_DIR/probe/<GSE>/ to test that values are integers.
# Per candidate: GEO series + GSM full records (species, technology, sorting, sample sources,
# patient field), supplementary file names (raw counts vs normalised), paper methods sentences
# (treatment, fresh/frozen, sorting) from PMC, and PubMed overlap with project sources.
# Usage: detailed_check_candidates_v1.py CURATED.tsv KNOWN_ACCESSIONS.txt OUT_DIR [PRIORITY=A]
import sys, os, re, json, time, gzip, io, urllib.request, urllib.parse, collections
import pandas as pd

cur, known_f, out = sys.argv[1:4]
prio = sys.argv[4] if len(sys.argv) > 4 else 'A'
os.makedirs(out, exist_ok=True)
known = set(open(known_f).read().split())
EU = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def get(url, binary=False, tries=4, limit=None):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'pancancer-audit'})
            with urllib.request.urlopen(req, timeout=180) as r:
                b = r.read(limit) if limit else r.read()
                return b if binary else b.decode('utf-8', 'replace').replace('\r', '')
        except Exception:
            time.sleep(2 ** i)
    return b'' if binary else ''


def soft(acc, targ):
    return get(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ={targ}&form=text&view=full')


def fields(text, prefix):
    return [l.split('=', 1)[1].strip() for l in text.splitlines() if l.startswith(prefix)]


TECH = [('10x', r'10x|chromium|cell ?ranger'), ('Smart-seq', r'smart-?seq'), ('Drop-seq', r'drop-?seq'),
        ('BD Rhapsody', r'rhapsody'), ('Singleron', r'singleron|gexscope'), ('inDrop', r'indrop'),
        ('Seq-Well', r'seq-?well'), ('MGI/DNBelab', r'dnbelab|c4 '), ('snRNA', r'single[- ]nucle|snrna|nuclei isolation')]
SORT = r'CD45\+|CD45 ?pos|CD45-negative|CD45-|FACS|sorted|MACS|magnetic|enrich|depleti'
NORM = r'tpm|fpkm|rpkm|norm|scaled|log|cpm|data\.txt'
RAW = r'matrix\.mtx|\.h5$|\.h5\.|filtered_feature|raw_feature|barcodes|features|genes\.tsv|count|umi|\.rds|\.h5ad|\.loom|_dge|expression_matrix'


def tobytes(sv):
    m = re.match(r'([\d.]+)([KMG]?)', sv or '')
    return float(m.group(1)) * {'': 1, 'K': 1e3, 'M': 1e6, 'G': 1e9}[m.group(2)] if m else None


def listing(url):
    html = get(url)
    return re.findall(r'href="([^"?/][^"]*)"', html)


def sizes(url):
    html = get(url)
    return {m.group(1): m.group(2) for m in re.finditer(r'href="([^"]+)">[^<]+</a>\s+\S+ \S+\s+(\S+)', html)}


def pmc_sentences(pmid):
    j = get(f'https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/?ids={pmid}&format=json&tool=pancancer_audit&email=none@example.org')
    try:
        pmc = json.loads(j)['records'][0].get('pmcid')
    except Exception:
        pmc = None
    if not pmc:
        return None, {}
    t = get(EU + f'efetch.fcgi?db=pmc&id={pmc[3:]}&retmode=xml&tool=pancancer_audit')
    t = re.sub(r'<ref-list>.*?</ref-list>', '', t, flags=re.S)
    t = re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', ' ', t))
    sents = re.split(r'(?<=[.;])\s+(?=[A-Z])', t)
    pats = {'treatment': r'treatment[- ]na[iï]ve|untreated|without (any )?(prior|previous|preoperative|neoadjuvant)|no (prior|previous|preoperative|neoadjuvant)|had not received|neoadjuvant',
            'fresh_frozen': r'\bfresh(ly)?\b|cryopreserv|frozen|thaw',
            'sorting': r'CD45\+|FACS|sorted|MACS|magnetic bead|enrich',
            'patients': r'\b\d+ (patients|individuals|donors|cases)\b'}
    ctx = r'patient|sample|specimen|tumou?r|biops|resect|tissue|cells were|dissociat|suspension'
    return pmc, {k: [s[:300] for s in sents if re.search(p, s, re.I) and re.search(ctx, s, re.I)][:5] for k, p in pats.items()}


def integer_probe(gse, files):
    """download the smallest likely-raw count file under 60 MB and test integrality"""
    cand = sorted((f for f in files if re.search(r'matrix\.mtx|count|umi|\.csv|\.tsv|\.txt', f['name'], re.I)
                   and not re.search(NORM + r'|barcodes|features|genes\.tsv|meta|annot|cluster|clinical', f['name'], re.I)
                   and f['bytes'] and f['bytes'] < 60e6), key=lambda f: f['bytes'])
    if not cand:
        return {'probe': 'no_small_count_file'}
    f = cand[0]
    raw = get(f['url'], binary=True)
    d = os.path.join(out, 'probe', gse); os.makedirs(d, exist_ok=True)
    open(os.path.join(d, f['name']), 'wb').write(raw)
    try:
        txt = gzip.decompress(raw).decode('utf-8', 'replace') if f['name'].endswith('.gz') else raw.decode('utf-8', 'replace')
    except Exception as e:
        return {'probe': 'unreadable', 'file': f['name'], 'error': repr(e)[:100]}
    vals = []
    for line in txt.splitlines()[:20000]:
        if line.startswith('%'):
            continue
        parts = re.split(r'[,\t ]+', line.strip())
        for p in (parts[2:3] if f['name'].lower().find('mtx') >= 0 else parts[1:]):
            try: vals.append(float(p.strip('"')))
            except ValueError: pass
    if not vals:
        return {'probe': 'no_numeric_values', 'file': f['name']}
    nonint = sum(abs(v - round(v)) > 1e-6 for v in vals)
    return {'probe': 'integer' if nonint == 0 else 'NON_INTEGER', 'file': f['name'], 'values_checked': len(vals),
            'non_integer': nonint, 'max': max(vals)}


c = pd.read_csv(cur, sep='\t')
accs = list(dict.fromkeys(c.loc[c.priority == prio, 'accession']))
targets = c[c.priority == prio].groupby('accession')['target'].agg(lambda s: '; '.join(dict.fromkeys(s)))
known_pm = {}
results = []
for acc in accs:
    s = soft(acc, 'self'); time.sleep(1)
    g = soft(acc, 'gsm'); time.sleep(1)
    pm = fields(s, '!Series_pubmed_id')
    summary = ' '.join(fields(s, '!Series_summary') + fields(s, '!Series_overall_design'))
    rel = fields(s, '!Series_relation')
    samples = []
    curd = None
    for line in g.splitlines():
        if line.startswith('^SAMPLE'):
            curd = {'gsm': line.split('=')[1].strip(), 'chars': {}, 'text': [], 'supp': []}; samples.append(curd)
        elif curd is None:
            continue
        elif line.startswith('!Sample_organism_ch1'): curd['org'] = line.split('=', 1)[1].strip()
        elif line.startswith('!Sample_library_strategy'): curd['lib'] = line.split('=', 1)[1].strip()
        elif line.startswith('!Sample_source_name_ch1'): curd['src'] = line.split('=', 1)[1].strip()
        elif line.startswith('!Sample_title'): curd['title'] = line.split('=', 1)[1].strip()
        elif line.startswith('!Sample_characteristics_ch1') and ':' in line:
            k, v = line.split('=', 1)[1].split(':', 1); curd['chars'].setdefault(k.strip().lower(), v.strip())
        elif line.startswith(('!Sample_extract_protocol', '!Sample_data_processing', '!Sample_description', '!Sample_treatment_protocol')):
            curd['text'].append(line.split('=', 1)[1])
        elif line.startswith('!Sample_supplementary_file') and 'NONE' not in line:
            curd['supp'].append(line.split('=', 1)[1].strip())
    alltext = ' '.join(summary.split()) + ' ' + ' '.join(' '.join(x['text']) for x in samples)
    tech = [n for n, rx in TECH if re.search(rx, alltext, re.I)]
    n = len(samples) or 1
    blob = [' '.join([x.get('title', ''), x.get('src', '')] + list(x['chars'].values())) for x in samples]
    pkeys = [k for k in collections.Counter(k for x in samples for k in x['chars']) if re.search(r'patient|donor|individual|subject|case', k)]
    patients = len({x['chars'].get(pkeys[0]) for x in samples}) if pkeys else None
    # supplementary files: series-level listing plus sample-level files
    gnum = re.sub(r'\d{3}$', 'nnn', acc)
    surl = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{gnum}/{acc}/suppl/'
    files = []
    sz = sizes(surl)
    for name in listing(surl):
        if name.startswith(('GSE', 'GSM')) or '.' in name:
            files.append({'name': name, 'url': surl + name, 'size': sz.get(name, ''), 'level': 'series'})
    for x in samples:
        for u in x['supp']:
            files.append({'name': u.rsplit('/', 1)[-1], 'url': u.replace('ftp://', 'https://'), 'size': '', 'level': 'sample'})
    for f in files: f['bytes'] = tobytes(f['size'])
    names = ' '.join(f['name'] for f in files)
    raw_files = [f['name'] for f in files if re.search(RAW, f['name'], re.I) and not re.search(NORM, f['name'], re.I)]
    norm_files = [f['name'] for f in files if re.search(NORM, f['name'], re.I)]
    # sample-level files lack sizes; probe with series files first
    probe = integer_probe(acc, [f for f in files if f['level'] == 'series'])
    if probe.get('probe') == 'no_small_count_file' and samples and samples[0]['supp']:
        # most 10x series keep matrices at sample level: list the first sample's folder with sizes
        sdir = samples[0]['supp'][0].replace('ftp://', 'https://').rsplit('/', 1)[0] + '/'
        ssz = sizes(sdir)
        probe = integer_probe(acc, [{'name': k, 'url': sdir + k, 'bytes': tobytes(v)} for k, v in ssz.items()])
        probe['probe_level'] = 'sample'
    time.sleep(1)
    pmc, sent = (None, {})
    if pm:
        pmc, sent = pmc_sentences(pm[0]); time.sleep(1)
    r = {'accession': acc, 'target': targets[acc], 'title': fields(s, '!Series_title')[0] if fields(s, '!Series_title') else '',
         'pubmed': ','.join(pm), 'pmc': pmc or '', 'relations': '; '.join(rel)[:300],
         'n_gsm': len(samples), 'human_frac': round(sum(x.get('org') == 'Homo sapiens' for x in samples) / n, 2),
         'library_strategies': '; '.join(f'{k}:{v}' for k, v in collections.Counter(x.get('lib', '') for x in samples).items()),
         'technology': ','.join(tech) or 'not_stated',
         'sorted_frac': round(sum(bool(re.search(SORT, b, re.I)) for b in blob) / n, 2),
         'sorting_words_in_protocol': bool(re.search(SORT, alltext, re.I)),
         'blood_frac': round(sum(bool(re.search(r'PBMC|peripheral blood|\bblood\b', b, re.I)) for b in blob) / n, 2),
         'patient_field': pkeys[0] if pkeys else '', 'patients_from_field': patients,
         'sample_sources': '; '.join(f'{k} x{v}' for k, v in collections.Counter(x.get('src', '') for x in samples).most_common(6)),
         'characteristic_keys': '; '.join(sorted({k for x in samples for k in x['chars']}))[:400],
         'n_files': len(files), 'raw_like_files': len(raw_files), 'normalised_files': len(norm_files),
         'example_files': '; '.join(sorted({re.sub(r'GSM\d+_', '', f['name']) for f in files})[:6])[:300],
         **{f'probe_{k}': v for k, v in probe.items()},
         'paper_treatment': ' || '.join(sent.get('treatment', []))[:600],
         'paper_fresh_frozen': ' || '.join(sent.get('fresh_frozen', []))[:600],
         'paper_sorting': ' || '.join(sent.get('sorting', []))[:600],
         'paper_patients': ' || '.join(sent.get('patients', []))[:400]}
    results.append(r)
    json.dump({'series': r, 'samples': [{k: (list(v) if isinstance(v, set) else v) for k, v in x.items() if k != 'text'} for x in samples], 'files': files},
              open(os.path.join(out, f'{acc}.json'), 'w'), indent=1)
    print(acc, r['technology'], 'raw', r['raw_like_files'], 'probe', r.get('probe_probe'), 'sorted', r['sorted_frac'], 'pts', patients, flush=True)
df = pd.DataFrame(results)
df.to_csv(os.path.join(out, 'detailed_check_v1.tsv'), sep='\t', index=False)
