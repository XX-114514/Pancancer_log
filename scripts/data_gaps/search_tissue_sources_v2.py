# Candidate search v2: same as search_tissue_sources_v1.py but the targets (analysis_main -> GEO term,
# CELLxGENE disease regex) are read from a JSON file, for the TCGA-coverage expansion (2026-11-01).
# Usage: search_tissue_sources_v2.py TARGETS.json KNOWN_ACCESSIONS.txt MENTIONED_ACCESSIONS.txt OUT_DIR
import sys, os, json, time, re, datetime, urllib.request, urllib.parse

known = set(open(sys.argv[2]).read().split())
mentioned = set(open(sys.argv[3]).read().split())
out = sys.argv[4]
os.makedirs(out, exist_ok=True)
EU = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
SC = '("single cell"[All Fields] OR "single-cell"[All Fields] OR scRNA-seq[All Fields] OR "scRNA seq"[All Fields])'
BASE = '"Homo sapiens"[Organism] AND gse[Entry Type] AND ' + SC
TARGETS = {k: tuple(v) for k, v in json.load(open(sys.argv[1])).items()}
EXCLUDE_TITLE = r'\b(mouse|murine|organoid|cell line|PDX|xenograft|CAR-?T|PBMC|blood only|in vitro)\b'


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except Exception:
            time.sleep(2 ** i)
    raise RuntimeError(url)


log, rows = [], []
for target, (term, _) in TARGETS.items():
    q = f'{BASE} AND {term}'
    ids = json.loads(get(EU + 'esearch.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'term': q, 'retmax': 500, 'retmode': 'json'})))['esearchresult']['idlist']
    log.append({'repository': 'GEO', 'target': target, 'query': q, 'hits': len(ids)})
    time.sleep(0.4)
    for s in range(0, len(ids), 100):
        summ = json.loads(get(EU + 'esummary.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'id': ','.join(ids[s:s + 100]), 'retmode': 'json'})))['result']
        time.sleep(0.4)
        for uid in summ.get('uids', []):
            d = summ[uid]
            acc = d.get('accession', '')
            if 'high throughput sequencing' not in d.get('gdstype', ''):
                continue
            pm = ','.join(str(p) for p in d.get('pubmedids', []))
            flags = []
            if acc in known: flags.append('KNOWN_SOURCE')
            elif acc in mentioned: flags.append('in_project_records')
            if re.search(EXCLUDE_TITLE, d.get('title', '') + ' ' + d.get('summary', '')[:400], re.I): flags.append('title_suggests_nontissue')
            if d.get('taxon', '') != 'Homo sapiens': flags.append('multi_species:' + d.get('taxon', ''))
            rows.append({'repository': 'GEO', 'target': target, 'accession': acc, 'title': d.get('title', ''),
                         'n_samples': d.get('n_samples', ''), 'pubmed': pm, 'date': d.get('pdat', ''),
                         'supp': d.get('suppfile', ''), 'platform': d.get('gpl', ''), 'flags': ';'.join(flags)})

# CELLxGENE curation API: all public collections with their datasets' disease/assay/suspension
cols = json.loads(get('https://api.cellxgene.cziscience.com/curation/v1/collections'))
log.append({'repository': 'CELLxGENE', 'target': 'all', 'query': 'curation/v1/collections', 'hits': len(cols)})
for target, (_, rx) in TARGETS.items():
    for c in cols:
        ds = c.get('datasets', [])
        hit = [d for d in ds if any(re.search(rx, x.get('label', ''), re.I) for x in d.get('disease', []))]
        if not hit:
            continue
        assays = sorted({a['label'] for d in hit for a in d.get('assay', [])})
        susp = sorted({s for d in hit for s in d.get('suspension_type', [])})
        tissues = sorted({t['label'] for d in hit for t in d.get('tissue', [])})
        cid = 'CXG_' + c['collection_id'].replace('-', '_')
        flags = ['KNOWN_SOURCE'] if any(('CXG_' + d.get('dataset_id', '').replace('-', '_')) in known for d in hit) else []
        rows.append({'repository': 'CELLxGENE', 'target': target, 'accession': cid, 'title': c.get('name', ''),
                     'n_samples': sum(d.get('cell_count', 0) for d in hit), 'pubmed': c.get('doi', '') or '',
                     'date': (c.get('published_at') or '')[:10], 'supp': '|'.join(assays),
                     'platform': '|'.join(susp) + ' / ' + '|'.join(tissues)[:150], 'flags': ';'.join(flags)})

# same publication as a known source -> flag (GEO pubmed ids of known sources are looked up once)
known_geo = sorted(a for a in known if a.startswith('GSE'))
kp = {}
for s in range(0, len(known_geo), 20):
    term = '(' + ' OR '.join(f'{a}[Accession]' for a in known_geo[s:s + 20]) + ') AND gse[Entry Type]'
    ids = json.loads(get(EU + 'esearch.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'term': term, 'retmax': 100, 'retmode': 'json'})))['esearchresult']['idlist']
    time.sleep(0.4)
    if ids:
        summ = json.loads(get(EU + 'esummary.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'id': ','.join(ids), 'retmode': 'json'})))['result']
        for uid in summ.get('uids', []):
            for p in summ[uid].get('pubmedids', []):
                kp[str(p)] = summ[uid].get('accession')
        time.sleep(0.4)
for r in rows:
    shared = [kp[p] for p in r['pubmed'].split(',') if p in kp and kp[p] != r['accession']]
    if shared:
        r['flags'] = ';'.join(filter(None, [r['flags'], 'same_PMID_as:' + ','.join(sorted(set(shared)))]))

stamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
json.dump({'run_UTC': stamp, 'queries': log}, open(os.path.join(out, 'search_log.json'), 'w'), indent=1)
cols_out = ['repository', 'target', 'accession', 'title', 'n_samples', 'pubmed', 'date', 'supp', 'platform', 'flags']
with open(os.path.join(out, 'candidates.tsv'), 'w') as f:
    f.write('\t'.join(cols_out) + '\n')
    for r in rows:
        f.write('\t'.join(str(r[c]).replace('\t', ' ').replace('\n', ' ') for c in cols_out) + '\n')
print(stamp, len(rows), 'candidate rows')
for e in log: print(e['repository'], e['target'], e['hits'])
