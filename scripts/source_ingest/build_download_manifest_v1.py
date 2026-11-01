# Build the download manifest for accepted supplementary sources: series-level supplementary files
# from GEO, filtered to raw counts + cell/sample metadata (no TCR/BCR, spatial images, normalised
# matrices or per-sample unfiltered 10x matrices). Sizes come from HTTP HEAD.
# Usage: build_download_manifest_v1.py OUT.tsv
import sys, re, time, urllib.request

KEEP = {  # accession -> regex of series-level files to fetch (decided from the detailed check)
    'GSE163203': r'_RAW\.tar$', 'GSE191288': r'_RAW\.tar$', 'GSE162454': r'_RAW\.tar$',
    'GSE299023': r'_RAW\.tar$', 'GSE205013': r'_RAW\.tar$', 'GSE242230': r'_RAW\.tar$|annotations',
    'GSE149614': r'count\.txt|metadata', 'GSE125449': r'Set[12]_', 'GSE189903': r'(barcodes|features|genes|matrix|[Ii]nfo|meta)[^/]*$',
    'GSE277524': r'_RAW\.tar$', 'GSE160269': r'CD45(neg|pos)_(UMIs|cells)', 'GSE203067': r'(barcodes|features|matrix|README)',
    'GSE215120': r'_RAW\.tar$',
}
SKIP = r'TCR|BCR|contig|normali[sz]ed|raw_feature_bc|spatial|\.jpg|\.png'


def get(url, method='GET'):
    req = urllib.request.Request(url, method=method, headers={'User-Agent': 'pancancer-ingest'})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r if method == 'HEAD' else r.read().decode('utf-8', 'replace')


rows = []
for acc, rx in KEEP.items():
    t = get(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief').replace('\r', '')
    files = [l.split('=', 1)[1].strip() for l in t.splitlines() if l.startswith('!Series_supplementary_file')]
    for f in files:
        name = f.rsplit('/', 1)[-1]
        if not re.search(rx, name) or re.search(SKIP, name, re.I):
            continue
        url = f.replace('ftp://', 'https://')
        try:
            size = int(get(url, 'HEAD').headers.get('Content-Length') or 0)
        except Exception:
            size = 0
        rows.append((acc, name, url, size))
        time.sleep(0.4)
    time.sleep(1)
with open(sys.argv[1], 'w') as f:
    f.write('source\tfile\turl\tbytes\n')
    for r in rows:
        f.write('\t'.join(map(str, r)) + '\n')
tot = {}
for a, _, _, s in rows:
    tot[a] = tot.get(a, 0) + s
for a in KEEP:
    print(a, sum(1 for r in rows if r[0] == a), 'files', round(tot.get(a, 0) / 1e9, 2), 'GB')
print('total', round(sum(tot.values()) / 1e9, 2), 'GB')
