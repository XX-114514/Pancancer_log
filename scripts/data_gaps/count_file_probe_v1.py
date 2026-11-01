# Integer-count probe for candidate GEO series (follow-up to detailed_check_candidates_v1.py,
# whose directory-listing probe sometimes picked filelist.txt). File names come from the GEO
# !Series_/!Sample_supplementary_file fields; only the first few MB of one count file per series
# are streamed and decompressed, so nothing large is downloaded. 10x .h5 and tar archives are
# reported by format, not opened.
# Usage: count_file_probe_v1.py DETAILED_DIR ACC [ACC ...]
import sys, os, re, json, zlib, time, urllib.request

ddir = sys.argv[1]
SKIP = r'filelist|readme|barcodes|features|genes\.tsv|meta|annot|info\.txt|samples\.txt|cells\.txt|cluster|clinical|tcr|bcr|contig|\.jpg|\.png|spatial|positions|scalefactors'
NORM = r'tpm|fpkm|rpkm|normali[sz]ed|scaled|logcounts|cpm'


def stream(url, nbytes=8_000_000):
    req = urllib.request.Request(url.replace('ftp://', 'https://'), headers={'User-Agent': 'pancancer-audit'})
    with urllib.request.urlopen(req, timeout=180) as r:
        size = r.headers.get('Content-Length')
        raw = r.read(nbytes)
    if url.endswith('.gz'):
        raw = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(raw)
    return raw.decode('utf-8', 'replace'), int(size) if size else None


def files_of(acc):
    t = urllib.request.urlopen(urllib.request.Request(
        f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=all&form=text&view=brief',
        headers={'User-Agent': 'pancancer-audit'}), timeout=180).read().decode('utf-8', 'replace')
    return [l.split('=', 1)[1].strip() for l in t.replace('\r', '').splitlines()
            if re.match(r'!(Series|Sample)_supplementary_file', l) and 'NONE' not in l]


res = {}
for acc in sys.argv[2:]:
    fs = files_of(acc); time.sleep(1)
    count_like = [f for f in fs if not re.search(SKIP, f.rsplit('/', 1)[-1], re.I)]
    norm_only = all(re.search(NORM, f, re.I) for f in count_like) if count_like else False
    pick = next((f for f in count_like if re.search(r'matrix\.mtx', f, re.I) and not re.search(NORM, f, re.I)), None) \
        or next((f for f in count_like if re.search(r'count|umi|dge|_mtx|exprs|expression|processed', f, re.I)
                 and not re.search(NORM, f, re.I) and re.search(r'\.(txt|csv|tsv)(\.gz)?$', f, re.I)), None)
    r = {'n_supp_files': len(fs), 'n_count_like': len(count_like), 'normalised_only': norm_only,
         'h5_files': sum(f.endswith('.h5') for f in fs), 'archives': sum(bool(re.search(r'\.tar(\.gz)?$', f)) for f in fs)}
    if pick:
        try:
            txt, size = stream(pick)
            vals = []
            mtx = 'mtx' in pick.lower()
            header_done = False
            for line in txt.splitlines()[:-1]:
                if line.startswith('%'):
                    continue
                parts = re.split(r'[,\t ]+', line.strip())
                if mtx and not header_done:
                    header_done = True; continue
                for p in (parts[2:3] if mtx else parts[1:]):
                    try: vals.append(float(p.strip('"')))
                    except ValueError: pass
                if len(vals) > 200000: break
            nonint = sum(abs(v - round(v)) > 1e-6 for v in vals)
            r.update(probe_file=pick.rsplit('/', 1)[-1], file_bytes=size, values_checked=len(vals), non_integer=nonint,
                     max_value=max(vals) if vals else None,
                     verdict='integer_counts' if vals and nonint == 0 else ('NON_INTEGER' if vals else 'no_values_parsed'))
        except Exception as e:
            r.update(probe_file=pick.rsplit('/', 1)[-1], verdict='probe_error', error=repr(e)[:120])
    elif r['h5_files']:
        r['verdict'] = 'h5_matrices (10x h5 stores integer counts by format; not opened)'
    elif r['archives']:
        r['verdict'] = 'archives_only (open after download)'
    else:
        r['verdict'] = 'no_count_file'
    res[acc] = r
    print(acc, json.dumps(r), flush=True)
    time.sleep(1)
json.dump(res, open(os.path.join(ddir, 'count_probe_v1.json'), 'w'), indent=1)
