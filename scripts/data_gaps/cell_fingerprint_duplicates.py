# Read-only screen for identical cells across sources (gap 3: cross-study aliases / reused samples).
# Each cell gets two gene-order-independent fingerprints of its raw count vector:
#   full: hash of the sorted multiset of non-zero counts; top50: hash of its 50 largest counts.
# Cells sharing a fingerprint but coming from different sources are candidate duplicates.
# top50 tolerates the symbol summation in the core store, which can alter small counts.
# Usage: cell_fingerprint_duplicates.py SOURCES_FILE OUT_DIR STORE=PATH:SOURCE_FIELD:MATRIX ...
#   MATRIX is X or layers/<name>; matrices must be CSR with integer-valued counts.
import sys, os, hashlib, json
import numpy as np, pandas as pd, zarr
from anndata.io import read_elem

sources = set(open(sys.argv[1]).read().split())
out = sys.argv[2]; os.makedirs(out, exist_ok=True)
MIN_TOTAL, MIN_NNZ, CHUNK = 500, 200, 20000
parts = []
for spec in sys.argv[3:]:
    name, rest = spec.split('=', 1)
    path, field, mat = rest.split(':')
    g = zarr.open(path, 'r')
    src = np.asarray(read_elem(g['obs'][field])).astype(str)
    keep = np.isin(src, list(sources))
    m = g[mat]
    assert m.attrs.get('encoding-type') == 'csr_matrix', (name, m.attrs.asdict())
    indptr, data = m['indptr'][:], m['data']
    rows = np.flatnonzero(keep)
    print(name, 'rows', len(rows), 'of', len(src), flush=True)
    recs = []
    # contiguous blocks of kept rows, read in chunks
    for s in range(0, len(rows), CHUNK):
        r = rows[s:s + CHUNK]
        a, b = indptr[r[0]], indptr[r[-1] + 1]
        blk = data[a:b]
        if not np.allclose(blk[:1000], np.round(blk[:1000])):
            raise SystemExit(f'{name}: {mat} is not integer counts')
        blk = np.rint(blk).astype(np.int64)
        for i in r:
            v = np.sort(blk[indptr[i] - a:indptr[i + 1] - a])
            tot, nnz = int(v.sum()), len(v)
            if tot < MIN_TOTAL or nnz < MIN_NNZ:
                continue
            recs.append((name, src[i], int(i), tot, nnz,
                         hashlib.blake2b(v.tobytes(), digest_size=12).hexdigest(),
                         hashlib.blake2b(v[-50:].tobytes(), digest_size=12).hexdigest()))
        print(name, s + len(r), '/', len(rows), flush=True)
    parts.append(pd.DataFrame(recs, columns=['store', 'source', 'row', 'total', 'nnz', 'fp_full', 'fp_top50']))
df = pd.concat(parts, ignore_index=True)
df.to_csv(os.path.join(out, 'cell_fingerprints.tsv.gz'), sep='\t', index=False)
summary = {'n_cells_fingerprinted': len(df), 'min_total': MIN_TOTAL, 'min_nnz': MIN_NNZ}
for fp in ('fp_full', 'fp_top50'):
    multi = df.groupby(fp)['source'].nunique()
    hit = df[df[fp].isin(multi[multi > 1].index)]
    pairs = (hit.groupby(fp)['source'].apply(lambda s: '|'.join(sorted(set(s))))
                .value_counts().rename('n_shared_fingerprints'))
    pairs.to_csv(os.path.join(out, f'cross_source_pairs_{fp}.tsv'), sep='\t')
    within = df.groupby(['source', fp]).size()
    summary[fp] = {'cells_in_cross_source_groups': int(len(hit)),
                   'source_pairs': int(len(pairs)),
                   'within_source_duplicate_cells': int((within[within > 1]).sum())}
json.dump(summary, open(os.path.join(out, 'summary.json'), 'w'), indent=1)
print(json.dumps(summary, indent=1))
