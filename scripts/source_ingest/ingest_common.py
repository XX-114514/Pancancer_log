# Shared readers / builders for supplementary-source ingest (factored out of ingest_supplement_v1.py
# unchanged, except that build() takes the output directory). Used by ingest_supplement_v2.py.
import sys, os, re, io, tarfile, gzip, json
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io, anndata as ad, h5py



def collapse(X, genes):
    """sum columns with the same symbol; X cells x genes"""
    genes = np.asarray(genes).astype(str)
    u, inv = np.unique(genes, return_inverse=True)
    if len(u) == len(genes):
        order = np.argsort(genes)
        return X.tocsc()[:, order].tocsr(), genes[order]
    agg = sp.csr_matrix((np.ones(len(genes)), (np.arange(len(genes)), inv)), shape=(len(genes), len(u)))
    return (X @ agg).tocsr(), u


def as_counts(X):
    X = sp.csr_matrix(X)
    if X.nnz and not np.allclose(X.data[:200000], np.round(X.data[:200000])):
        raise ValueError('non-integer counts')
    X.data = np.rint(X.data)
    return X.astype(np.int32)


def tar_members(path):
    t = tarfile.open(path)
    return t, {m.name: m for m in t.getmembers()}


def read_mm(fobj):
    """MatrixMarket coordinate reader tolerant of blank '%' comment lines (Singleron); returns cells x genes"""
    head = None
    while True:
        line = fobj.readline().decode().strip()
        if line.startswith('%') or not line:
            continue
        head = [int(x) for x in line.split()]; break
    nr, nc, nnz = head
    t = pd.read_csv(fobj, sep=r'\s+', header=None, dtype={0: np.int32, 1: np.int32, 2: np.float64}, engine='c')
    assert len(t) == nnz, (len(t), nnz)
    return sp.csr_matrix((t[2].values, (t[1].values - 1, t[0].values - 1)), shape=(nc, nr))


def floor_call(X, bc):
    """cell calling for unfiltered matrices: keep barcodes passing the QC floor (>=200 genes, >=500 UMI)"""
    X = sp.csr_matrix(X)
    ok = (np.diff(X.indptr) >= 200) & (np.asarray(X.sum(1)).ravel() >= 500)
    return X[ok], bc[ok]


def read_mtx_triplet(open_fn, prefix):
    # file objects from tar members / open() are not decompressed automatically
    gz = lambda name: gzip.open(open_fn(name))
    mtx = read_mm(gz(prefix + 'matrix.mtx.gz'))
    bc = pd.read_csv(gz(prefix + 'barcodes.tsv.gz'), header=None, sep='\t')[0].astype(str).values
    ft = None
    for fname in ('features.tsv.gz', 'genes.tsv.gz'):
        try:
            ft = pd.read_csv(gz(prefix + fname), header=None, sep='\t'); break
        except (KeyError, FileNotFoundError):
            continue
    if ft.shape[1] >= 3:
        ft = ft[ft[2].astype(str).str.contains('Gene Expression')] if (ft[2] != 'Gene Expression').any() else ft
        mtx = mtx[:, ft.index.values]
    sym = ft[1] if ft.shape[1] >= 2 else ft[0]
    return mtx, bc, sym.astype(str).values


def read_10x_h5(fobj):
    with h5py.File(fobj, 'r') as h:
        g = h['matrix']
        X = sp.csc_matrix((g['data'][:], g['indices'][:], g['indptr'][:]), shape=tuple(g['shape'][:])).T.tocsr()
        bc = g['barcodes'][:].astype(str)
        f = g['features']
        keep = f['feature_type'][:].astype(str) == 'Gene Expression'
        return X[:, keep], bc, f['name'][:].astype(str)[keep]


def dense_to_sparse(path_or_buf, sep='\t', index_col=0, chunksize=1000):
    rows, genes = [], []
    cells = None
    for ch in pd.read_csv(path_or_buf, sep=sep, index_col=index_col, chunksize=chunksize):
        cells = ch.columns.astype(str).values
        genes.extend(ch.index.astype(str))
        rows.append(sp.csr_matrix(ch.values))
    M = sp.vstack(rows).T.tocsr()
    return M, cells, np.asarray(genes)


def build(src, parts, cancer_main, extra_uns, out):
    """parts: list of (X, barcodes, genes, obs-dict)"""
    allg = sorted(set().union(*[set(p[2]) for p in parts]))
    gi = {g: i for i, g in enumerate(allg)}
    mats, obs = [], []
    for X, bc, genes, meta in parts:
        X, genes = collapse(as_counts(X), genes)
        cols = np.array([gi[g] for g in genes])
        X = sp.csr_matrix((X.data, cols[X.indices], X.indptr), shape=(X.shape[0], len(allg)))
        mats.append(X)
        o = pd.DataFrame({'origin_barcode': bc})
        for k, v in meta.items():
            o[k] = v
        obs.append(o)
    X = sp.vstack(mats).tocsr()
    o = pd.concat(obs, ignore_index=True)
    o['source'] = src
    o['cell_id'] = src + '|' + o['sample_id'].astype(str) + '|' + o['origin_barcode'].astype(str)
    if o.cell_id.duplicated().any():
        raise ValueError(f'{src}: duplicated cell ids')
    o.index = o.cell_id.values
    if 'analysis_main' not in o:
        o['analysis_main'] = cancer_main
    a = ad.AnnData(X=X, obs=o, var=pd.DataFrame(index=allg))
    a.uns['ingest'] = json.dumps(extra_uns)
    a.write_h5ad(os.path.join(out, f'{src}.h5ad'), compression='gzip')
    s = o.groupby('sample_id').size()
    print(json.dumps({'source': src, 'cells': int(a.n_obs), 'genes': int(a.n_vars), 'samples': int(len(s)),
                      'patients': int(o.patient_id.nunique())}), flush=True)


def common(**kw):
    base = {'treatment_status': 'unknown', 'treatment_detail': '', 'fresh_frozen': 'unknown',
            'capture_scope': 'whole_tissue', 'suspension': 'cell', 'cell_calling': 'author_filtered_matrix'}
    base.update(kw)
    return base


