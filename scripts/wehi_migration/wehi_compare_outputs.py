#!/usr/bin/env python3
"""Compare a WEHI replay output tree with the migrated old-server output tree, file by file.

Usage: wehi_compare_outputs.py ORIG_DIR REPLAY_DIR OUT_TSV [--rtol R] [--atol A]
TSV and TSV.GZ files are compared numerically (same shape, labels and non-numeric cells; numeric
cells within rtol/atol; gzip headers carry a timestamp); NPY/NPZ and H5/HDF5 array by array (same
names, shapes, dtypes; numeric values within rtol/atol, others exact); JSON files key by key
(numbers within rtol/atol, listing the keys that differ: run reports carry timestamps, paths and
SHA256 of non-identical files); other files by SHA256. RDS differences are reported, not judged:
RDS embeds the R session.
"""
import argparse, hashlib, json, os, sys
import numpy as np, pandas as pd

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def _num_close(x, y, rtol, atol):
    return isinstance(x, (int, float)) and isinstance(y, (int, float)) and not isinstance(x, bool) \
        and not isinstance(y, bool) and np.isclose(x, y, rtol=rtol, atol=atol)

def diff_keys(x, y, path='', rtol=0.0, atol=0.0):
    if isinstance(x, dict) and isinstance(y, dict):
        out = []
        for k in sorted(set(x) | set(y)):
            out += diff_keys(x.get(k), y.get(k), f'{path}.{k}' if path else str(k), rtol, atol)
        return out
    if isinstance(x, list) and isinstance(y, list) and len(x) == len(y):
        out = []
        for i, (u, v) in enumerate(zip(x, y)):
            out += diff_keys(u, v, f'{path}[{i}]', rtol, atol)
        return out
    return [] if (x == y or _num_close(x, y, rtol, atol)) else [path or '<root>']

def cmp_arrays(xs, ys, rtol, atol):
    """xs, ys: dict name -> ndarray."""
    if set(xs) != set(ys):
        return 'DIFF_ARRAY_NAMES', f'{sorted(set(xs) ^ set(ys))[:6]}'
    worst = 0.0
    for k in sorted(xs):
        x, y = xs[k], ys[k]
        if x.shape != y.shape or x.dtype != y.dtype:
            return 'DIFF_SHAPE', f'{k}: {x.shape}/{x.dtype} vs {y.shape}/{y.dtype}'
        if np.issubdtype(x.dtype, np.number) and not np.issubdtype(x.dtype, np.bool_):
            xf, yf = x.astype(float), y.astype(float)
            if not np.array_equal(np.isnan(xf), np.isnan(yf)):
                return 'DIFF_NA', f'{k}: NaN pattern differs'
            m = np.isfinite(xf)
            if m.any():
                d = float(np.abs(xf[m] - yf[m]).max()); worst = max(worst, d)
                if not np.allclose(xf[m], yf[m], rtol=rtol, atol=atol):
                    return 'DIFF_VALUES', f'{k}: max_abs_diff={d:.3g}'
        elif not np.array_equal(x, y):
            return 'DIFF_TEXT', f'{k}: non-numeric values differ'
    return ('IDENTICAL', 'array values equal') if worst == 0 else ('WITHIN_TOL', f'max_abs_diff={worst:.3g}')

def load_np(p):
    if p.endswith('.npy'):
        return {'': np.load(p, allow_pickle=False)}
    with np.load(p, allow_pickle=False) as f:
        return {k: f[k] for k in f.files}

def load_h5(p):
    import h5py
    out = {}
    with h5py.File(p, 'r') as f:
        f.visititems(lambda name, obj: out.__setitem__(name, obj[()]) if isinstance(obj, h5py.Dataset) else None)
    return {k: np.asarray(v) for k, v in out.items()}

def cmp_json(a, b, rtol=0.0, atol=0.0):
    d = diff_keys(json.load(open(a)), json.load(open(b)), '', rtol, atol)
    return ('IDENTICAL', 'json values equal') if not d else ('DIFF_JSON_KEYS', ';'.join(d[:12]) + (f';+{len(d)-12}' if len(d) > 12 else ''))

def cmp_tsv(a, b, rtol, atol):
    x = pd.read_csv(a, sep='\t', dtype=str, keep_default_na=False)
    y = pd.read_csv(b, sep='\t', dtype=str, keep_default_na=False)
    if x.shape != y.shape or list(x.columns) != list(y.columns):
        return 'DIFF_SHAPE', f'{x.shape} vs {y.shape}'
    xn = x.apply(pd.to_numeric, errors='coerce'); yn = y.apply(pd.to_numeric, errors='coerce')
    num = xn.notna() & yn.notna()
    if not ((x == y) | num).all().all():
        return 'DIFF_TEXT', 'non-numeric cells differ'
    if not (xn.isna() == yn.isna()).all().all():
        return 'DIFF_NA', 'numeric/NA pattern differs'
    xv, yv = xn.to_numpy(float)[num.to_numpy()], yn.to_numpy(float)[num.to_numpy()]
    if xv.size == 0:
        return 'IDENTICAL', 'no numeric cells'
    d = np.abs(xv - yv); m = float(d.max())
    if m == 0:
        return 'IDENTICAL', 'max_abs_diff=0'
    ok = np.allclose(xv, yv, rtol=rtol, atol=atol)
    return ('WITHIN_TOL' if ok else 'DIFF_VALUES'), f'max_abs_diff={m:.3g}'

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('orig'); ap.add_argument('replay'); ap.add_argument('out')
    ap.add_argument('--rtol', type=float, default=1e-9); ap.add_argument('--atol', type=float, default=1e-12)
    a = ap.parse_args()
    rows = []
    for root, _, files in os.walk(a.orig):
        for fn in sorted(files):
            p = os.path.join(root, fn); rel = os.path.relpath(p, a.orig); q = os.path.join(a.replay, rel)
            if not os.path.exists(q):
                rows.append((rel, 'MISSING_IN_REPLAY', '')); continue
            if sha(p) == sha(q):
                rows.append((rel, 'IDENTICAL', 'sha256')); continue
            if fn.endswith(('.tsv', '.tsv.gz')):
                rows.append((rel,) + cmp_tsv(p, q, a.rtol, a.atol))
            elif fn.endswith('.json'):
                rows.append((rel,) + cmp_json(p, q, a.rtol, a.atol))
            elif fn.endswith(('.npy', '.npz')):
                rows.append((rel,) + cmp_arrays(load_np(p), load_np(q), a.rtol, a.atol))
            elif fn.endswith(('.h5', '.hdf5')):
                rows.append((rel,) + cmp_arrays(load_h5(p), load_h5(q), a.rtol, a.atol))
            else:
                rows.append((rel, 'DIFF_BYTES', 'not judged (json/rds/other)'))
    pd.DataFrame(rows, columns=['file', 'result', 'detail']).to_csv(a.out, sep='\t', index=False)
    s = pd.Series([r[1] for r in rows]).value_counts()
    print(s.to_string())
    tsv_bad = [r for r in rows if r[0].endswith(('.tsv', '.tsv.gz', '.npy', '.npz', '.h5', '.hdf5')) and r[1] not in ('IDENTICAL', 'WITHIN_TOL')]
    sys.exit(1 if tsv_bad else 0)

if __name__ == '__main__':
    main()
