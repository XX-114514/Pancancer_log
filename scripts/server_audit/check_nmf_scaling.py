"""Read existing W matrices and audit factor-scale sensitivity; no NMF is fitted.

Usage: python check_nmf_scaling.py RUN REFERENCE_W OUTPUT_TSV
OUTPUT_TSV must not exist. Canonicalisation makes every nonzero W column sum to
one before applying the original row-membership cosine and Hungarian matching.
This is a sensitivity audit, not a replacement factorisation or validation set.
"""
import itertools
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


def member(w, canonical=False):
    x = w.to_numpy(float).copy()
    if canonical:
        x /= np.maximum(x.sum(axis=0), 1e-300)
    x /= np.maximum(x.sum(axis=1, keepdims=True), 1e-300)
    return x / np.maximum(np.linalg.norm(x, axis=0), 1e-300)


def compare(a, b, canonical):
    c = member(a, canonical).T @ member(b, canonical)
    i, j = linear_sum_assignment(-c)
    return c[i, j].mean(), tuple(j)


def main():
    run, refpath, output = map(Path, sys.argv[1:])
    ref = pd.read_csv(refpath, sep='\t', index_col=0)
    records = []
    for method in sorted((run / 'methods').iterdir()):
        files = sorted(method.glob('seed_*/W.tsv'))
        ws = [pd.read_csv(p, sep='\t', index_col=0).loc[ref.index] for p in files]
        if not ws:
            continue
        old = [compare(ref, w, False) for w in ws]
        new = [compare(ref, w, True) for w in ws]
        oldstab = [compare(a, b, False)[0] for a, b in itertools.combinations(ws, 2)]
        newstab = [compare(a, b, True)[0] for a, b in itertools.combinations(ws, 2)]
        # Equivalent gauge: W*D, H/D represents exactly the same fitted matrix.
        scaled = ws[0] * np.logspace(-2, 2, ws[0].shape[1])
        records.append(dict(method=method.name, seeds=len(ws),
            original_ref_cos=np.mean([r[0] for r in old]),
            canonical_ref_cos=np.mean([r[0] for r in new]),
            changed_assignments=sum(a[1] != b[1] for a, b in zip(old, new)),
            original_seed_cos=np.mean(oldstab), canonical_seed_cos=np.mean(newstab),
            equivalent_rescaling_original_cos=compare(ws[0], scaled, False)[0],
            equivalent_rescaling_canonical_cos=compare(ws[0], scaled, True)[0]))
    with output.open('x') as f:
        pd.DataFrame(records).to_csv(f, sep='\t', index=False)
    print(pd.DataFrame(records).round(4).to_string(index=False))


if __name__ == '__main__':
    main()
