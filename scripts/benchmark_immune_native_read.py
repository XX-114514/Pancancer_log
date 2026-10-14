#!/usr/bin/env python3
"""Measure actual per-library reads against a bounded native CSR source cache."""
import argparse,hashlib,json,resource,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy import sparse
from immune_atlas_io import Source


def benchmark(a):
    config=json.loads(a.contract.read_text())
    record=next(r for r in config['records'] if r['dataset']==a.dataset)
    source=Source(record)
    rows=source.rows
    if len(rows)!=len(np.unique(rows)) or not np.array_equal(rows,np.arange(rows[0],rows[-1]+1)):
        raise ValueError('This diagnostic requires an actual contiguous source selector')
    group=source.handle[record['counts_key']]
    ptr=group['indptr'];nnz=int(ptr[int(rows[-1])+1])-int(ptr[int(rows[0])])
    estimate=nnz*(group['data'].dtype.itemsize+group['indices'].dtype.itemsize)+(len(rows)+1)*8
    if 3*estimate>a.peak_budget_bytes:raise ValueError('Diagnostic source cache exceeds memory budget')
    library=source.obs.library_id.iloc[rows].astype(str).to_numpy()
    positions=np.flatnonzero(library==library[0])
    print(json.dumps(dict(phase='source_cache_read',n_cells=len(rows),n_library_cells=len(positions),estimated_CSR_bytes=estimate)),flush=True)
    started=time.perf_counter()
    cache=sparse.csr_matrix(source.matrix[slice(int(rows[0]),int(rows[-1])+1),:])[:,source.columns].tocsr()
    cache_seconds=time.perf_counter()-started
    if not np.isfinite(cache.data).all() or (cache.data<0).any() or (cache.data!=np.rint(cache.data)).any():raise ValueError('Invalid cached native counts')
    print(json.dumps(dict(phase='original_library_read',cache_seconds=cache_seconds)),flush=True)
    started=time.perf_counter();original=source.raw(positions);original_seconds=time.perf_counter()-started
    started=time.perf_counter();cached=cache[positions].tocsr();cached_seconds=time.perf_counter()-started
    diff=original.astype(np.int64)-cached.astype(np.int64);diff.eliminate_zeros()
    if diff.nnz:raise ValueError('Native count parity differs')
    digest=hashlib.sha256()
    for array in [cached.indptr,cached.indices,cached.data]:digest.update(array.tobytes())
    report=dict(observed_utc=datetime.now(timezone.utc).isoformat(),dataset=a.dataset,
        n_source_cells=len(rows),n_source_genes=cache.shape[1],n_library_cells=len(positions),
        source_cache_seconds=cache_seconds,original_library_read_seconds=original_seconds,
        cached_library_read_seconds=cached_seconds,all_library_counts_exact=True,
        native_raw_counts_not_corrected=True,source_cache_estimated_bytes=estimate,
        maximum_process_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        native_library_count_digest=digest.hexdigest(),biological_inference_not_repeated=True)
    with a.output.open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--contract',type=Path,required=True);p.add_argument('--dataset',required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--peak-budget-bytes',type=int,default=40000000000)
    benchmark(p.parse_args())
