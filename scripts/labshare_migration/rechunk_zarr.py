#!/usr/bin/env python3
"""Rewrite a zarr v2 store (e.g. an AnnData .zarr) with larger chunks, to cut its file count.

Usage: rechunk_zarr.py SRC DST [--target-mib 16] [--min-chunks 64] [--jobs N] --report TSV
Every array with more than --min-chunks chunks (object dtype excepted) gets chunks of about
--target-mib MiB uncompressed along axis 0 (never smaller than before; other axes whole).
Dtype, compressor, filters, fill value, order and attributes are kept. All other arrays and
group/attribute files are copied byte for byte. Metadata are consolidated as in the source.
Then every rewritten array is compared value by value with the source (NaN == NaN), and every
copied file by SHA256. Exits non-zero on any difference. SRC is only read.
"""
import argparse
import hashlib
import json
import multiprocessing
import os
import shutil
import sys

import numpy as np
import zarr
from numcodecs import get_codec

BLOCK_CHUNKS = 8  # new chunks per copy/compare task


def arrays_and_groups(root):
    arrays, others = [], []
    for d, dirs, names in os.walk(root):
        rel = os.path.relpath(d, root)
        if '.zarray' in names:
            arrays.append(rel)
            dirs[:] = []
            continue
        others += [os.path.normpath(os.path.join(rel, n)) for n in names]
    return sorted(arrays), sorted(others)


def plan(meta, target, min_chunks):
    shape, chunks = meta['shape'], meta['chunks']
    if not shape or meta['dtype'] in ('|O', '<O') or meta['dtype'].startswith('|O'):
        return None
    n = 1
    for s, c in zip(shape, chunks):
        n *= -(-s // c)
    if n <= min_chunks:
        return None
    itemsize = np.dtype(meta['dtype']).itemsize
    other = 1
    for s in shape[1:]:
        other *= s
    rows = max(chunks[0], target // (itemsize * other))
    return [min(rows, shape[0])] + list(shape[1:])


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def copy_block(task):
    src, dst, a, b = task
    s = zarr.open_array(src, mode='r')
    d = zarr.open_array(dst, mode='r+')
    d[a:b] = s[a:b]
    return b - a


def check_block(task):
    src, dst, a, b = task
    x = zarr.open_array(src, mode='r')[a:b]
    y = zarr.open_array(dst, mode='r')[a:b]
    eq_nan = x.dtype.kind in 'fc'
    return (src, a, b, bool(np.array_equal(x, y, equal_nan=eq_nan)))


def check_file(task):
    s, d = task
    return (s, sha(s) == sha(d))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src')
    ap.add_argument('dst')
    ap.add_argument('--target-mib', type=float, default=16)
    ap.add_argument('--min-chunks', type=int, default=64)
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--report', required=True)
    a = ap.parse_args()
    if os.path.exists(a.dst):
        sys.exit(f'exists: {a.dst}')
    target = int(a.target_mib * 1024 * 1024)
    arrays, others = arrays_and_groups(a.src)
    copied_files, tasks, rows = [], [], []
    for rel in others:
        s, d = os.path.join(a.src, rel), os.path.join(a.dst, rel)
        os.makedirs(os.path.dirname(d), exist_ok=True)
        if os.path.basename(rel) == '.zmetadata':
            continue
        shutil.copy2(s, d)
        copied_files.append((s, d))
    for rel in arrays:
        s, d = os.path.join(a.src, rel), os.path.join(a.dst, rel)
        meta = json.load(open(os.path.join(s, '.zarray')))
        new = plan(meta, target, a.min_chunks)
        if new is None:
            shutil.copytree(s, d)
            for n in os.listdir(s):
                copied_files.append((os.path.join(s, n), os.path.join(d, n)))
            rows.append((rel, meta['chunks'], meta['chunks'], 'copied'))
            continue
        comp = get_codec(meta['compressor']) if meta['compressor'] else None
        filt = [get_codec(f) for f in meta['filters']] if meta['filters'] else None
        extra = {'dimension_separator': meta['dimension_separator']} if 'dimension_separator' in meta else {}
        zarr.open_array(d, mode='w', shape=meta['shape'], chunks=new, dtype=meta['dtype'],
                        compressor=comp, filters=filt, fill_value=meta['fill_value'],
                        order=meta['order'], **extra)
        if os.path.exists(os.path.join(s, '.zattrs')):
            shutil.copy2(os.path.join(s, '.zattrs'), os.path.join(d, '.zattrs'))
            copied_files.append((os.path.join(s, '.zattrs'), os.path.join(d, '.zattrs')))
        step = new[0] * BLOCK_CHUNKS
        tasks += [(s, d, i, min(i + step, meta['shape'][0])) for i in range(0, meta['shape'][0], step)]
        rows.append((rel, meta['chunks'], new, 'rechunked'))
    with multiprocessing.Pool(a.jobs) as pool:
        for _ in pool.imap_unordered(copy_block, tasks):
            pass
        print(f'written {len(tasks)} blocks', flush=True)
        if os.path.exists(os.path.join(a.src, '.zmetadata')):
            zarr.consolidate_metadata(zarr.DirectoryStore(a.dst))
        bad = [r for r in pool.imap_unordered(check_block, tasks) if not r[3]]
        bad += [(s, None, None, False) for s, ok in pool.imap_unordered(check_file, copied_files) if not ok]
    # consolidated metadata must match the source apart from the chunk sizes
    if os.path.exists(os.path.join(a.src, '.zmetadata')):
        ms = json.load(open(os.path.join(a.src, '.zmetadata')))['metadata']
        md = json.load(open(os.path.join(a.dst, '.zmetadata')))['metadata']
        for k in ms:
            x, y = dict(ms[k]), dict(md.get(k, {}))
            if k.endswith('.zarray'):
                x.pop('chunks', None)
                y.pop('chunks', None)
            if x != y:
                bad.append((k, None, None, False))
        bad += [(k, None, None, False) for k in md if k not in ms]
    n_src = sum(len(f) for _, _, f in os.walk(a.src))
    n_dst = sum(len(f) for _, _, f in os.walk(a.dst))
    with open(a.report, 'w') as out:
        out.write('array\told_chunks\tnew_chunks\taction\n')
        for r in rows:
            out.write('\t'.join(map(str, r)) + '\n')
        for b in bad:
            out.write(f'MISMATCH\t{b[0]}\t{b[1]}\t{b[2]}\n')
    print(f'files_src {n_src} files_dst {n_dst} arrays {len(arrays)} '
          f'rechunked {sum(r[3] == "rechunked" for r in rows)} blocks {len(tasks)} mismatches {len(bad)}')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
