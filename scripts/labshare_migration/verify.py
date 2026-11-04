#!/usr/bin/env python3
"""Compare a copied tree with its source: same file set, same SHA256, same symlink targets.

Usage: verify.py SRC DST [--exclude 'a:b/c'] [--jobs N] --out PREFIX
Excludes are ':'-separated paths relative to SRC (skipped on both sides).
Writes PREFIX.sha256 (sha256  relpath of every regular file, source side) and PREFIX.verify.tsv
(problems only), and exits non-zero on any missing, extra or changed entry.
"""
import argparse
import hashlib
import multiprocessing
import os
import sys


def sha256(path):
    h = hashlib.sha256()
    try:
        with open(path, 'rb') as f:
            for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
                h.update(block)
    except OSError as e:
        return 'ERROR:' + e.__class__.__name__
    return h.hexdigest()


def pair(task):
    rel, src, dst = task
    return rel, sha256(os.path.join(src, rel)), sha256(os.path.join(dst, rel))


def walk(root, excludes):
    files, links = set(), {}
    for d, dirs, names in os.walk(root):
        rd = os.path.relpath(d, root)
        rd = '' if rd == '.' else rd
        keep = []
        for n in dirs:
            rel = os.path.join(rd, n)
            if rel in excludes:
                continue
            if os.path.islink(os.path.join(d, n)):
                links[rel] = os.readlink(os.path.join(d, n))
            else:
                keep.append(n)
        dirs[:] = keep
        for n in names:
            rel = os.path.join(rd, n)
            if rel in excludes:
                continue
            p = os.path.join(d, n)
            if os.path.islink(p):
                links[rel] = os.readlink(p)
            else:
                files.add(rel)
    return files, links


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src')
    ap.add_argument('dst')
    ap.add_argument('--exclude', default='')
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    excludes = {e.strip('/') for e in a.exclude.split(':') if e}
    sf, sl = walk(a.src, excludes)
    df, dl = walk(a.dst, excludes)
    problems = [('missing_in_dst', r) for r in sorted(sf - df)]
    problems += [('extra_in_dst', r) for r in sorted(df - sf)]
    for r in sorted(set(sl) | set(dl)):
        if sl.get(r) != dl.get(r):
            problems.append(('symlink_differs', r))
    common = sorted(sf & df)
    n_ok = 0
    with open(a.out + '.sha256', 'w', encoding='utf-8', errors='surrogateescape') as out, \
            multiprocessing.Pool(a.jobs) as pool:
        for rel, hs, hd in pool.imap(pair, [(r, a.src, a.dst) for r in common], chunksize=64):
            out.write(f'{hs}  {rel}\n')
            if hs != hd or hs.startswith('ERROR'):
                problems.append(('sha_differs', rel))
            else:
                n_ok += 1
    with open(a.out + '.verify.tsv', 'w', encoding='utf-8', errors='surrogateescape') as out:
        out.write('problem\tpath\n')
        for p, r in problems:
            out.write(f'{p}\t{r}\n')
    print(f'files_src {len(sf)} files_dst {len(df)} sha_ok {n_ok} symlinks {len(sl)} problems {len(problems)}')
    sys.exit(1 if problems else 0)


if __name__ == '__main__':
    main()
