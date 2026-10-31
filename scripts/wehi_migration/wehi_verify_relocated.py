#!/usr/bin/env python3
"""WEHI: re-check every migrated file at its new place after wehi_organize.py.

Reads each item's old-server per-file SHA256 list (migration/items/<item>/files.sha256, paths
relative to the directory the item was unpacked into), maps every path through
migration/relocation.tsv (longest old-path prefix wins) and recomputes SHA256 in parallel.

Usage: wehi_verify_relocated.py PC_ROOT DATA_ROOT [--jobs N]
Writes migration/verify_relocated.tsv and exits non-zero on any missing or changed file.
"""
import argparse
import hashlib
import multiprocessing
import os
import sys

UNPACK_DIR = {'m4a_inputs': 'm4', 'm4b_runs': 'm4', 'm4c_env': 'm4'}


def sha256(path):
    h = hashlib.sha256()
    try:
        with open(path, 'rb') as f:
            for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
                h.update(block)
    except OSError as e:
        return 'ERROR:' + e.__class__.__name__
    return h.hexdigest()


def check(task):
    want, path = task
    got = sha256(path)
    return want == got, path, got


def read_list(path):
    for line in open(path, encoding='utf-8', errors='surrogateescape'):
        line = line.rstrip('\n')
        escaped = line.startswith('\\')
        if escaped:
            line = line[1:]
        digest, rel = line[:64], line[66:]
        if escaped:
            rel = rel.replace('\\n', '\n').replace('\\\\', '\\')
        yield digest, rel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pc_root')
    ap.add_argument('data_root')
    ap.add_argument('--jobs', type=int, default=os.cpu_count())
    a = ap.parse_args()
    pc, data = os.path.abspath(a.pc_root), os.path.abspath(a.data_root)
    moves = []
    for line in open(os.path.join(pc, 'migration', 'relocation.tsv')):
        if line.startswith('old_path\t'):
            continue
        old, new = line.rstrip('\n').split('\t')
        moves.append((old, new))
    moves.sort(key=lambda m: len(m[0]), reverse=True)

    def relocate(old):
        for o, n in moves:
            if old == o or old.startswith(o + '/'):
                return n + old[len(o):]
        return None

    report = open(os.path.join(pc, 'migration', 'verify_relocated.tsv'), 'w')
    report.write('item\tn_files\tn_ok\tn_unmapped\tn_bad\n')
    failed = False
    items_dir = os.path.join(pc, 'migration', 'items')
    with multiprocessing.Pool(a.jobs) as pool:
        for item in sorted(os.listdir(items_dir)):
            lst = os.path.join(items_dir, item, 'files.sha256')
            if not os.path.isfile(lst):
                continue
            base = os.path.join(data, UNPACK_DIR.get(item, ''))
            tasks, unmapped = [], []
            for digest, rel in read_list(lst):
                new = relocate(os.path.normpath(os.path.join(base, rel)))
                (tasks.append((digest, new)) if new else unmapped.append(rel))
            bad = [r for r in pool.imap_unordered(check, tasks, chunksize=64) if not r[0]]
            report.write('%s\t%d\t%d\t%d\t%d\n' % (item, len(tasks) + len(unmapped),
                                                   len(tasks) - len(bad), len(unmapped), len(bad)))
            report.flush()
            for _, path, got in bad[:20]:
                print('%s BAD %s %s' % (item, got, path), file=sys.stderr)
            for rel in unmapped[:20]:
                print('%s UNMAPPED %s' % (item, rel), file=sys.stderr)
            failed |= bool(bad or unmapped)
    report.close()
    print(open(os.path.join(pc, 'migration', 'verify_relocated.tsv')).read())
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
