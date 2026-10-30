#!/usr/bin/env python3
"""Old server: resolve the run directories the current TME working layer depends on (M2).

Starting from seed run names, read every JSON under each run's config/, provenance/ and the top
level of outputs/, collect referenced run directories under PROJECT_V5/runs (and any other
absolute paths), and recurse. Writes a TSV of run name, depth reached, size in MB (if a size
table is given) and the run that first referenced it, plus a list of referenced paths outside
the runs tree. Read-only; Python 3.6 compatible.

Usage: old_m2_closure.py PROJECT_V5 OUT_PREFIX [--sizes du_runs.txt] SEED_RUN [SEED_RUN ...]
"""
import argparse
import collections
import os
import re
import sys

PATH_RE = re.compile(r'(/(?:home|data[0-9])/[^"\'\s,\]\}]+)')


def run_jsons(run_dir):
    for sub in ('config', 'provenance'):
        top = os.path.join(run_dir, sub)
        for root, _, files in os.walk(top):
            for f in files:
                if f.endswith('.json'):
                    yield os.path.join(root, f)
    out = os.path.join(run_dir, 'outputs')
    if os.path.isdir(out):
        for f in os.listdir(out):
            if f.endswith('.json'):
                yield os.path.join(out, f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('project_v5')
    ap.add_argument('out_prefix')
    ap.add_argument('--sizes')
    ap.add_argument('--max-json-bytes', type=int, default=64 * 1024 ** 2)
    ap.add_argument('seeds', nargs='+')
    a = ap.parse_args()
    runs_root = os.path.join(os.path.realpath(a.project_v5), 'runs')
    # Runs are reached through more than one mount alias (home symlink vs data disk path).
    run_re = re.compile(r'Project_v5/runs/([^/"\'\s]+)')

    sizes = {}
    if a.sizes:
        for line in open(a.sizes):
            parts = line.split()
            if len(parts) == 2:
                sizes[parts[1].rstrip('/')] = int(parts[0])

    parent = {s: '<seed>' for s in a.seeds}
    depth = {s: 0 for s in a.seeds}
    external = collections.Counter()
    queue = collections.deque(a.seeds)
    while queue:
        run = queue.popleft()
        run_dir = os.path.join(runs_root, run)
        if not os.path.isdir(run_dir):
            print('missing run dir: %s' % run, file=sys.stderr)
            continue
        for path in run_jsons(run_dir):
            try:
                if os.path.getsize(path) > a.max_json_bytes:
                    print('skipped large json: %s' % path, file=sys.stderr)
                    continue
                text = open(path, encoding='utf-8', errors='replace').read()
            except OSError as e:
                print('unreadable: %s (%s)' % (path, e), file=sys.stderr)
                continue
            for m in PATH_RE.finditer(text):
                p = m.group(1)
                r = run_re.search(p)
                if r:
                    ref = r.group(1)
                    if ref not in parent:
                        parent[ref] = run
                        depth[ref] = depth[run] + 1
                        queue.append(ref)
                else:
                    external[p] += 1

    with open(a.out_prefix + '.runs.tsv', 'w') as f:
        f.write('run\tdepth\tsize_MB\tfirst_referenced_by\n')
        for run in sorted(parent, key=lambda r: (depth[r], r)):
            f.write('%s\t%d\t%s\t%s\n' % (run, depth[run], sizes.get(run, ''), parent[run]))
    with open(a.out_prefix + '.external.tsv', 'w') as f:
        for p, n in external.most_common():
            f.write('%d\t%s\n' % (n, p))
    total = sum(sizes.get(r, 0) for r in parent)
    print('runs %d, known size %d MB, external paths %d' % (len(parent), total, len(external)))


if __name__ == '__main__':
    main()
