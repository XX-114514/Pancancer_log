#!/usr/bin/env python3
"""Rewrite an absolute path prefix inside a copied tree of conda envs / R libraries.

Usage: relocate_prefix.py ROOT OLD NEW [--dry-run] --report TSV
NEW must not be longer than OLD. Same rules as conda's own prefix replacement:
- text files (no NUL byte anywhere): every OLD is replaced by NEW;
- binary files: each NUL-terminated string that contains OLD is rewritten in place and padded
  with NULs to its original length, so offsets in ELF and other binaries stay valid;
- symlinks whose target contains OLD are re-pointed;
- .pyc files are left alone (Python only checks the source mtime; the stale path only shows in
  tracebacks), as are files without OLD.
Permissions and mtimes are kept. The report lists every changed file and how.
"""
import argparse
import os
import re
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root')
    ap.add_argument('old')
    ap.add_argument('new')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--report', required=True)
    a = ap.parse_args()
    old, new = a.old.encode(), a.new.encode()
    if len(new) > len(old):
        sys.exit('NEW is longer than OLD; binary replacement impossible')
    pad = len(old) - len(new)
    # a NUL-terminated string starting at OLD; every OLD inside it is replaced, NULs pad it back
    binpat = re.compile(re.escape(old) + rb'[^\0]*\0')
    counts = {'text': 0, 'binary': 0, 'binary_unterminated_left': 0, 'symlink': 0, 'pyc_skipped': 0}
    with open(a.report, 'w') as rep:
        rep.write('kind\tpath\n')
        for d, dirs, names in os.walk(a.root):
            for n in dirs + names:
                p = os.path.join(d, n)
                if os.path.islink(p):
                    t = os.readlink(p)
                    if a.old in t:
                        counts['symlink'] += 1
                        rep.write(f'symlink\t{p}\n')
                        if not a.dry_run:
                            os.unlink(p)
                            os.symlink(t.replace(a.old, a.new), p)
                    continue
                if n in dirs or not os.path.isfile(p):
                    continue
                with open(p, 'rb') as f:
                    data = f.read()
                if old not in data:
                    continue
                if n.endswith('.pyc'):
                    counts['pyc_skipped'] += 1
                    continue
                if b'\0' in data:
                    kind = 'binary'
                    out = binpat.sub(lambda m: m.group(0).replace(old, new) + b'\0' * (m.group(0).count(old) * pad), data)
                    assert len(out) == len(data)
                    if old in out:
                        kind = 'binary_unterminated_left'
                else:
                    kind = 'text'
                    out = data.replace(old, new)
                counts[kind] += 1
                rep.write(f'{kind}\t{p}\n')
                if not a.dry_run:
                    st = os.stat(p)
                    if not os.access(p, os.W_OK):
                        os.chmod(p, st.st_mode | 0o200)
                    with open(p, 'r+b') as f:
                        f.write(out)
                        f.truncate()
                    os.chmod(p, st.st_mode)
                    os.utime(p, ns=(st.st_atime_ns, st.st_mtime_ns))
    print(' '.join(f'{k} {v}' for k, v in counts.items()))


if __name__ == '__main__':
    main()
