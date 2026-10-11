"""Manual, frozen-whitelist cleanup; default is read-only preview.

Run on the reviewed server with its private manifest. --apply is interactive.
GPU duplicates become symlinks to the verified retained copies, preserving
old paths. They then share the retained file and must not be modified in place.
This tool never deletes directories, installed environments or source datasets.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import sys


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8388608), b''):
            digest.update(block)
    return digest.hexdigest()


def verify(row, roots):
    path = Path(row['path'])
    assert not path.is_symlink() and path.is_file(), 'Missing or symlink: ' + str(path)
    resolved = path.resolve()
    assert any(resolved.is_relative_to(root) for root in roots), 'Outside whitelist root'
    stat = path.stat()
    assert (stat.st_dev, stat.st_ino, stat.st_size, stat.st_nlink) == (
        row['dev'], row['inode'], row['bytes'], 1), 'File changed or hard linked: ' + str(path)
    assert sha(path) == row['sha256'], 'SHA changed: ' + str(path)
    if row.get('extracted_directory'):
        assert Path(row['extracted_directory']).is_dir(), 'Unpacked package absent'
    if row.get('keep'):
        kept = Path(row['keep'])
        assert kept.is_file() and not kept.is_symlink()
        assert kept.resolve() != resolved and sha(kept) == row['keep_sha256']
    if row['group'] == 'package_archives':
        assert path.name.endswith(('.conda', '.tar.bz2'))
    if row['group'] == 'failed_chunks':
        assert path.name.endswith('.bin') and 'range_chunk_attempts_v2' in path.parts
    return path


def check_processes(paths, groups):
    assert Path('/proc').exists(), 'Run on the reviewed Linux server'
    blocked = []
    for process in Path('/proc').iterdir():
        if not process.name.isdigit():
            continue
        try:
            if process.stat().st_uid != os.getuid():
                continue
            argv = (process / 'cmdline').read_bytes().split(b'\0')
            basenames = [Path(value.decode(errors='replace')).name for value in argv if value]
            if ('package_archives' in groups and
                    any(name in ['conda', 'mamba', 'micromamba', 'pip', 'pip3'] for name in basenames)):
                blocked.append((process.name, 'package manager is active'))
            if any(group.startswith('failed_') for group in groups):
                if any(name.startswith('fetch_server_priors') for name in basenames):
                    blocked.append((process.name, 'resource worker is active'))
            for descriptor in (process / 'fd').iterdir():
                try:
                    target = str(Path(os.readlink(descriptor)).resolve())
                    if target in paths:
                        blocked.append((process.name, 'candidate/retained file is open'))
                except FileNotFoundError:
                    pass
            maps = process / 'maps'
            if maps.exists():
                for line in maps.read_text().splitlines():
                    parts = line.split(None, 5)
                    if len(parts) == 6 and parts[5] in paths:
                        blocked.append((process.name, 'candidate/retained file is mapped'))
        except (FileNotFoundError, ProcessLookupError):
            pass
    assert not blocked, 'Abort: active process evidence: ' + repr(blocked)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--groups', nargs='+', choices=[
        'package_archives', 'failed_chunks', 'failed_tf_partial', 'gpu_duplicates'],
        default=['package_archives', 'failed_chunks', 'failed_tf_partial'])
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    manifest_path = Path(args.manifest)
    assert sha(manifest_path) == args.manifest_sha256, 'Manifest SHA mismatch'
    plan = json.loads(manifest_path.read_text())
    assert plan['schema'] == 'manual_cleanup_reviewed_plan_v1'
    roots = [Path(value).resolve() for value in plan['allowed_roots']]
    chosen = [row for row in plan['files'] if row['group'] in args.groups]
    assert chosen and len({row['path'] for row in chosen}) == len(chosen)
    for group in args.groups:
        rows = [row for row in chosen if row['group'] == group]
        size = sum(row['allocated_bytes'] for row in rows) / 1024**3
        print('{}: {} files, estimated {:.3f} GiB'.format(group, len(rows), size))
        for row in rows[:3]:
            print('  ' + row['path'])
    print('Checking frozen file SHA, retained copies and active processes...', flush=True)
    for row in chosen:
        verify(row, roots)
    paths = {str(Path(row['path']).resolve()) for row in chosen}
    paths.update(str(Path(row['keep']).resolve()) for row in chosen if row.get('keep'))
    check_processes(paths, args.groups)
    if not args.apply:
        print('DRY RUN: checks passed; no file deleted or replaced.')
        return
    assert sys.stdin.isatty(), '--apply requires interactive manual confirmation'
    for group in args.groups:
        if any(row['group'] == group for row in chosen):
            token = 'APPLY ' + group
            answer = input('To execute {}, type {}: '.format(group, token))
            assert answer == token, 'Cancelled; no deletion has started'
    check_processes(paths, args.groups)
    log = Path('/tmp') / ('pancancer_manual_cleanup_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S') + '.jsonl')
    with log.open('x') as handle:
        handle.write(json.dumps(dict(event='MANUAL_USER_APPLY_STARTED',
                                     UTC=datetime.now(timezone.utc).isoformat(),
                                     manifest_sha256=args.manifest_sha256, groups=args.groups)) + '\n')
        handle.flush()
        for row in chosen:
            path = verify(row, roots)
            handle.write(json.dumps(dict(event='PRE_OPERATION', row=row)) + '\n')
            handle.flush()
            if row['group'] == 'gpu_duplicates':
                temporary = path.with_name(path.name + '.manual_dedup_link_' + str(os.getpid()))
                assert not temporary.exists() and not temporary.is_symlink()
                temporary.symlink_to(os.path.relpath(row['keep'], path.parent))
                # Atomic replacement discards only the verified redundant copy.
                # The retained file is never renamed, written or deleted.
                os.replace(temporary, path)
                assert path.is_symlink() and sha(path) == row['keep_sha256']
                operation = 'REDUNDANT_COPY_REPLACED_WITH_SYMLINK'
            else:
                path.unlink()
                operation = 'WHITELIST_CACHE_FILE_REMOVED'
            handle.write(json.dumps(dict(event=operation, path=str(path),
                                         UTC=datetime.now(timezone.utc).isoformat())) + '\n')
            handle.flush()
    print('Finished manual operations. Execution receipt: ' + str(log))
    print('Run quota -sg to check actual released quota; observed estimates are not guarantees.')


if __name__ == '__main__':
    main()
