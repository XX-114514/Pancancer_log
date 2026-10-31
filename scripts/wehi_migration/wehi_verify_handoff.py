# Read-only: check the source-server handoff manifest and its dependencies on WEHI via PATH_MAP.
import json, sys, hashlib, datetime, os
from pathlib import Path
src_dir, path_map, out = map(Path, sys.argv[1:4])
pm = [l.rstrip('\n').split('\t') for l in open(path_map)][1:]
pm = sorted([(a, b) for a, b, *_ in pm if b.startswith('/')], key=lambda x: -len(x[0]))
def wehi(p):
    for a, b in pm:
        if p == a or p.startswith(a + '/'):
            return b + p[len(a):]
    return None
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 22), b''): h.update(blk)
    return h.hexdigest()
man = json.load(open(src_dir / 'Claude_sync_file_manifest.json'))
rec = json.load(open(src_dir / 'Claude_manifest_sha_dependency_recheck_20261031.json'))
srcsha = {i['path']: i.get('sha256') for i in rec['items']}
srcsize = {i['path']: i.get('size') for i in rec['items']}
for i in rec['items']:
    for d in i.get('deps', []):
        if d.get('resolved'):
            srcsha.setdefault(d['resolved'], d.get('sha256')); srcsize.setdefault(d['resolved'], d.get('size'))
targets = [(it['server_path'], 'item', it['group']) for it in man['items']]
seen = {t[0] for t in targets}
for i in rec['items']:
    for d in i.get('deps', []):
        p = d.get('resolved') or d['ref']
        if p not in seen: seen.add(p); targets.append((p, 'dep', Path(i['path']).name))
rows = []
for p, kind, grp in targets:
    w = wehi(p); r = {'source_path': p, 'kind': kind, 'group': grp, 'wehi_path': w,
                     'source_sha256': srcsha.get(p), 'source_size': srcsize.get(p)}
    if w is None: r['status'] = 'no_path_map'
    elif not os.path.exists(w): r['status'] = 'missing_on_WEHI'
    elif os.path.isdir(w): r['status'] = 'dir_present'
    else:
        r['wehi_size'] = os.path.getsize(w)
        if r['source_sha256'] in (None, 'skipped_large'):
            r['status'] = 'present_size_only' + ('' if r['source_size'] in (None, r['wehi_size']) else '_SIZE_MISMATCH')
        else:
            r['wehi_sha256'] = sha(w); r['status'] = 'sha_match' if r['wehi_sha256'] == r['source_sha256'] else 'SHA_MISMATCH'
    rows.append(r)
res = {'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'host': os.uname().nodename,
       'source_manifest_sha256': sha(src_dir / 'Claude_sync_file_manifest.json'),
       'source_recheck_sha256': sha(src_dir / 'Claude_manifest_sha_dependency_recheck_20261031.json'),
       'n_targets': len(rows), 'status_counts': {}, 'rows': rows}
for r in rows: res['status_counts'][r['status']] = res['status_counts'].get(r['status'], 0) + 1
json.dump(res, open(out, 'w'), indent=1, ensure_ascii=False)
print(json.dumps(res['status_counts']))
for r in rows:
    if r['status'] not in ('sha_match',): print(r['status'], r['kind'], r['source_path'])
