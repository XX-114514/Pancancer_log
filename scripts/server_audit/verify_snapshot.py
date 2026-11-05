"""Verify public-copy hashes and numerical claims without server data or a network.

Usage: python scripts/server_audit/verify_snapshot.py [SNAPSHOT_DIRECTORY]
This does not assert the validity of the underlying biological models.
"""
import csv
import hashlib
import json
import math
import sys
from pathlib import Path


def rows(p):
    with p.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def main():
    repo = Path(__file__).resolve().parents[2]
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else repo / 'inventories/server_analysis/2026-11-05'
    checked = 0
    for manifest in ('source_evidence.tsv', 'supplementary_evidence.tsv'):
        for r in rows(root / manifest):
            if len(r['published_sha256']) == 64:
                p = root / r['published_copy']
                require(p.is_file(), 'missing published copy: ' + str(p))
                require(hashlib.sha256(p.read_bytes()).hexdigest() == r['published_sha256'], 'changed published evidence: ' + str(p))
                checked += 1
    code = rows(root / 'imported_code.tsv')
    require(len(code) == 43, 'unexpected source code scope')
    for r in code:
        p = repo / r['published_path']
        require(hashlib.sha256(p.read_bytes()).hexdigest() == r['published_sha256'], 'changed imported code: ' + str(p))
    for p in root.rglob('*.tsv'):
        with p.open() as f:
            lines = list(csv.reader(f, delimiter='\t'))
        require(bool(lines) and all(len(r) == len(lines[0]) for r in lines), 'ragged TSV: ' + str(p))
    tasks = rows(root / 'nmf_tasks.tsv')
    require(len(tasks) == 45, 'expected 45 NMF tasks')
    require(len({(r['method'], r['seed']) for r in tasks}) == 45, 'duplicate NMF task')
    require(all(int(r['K']) == 9 for r in tasks), 'NMF controlled K changed')
    require(len({r['method'] for r in tasks}) == 9, 'expected nine methods')
    for name in ('major_v5_to_v7_transition.tsv', 'malignancy_v3_to_v7_transition.tsv'):
        p = root / 'tables/v5_v7' / name
        rr = rows(p)
        total = sum(int(v) for r in rr for k, v in r.items() if k != list(r)[0])
        require(total == 4676787, 'transition table universe mismatch: ' + name)
    pilot = json.loads((root / 'annotation_pilots.json').read_text())
    routed = pilot['20261103_marker_routed_annotation/GSE138709_v2/summary.json']
    require(sum(routed['route_counts'].values()) == routed['source_cells'] == 11381, 'route count mismatch')
    require(routed['accepted_fine_identity'] is False, 'pilot became final identity')
    for br in ('primary_strict', 'legacy_all_display'):
        cfg = json.loads((root / 'tables/v7_cm' / ('covarnet_' + br + '.json')).read_text())
        require(cfg['rank_override'] == 9, 'K9 override no longer explicit')
        rank = rows(root / 'tables/v7_cm' / br / 'rank_survey_v1.tsv')
        require([int(r['rank']) for r in rank] == list(range(2, 21)), 'incomplete rank survey')
    complete = rows(root / 'spatial_completion.tsv')
    require(all(int(r['summaries_present']) == 29 for r in complete), 'spatial summary missing')
    scale = rows(root / 'nmf_scale_sensitivity.tsv')
    require(all(abs(float(r['equivalent_rescaling_canonical_cos']) - 1) < 1e-10 for r in scale), 'canonical gauge is not invariant')
    require(any(float(r['equivalent_rescaling_original_cos']) < .9 for r in scale), 'scale counterexample not reproduced')
    # Independent exact binomial arithmetic for the reported directional counts.
    ps = {k: sum(math.comb(18, i) for i in range(k, 19)) / 2**18 for k in (9, 11, 15, 16)}
    require(abs(ps[15] - .0037689208984375) < 1e-14, 'sign calculation mismatch')
    require(abs(ps[16] - .0006561279296875) < 1e-14, 'sign calculation mismatch')
    print(json.dumps(dict(status='PASS', published_evidence_hashes=checked,
        imported_code_hashes=len(code), nmf_tasks=len(tasks), transition_cells=4676787,
        exact_sign_p=ps, biological_validation_claimed=False), indent=2))


if __name__ == '__main__':
    main()
