"""Finite official-source acquisition on the research server, then verified ready.

Use a process-local direct opener: the earlier temporary loopback proxy may be
stale. No persistent proxy settings or SSH configuration are changed.
"""
from pathlib import Path
import json, os, sys, urllib.request, traceback, hashlib
from datetime import datetime, timezone
import anndata as ad
from selective_public_zip_v3 import download_selected
from native_coordinate_contract_v1 import native_micron_coordinates

root = Path(sys.argv[1]).resolve()
config = json.loads((root / 'config/cohort_PRIVATE.json').read_text())
source = json.loads((root / 'config/direct_acquisition_PRIVATE.json').read_text())
out = root / 'acquisition_direct_v1'
out.mkdir(exist_ok=False)
urllib.request.install_opener(urllib.request.build_opener(urllib.request.ProxyHandler({})))


def event(**data):
    with (out / 'events_PRIVATE.jsonl').open('a') as f:
        f.write(json.dumps({'UTC': datetime.now(timezone.utc).isoformat(), **data}) + '\n')


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def checked_link(path, target, expected):
    assert path.stat().st_size == expected['bytes']
    assert sha256(path) == expected['sha256']
    if target.exists():
        assert target.stat().st_size == expected['bytes'] and sha256(target) == expected['sha256']
    else:
        os.link(path, target)


try:
    event(stage='mapping', status='OFFICIAL_PUBLIC_DIRECT_RANGE_STARTED')
    expected_maps = config['mapping_receipt']['members']
    receipt = download_selected(source['mapping_url'], [v['member'] for v in expected_maps],
                                out / 'mapping', source['mapping_budget_bytes'])
    by_member = {v['member']: v for v in receipt['members']}
    for expected in expected_maps:
        accepted = by_member[expected['member']]
        assert accepted['member_CRC_validated']
        checked_link(out / 'mapping' / accepted['local_name'],
                     root / 'inputs_direct_v1/mapping' / expected['local_name'], expected)
    with (root / 'inputs_direct_v1/mapping/complete_PRIVATE.json').open('x') as f:
        json.dump(config['mapping_receipt'], f, indent=2)
    event(stage='mapping', status='ALL18_SERVER_CRC_SIZE_AND_INDEPENDENT_SHA_ACCEPTED')
    for case in config['cases']:
        unit = case['unit']
        try:
            event(stage=unit, status='OFFICIAL_PUBLIC_NATIVE_RANGE_STARTED')
            receipt = download_selected(source['native_url'], [case['expected_member']],
                                        out / unit, source['native_budget_bytes'][unit])
            accepted = receipt['members'][0]
            assert accepted['member_CRC_validated']
            expected = case['receipt']['members'][0]
            path = out / unit / accepted['local_name']
            target = root / 'inputs_direct_v1' / unit
            checked_link(path, target / expected['local_name'], expected)
            data = ad.read_h5ad(path, backed='r')
            _, coordinates = native_micron_coordinates(data)
            data.file.close()
            with (target / 'complete_PRIVATE.json').open('x') as f:
                json.dump(case['receipt'], f, indent=2)
            with (target / 'ready_PRIVATE.json').open('x') as f:
                json.dump({'UTC': datetime.now(timezone.utc).isoformat(),
                           'member_CRC_and_server_size_SHA_match': True,
                           'coordinate_contract': coordinates}, f, indent=2)
            event(stage=unit, status='NATIVE_AND_COORDINATE_CONTRACT_ACCEPTED_READY')
        except Exception as error:
            (out / (unit + '_failure_PRIVATE.txt')).write_text(traceback.format_exc())
            event(stage=unit, status='FAILED_PARTIAL_RETAINED', exception_type=type(error).__name__)
    event(stage='all6', status='FINITE_SOURCE_LIST_ENDED_CHECK_INDIVIDUAL_READY_RECEIPTS')
except Exception as error:
    (out / 'mapping_failure_PRIVATE.txt').write_text(traceback.format_exc())
    event(stage='mapping', status='FAILED_PARTIAL_RETAINED', exception_type=type(error).__name__)
    raise
