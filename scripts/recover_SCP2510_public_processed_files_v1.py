"""Read public SCP files through the same anonymous GCS scope as its viewer.

The page-provided short-lived public read token stays in memory. No token or
signed download URL is written to receipts. Matrices remain outside public Git.
"""
import argparse
import base64
import hashlib
import json
import re
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path


def main(a):
    inventory = json.loads(a.inventory.read_text())
    explore = json.loads(a.explore.read_text())
    assert inventory['public'] is True and explore['isPublic'] is True
    page = a.page.read_text()
    match = re.search(r'window\.SCP\.readOnlyTokenObject\s*=\s*', page)
    assert match, 'Public viewer read scope was not provided'
    token = json.JSONDecoder().raw_decode(page[match.end():])[0]['access_token']
    wanted = [f for f in inventory['study_files'] if f['name'].startswith('spatial_Pt')
              or f['name'] == 'single_cell_portal/metadata.txt.gz'
              or f['name'].startswith('single_cell_portal/counts/')]
    assert len(wanted) == 8
    assert sum(f['upload_file_size'] for f in wanted) < 256 * 1024**2
    assert not a.output.exists()
    a.output.mkdir(mode=0o700)

    def fetch(f):
        location = f['bucket_location']
        assert location.startswith('single_cell_portal/') and '..' not in location
        url = ('https://storage.googleapis.com/download/storage/v1/b/'
               + urllib.parse.quote(explore['bucketId'], safe='') + '/o/'
               + urllib.parse.quote(location, safe='') + '?alt=media')
        destination = a.output / location
        destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        expected = int(f['upload_file_size'])
        request = urllib.request.Request(url, headers={
            'Authorization': 'Bearer ' + token, 'Accept-Encoding': 'identity',
            'User-Agent': 'CM02-public-processed-data-recovery/1'})
        md5 = hashlib.md5()
        sha = hashlib.sha256()
        written = 0
        with urllib.request.urlopen(request, timeout=90) as response, destination.open('xb') as sink:
            assert response.status == 200
            assert int(response.headers['Content-Length']) == expected
            official_hash = response.headers.get('x-goog-hash', '')
            for chunk in iter(lambda: response.read(1024**2), b''):
                written += len(chunk)
                assert written <= expected
                md5.update(chunk)
                sha.update(chunk)
                sink.write(chunk)
        assert written == expected
        official_md5 = next((x.strip().split('=', 1)[1] for x in official_hash.split(',')
                             if x.strip().startswith('md5=')), None)
        assert official_md5 and base64.b64decode(official_md5).hex() == md5.hexdigest()
        return dict(name=f['name'], relative_path=location, bytes=written,
                    sha256=sha.hexdigest(), official_md5=md5.hexdigest(),
                    source='SCP2510 public viewer file inventory',
                    MD5_and_length_verified=True)

    receipts = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        for future in as_completed([pool.submit(fetch, f) for f in wanted]):
            row = future.result()
            receipts.append(row)
            print(json.dumps(row), flush=True)
    receipt = dict(status='ALL8_PUBLIC_PROCESSED_FILES_OFFICIAL_MD5_LENGTH_ACCEPTED',
                   finished_utc=datetime.now(timezone.utc).isoformat(),
                   source='https://singlecell.broadinstitute.org/single_cell/study/SCP2510',
                   total_bytes=sum(x['bytes'] for x in receipts), files=receipts,
                   credentials_or_signed_urls_not_recorded=True,
                   numeric_segmentation_and_clinical_qualification_pending=True)
    (a.output / 'download_technical_public.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for field in ['inventory', 'explore', 'page', 'output']:
        p.add_argument('--' + field, type=Path, required=True)
    main(p.parse_args())
