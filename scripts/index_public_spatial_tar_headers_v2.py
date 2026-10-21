"""Index a public spatial archive by HTTP ranges, without member payloads.

Only tar/PAX/long-name metadata is read. Never extracts archive paths to disk.
Both a manual checksum/number parser and stdlib TarInfo validate every header.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import re
import tarfile
import urllib.parse
import urllib.request


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def number(value):
    if value[0] in (0x80, 0xff):
        return int.from_bytes(value[1:], 'big', signed=value[0] == 0xff)
    return int(value.strip(b'\x00 ') or b'0', 8)


def role(name):
    low = name.lower()
    for text, label in [('cell_feature_matrix', 'cell_feature_matrix'), ('gene_panel', 'gene_panel'), ('features.tsv', 'feature_axis'), ('experiment.xenium', 'experiment_metadata'), ('cells.csv', 'cell_coordinates_and_metrics'), ('cells.parquet', 'cell_coordinates_and_metrics'), ('cell_boundaries', 'cell_boundaries'), ('nucleus_boundaries', 'nucleus_boundaries'), ('transcripts', 'transcript_coordinates'), ('morphology', 'morphology_image'), ('metrics', 'metrics'), ('analysis', 'analysis_bundle')]:
        if text in low:
            return label
    return 'other_or_nested_archive'


def main():
    url = os.environ['PANCANCER_ARCHIVE_URL']
    assert urllib.parse.urlparse(url).scheme == 'https' and urllib.parse.urlparse(url).hostname == 'ftp.ncbi.nlm.nih.gov'
    out = Path(os.environ['PANCANCER_ARCHIVE_OUTPUT'])
    assert not out.exists()
    out.mkdir(mode=0o700)
    save = lambda n, v: (out / n).write_text(json.dumps(v, indent=2) + '\n')
    save('plan_before_tar_index.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), url=url, max_members=2000, max_HTTP_requests=4000, max_metadata_bytes=4*1024**2, regular_file_payloads_never_read=True, no_archive_paths_extracted=True, exact_ranges_and_source_ETag_or_Last_Modified_pinned=True, every_header_manual_and_stdlib_checked=True))
    etag, modified, total, transferred, requests = None, None, None, 0, 0
    def get(start, size):
        nonlocal etag, modified, total, transferred, requests
        assert size <= 65536 and transferred + size <= 4*1024**2 and requests < 4000
        headers = {'Range': f'bytes={start}-{start+size-1}', 'Accept-Encoding': 'identity'}
        if etag:
            headers['If-Match'] = etag
        elif modified:
            headers['If-Unmodified-Since'] = modified
        request = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(request, timeout=25) as response:
            assert response.status == 206 and urllib.parse.urlparse(response.geturl()).hostname == 'ftp.ncbi.nlm.nih.gov'
            match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', response.headers.get('Content-Range', ''))
            assert match and (int(match[1]), int(match[2])) == (start, start+size-1)
            if total is None:
                total, etag = int(match[3]), response.headers.get('ETag')
                modified = response.headers.get('Last-Modified')
                assert etag or modified
            assert total == int(match[3]) and etag == response.headers.get('ETag') and modified == response.headers.get('Last-Modified')
            raw = response.read(size+1)
        assert len(raw) == size
        transferred += len(raw)
        requests += 1
        return raw
    members, header_hashes = [], []
    offset, longname, pax, complete = 0, None, {}, False
    while len(members) < 2000:
        b = get(offset, 512)
        if b == b'\x00'*512:
            assert get(offset+512, 512) == b
            complete = True
            break
        assert sum(b[:148])+8*32+sum(b[156:]) == number(b[148:156])
        standard = tarfile.TarInfo.frombuf(b, 'utf-8', 'surrogateescape')
        size = number(b[124:136])
        assert size == standard.size and size >= 0
        name = b[:100].split(b'\x00')[0].decode('utf-8', 'surrogateescape')
        prefix = b[345:500].split(b'\x00')[0].decode('utf-8', 'surrogateescape')
        if prefix and b[257:265] == b'ustar\x0000':
            name = prefix+'/'+name
        if standard.isdir():
            name = name.rstrip('/')
        assert name == standard.name
        header_hashes.append(dict(offset=offset, SHA256=hashlib.sha256(b).hexdigest()))
        kind = b[156:157]
        if kind in [b'x', b'g', b'L', b'K']:
            assert size <= 65536
            metadata = get(offset+512, size)
            if kind == b'L':
                longname = metadata.rstrip(b'\x00').decode('utf-8', 'surrogateescape')
            elif kind in [b'x', b'g']:
                cursor = 0
                while cursor < len(metadata):
                    end = metadata.index(b' ', cursor)
                    length = int(metadata[cursor:end])
                    record = metadata[end+1:cursor+length-1]
                    key, value = record.split(b'=', 1)
                    if key in [b'path', b'size']:
                        pax[key.decode()] = value.decode('utf-8', 'surrogateescape')
                    cursor += length
                assert cursor == size
        else:
            name = pax.get('path', longname or name)
            size = int(pax.get('size', size))
            assert size >= 0 and offset+512+size <= total
            members.append(dict(name=name, type=kind.decode('ascii') if kind else '', header_offset=offset, payload_offset=offset+512, bytes=size, role=role(name), header_SHA256=hashlib.sha256(b).hexdigest()))
            longname, pax = None, {}
            temporary = out / 'private_index_checkpoint.tmp'
            temporary.write_text(json.dumps(dict(scanned_members=len(members), HTTP_requests=requests, metadata_bytes=transferred, next_offset=offset+512+((size+511)//512)*512), indent=2)+'\n')
            temporary.replace(out / 'private_index_checkpoint.json')
        offset += 512+((size+511)//512)*512
    assert complete
    assert hashlib.sha256(get(0, 512)).hexdigest() == header_hashes[0]['SHA256']
    save('private_complete_tar_member_index.json', dict(url=url, source_ETag=etag, source_Last_Modified=modified, archive_bytes=total, metadata_bytes_read=transferred, HTTP_requests=requests, all_regular_file_payloads_skipped=True, complete_two_zero_block_terminator=True, members=members, headers=header_hashes))
    roles = {k: sum(m['role'] == k for m in members) for k in sorted({m['role'] for m in members})}
    save('private_terminal.json', dict(status='COMPLETE_PUBLIC_TAR_HEADER_INDEX_NO_REGULAR_MEMBER_PAYLOAD', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), member_index_SHA256=sha(out / 'private_complete_tar_member_index.json'), archive_bytes=total, indexed_members=len(members), HTTP_requests=requests, metadata_bytes_read=transferred, member_roles=roles, no_expression_cell_coordinates_TCR_or_clinical_payload_read=True, all_archive_paths_data_only=True))
    print(json.dumps(dict(status='COMPLETE_PUBLIC_TAR_HEADER_INDEX_NO_REGULAR_MEMBER_PAYLOAD', indexed_members=len(members), metadata_bytes_read=transferred, member_roles=roles)))


if __name__ == '__main__':
    main()
