"""Bounded local transport pilot of the same selected official ZIP member.

No decompressed cell rows, coordinates or clinical data are parsed/persisted.
Default TLS and official same-origin HTTP ranges remain mandatory. This does
not restart the remote process and is not real-data scientific acceptance.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import time
import traceback
import zipfile


def sha(blob):
    return hashlib.sha256(blob).hexdigest()


def load(path, expected, name):
    assert sha(path.read_bytes()) == expected
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    out = Path(os.environ['PANCANCER_LOCAL_TRANSPORT_PILOT_OUTPUT'])
    assert not out.exists()
    out.mkdir(mode=0o700, parents=True)
    def save(name, value):
        (out/name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    save('plan_before_local_pilot.json', dict(started_UTC=datetime.now(timezone.utc).isoformat(),
        source_SHA256=sha(Path(__file__).read_bytes()), selected_compressed_prefix_bytes=64*1024,
        strict_default_TLS_and_official_same_origin=True, byte_budget_including_retries=8*1024**2,
        range_start_deadline_seconds=180, per_request_deadline_seconds=45, maximum_transient_tries=3,
        no_whole_member_CRC_or_science_claim=True, no_cell_rows_or_individual_arrays_written=True))
    source = load(Path(__file__).parent/'preflight_melanoma_IMC_public_archive_v2.py',
        '9fb200d98e321584fde8bca6cb92f72151ac33ec00d64e04d93e7e2b33e870e2', 'local_source')
    transport = load(Path(__file__).parent/'melanoma_IMC_local_bounded_transport_v2.py',
        'b16f8537602a75529a9df48eb7b1dd7f26d73a04e4146178dbc6c03dec6c4d21', 'local_ranges')
    record = json.loads(source.fetch(source.BASE+'api/records/6838169', 2*1024**2))
    item = next(x for x in record['files'] if x['key'] == 'Melanoma_code_data.zip')
    assert item['size'] == 668227028 and item['checksum'] == 'md5:da74b7228763b76dcd2d9830c1a10cd9'
    class PilotRange(transport.RangeFile):
        def _read_chunk(self, n):
            assert time.monotonic()-self.started <= 180
            assert self.requested+n*3 <= 8*1024**2
            return super()._read_chunk(n)
    rf = PilotRange(source, item['links']['self'], item['size'],
        lambda value: save('private_pilot_progress.json', value))
    def signature(entries):
        return sha(json.dumps([(x.filename, x.file_size, x.compress_size, x.CRC) for x in entries],
            separators=(',', ':')).encode())
    with zipfile.ZipFile(rf) as archive:
        directory = signature(archive.infolist())
        assert directory == '26339da50b39ffdf7beb641a67f50523c867f2d3a86b24f42175f3a16dfe8992'
        info = archive.getinfo('data/input/sc_data.csv')
        assert info.file_size == 172707928 and info.compress_size == 51454934
        with archive.open(info):
            body_start = time.monotonic()
            compressed_prefix = rf.read(64*1024)
            body_seconds = time.monotonic()-body_start
            assert len(compressed_prefix) == 64*1024
    with zipfile.ZipFile(rf) as replay:
        assert signature(replay.infolist()) == directory
    final = json.loads(source.fetch(source.BASE+'api/records/6838169', 2*1024**2))
    again = next(x for x in final['files'] if x['key'] == item['key'])
    assert all(item.get(k) == again.get(k) for k in ['id', 'size', 'checksum', 'links'])
    result = dict(status='COMPLETE_LOCAL_STRICT_COMPRESSED_PREFIX_TRANSPORT_PILOT_ONLY',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(Path(__file__).read_bytes()),
        directory_SHA256=directory, compressed_prefix_bytes=len(compressed_prefix), body_seconds=body_seconds,
        selected_payload_byte_rate=len(compressed_prefix)/body_seconds,
        successful_range_bytes=rf.transferred, requested_bytes_including_retries=rf.requested,
        strict_default_TLS=True, published_file_identity_and_directory_replayed=True,
        no_complete_member_CRC_or_scientific_acceptance=True, no_uncompressed_rows_parsed_or_written=True)
    rf.close()
    save('private_local_pilot_complete.json', result)
    print(json.dumps(result))


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        out = Path(os.environ['PANCANCER_LOCAL_TRANSPORT_PILOT_OUTPUT'])
        out.mkdir(mode=0o700, parents=True, exist_ok=True)
        result = dict(status='LOCAL_PILOT_FAILED_NOT_ACCEPTED',
            observed_UTC=datetime.now(timezone.utc).isoformat(), error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name, line=x.lineno, function=x.name)
                for x in traceback.extract_tb(error.__traceback__)], no_raw_payload_or_record_values=True)
        (out/'private_local_pilot_failed.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(result))
        raise SystemExit(1)
