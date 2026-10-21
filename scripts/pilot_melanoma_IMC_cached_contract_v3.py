"""Same fixed source version, avoiding an optional metadata API dependency.

Only bounded compressed prefix bytes stay in RAM. Two official same-origin
file routes are eligible; the second is tried only after a transient network
failure. TLS/HTTP status/identity failures do not trigger route fallback.
The accepted immutable directory digest anchors the cached source contract.
Full scientific input still requires the previously accepted full-member SHA.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import http.client
import importlib.util
import json
import os
import time
import traceback
import zipfile


def sha(value):
    return hashlib.sha256(value).hexdigest()


def main():
    out = Path(os.environ['PANCANCER_LOCAL_TRANSPORT_PILOT_OUTPUT'])
    assert not out.exists()
    out.mkdir(mode=0o700, parents=True)
    def save(name, value):
        (out/name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    path = Path(__file__).parent/'melanoma_IMC_coalesced_ranges_v2.py'
    assert sha(path.read_bytes()) == '65547a1952c7c209c4179836bc7e05074143d19530f6e800bfdcb7c4eb90a500'
    spec = importlib.util.spec_from_file_location('cached_contract_ranges', path)
    transport = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(transport)
    urls = ['https://zenodo.org/api/records/6838169/files/Melanoma_code_data.zip/content',
        'https://zenodo.org/records/6838169/files/Melanoma_code_data.zip?download=1']
    directory_sha = '26339da50b39ffdf7beb641a67f50523c867f2d3a86b24f42175f3a16dfe8992'
    save('plan_before_cached_contract_pilot.json', dict(started_UTC=datetime.now(timezone.utc).isoformat(),
        source_SHA256=sha(Path(__file__).read_bytes()), source_record=6838169, source_version='v2',
        source_published_zip_size=668227028, source_published_zip_MD5='da74b7228763b76dcd2d9830c1a10cd9',
        required_directory_SHA256=directory_sha,
        future_full_member_required_SHA256='45fdf21775bf38fc3cf869ed60d45cb51caac0b1076b2a686572961082589c8e',
        maximum_official_route_trials=2, per_route_range_budget_including_retries=8*1024**2,
        cached_source_contract_not_fresh_API_replay=True,
        coalesced_block_bytes=65536,maximum_RAM_cache_blocks=8,fresh_wire_directory_replay=True,
        prefix_bytes=65536, strict_default_TLS_HTTP11_context=True,
        no_whole_archive_or_scientific_acceptance=True, no_decompressed_rows_parsed_or_written=True))
    def signature(entries):
        return sha(json.dumps([(x.filename, x.file_size, x.compress_size, x.CRC) for x in entries],
            separators=(',', ':')).encode())
    class PilotRange(transport.RangeFile):
        def _fetch_block(self, start):
            n = min(self.block_bytes, self.size-start)
            assert time.monotonic()-self.started <= 180
            assert self.requested+n*3 <= 8*1024**2
            return super()._fetch_block(start)
    for index, url in enumerate(urls, 1):
        rf = PilotRange(None, url, 668227028,
            lambda value: save('private_route_'+str(index)+'_progress.json', value))
        try:
            with zipfile.ZipFile(rf) as archive:
                assert signature(archive.infolist()) == directory_sha
                info = archive.getinfo('data/input/sc_data.csv')
                assert info.file_size == 172707928 and info.compress_size == 51454934
                with archive.open(info):
                    started = time.monotonic()
                    prefix = rf.read(65536)
                    body_seconds = time.monotonic()-started
                    assert len(prefix) == 65536
            rf.clear_cache()
            with zipfile.ZipFile(rf) as replay:
                assert signature(replay.infolist()) == directory_sha
            result = dict(status='COMPLETE_CACHED_CONTRACT_LOCAL_PREFIX_PILOT_ONLY',
                completed_UTC=datetime.now(timezone.utc).isoformat(), producer_SHA256=sha(Path(__file__).read_bytes()),
                successful_official_route=index, compressed_prefix_bytes=len(prefix), body_seconds=body_seconds,
                successful_range_bytes=rf.transferred, requested_range_bytes=rf.requested,
                directory_replayed_and_equal_prior_accepted_digest=True,
                fresh_metadata_API_not_selected=True, full_member_CRC_and_SHA_still_required=True,
                no_individual_rows_or_scientific_effects_written=True)
            save('private_cached_contract_complete.json', result)
            print(json.dumps(result))
            return
        except Exception as error:
            transient = isinstance(error, (TimeoutError, ConnectionResetError, BrokenPipeError,
                http.client.IncompleteRead, http.client.RemoteDisconnected))
            result = dict(status='ROUTE_PILOT_NOT_ACCEPTED', observed_UTC=datetime.now(timezone.utc).isoformat(),
                route=index, error_class=type(error).__name__, transient_network_failure=transient,
                frames=[dict(file=Path(x.filename).name, line=x.lineno, function=x.name)
                    for x in traceback.extract_tb(error.__traceback__)],
                no_raw_payload_or_record_values=True)
            save('private_route_'+str(index)+'_failed.json', result)
            if not transient:
                raise
        finally:
            rf.close()
    save('private_cached_contract_failed.json', dict(status='ALL_BOUNDED_CACHED_CONTRACT_ROUTES_NOT_ACCEPTED',
        observed_UTC=datetime.now(timezone.utc).isoformat(), no_scientific_acceptance=True))
    raise SystemExit(1)


if __name__ == '__main__':
    main()
