"""Restore the previously accepted request-only direct GEO route.

The SSH shell's inherited Codex loopback proxy is not applied to these two
public reference requests. Original URLs, default TLS, complete SHA/length
checks and bounded cache reader remain unchanged. No OS/Clash/Codex edits.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import urllib.request
from datetime import datetime,timezone


def with_direct_requests(callback):
    original=urllib.request.urlopen
    opener=urllib.request.build_opener(urllib.request.ProxyHandler({}))
    try:
        urllib.request.urlopen=opener.open
        return callback()
    finally:urllib.request.urlopen=original


def main():
    file=Path(__file__).with_name('qualify_CM02_RNA_programme_cached_sources_v2.py')
    assert hashlib.sha256(file.read_bytes()).hexdigest()=='ddf0e74118f569401bb219a21ed4f21be35905153a57be50e4d4e0509b35fb8f'
    spec=importlib.util.spec_from_file_location('same_full_source_cache',file)
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    try:with_direct_requests(core.main)
    finally:
        proof=Path(__file__).parent.parent/'private_request_route_proof_v3.json'
        proof.write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),
            route='Existing successful workflow: direct GEO in this own process only',
            original_urls_TLS_SHA_gates_unchanged=True,inherited_proxy_settings_not_modified=True,
            request_function_restored_after_own_step=True,
            wrapper_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+chr(10))


if __name__=='__main__':main()
