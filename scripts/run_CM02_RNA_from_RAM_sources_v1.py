"""Run the fixed RNA pipeline with exact public original source bytes in RAM.

The frozen manifest supplies code hashes only. Public RDS source bytes are
accepted through the existing whitelist/SHA receiver and never extracted.
The callback may only request those original URLs, via a finally-restored
process-local urlopen wrapper. This does not change network configuration.
"""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import sys
import traceback


def main():
    root=Path(__file__).parent.parent;cache={}
    try:
        manifest=json.loads((root/'frozen_code_manifest.json').read_text())
        for name,expected in manifest['files'].items():
            assert name==Path(name).name
            assert hashlib.sha256((root/'scripts'/name).read_bytes()).hexdigest()==expected
        import qualify_CM02_RNA_programme_RAM_public_sources_v4 as receiver
        cache,proof=receiver.receive_exact_sources(sys.stdin.buffer)
        (root/'private_original_public_source_RAM_receipt.json').write_text(json.dumps(dict(
            status='COMPLETE_EXACT_PUBLIC_SOURCE_BYTES_IN_RAM',completed_UTC=datetime.now(timezone.utc).isoformat(),
            sources=proof,no_original_source_files_or_local_derived_records_saved=True),indent=2)+'\n')
        import analyze_CM02_RNA_multigene_in_RAM_v1 as pipeline
        receiver.with_cached_sources(cache,pipeline.main)
    except BaseException as error:
        (root/'private_driver_failure.json').write_text(json.dumps(dict(
            status='RAM_RNA_DRIVER_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),
            error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],
            no_full_error_or_individual_values_written=True),indent=2)+'\n')
        raise SystemExit(1)
    finally:cache.clear()


if __name__=='__main__':main()
