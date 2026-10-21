"""Small published metadata, redacted date grammar only, no actual date values."""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import os
import re
import zipfile

DIRECTORY_SHA = 'dff220e35d38da2b60c42950190d71b2259f12aae3c21ad4bc6c05ad5c5d2893'
CLINICAL_SHA = 'ac4a02e4d0e23504aa980b21585de74c4d7390d2e6e0ec2399303c5c29ca5b25'


def grammar(rows):
    output = {}
    for field in ['Date_surgery', 'Start_treatment']:
        counts = Counter()
        for row in rows:
            text = row[field].strip()
            if not text or text.lower() in {'na', 'nan', 'n/a', 'none', 'null'}:
                shape = 'MISSING'
            else:
                shape = ''.join('D' if c.isdigit() else 'L' if c.isalpha() else c for c in text)
            counts[shape] += 1
        output[field] = dict(sorted(counts.items()))
    return output


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--scripts-root', type=Path, required=True); args = parser.parse_args()
    assert not args.output.exists()
    path = args.scripts_root/'melanoma_IMC_existing_clash_transport_v1.py'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == '546b491e20e6797e5e255bfa226fa568a83dca08eb01334a27128fb0587755ac'
    spec = importlib.util.spec_from_file_location('bounded_grammar_transport', path)
    transport = importlib.util.module_from_spec(spec); spec.loader.exec_module(transport)
    transport.PROXY_PORT = int(os.environ['PANCANCER_EXISTING_CLASH_PORT'])
    rf = transport.RangeFile(None, 'https://zenodo.org/api/records/5994136/files/data_for_analysis.zip/content', 7123744479)
    rf.max_chunk = 64*1024
    def directory(archive):
        return hashlib.sha256(json.dumps([(x.filename, x.file_size, x.compress_size, x.CRC)
            for x in archive.infolist()], separators=(',', ':')).encode()).hexdigest()
    try:
        with zipfile.ZipFile(rf) as archive:
            assert directory(archive) == DIRECTORY_SHA
            raw = archive.read('data_for_analysis_zenodo/protein/clinical_data_protein.csv')
            assert hashlib.sha256(raw).hexdigest() == CLINICAL_SHA
            # Select only the two date strings; no other individual values retained.
            selected = [{key: row[key] for key in ['Date_surgery', 'Start_treatment']}
                        for row in csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))]
            result = grammar(selected)
        with zipfile.ZipFile(rf) as replay:
            assert directory(replay) == DIRECTORY_SHA
        args.output.write_text(json.dumps(dict(status='COMPLETE_REDACTED_DATE_GRAMMAR_ONLY',
            completed_UTC=datetime.now(timezone.utc).isoformat(), source_date_shape_counts=result,
            actual_dates_patient_rows_outcomes_not_returned=True,
            no_calendar_convention_or_clinical_phase_inferred=True,
            clinical_member_SHA256=CLINICAL_SHA, directory_SHA256=DIRECTORY_SHA,
            full_small_member_CRC_checked=True, wire_bytes=rf.transferred,
            producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2)+'\n')
    finally:
        rf.close()


if __name__ == '__main__':main()
