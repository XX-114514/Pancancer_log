"""Independent standard-library CSV/record linkage of clinical-source summaries."""
import base64
import csv
import hashlib
import io
import json
import math
import os
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    payload = json.load(sys.stdin)
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    production = Path(os.environ['PANCANCER_FOLLOWUP_PRODUCTION'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    terminal_path = production / 'private_terminal.json'
    terminal = json.loads(terminal_path.read_text())
    assert terminal['status'] == 'COMPLETE_ORIGINAL_AUTHOR_FOLLOWUP_SOURCE_QUALIFICATION_PENDING_AUDIT'
    for name, proof in terminal['files'].items():
        p = production / name
        assert p.stat().st_size == proof['bytes'] and sha(p) == proof['sha256']
    plan = json.loads((production / 'plan_before_record_linkage_and_followup_counts.json').read_text())
    raw = base64.b64decode(payload['public_csv_base64'], validate=True)
    assert len(raw) == 146428
    assert hashlib.sha256(raw).hexdigest() == plan['public_source_sha256'] == payload['public_csv_sha256']
    assert hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest() == plan['public_source_git_blob']
    fields = plan['selected_fields']
    reader = csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))
    assert len(reader.fieldnames) == len(set(reader.fieldnames))
    selected = [{k: row[k] for k in fields} for row in reader]
    with (production / 'private_selected_clinical_metadata.tsv').open() as h:
        actual = list(csv.DictReader(h, delimiter='\t'))
    assert actual == selected
    native_path = root / 'native_sample_clinical_registry.tsv'
    assert sha(native_path) == plan['native_registry_sha256']
    with native_path.open() as h:
        native = list(csv.DictReader(h, delimiter='\t'))
    missing = set(plan['missing_literals'])
    def valid(value):
        return value.strip().casefold() not in missing
    def numeric(value):
        if not valid(value):
            return None
        try:
            x = float(value)
        except ValueError:
            return None
        return x if math.isfinite(x) else None
    counts = Counter(r['SampleID'] for r in selected if valid(r['SampleID']))
    nc = Counter(r['sampleID'] for r in native if valid(r['sampleID']))
    unique = {r['SampleID']: r for r in selected if valid(r['SampleID']) and counts[r['SampleID']] == 1}
    nv = {r['sampleID']: r for r in native if valid(r['sampleID']) and nc[r['sampleID']] == 1}
    common = sorted(set(unique) & set(nv))
    shared = []
    for left, right in [('Pathological Response', 'pathological_response'),
                        ('Pre-treatment Staging', 'pre_treatment_staging'), ('histology', 'cancer_type')]:
        pairs = [(unique[k][left], nv[k][right]) for k in common if valid(unique[k][left]) and valid(nv[k][right])]
        shared.append(dict(source_field=left, native_field=right, comparable_records=len(pairs),
                           exact_literal_matches=sum(a == b for a, b in pairs),
                           strip_casefold_matches=sum(a.strip().casefold() == b.strip().casefold() for a, b in pairs)))
    time_rows = []
    for field in ['RFS_months', 'RFS_months_new', 'OS_months']:
        s = [r[field] for r in selected]
        x = [numeric(v) for v in s]
        finite = [v for v in x if v is not None]
        time_rows.append(dict(field=field, rows=len(s), missing=sum(not valid(v) for v in s),
                              nonnumeric_nonmissing=sum(valid(v) and n is None for v, n in zip(s, x)),
                              finite=len(finite), positive=sum(v > 0 for v in finite),
                              nonpositive=sum(v <= 0 for v in finite),
                              minimum=min(finite) if finite else None, maximum=max(finite) if finite else None))
    statuses = []
    for field in ['RFS_status', 'RFS_status_new', 'OS_status(surgery_to_death)']:
        for value, n in sorted(Counter(r[field] for r in selected).items()):
            statuses.append(dict(field=field, source_literal=value, records=n))
    pairs = [(numeric(r['RFS_months']), numeric(r['RFS_months_new'])) for r in selected]
    both = [(a, b) for a, b in pairs if a is not None and b is not None]
    comparisons = dict(time_both_numeric=len(both), time_equal=sum(a == b for a, b in both),
                       time_new_less=sum(b < a for a, b in both), time_new_greater=sum(b > a for a, b in both),
                       status_literal_equal=sum(r['RFS_status'] == r['RFS_status_new'] for r in selected),
                       time_original_only_numeric=sum(a is not None and b is None for a, b in pairs),
                       time_new_only_numeric=sum(a is None and b is not None for a, b in pairs))
    summary = dict(source_rows=len(selected), source_unique_valid_keys=len(counts),
                   source_rows_missing_key=sum(not valid(r['SampleID']) for r in selected),
                   source_repeated_key_tokens=sum(n > 1 for n in counts.values()), native_rows=len(native),
                   native_unique_valid_keys=len(nc), unique_exact_key_overlap=len(common),
                   source_unique_keys_absent_native=len(set(unique) - set(nv)),
                   native_unique_keys_absent_source=len(set(nv) - set(unique)), shared_field_comparisons=shared,
                   time_fields=time_rows, original_status_distributions=statuses, original_vs_new_RFS=comparisons,
                   author_patient_identity_bridge_accepted=False, endpoint_semantics_or_time_origin_accepted=False,
                   effects_computed=False)
    expected = json.loads((production / 'private_source_only_linkage_and_followup_summary.json').read_text())
    # Comparisons are exact integers/categories; float extrema only concern time units.
    assert summary == expected
    assert sha(native_path) == plan['native_registry_sha256']
    out.mkdir(mode=0o700)
    result = dict(status='ACCEPTED_ALL_SELECTED_CLINICAL_RECORDS_KEYS_AND_SOURCE_AGGREGATES',
                  completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                  producer_terminal_sha256=sha(terminal_path), source_csv_SHA256=hashlib.sha256(raw).hexdigest(),
                  independently_retrieved_same_source_csv=True, selected_records=len(selected),
                  selected_fields=len(fields), all_selected_record_values_and_aggregates_equal=True,
                  source_patient_truth_or_endpoint_definition_verified=False,
                  new_expression_or_clinical_effects=False)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
