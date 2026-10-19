"""RAM-only public clinical-source qualification; no expression or effects."""
import base64
import csv
import hashlib
import io
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

COMMIT = '63ba67130a8e4c11e0bb7f5fac31b12ffcb2fab4'
BLOB = '10b5378f82c6fa94eb1626f4cb795f1c37d9a8fb'
SOURCE_BYTES = 146428
FIELDS = ['SampleID', 'Pathological Response', 'PRR', 'histology',
          'Pre-treatment TNM', 'Pre-treatment Staging', 'Surgery TNM',
          'RFS_months', 'RFS_status', 'RFS_months_new', 'RFS_status_new',
          'OS_months', 'OS_status(surgery_to_death)', 'pCR', 'MPR']
TIME = ['RFS_months', 'RFS_months_new', 'OS_months']
STATUS = ['RFS_status', 'RFS_status_new', 'OS_status(surgery_to_death)']
MISSING = {'', 'na', 'n/a', 'nan', 'none', 'null', 'unknown', 'not_reviewed'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def normalized(series):
    return series.str.strip().str.casefold()


def available(series):
    return ~normalized(series).isin(MISSING)


def summarize(table, native):
    key = table.SampleID
    valid = available(key)
    counts = key[valid].value_counts()
    nk = native.sampleID
    native_valid = available(nk)
    unique = table[valid & ~key.duplicated(keep=False)].set_index('SampleID')
    nv = native[native_valid & ~nk.duplicated(keep=False)].set_index('sampleID')
    common = unique.index.intersection(nv.index, sort=False)
    shared = []
    for left, right in [('Pathological Response', 'pathological_response'),
                        ('Pre-treatment Staging', 'pre_treatment_staging'),
                        ('histology', 'cancer_type')]:
        a = unique.loc[common, left]
        b = nv.loc[common, right]
        both = available(a) & available(b)
        shared.append(dict(source_field=left, native_field=right,
                           comparable_records=int(both.sum()),
                           exact_literal_matches=int((a[both] == b[both]).sum()),
                           strip_casefold_matches=int((normalized(a[both]) == normalized(b[both])).sum())))
    time_rows = []
    for field in TIME:
        s = table[field]
        x = pd.to_numeric(s.where(available(s)), errors='coerce')
        finite = np.isfinite(x.to_numpy(float))
        time_rows.append(dict(field=field, rows=len(s), missing=int((~available(s)).sum()),
                              nonnumeric_nonmissing=int((available(s) & x.isna()).sum()),
                              finite=int(finite.sum()), positive=int((x > 0).sum()),
                              nonpositive=int((x <= 0).sum()),
                              minimum=float(x[finite].min()) if finite.any() else None,
                              maximum=float(x[finite].max()) if finite.any() else None))
    status_rows = []
    for field in STATUS:
        # Only aggregate original categorical literals; never disclose case keys.
        for value, n in table[field].value_counts(dropna=False).sort_index().items():
            status_rows.append(dict(field=field, source_literal=str(value), records=int(n)))
    a = pd.to_numeric(table.RFS_months.where(available(table.RFS_months)), errors='coerce')
    b = pd.to_numeric(table.RFS_months_new.where(available(table.RFS_months_new)), errors='coerce')
    both = np.isfinite(a) & np.isfinite(b)
    comparisons = dict(time_both_numeric=int(both.sum()), time_equal=int((a[both] == b[both]).sum()),
                       time_new_less=int((b[both] < a[both]).sum()), time_new_greater=int((b[both] > a[both]).sum()),
                       status_literal_equal=int(table.RFS_status.eq(table.RFS_status_new).sum()),
                       time_original_only_numeric=int((np.isfinite(a) & ~np.isfinite(b)).sum()),
                       time_new_only_numeric=int((~np.isfinite(a) & np.isfinite(b)).sum()))
    return dict(source_rows=len(table), source_unique_valid_keys=len(counts),
                source_rows_missing_key=int((~valid).sum()), source_repeated_key_tokens=int((counts > 1).sum()),
                native_rows=len(native), native_unique_valid_keys=int(nk[native_valid].nunique()),
                unique_exact_key_overlap=len(common), source_unique_keys_absent_native=int((~unique.index.isin(nv.index)).sum()),
                native_unique_keys_absent_source=int((~nv.index.isin(unique.index)).sum()),
                shared_field_comparisons=shared, time_fields=time_rows, original_status_distributions=status_rows,
                original_vs_new_RFS=comparisons, author_patient_identity_bridge_accepted=False,
                endpoint_semantics_or_time_origin_accepted=False, effects_computed=False)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    payload = json.load(sys.stdin)
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    raw = base64.b64decode(payload['public_csv_base64'], validate=True)
    assert len(raw) == SOURCE_BYTES
    assert hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest() == BLOB
    source_sha = hashlib.sha256(raw).hexdigest()
    assert source_sha == payload['public_csv_sha256'] and payload['commit'] == COMMIT
    header = next(csv.reader(io.StringIO(raw.decode('utf-8-sig'))))
    assert len(header) == len(set(header)) and set(FIELDS).issubset(header)
    table = pd.read_csv(io.BytesIO(raw), usecols=FIELDS, dtype=str, keep_default_na=False)[FIELDS]
    native_path = root / 'native_sample_clinical_registry.tsv'
    before = sha(native_path)
    native = pd.read_csv(native_path, sep='\t', dtype=str, keep_default_na=False)
    assert native.sampleID.is_unique
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                public_repository_commit=COMMIT, public_source_git_blob=BLOB, public_source_sha256=source_sha,
                native_registry_sha256=before, selected_fields=FIELDS, missing_literals=sorted(MISSING),
                public_source_body_RAM_only=True, learned_fractions_or_programmes_not_parsed=True,
                clinical_effects_not_computed=True, native_registry_readonly=True,
                patient_independence_not_assumed=True, output_cap_bytes=2 * 1024 ** 2)
    (out / 'plan_before_record_linkage_and_followup_counts.json').write_text(json.dumps(plan, indent=2) + '\n')
    table.to_csv(out / 'private_selected_clinical_metadata.tsv', sep='\t', index=False)
    summary = summarize(table, native)
    assert sha(native_path) == before
    (out / 'private_source_only_linkage_and_followup_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    files = {p.name: dict(bytes=p.stat().st_size, sha256=sha(p)) for p in out.iterdir() if p.is_file()}
    assert sum(v['bytes'] for v in files.values()) < 2 * 1024 ** 2
    terminal = dict(status='COMPLETE_ORIGINAL_AUTHOR_FOLLOWUP_SOURCE_QUALIFICATION_PENDING_AUDIT',
                    completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                    native_registry_unchanged=True, effects_computed=False, files=files)
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2) + '\n')
    print(json.dumps(dict(status=terminal['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
