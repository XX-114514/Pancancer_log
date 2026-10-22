"""Explicitly authorized clinical-column qualification, local RAM only.

No records, selected clinical tables, omics columns or individual hashes are
saved, uploaded or printed. An optional reference is cohort-level digests only.
The old remote producer main (which wrote selected rows) is never called.
"""
from pathlib import Path
from datetime import datetime, timezone
from urllib.request import Request, build_opener, ProxyHandler
from urllib.parse import urlparse
import hashlib
import json
import math
import os
import pandas as pd
import NSCLC_selected_clinical_RAM_v2 as selected


def fingerprint(table, key, columns):
    """Whole-cohort digests; no per-record hashes or IDs in the return value."""
    valid = ~table[key].str.strip().str.casefold().isin(selected.MISSING)
    unique = valid & ~table[key].duplicated(keep=False)
    work = table.loc[unique].sort_values(key)
    digest = lambda value: hashlib.sha256(json.dumps(value, separators=(',', ':'), ensure_ascii=True).encode()).hexdigest()
    keys = work[key].tolist()
    fields = []
    for name in columns:
        values = work[name].tolist()
        fields.append(dict(field=name,
            literal_key_value_digest=digest(list(zip(keys, values))),
            strip_casefold_value_digest=digest([(k, v.strip().casefold()) for k, v in zip(keys, values)])))
    return dict(rows=len(table), valid_unique_key_rows=int(unique.sum()),
        key_set_digest=digest(keys), fields=fields,
        digest_covers_entire_cohort_not_individual_tokens=True)


def clinical_only(table):
    empty = pd.DataFrame({k: pd.Series(dtype=str) for k in
        ['sampleID', 'pathological_response', 'cancer_type', 'pre_treatment_staging']})
    answer = selected.qualify(table, empty)
    for name in ['native_rows', 'native_unique_valid_keys', 'unique_exact_key_overlap',
                 'source_unique_keys_absent_native', 'native_unique_keys_absent_source',
                 'shared_field_comparisons']:
        answer.pop(name)
    # A zero-row sentinel reuses source QC only; it is never evidence of linkage.
    answer['no_individual_reference_records_read_or_used'] = True
    pairs = []
    def number(value):
        if value.strip().casefold() in selected.MISSING:
            return None
        try:
            x = float(value)
            return x if math.isfinite(x) else None
        except ValueError:
            return None
    for time_field, status_field in zip(selected.original.TIME, selected.original.STATUS):
        ts, ss = table[time_field].tolist(), table[status_field].tolist()
        tx, sx = list(map(number, ts)), list(map(number, ss))
        expected = next(x for x in answer['time_fields'] if x['field'] == time_field)
        assert expected['finite'] == sum(x is not None for x in tx)
        assert expected['positive'] == sum(x is not None and x > 0 for x in tx)
        assert expected['nonpositive'] == sum(x is not None and x <= 0 for x in tx)
        expected.pop('minimum'); expected.pop('maximum')
        pairs.append(dict(time_field=time_field, status_field=status_field,
            positive_finite_time_and_numeric_binary_status_records=sum(
                x is not None and x > 0 and s in [0, 1] for x, s in zip(tx, sx)),
            numeric_binary_status_records=sum(s in [0, 1] for s in sx),
            missing_or_other_status_records=sum(s not in [0, 1] for s in sx),
            numeric_eligibility_not_endpoint_semantics=True))
    answer['time_status_joint_qualification'] = pairs
    answer['OS_header_mentions_surgery_to_death_not_censoring_definition'] = True
    answer['RFS_new_transformation_and_event_definition_not_established'] = True
    answer['no_effects_candidate_selection_or_survival_models'] = True
    answer['clinical_source_fingerprint'] = fingerprint(table, 'SampleID',
        ['Pathological Response', 'histology', 'Pre-treatment Staging'])
    return answer


def synthetic():
    selected.synthetic()
    table = pd.DataFrame([{k: '' for k in selected.FIELDS} for _ in range(3)])
    table['SampleID'] = ['FAKE_A', 'FAKE_B', 'FAKE_C']
    table['RFS_months'] = ['12', '0', 'bad']
    table['RFS_status'] = ['1', '0', 'NA']
    table['RFS_months_new'] = ['12', '6', '']
    table['RFS_status_new'] = ['1', '2', '1']
    table['OS_months'] = ['1', '4', '-1']
    table['OS_status(surgery_to_death)'] = ['1', '0', '1']
    answer = clinical_only(table)
    assert [x['positive_finite_time_and_numeric_binary_status_records'] for x in
            answer['time_status_joint_qualification']] == [1, 1, 2]
    assert 'FAKE_A' not in json.dumps(answer)
    reversed_table = table.iloc[::-1]
    assert fingerprint(table, 'SampleID', ['histology']) == fingerprint(reversed_table, 'SampleID', ['histology'])


def main():
    root = Path(os.environ['PANCANCER_CLINICAL_SCOPE_TASK'])
    frozen = root/'private_frozen_NSCLC_clinical_columns_local_RAM_v3.json'
    plan = json.loads(frozen.read_text())
    assert plan['explicit_human_clinical_only_permission_received']
    for name, digest in plan['code_SHA256'].items():
        assert hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest() == digest
    url = plan['source_url']
    assert urlparse(url).hostname == 'raw.githubusercontent.com'
    opener = build_opener(ProxyHandler({'https': 'http://127.0.0.1:7897'}))
    with opener.open(Request(url), timeout=30) as response:
        assert response.status == 200 and response.url == url
        assert 'text/html' not in response.headers.get('Content-Type', '').casefold()
        raw = response.read(plan['expected_bytes']+1)
    assert len(raw) == plan['expected_bytes']
    assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest() == plan['expected_git_blob_SHA1']
    source_sha = hashlib.sha256(raw).hexdigest()
    table = selected.parse_selected(raw)
    del raw
    answer = clinical_only(table)
    del table
    answer.update(status='COMPLETE_LOCAL_RAM_CLINICAL_COLUMNS_ONLY_PENDING_REFERENCE_DIGEST',
        observed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=source_sha,
        pinned_public_commit=plan['public_commit'], selected_fields=selected.FIELDS,
        independent_selected_value_and_numeric_checks_passed=True,
        no_individual_rows_files_upload_or_display=True,
        frozen_contract_SHA256=hashlib.sha256(frozen.read_bytes()).hexdigest())
    out = root/'private_NSCLC_local_clinical_column_qualification_v3.json'
    assert not out.exists()
    out.write_text(json.dumps(answer, indent=2)+'\n')
    print(json.dumps(dict(status=answer['status'], aggregate_receipt_SHA256=hashlib.sha256(out.read_bytes()).hexdigest(),
        no_individual_records_emitted=True)))


if __name__ == '__main__':
    main()
