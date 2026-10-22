"""Exact clinical material linkage using two incoming public sources in RAM.

Only four clinical columns are selected from the same previously accepted
official metadata object. No cell IDs, expression, TCR, rows or individual
hashes leave local RAM. The source-side reference contributes whole-cohort
digests only; no native records are transferred from the server.
"""
from pathlib import Path
from datetime import datetime, timezone
from urllib.request import Request, build_opener, ProxyHandler
import hashlib
import io
import json
import os
import pandas as pd
import qualify_NSCLC_clinical_columns_local_RAM_v3 as qualification
import run_NSCLC_TCR_local_RNA_bridge_v1 as public_input


def main():
    root = Path(os.environ['PANCANCER_CLINICAL_SCOPE_TASK'])
    frozen = root/'private_frozen_NSCLC_clinical_public_metadata_join_v1.json'
    plan = json.loads(frozen.read_text())
    for name, digest in plan['code_SHA256'].items():
        assert hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest() == digest
    assert plan['no_effects_or_individual_record_exports']
    opener = build_opener(ProxyHandler({'https': 'http://127.0.0.1:7897'}))
    fields = ['sampleID', 'pathological_response', 'cancer_type', 'pre_treatment_staging']
    buffer, receipt = public_input.fetch(plan['official_original_metadata'], opener)
    cases, ambiguous, rows = {}, set(), 0
    for chunk in pd.read_csv(buffer, compression='gzip', usecols=fields, dtype=str,
                             keep_default_na=False, chunksize=50000):
        rows += len(chunk)
        assert rows <= plan['expected_original_metadata_rows']
        for row in chunk[fields].drop_duplicates().itertuples(index=False, name=None):
            key, values = row[0], row[1:]
            if key in cases and cases[key] != values:
                ambiguous.add(key)
            else:
                cases[key] = values
    buffer.close()
    del buffer, chunk
    assert rows == plan['expected_original_metadata_rows'] and not ambiguous
    native = pd.DataFrame([(k, *v) for k, v in cases.items()], columns=fields)
    del cases
    reference_path = root/'private_NSCLC_native_clinical_reference_cohort_digest_v1.json'
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == plan['source_reference_digest_receipt_SHA256']
    reference = json.loads(reference_path.read_text())
    actual = qualification.fingerprint(native, fields[0], fields[1:])
    key_equal = all(actual[k] == reference[k] for k in ['rows', 'valid_unique_key_rows', 'key_set_digest'])
    literal_equal = actual['fields'] == reference['fields']
    assert key_equal and literal_equal, 'PUBLIC_METADATA_NOT_IDENTICAL_TO_NATIVE_REFERENCE'
    with opener.open(Request(plan['clinical_source_url']), timeout=30) as response:
        assert response.status == 200 and response.url == plan['clinical_source_url']
        raw = response.read(plan['expected_clinical_bytes']+1)
    assert len(raw) == plan['expected_clinical_bytes']
    assert hashlib.sha256(raw).hexdigest() == plan['expected_clinical_SHA256']
    table = qualification.selected.parse_selected(raw)
    del raw
    q = qualification.selected.qualify(table, native)
    names = ['source_rows', 'source_unique_valid_keys', 'source_rows_missing_key',
        'source_repeated_key_tokens', 'native_rows', 'native_unique_valid_keys',
        'unique_exact_key_overlap', 'source_unique_keys_absent_native',
        'native_unique_keys_absent_source', 'shared_field_comparisons']
    answer = {k: q[k] for k in names}
    answer.update(status='COMPLETE_LOCAL_EXACT_CLINICAL_MATERIAL_LINKAGE_NO_EFFECTS',
        observed_UTC=datetime.now(timezone.utc).isoformat(),
        original_public_metadata_receipt=receipt,
        original_public_metadata_clinical_columns=fields,
        entire_reconstructed_reference_keys_and_clinical_values_identical_to_server_reference=True,
        server_individual_records_not_transferred=True,
        source_reference_registry_unchanged=reference['native_registry_unchanged'],
        independent_reviewed_patient_identity_not_established=True,
        RFS_event_and_time_origin_not_established=True,
        no_clinical_effects_or_survival_models_computed=True,
        frozen_contract_SHA256=hashlib.sha256(frozen.read_bytes()).hexdigest(),
        all_source_individual_rows_and_arrays_local_RAM_only=True)
    del table, native, q
    out = root/'private_NSCLC_clinical_material_exact_linkage_v1.json'
    assert not out.exists()
    out.write_text(json.dumps(answer, indent=2)+'\n')
    print(json.dumps(dict(status=answer['status'], aggregate_receipt_SHA256=hashlib.sha256(out.read_bytes()).hexdigest(),
        no_individual_data_displayed_saved_or_uploaded=True)))


if __name__ == '__main__':
    main()
