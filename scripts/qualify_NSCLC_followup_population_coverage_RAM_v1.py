"""Clinical coverage/missingness only, immutable public inputs in local RAM.

No workbook, source CSV, individual records, expression or programme arrays
are saved. Group counts are qualification summaries, never survival effects.
"""
from collections import Counter
from pathlib import Path
from datetime import datetime, timezone
from urllib.request import Request, build_opener, ProxyHandler
import hashlib
import json
import math
import os
import re
import resource
import pandas as pd
import qualify_NSCLC_clinical_columns_local_RAM_v3 as qualification
import run_NSCLC_TCR_local_RNA_bridge_v1 as public_input

AXES = {
    'histology': ['LUAD', 'LUSC', 'OTHER_OR_UNKNOWN'],
    'pathological_response': ['pCR', 'MPR', 'non-MPR', 'OTHER_OR_UNKNOWN'],
    'pretreatment_stage_prefix': ['I', 'II', 'III', 'IV', 'OTHER_OR_UNKNOWN'],
    'strict_original_ICI_role': ['EXPLICIT_ICI_AUTHOR_CASE', 'NO_ICI', 'AMBIGUOUS_OR_UNKNOWN_ICI']}


def labels(native, pattern):
    def stage(value):
        return next((prefix for prefix in ['IV', 'III', 'II', 'I'] if value.startswith(prefix)), 'OTHER_OR_UNKNOWN')
    def role(value):
        return 'EXPLICIT_ICI_AUTHOR_CASE' if value == 'Yes' or pattern.fullmatch(value) else 'NO_ICI' if value == 'No' else 'AMBIGUOUS_OR_UNKNOWN_ICI'
    return pd.DataFrame({
        'histology': native.cancer_type.map(lambda v: v if v in ['LUAD', 'LUSC'] else 'OTHER_OR_UNKNOWN'),
        'pathological_response': native.pathological_response.map(lambda v: v if v in ['pCR', 'MPR', 'non-MPR'] else 'OTHER_OR_UNKNOWN'),
        'pretreatment_stage_prefix': native.pre_treatment_staging.map(stage),
        'strict_original_ICI_role': native.anti_PD1.map(role)}, index=native.index)


def coverage(native, clinical, pattern):
    assert native.sampleID.is_unique and clinical.SampleID.is_unique
    assert set(clinical.SampleID).issubset(set(native.sampleID))
    clinical = clinical.set_index('SampleID', drop=False)
    design = labels(native, pattern)
    present = native.sampleID.isin(clinical.index)
    metrics = {'in_followup_table': present.astype(int)}
    for name, time_field, status_field in [('RFS_old', 'RFS_months', 'RFS_status'),
            ('RFS_new', 'RFS_months_new', 'RFS_status_new'),
            ('OS', 'OS_months', 'OS_status(surgery_to_death)')]:
        times = pd.to_numeric(native.sampleID.map(clinical[time_field]), errors='coerce')
        status = pd.to_numeric(native.sampleID.map(clinical[status_field]), errors='coerce')
        metrics[name+'_numeric_complete'] = (times.map(lambda v: pd.notna(v) and math.isfinite(v) and v > 0) & status.isin([0, 1])).astype(int)
        metrics[name+'_code1'] = status.eq(1).astype(int)
    result = []
    for axis, groups in AXES.items():
        oracle = Counter()
        for i in range(len(native)):
            for metric, values in metrics.items():
                oracle[(design[axis].iloc[i], metric)] += int(values.iloc[i])
            oracle[(design[axis].iloc[i], 'native_records')] += 1
        frame = pd.DataFrame({'group': design[axis], **metrics})
        grouped = frame.groupby('group', dropna=False).sum(numeric_only=True)
        group_sizes = frame.groupby('group', dropna=False).size()
        table = []
        for group in groups:
            n = int(group_sizes.get(group, 0))
            row = dict(group=group, native_records=n)
            assert n == oracle[(group, 'native_records')]
            for metric in metrics:
                row[metric] = int(grouped.loc[group, metric]) if group in grouped.index else 0
                assert row[metric] == oracle[(group, metric)]
            row['not_in_table'] = n-row['in_followup_table']
            row['coverage_fraction'] = row['in_followup_table']/n if n else None
            assert 0 <= row['in_followup_table'] <= n
            table.append(row)
        assert sum(v['native_records'] for v in table) == len(native)
        assert sum(v['in_followup_table'] for v in table) == len(clinical)
        for metric, values in metrics.items():
            assert sum(v[metric] for v in table) == int(values.sum())
        result.append(dict(axis=axis, categories=table))
    return dict(native_records=len(native), followup_table_records=len(clinical),
        native_records_not_in_this_table=len(native)-len(clinical),
        overall_coverage_fraction=len(clinical)/len(native), axes=result,
        all_categorical_counts_independently_recomputed_by_counter=True,
        absent_from_table_reason_not_inferred=True, no_hypothesis_tests_or_survival_effects=True,
        code1_not_yet_adopted_as_specific_clinical_event=True,
        record_units_not_reviewed_independent_patients=True)


def main():
    root = Path(os.environ['PANCANCER_CLINICAL_SCOPE_TASK'])
    frozen = root/'private_frozen_NSCLC_followup_population_coverage_RAM_v1.json'
    plan = json.loads(frozen.read_text())
    for name, digest in plan['code_SHA256'].items():
        assert hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest() == digest
    assert plan['clinical_qualification_only_no_effects_or_individual_files']
    opener = build_opener(ProxyHandler({'https': 'http://127.0.0.1:7897'}))
    fields = ['sampleID', 'pathological_response', 'cancer_type', 'pre_treatment_staging', 'anti-PD1_therapy']
    buffer, receipt = public_input.fetch(plan['official_original_metadata'], opener)
    cases, ambiguous, nrows = {}, set(), 0
    for chunk in pd.read_csv(buffer, compression='gzip', usecols=fields, dtype=str, keep_default_na=False, chunksize=50000):
        nrows += len(chunk)
        assert nrows <= plan['expected_original_metadata_rows']
        for row in chunk[fields].drop_duplicates().itertuples(index=False, name=None):
            key, values = row[0], row[1:]
            if key in cases and cases[key] != values:
                ambiguous.add(key)
            else:
                cases[key] = values
    buffer.close()
    del buffer, chunk
    assert nrows == plan['expected_original_metadata_rows'] and not ambiguous
    native = pd.DataFrame([(k, *v) for k, v in cases.items()], columns=fields).rename(columns={'anti-PD1_therapy': 'anti_PD1'})
    del cases
    reference_path = root/'private_NSCLC_native_clinical_reference_cohort_digest_v1.json'
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == plan['native_reference_digest_receipt_SHA256']
    reference = json.loads(reference_path.read_text())
    actual = qualification.fingerprint(native, fields[0], fields[1:4])
    assert all(actual[k] == reference[k] for k in ['rows', 'valid_unique_key_rows', 'key_set_digest', 'fields'])
    with opener.open(Request(plan['clinical_source_url']), timeout=30) as response:
        assert response.status == 200 and response.url == plan['clinical_source_url']
        raw = response.read(plan['expected_clinical_bytes']+1)
    assert len(raw) == plan['expected_clinical_bytes'] and hashlib.sha256(raw).hexdigest() == plan['expected_clinical_SHA256']
    clinical = qualification.selected.parse_selected(raw)
    del raw
    linked = qualification.selected.qualify(clinical, native)
    assert linked['unique_exact_key_overlap'] == len(clinical) and linked['source_unique_keys_absent_native'] == 0
    assert all(v['exact_literal_matches'] == v['comparable_records'] == len(clinical) for v in linked['shared_field_comparisons'])
    pattern = re.compile(plan['strict_original_therapy_regex'])
    design = labels(native, pattern)
    assert design.strict_original_ICI_role.value_counts().to_dict() == plan['expected_original_therapy_role_counts']
    answer = coverage(native, clinical, pattern)
    del native, clinical, linked, design
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss < plan['local_peak_RSS_cap_bytes']
    answer.update(status='COMPLETE_LOCAL_CLINICAL_FOLLOWUP_POPULATION_COVERAGE_NO_EFFECTS',
        observed_UTC=datetime.now(timezone.utc).isoformat(),
        source='GSE243013', original_public_metadata_receipt=receipt,
        frozen_contract_SHA256=hashlib.sha256(frozen.read_bytes()).hexdigest(),
        same_original_public_metadata_and_native_cohort_digests=True,
        fixed_clinical_source_SHA256=plan['expected_clinical_SHA256'],
        original_therapy_role_grammar_and_source_counts_unchanged=True,
        no_cell_ID_expression_TCR_or_programme_fields_selected=True,
        no_individual_rows_or_hash_tokens_saved_uploaded_or_displayed=True,
        no_new_survival_or_communication_effect_family=True)
    out = root/'private_NSCLC_followup_population_coverage_v1.json'
    assert not out.exists()
    out.write_text(json.dumps(answer, indent=2)+'\n')
    print(json.dumps(dict(status=answer['status'], aggregate_receipt_SHA256=hashlib.sha256(out.read_bytes()).hexdigest(),
        no_individual_records_emitted=True)))


if __name__ == '__main__':
    main()
