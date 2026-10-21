"""Source-only coverage of already selected imaging stage/site covariates.

No outcome values or expression are parsed. Records and category labels stay
in RAM. Report covariate consistency/levels, never patient or ROI rows.
"""
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import importlib.util
import io
import json
import os
import zipfile


def scope(rows):
    assert rows and all('patient' in row for row in rows)
    units = {row['patient'] for row in rows}
    assert all(units)
    result = []
    for field, role in [('Site', 'site'), ('TNMstage', 'stage')]:
        assert all(field in row for row in rows)
        values = defaultdict(set)
        for row in rows:
            value = row[field].strip()
            if value.casefold() not in {'', 'na', 'n/a', 'unknown', 'nan', 'none'}:
                values[row['patient']].add(value)
        conflicts = sum(len(values[u]) > 1 for u in units)
        missing = sum(len(values[u]) == 0 for u in units)
        levels = set().union(*(values[u] for u in units))
        result.append(dict(role=role, explicit_field=field, source_case_units=len(units),
            cases_with_conflicting_ROI_values=conflicts, cases_with_no_recognized_value=missing,
            distinct_nonmissing_categories=len(levels),
            full_case_constant_covariate=bool(len(levels) == 1 and not missing and not conflicts),
            necessary_dummy_columns_if_complete_consistent=max(0, len(levels)-1) if not missing and not conflicts else None,
            category_values_not_exported=True, no_effect_or_outcome_test=True))
    return dict(source_case_units=len(units), covariate_scope=result,
        source_ROI_units=len(rows), no_outcome_values_read=True)


def main():
    helper = Path(__file__).parent/'review_melanoma_IMC_coordinate_semantics_v1.py'
    assert hashlib.sha256(helper.read_bytes()).hexdigest() == '1e35e7f09cb5f792f1aa0a596ad9bbba7a04bf13d1616ded1e8756bfefbfc2cf'
    spec = importlib.util.spec_from_file_location('own_generic_review', helper)
    review = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(review)
    source = review.load(Path(os.environ['PANCANCER_IMC_TRANSPORT_HELPER']),
        '9fb200d98e321584fde8bca6cb92f72151ac33ec00d64e04d93e7e2b33e870e2', 'imaging_source')
    transport = review.load(Path(os.environ['PANCANCER_IMC_BOUNDED_TRANSPORT']),
        '709f9b5ced875fea663557bad936f32d004e12d5daf291048b7de476064b0d68', 'imaging_transport')
    record = json.loads(source.fetch(source.BASE+'api/records/6838169', 2*1024**2))
    item = next(x for x in record['files'] if x['key'] == 'Melanoma_code_data.zip')
    assert item['size'] == 668227028 and item['checksum'] == 'md5:da74b7228763b76dcd2d9830c1a10cd9'
    rf = transport.RangeFile(source, item['links']['self'], item['size'])
    with zipfile.ZipFile(rf) as archive:
        blob = archive.read('data/input/hc_result/hc_result_sample_anno_invasive.csv')
        assert review.sha(blob) == '59c9f04e4e12cc7484f1cc262b96b77c7a5eae4efff7d2510fbc23921b6f55dc'
        rows = [{field: row[field] for field in ['patient', 'Site', 'TNMstage']}
            for row in csv.DictReader(io.StringIO(blob.decode('utf-8-sig')))]
        result = scope(rows)
    again = next(x for x in json.loads(source.fetch(source.BASE+'api/records/6838169', 2*1024**2))['files']
        if x['key'] == item['key'])
    assert all(item.get(k) == again.get(k) for k in ['id', 'size', 'checksum', 'links'])
    rf.close()
    print(json.dumps(dict(result, status='COMPLETE_IMAGING_STAGE_SITE_COVARIATE_SCOPE_ONLY',
        observed_UTC=datetime.now(timezone.utc).isoformat(),
        producer_SHA256=review.sha(Path(__file__).read_bytes()), complete_input_member_SHA256=review.sha(blob),
        no_clinical_or_spatial_effect=True, no_individual_or_category_values_written=True), indent=2))


if __name__ == '__main__':
    main()
