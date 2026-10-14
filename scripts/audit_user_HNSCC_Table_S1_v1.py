"""Read original author Table S1 privately; audit keys and response consistency.

No clinical testing, inferred response overwrite, or patient-row export to stdout.
"""
import argparse
import csv
import hashlib
import json
import re
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import xml.etree.ElementTree as ET

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def read_table(path):
    with zipfile.ZipFile(path) as archive:
        strings = [''.join(s.itertext()) for s in ET.fromstring(
            archive.read('xl/sharedStrings.xml')).findall('m:si', NS)]
        wb = ET.fromstring(archive.read('xl/workbook.xml'))
        sheets = wb.findall('m:sheets/m:sheet', NS)
        assert len(sheets) == 1 and sheets[0].attrib['name'] == 'S1'
        rels = ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))
        rid = sheets[0].attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
        target = next(r.attrib['Target'] for r in rels if r.attrib['Id'] == rid)
        target = target.lstrip('/') if target.startswith('/') else 'xl/' + target
        rows = []
        for row in ET.fromstring(archive.read(target)).findall('m:sheetData/m:row', NS):
            values = {}
            for cell in row.findall('m:c', NS):
                value = cell.find('m:v', NS)
                kind = cell.attrib.get('t')
                text = strings[int(value.text)] if kind == 's' and value is not None else (
                    ''.join(cell.find('m:is', NS).itertext()) if kind == 'inlineStr' else
                    value.text if value is not None else '')
                values[re.sub(r'\d', '', cell.attrib['r'])] = text
            rows.append(values)
        assert 'Clinical metadata' in rows[0]['A']
        headers = rows[1]
        return [{name: row.get(col, '') for col, name in headers.items()} for row in rows[2:]]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workbook', type=Path, required=True)
    parser.add_argument('--capture-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    data = read_table(args.workbook)
    assert len(data) == 57
    cross = list(csv.DictReader((args.capture_root / 'private_GEO_patient_title_crosscheck.tsv').open(), delimiter='\t'))
    counts = {r['library_id']: r for r in csv.DictReader((args.capture_root / 'all57_private_cell_capture.tsv').open(), delimiter='\t')}
    lookup = {(r['patient.tag'], r['time']): r for r in data}
    assert len(lookup) == len(data)
    phases = {'pre': 'Pre', 'post': 'Post', 'post2': 'Recurrence'}
    joined = []
    for row in cross:
        key = (row['patient_title_token'], phases[row['title_timepoint']])
        assert key in lookup, 'Exact GEO patient/time key missing'
        table = lookup[key]
        assert table['neo_type'] == counts[row['library_id']]['treatment']
        joined.append(dict(row, **table, pathology_threshold_response=int(float(table['Tumor.regression']) >= 50)))
    assert len(joined) == len(data) and len({(r['patient.tag'], r['time']) for r in joined}) == len(data)
    groups = defaultdict(list)
    for row in joined:
        groups[row['patient.tag']].append(row)
    conflicts = {p for p, rows in groups.items() if len({r['Partial.pathologic.response'] for r in rows}) > 1}
    regression_conflicts = {p for p, rows in groups.items() if len({r['Tumor.regression'] for r in rows}) > 1}
    mismatch = [r for r in joined if int(float(r['Partial.pathologic.response'])) != r['pathology_threshold_response']]
    summary = dict(
        status='AUTHOR_TABLE_EXACT_KEYS_MATCHED_RESPONSE_DISCREPANCIES_QUARANTINED',
        observed_UTC=datetime.now(timezone.utc).isoformat(),
        sha256=hashlib.sha256(args.workbook.read_bytes()).hexdigest(),
        exact_matched_libraries=len(joined), source_patient_tokens=len(groups),
        time_counts=dict(Counter(r['time'] for r in joined)),
        treatment_matches=len(joined), response_threshold_mismatch_rows=len(mismatch),
        mismatch_by_time=dict(Counter(r['time'] for r in mismatch)),
        patients_with_conflicting_binary_response=len(conflicts),
        patients_with_conflicting_regression=len(regression_conflicts),
        response_field_adopted=False, clinical_tests_executed=False,
        clinical_analysis_gate='HOLD_FOR_AUTHOR_BINARY_VS_PATHOLOGY_DISCREPANCY_RECONCILIATION',
        recurrence_mapping='Exact patient + GEO post2 matched to author Recurrence; separate from primary Post',
        doi='10.1016/j.xcrm.2025.102408')
    args.output.mkdir()
    for name, records in [('all57_private_author_table_join.tsv', joined), ('private_response_discrepancies.tsv', mismatch)]:
        with (args.output / name).open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(joined[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(records)
    (args.output / 'aggregate_audit.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
