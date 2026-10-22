"""Clinical-column-only local-RAM parser and qualification, never downloads.

The specifically denied source must not be fetched without new explicit human
authorization. This module only prepares a reviewable parser on synthetic data.
No individual records may be written or returned to chat/public outputs.
"""
from collections import Counter
import csv
import io
import math
import numpy as np
import pandas as pd
import qualify_NSCLC_original_author_followup_v1 as original

FIELDS = original.FIELDS
MISSING = original.MISSING


def select_csv_fields(text, fields):
    """Independent quoted-CSV lexer, retaining selected field text only."""
    header_end = text.find('\n')
    assert header_end >= 0
    header = next(csv.reader([text[:header_end].rstrip('\r')]))
    assert len(header) == len(set(header)) and set(fields).issubset(header)
    indexes = {header.index(name): name for name in fields}
    rows, row, value = [], {}, []
    col, quoted, just_closed, i = 0, False, False, header_end+1
    while i < len(text):
        char = text[i]
        if quoted:
            if char == '"':
                if i+1 < len(text) and text[i+1] == '"':
                    if col in indexes:
                        value.append('"')
                    i += 2
                    continue
                quoted = False
                just_closed = True
            elif col in indexes:
                value.append(char)
        elif char == '"':
            assert not value and not just_closed
            quoted = True
        elif char in ',\r\n':
            if col in indexes:
                row[indexes[col]] = ''.join(value)
            value, just_closed = [], False
            col += 1
            if char != ',':
                assert col == len(header) and set(row) == set(fields)
                rows.append(row)
                row, col = {}, 0
                if char == '\r' and i+1 < len(text) and text[i+1] == '\n':
                    i += 1
        else:
            assert not just_closed
            if col in indexes:
                value.append(char)
        i += 1
    assert not quoted
    if col or value or row:
        if col in indexes:
            row[indexes[col]] = ''.join(value)
        assert col+1 == len(header) and set(row) == set(fields)
        rows.append(row)
    return rows


def parse_selected(raw):
    text = raw.decode('utf-8-sig')
    header = next(csv.reader([text.split('\n', 1)[0].rstrip('\r')]))
    assert len(header) == len(set(header)) and set(FIELDS).issubset(header)
    table = pd.read_csv(io.StringIO(text), usecols=FIELDS, dtype=str, keep_default_na=False, nrows=2001)[FIELDS]
    assert 0 < len(table) <= 2000
    independent = select_csv_fields(text, FIELDS)
    assert independent == table.to_dict('records')
    return table


def safe_category_counts(table):
    rows = []
    allowed = {'', '0', '1', '2', '0.0', '1.0', '2.0', 'NA', 'NaN', 'nan', 'N/A', 'None', 'unknown'}
    for field in original.STATUS:
        counts = Counter(table[field])
        rows.append(dict(field=field, expected_code_counts={k: counts[k] for k in sorted(allowed) if counts[k]},
                         other_literal_category_count=sum(k not in allowed for k in counts),
                         other_literal_rows=sum(n for k, n in counts.items() if k not in allowed)))
    return rows


def qualify(table, native):
    """Return source aggregates only. Individual frames must stay local RAM."""
    assert native.sampleID.is_unique
    q = original.summarize(table, native)
    # The original implementation reports unrestricted status strings. Suppress
    # those values before any output; retain finite expected clinical codes only.
    q.pop('original_status_distributions')
    q['status_code_qualification'] = safe_category_counts(table)
    q['individual_patient_independence_not_established_by_unique_keys'] = True
    q['source_origin_event_definition_not_established_by_numeric_codes'] = True
    q['selected_values_independently_checked_by_nonretaining_CSV_lexer'] = True
    return q


def synthetic():
    fields = ['UNSELECTED_PROGRAM']+FIELDS+['UNSELECTED_FREE_TEXT']
    rows = []
    for i in range(3):
        row = {k: '1' for k in FIELDS}
        row.update(SampleID='SYNTHETIC_'+str(i), **{'Pathological Response': 'MPR', 'histology': 'LUAD', 'Pre-treatment Staging': 'II'})
        row['UNSELECTED_PROGRAM'] = 'DO_NOT_RETAIN,PROGRAM'
        row['UNSELECTED_FREE_TEXT'] = 'DO_NOT_RETAIN\n"FREE_TEXT"'
        rows.append(row)
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator='\r\n')
    writer.writeheader()
    writer.writerows(rows)
    table = parse_selected(output.getvalue().encode())
    assert 'DO_NOT_RETAIN' not in repr(table.to_dict('records'))
    assert select_csv_fields(output.getvalue().rstrip('\r\n'), FIELDS) == table.to_dict('records')
    native = pd.DataFrame(dict(sampleID=[row['SampleID'] for row in rows], pathological_response=['MPR']*3,
                               cancer_type=['LUAD']*3, pre_treatment_staging=['II']*3))
    table.loc[0, 'RFS_status'] = 'DO_NOT_EXPORT_UNKNOWN_LITERAL'
    summary = qualify(table, native)
    assert 'DO_NOT_EXPORT_UNKNOWN_LITERAL' not in repr(summary)
    assert summary['unique_exact_key_overlap'] == 3 and not summary['effects_computed']
    table.loc[1, 'SampleID'] = table.loc[0, 'SampleID']
    summary = qualify(table, native)
    assert summary['source_repeated_key_tokens'] == 1 and summary['unique_exact_key_overlap'] == 1


if __name__ == '__main__':
    synthetic()
    print('PASS_SYNTHETIC_SELECTED_CLINICAL_FIELDS_AND_UNKNOWN_LITERAL_EXCLUSION_NO_DOWNLOAD')
