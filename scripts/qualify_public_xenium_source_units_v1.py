"""RAM-only public source-unit metadata; exports no individual/sample records."""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import io
import json
import os
import re
import urllib.request


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    out.mkdir(mode=0o700)
    save = lambda n, v: (out / n).write_text(json.dumps(v, indent=2) + '\n')
    url = os.environ['PANCANCER_GEO_SOURCE_URL']
    assert url == 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE300nnn/GSE300147/soft/GSE300147_family.soft.gz'
    save('plan_before_source_unit_metadata.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), source=url, metadata_bytes_RAM_only=True, selected_fields=['Sample accession', 'Sample_title for explicit source-case and run grouping', 'diagnosis/tissue/treatment/timepoint characteristics only'], no_outcome_values_sequences_or_patient_table=True, grouping='exact source title pattern only, no cross-study identity or followup joins', output='GSE-level source-unit counts and metadata schema only', index_SHA256=os.environ['PANCANCER_TAR_INDEX_SHA256']))
    with urllib.request.urlopen(url, timeout=25) as response:
        raw = response.read(2*1024**2+1)
    assert len(raw) <= 2*1024**2
    compressed_hash = hashlib.sha256(raw).hexdigest()
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as f:
        text = f.read(16*1024**2+1)
    assert len(text) <= 16*1024**2
    sample, records, keys = None, {}, set()
    selected = {'diagnosis','disease','histology','tissue','tissue_type','tumor_type','treatment','time_point','timepoint','sample_type'}
    for line in text.decode().splitlines():
        if line.startswith('^SAMPLE = '):
            sample = line.split(' = ', 1)[1]
            assert sample not in records
            records[sample] = dict(title=None, selected_characteristics={})
        elif line.startswith('^'):
            sample = None
        elif sample and line.startswith('!Sample_title = '):
            records[sample]['title'] = line.split(' = ', 1)[1]
        elif sample and line.startswith('!Sample_characteristics_ch1 = '):
            value = line.split(' = ', 1)[1]
            if ':' in value:
                key, value = value.split(':', 1)
                key = key.strip().casefold().replace(' ', '_')
                keys.add(key)
                if key in selected:
                    records[sample]['selected_characteristics'].setdefault(key, []).append(value.strip())
    index_path = Path(os.environ['PANCANCER_TAR_INDEX'])
    assert sha(index_path) == os.environ['PANCANCER_TAR_INDEX_SHA256']
    index = json.loads(index_path.read_text())
    samples = {re.match(r'GSM[0-9]+', r['name'])[0] for r in index['members'] if r['role'] == 'cell_feature_matrix'}
    assert samples == set(records) and len(samples) == 18
    cases, case_runs, diagnosis = Counter(), {}, Counter()
    unparsed = 0
    diagnostic_non_HNSCC = 0
    for key, r in records.items():
        m = re.fullmatch(r'Patient ([0-9]+), run ([0-9]+)', r['title'] or '')
        if m:
            cases[m[1]] += 1
            assert (m[1], m[2]) not in case_runs
            case_runs[m[1], m[2]] = key
        else:
            unparsed += 1
        for field in ['diagnosis','disease','histology','tumor_type']:
            for value in r['selected_characteristics'].get(field, []):
                diagnosis[value] += 1
        if any('ameloblastoma' in value.casefold() for field in ['diagnosis','disease','histology','tumor_type'] for value in r['selected_characteristics'].get(field, [])):
            diagnostic_non_HNSCC += 1
    result = dict(source='GSE300147', author_source_cases_not_cross_study_patient_review=True, GEO_samples=len(records), archive_samples=len(samples), exact_sample_key_set_match=True, title_pattern_unparsed_samples=unparsed, explicit_source_case_groups=len(cases), distribution_of_samples_per_source_case={str(k): v for k, v in Counter(cases.values()).items()}, diagnosis_value_counts=dict(diagnosis), explicit_ameloblastoma_sample_count=diagnostic_non_HNSCC, selected_characteristic_schema=sorted(keys & selected), no_new_clinical_endpoint_or_ICI_phase_acceptance=True, raw_SOFT_metadata_and_case_map_not_saved=True, no_sample_or_case_rows_exported=True)
    save('private_source_unit_summary.json', result)
    save('private_terminal.json', dict(status='COMPLETE_PUBLIC_GEO_SOURCE_UNIT_AND_ARCHIVE_KEY_QUALIFICATION', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), compressed_SOFT_SHA256=compressed_hash, source_summary_SHA256=sha(out / 'private_source_unit_summary.json'), source_case_units_not_new_independent_patient_truth=True, no_individual_rows_or_raw_metadata_saved=True))
    print(json.dumps(dict(status='COMPLETE_PUBLIC_GEO_SOURCE_UNIT_AND_ARCHIVE_KEY_QUALIFICATION', source_summary=result)))


if __name__ == '__main__':
    main()
