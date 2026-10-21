"""Independent MINiML source-case/run audit, with source-stage aggregation.

No sample/case rows, outcomes or raw metadata are saved. Technical sections
never become additional patients; non-HNSCC source names are explicitly counted.
"""
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import io
import json
import os
import re
import tarfile
import urllib.request
import xml.etree.ElementTree as ET


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    parent = Path(os.environ['PANCANCER_SOURCE_UNIT_PARENT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    terminal = json.loads((parent / 'private_terminal.json').read_text())
    assert sha(parent / 'private_source_unit_summary.json') == terminal['source_summary_SHA256']
    previous = json.loads((parent / 'private_source_unit_summary.json').read_text())
    out.mkdir(mode=0o700)
    save = lambda n, v: (out / n).write_text(json.dumps(v, indent=2) + '\n')
    url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE300nnn/GSE300147/miniml/GSE300147_family.xml.tgz'
    save('plan_before_MINiML_source_units_audit.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), parent_terminal_SHA256=sha(parent / 'private_terminal.json'), source=url, selected_fields=['Sample accession and explicit title pattern', 'Channel Source', 'recurrence and p16 status only'], no_outcome_values_or_cross_study_patient_identity=True, output='source aggregate diagnosis, repeated-section and stage counts only'))
    with urllib.request.urlopen(url, timeout=25) as r:
        raw = r.read(2*1024**2+1)
    assert len(raw) <= 2*1024**2
    with tarfile.open(fileobj=io.BytesIO(raw), mode='r:gz') as archive:
        document = None
        for entry in archive:
            if entry.isfile() and entry.name.endswith('.xml'):
                assert entry.size <= 16*1024**2
                document = archive.extractfile(entry).read(16*1024**2+1)
                break
    assert document is not None and len(document) <= 16*1024**2
    root = ET.fromstring(document)
    records = []
    for e in root.findall('.//{*}Sample'):
        source = e.findtext('{*}Channel/{*}Source') or ''
        title = e.findtext('{*}Title') or ''
        m = re.fullmatch(r'Patient ([0-9]+), run ([0-9]+)', title)
        assert m
        fields = {v.attrib.get('tag', '').strip().casefold().replace(' ', '_'): (v.text or '').strip() for v in e.findall('{*}Channel/{*}Characteristics') if v.attrib.get('tag', '').strip().casefold().replace(' ', '_') in ['recurrence', 'p16_status']}
        name = source.casefold()
        diagnosis = 'ameloblastoma' if 'ameloblastoma' in name else ('HNSCC_or_squamous' if 'hnscc' in name or 'squamous' in name else 'other_or_unresolved')
        recurrence = fields.get('recurrence', '').casefold()
        stage = 'primary' if recurrence == 'primary' else ('recurrent' if recurrence in ['recurrent','recurrence','locoregional recurrence','local recurrence'] else 'other_or_unknown')
        p16 = fields.get('p16_status', '').casefold()
        p16 = p16 if p16 in ['positive','negative'] else 'unknown_or_not_applicable'
        records.append(dict(accession=e.attrib['iid'], source=source, case=m[1], run=m[2], diagnosis=diagnosis, recurrence=stage, p16=p16))
    index_path = Path(os.environ['PANCANCER_TAR_INDEX'])
    assert sha(index_path) == os.environ['PANCANCER_TAR_INDEX_SHA256']
    index = json.loads(index_path.read_text())
    archived = {re.match(r'GSM[0-9]+', e['name'])[0] for e in index['members'] if e['role'] == 'cell_feature_matrix'}
    assert archived == {r['accession'] for r in records} and len(records) == len(archived) == 18
    assert len({(r['case'], r['run']) for r in records}) == len(records)
    groups = defaultdict(list)
    for r in records:
        groups[r['case']].append(r)
    source_names = Counter(r['source'] for r in records)
    assert dict(source_names) == previous['source_name_categories']
    assert len(groups) == previous['explicit_source_case_groups']
    assert {str(k): v for k, v in Counter(len(v) for v in groups.values()).items()} == previous['distribution_of_samples_per_source_case']
    case_diagnosis, sample_diagnosis, case_stages, case_p16 = Counter(), Counter(), Counter(), Counter()
    for values in groups.values():
        ds = {r['diagnosis'] for r in values}
        assert len(ds) == 1
        diagnosis = ds.pop()
        case_diagnosis[diagnosis] += 1
        sample_diagnosis[diagnosis] += len(values)
        if diagnosis == 'HNSCC_or_squamous':
            stages = {r['recurrence'] for r in values}
            p16s = {r['p16'] for r in values}
            case_stages[next(iter(stages)) if len(stages) == 1 else 'conflicting'] += 1
            case_p16[next(iter(p16s)) if len(p16s) == 1 else 'conflicting'] += 1
    assert dict(case_diagnosis) == previous['source_case_diagnosis_groups']
    result = dict(source='GSE300147', archive_and_both_metadata_sample_key_sets_exact=True, source_case_groups=len(groups), source_sample_count=len(records), source_case_diagnoses=dict(case_diagnosis), source_sample_diagnoses=dict(sample_diagnosis), HNSCC_source_case_recurrence=dict(case_stages), HNSCC_source_case_p16=dict(case_p16), source_case_groups_not_new_cross_study_patient_review=True, raw_metadata_or_case_sample_rows_not_saved=True, no_new_ICI_phase_or_clinical_endpoint_acceptance=True)
    save('private_independent_source_units_and_stage_summary.json', result)
    save('private_terminal.json', dict(status='ACCEPTED_MINiML_SOFT_ARCHIVE_SOURCE_UNIT_AND_STAGE_AGGREGATION', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), parent_terminal_SHA256=sha(parent / 'private_terminal.json'), MINiML_archive_SHA256=hashlib.sha256(raw).hexdigest(), independent_XML_parser_replayed_source_names_titles_case_runs_and_archive_keys=True, all_individual_records_RAM_only=True, no_new_spatial_or_clinical_inference=True, files={p.name: dict(SHA256=sha(p), bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='ACCEPTED_MINiML_SOFT_ARCHIVE_SOURCE_UNIT_AND_STAGE_AGGREGATION', source_summary=result)))


if __name__ == '__main__':
    main()
