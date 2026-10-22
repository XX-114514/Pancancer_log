"""Read generic published panel CSVs only; never cells or clinical observations."""
from pathlib import Path
import csv
import datetime
import hashlib
import io
import json
import os
import re
import urllib.request


ALIASES = {
    'CD3_family': ('CD3', 'CD3E'), 'CD4': ('CD4',),
    'CD8_family': ('CD8', 'CD8A'), 'FOXP3': ('FOXP3',),
    'CTLA4': ('CTLA4', 'CD152'), 'PD1': ('PD1', 'PDCD1', 'CD279'),
    'TIM3': ('TIM3', 'HAVCR2', 'CD366'), 'MKI67': ('MKI67', 'KI67'),
    'CCR7': ('CCR7', 'CD197'), 'CD127': ('CD127', 'IL7R'),
}


def fetch(url, maximum):
    with urllib.request.urlopen(url, timeout=30) as response:
        assert response.status == 200
        data = response.read(maximum + 1)
        assert len(data) <= maximum, 'METADATA_TOO_LARGE'
        return data, dict(url=url, final_url=response.url,
            content_type=response.headers.get('Content-Type'), bytes=len(data),
            SHA256=hashlib.sha256(data).hexdigest(),
            observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())


def normalized_words(labels):
    # Only exact published names or listed aliases; do not infer absent markers.
    return {re.sub('[^A-Z0-9]', '', token.upper()) for label in labels
            for token in re.findall(r'[A-Za-z][A-Za-z0-9-]*', label)}


def main():
    task = Path(os.environ['PANCANCER_PUBLISHED_PANEL_TASK']).resolve()
    repository = Path(__file__).resolve().parents[1]
    assert task != repository and repository not in task.parents
    config = json.loads((task / 'private_published_panel_config_v1.json').read_text())
    assert config['scope'] == 'GENERIC_PANEL_METADATA_ONLY_NO_CELL_OR_CLINICAL_DATA'
    raw, record_proof = fetch('https://zenodo.org/api/records/' + config['record'], 1024**2)
    record = json.loads(raw)
    assert str(record['id']) == config['record']
    assert record['metadata']['access_right'] == 'open'
    choices = [v for v in record['files'] if hashlib.sha256(v['key'].encode()).hexdigest()
               == config['channel_file_key_SHA256']]
    assert len(choices) == 1
    entry = choices[0]
    assert entry['size'] == config['channel_file_bytes']
    assert entry['checksum'] == 'md5:' + config['channel_file_MD5']
    channel_bytes, channel_proof = fetch(config['channel_file_URL'], 16384)
    assert len(channel_bytes) == entry['size']
    assert hashlib.md5(channel_bytes).hexdigest() == config['channel_file_MD5']
    channel_rows = list(csv.reader(io.StringIO(channel_bytes.decode('utf-8-sig'))))
    assert len(channel_rows) <= 100 and all(len(r) == 1 for r in channel_rows)
    channels = [r[0].strip() for r in channel_rows if r[0].strip()]
    protocol_bytes, protocol_proof = fetch(config['protocol_URL'], 65536)
    rows = list(csv.reader(io.StringIO(protocol_bytes.decode('utf-8-sig'))))
    assert 1 < len(rows) <= 150 and max(map(len, rows)) <= 30
    columns = [j for j, value in enumerate(rows[0]) if re.fullmatch(
        r'(biomarker|marker|target|antibody|antigen|protein)( name)?', value.strip(), re.I)]
    assert len(columns) == 1, 'AMBIGUOUS_OR_UNKNOWN_TARGET_COLUMN'
    targets = [r[columns[0]].strip() for r in rows[1:] if len(r) > columns[0] and r[columns[0]].strip()]
    channel_words, target_words = normalized_words(channels), normalized_words(targets)
    presence = {name: dict(channel_declared=bool(set(aliases) & channel_words),
                          protocol_declared=bool(set(aliases) & target_words))
                for name, aliases in ALIASES.items()}
    for row in presence.values():
        row['status'] = ('PHYSICALLY_DECLARED_FINAL_QC_AND_NUMERIC_MEASUREMENT_UNREVIEWED'
                         if row['channel_declared'] and row['protocol_declared'] else
                         'NOT_IN_EITHER_DECLARED_PHYSICAL_PANEL'
                         if not row['channel_declared'] and not row['protocol_declared'] else
                         'CHANNEL_PROTOCOL_LABEL_DISAGREEMENT_REQUIRES_REVIEW')
    states = {
        'regulatory_candidate': ('CD4', 'FOXP3', 'CTLA4'),
        'receptor_pair_candidate': ('CD8_family', 'PD1', 'TIM3'),
        'MKI67_candidate': ('CD8_family', 'MKI67'),
        'CCR7_CD127_candidate': ('CD8_family', 'CCR7', 'CD127'),
    }
    report = dict(status='PUBLISHED_GENERIC_PANEL_QUALIFICATION_ONLY',
        observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        record=record_proof, channel_list=channel_proof, protocol=protocol_proof,
        source_panel_metadata_receipt=config['previous_metadata_receipt_SHA256'],
        marker_presence=presence,
        necessary_state_marker_coverage={name: all(presence[x]['channel_declared'] and
                                                 presence[x]['protocol_declared'] for x in needed)
                                        for name, needed in states.items()},
        protocol_target_column=rows[0][columns[0]],
        no_cell_expression_coordinates_clinical_or_sequences_read=True,
        not_actual_filtered_assay_or_state_measurement_acceptance=True,
        no_positivity_thresholds_or_new_statistical_tests=True,
        reader_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    output = task / 'private_Monkman2024_declared_marker_coverage_v2.json'
    assert not output.exists()
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
