"""Read a small public GEO SOFT metadata file; leave atlas and outcomes alone."""
import gzip
import hashlib
import json
import os
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'GSE305078_official_sampling_sidecar_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305078/soft/GSE305078_family.soft.gz'
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), public_metadata_url=url,
                scope='Public GEO family metadata only; download <=2MiB compressed, no expression/sequence/prior files. Original source32 libraries only.',
                mapping='Exact original library_id to unique official sample title or accession, with no fuzzy matching or stage inferred from prefix.',
                limits='Source sampling phase can be supported by explicit overall-design statement. Stage from each official sample characteristics. No pCR/response, patient-ID verification or independent validation adopted. Authority remains read-only.')
    (out / 'plan_before_metadata_download.json').write_text(json.dumps(plan, indent=2) + '\n')
    try:
        request = urllib.request.Request(url, headers={'User-Agent': 'PanCancer-research-metadata-review/1.0'})
        with urllib.request.urlopen(request, timeout=35) as response:
            assert response.status == 200
            payload = response.read(2 * 1024 ** 2 + 1)
        assert len(payload) <= 2 * 1024 ** 2 and payload[:2] == b'\x1f\x8b'
        text = gzip.decompress(payload).decode('utf-8')
        assert len(text) < 32 * 1024 ** 2
        (out / 'private_original_GEO_family.soft.gz').write_bytes(payload)
        series = []; samples = []; current = None
        for line in text.splitlines():
            if line.startswith('^SAMPLE = '):
                if current is not None:
                    samples.append(current)
                current = {'GSM': line.split(' = ', 1)[1], 'characteristics': []}
            elif line.startswith('^'):
                if current is not None:
                    samples.append(current); current = None
            elif line.startswith('!Series_overall_design = '):
                series.append(line.split(' = ', 1)[1])
            elif current is not None and ' = ' in line:
                key, value = line.split(' = ', 1)
                if key == '!Sample_title':
                    current['title'] = value
                elif key == '!Sample_characteristics_ch1':
                    current['characteristics'].append(value)
        if current is not None:
            samples.append(current)
        assert len(samples) == 32 and len({r['GSM'] for r in samples}) == 32
        assert len({r['title'] for r in samples}) == 32
        design = '\n'.join(series)
        all_pre = 'all samples were collected before any treatment was administered' in design.lower()
        assert all_pre, 'Explicit all-samples pre-treatment statement absent; do not impute phase'
        lookup = {}
        for sample in samples:
            fields = {}
            for item in sample['characteristics']:
                if ':' in item:
                    k, v = item.split(':', 1)
                    fields[k.strip().lower()] = v.strip()
            sample['stage'] = fields.get('stage', '')
            sample['tissue'] = fields.get('tissue', '')
            for token in [sample['title'], sample['GSM']]:
                assert token not in lookup
                lookup[token] = sample
        libs = pd.read_csv(root / 'existing_immune_cohort_candidate_extension_preflight_20261016_v2/private_all_library_pool_input_coverage.tsv', sep='\t', dtype=str, keep_default_na=False)
        libs = libs[libs.dataset.eq('GSE305078')]
        assert len(libs) == 32 and libs.library_id.is_unique
        rows = []
        for lib in libs.itertuples():
            sample = lookup.get(lib.library_id)
            rows.append(dict(dataset=lib.dataset, library_id=lib.library_id, exact_official_key_matched=sample is not None,
                             official_accession=sample['GSM'] if sample else '', official_title=sample['title'] if sample else '',
                             official_stage=sample['stage'] if sample else '', official_tissue=sample['tissue'] if sample else '',
                             author_sampling_phase='pre-treatment' if sample else '',
                             sampling_phase_evidence='explicit_series_overall_design_all32' if sample else 'original_key_mapping_pending',
                             original_patient_identity_verified=lib.patient_identity_verified, patient_identity_newly_verified=False,
                             clinical_outcome_adopted=False, original_authority_changed=False))
        table = pd.DataFrame(rows)
        table.to_csv(out / 'private_exact_official_sampling_sidecar.tsv', sep='\t', index=False)
        # Parsed public source retained privately; no gene expression or clinical effect.
        (out / 'private_parsed_official_samples.json').write_text(json.dumps(dict(overall_design=design, samples=samples), indent=2) + '\n')
        status = 'EXACT_ALL32_OFFICIAL_SAMPLING_KEYS_AND_PRETREATMENT_PHASE_SUPPORTED_PATIENT_OUTCOME_REVIEW_PENDING' if table.exact_official_key_matched.all() else 'OFFICIAL_METADATA_PARSED_ORIGINAL_LIBRARY_KEY_MAPPING_PENDING_NO_PHASE_IMPUTATION'
        result = dict(status=status, finished_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__),
                      official_samples=len(samples), original_libraries=len(libs), exact_matches=int(table.exact_official_key_matched.sum()),
                      all_samples_pre_treatment_explicitly_stated=all_pre, original_authority_changed=False, patient_identity_or_outcomes_adopted=False,
                      files={p.name: sha(p) for p in out.iterdir() if p.is_file()})
    except Exception as error:
        result = dict(status='METADATA_INPUT_NOT_ACCEPTED_NO_AUTHORITY_OR_OUTCOME_CHANGE', finished_UTC=datetime.now(timezone.utc).isoformat(),
                      code_sha256=sha(__file__), error_type=type(error).__name__, error=str(error), original_authority_changed=False)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
