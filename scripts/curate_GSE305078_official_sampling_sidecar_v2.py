"""Exact dual-token mapping of cached source-library names to official GEO."""
import gzip
import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    old = root / 'GSE305078_official_sampling_sidecar_20261016_v1'
    terminal = json.loads((old / 'private_terminal.json').read_text())
    assert terminal['status'] == 'OFFICIAL_METADATA_PARSED_ORIGINAL_LIBRARY_KEY_MAPPING_PENDING_NO_PHASE_IMPUTATION'
    for name, digest in terminal['files'].items():
        assert sha(old / name) == digest
    public = json.loads((old / 'private_parsed_official_samples.json').read_text())
    assert 'all samples were collected before any treatment was administered' in public['overall_design'].lower()
    official = {r['GSM']: r for r in public['samples']}
    assert len(official) == 32
    source = root / 'existing_immune_cohort_candidate_extension_preflight_20261016_v2/private_all_library_pool_input_coverage.tsv'
    libs = pd.read_csv(source, sep='\t', dtype=str, keep_default_na=False)
    libs = libs[libs.dataset.eq('GSE305078')]
    rows = []
    raw_text = gzip.decompress((old / 'private_original_GEO_family.soft.gz').read_bytes()).decode()
    blocks = {b.splitlines()[0].strip(): b for b in raw_text.split('^SAMPLE = ')[1:]}
    for r in libs.itertuples():
        gsms = set(re.findall(r'(?<![A-Za-z0-9])GSM\d+(?![A-Za-z0-9])', r.library_id))
        titles = set(re.findall(r'(?<![A-Za-z0-9])(?:E|A)\d{3}(?![A-Za-z0-9])', r.library_id))
        assert len(gsms) == len(titles) == 1, 'Library filename does not have unique dual official tokens'
        gsm = next(iter(gsms)); title = next(iter(titles)); sample = official[gsm]
        assert sample['title'] == title
        # Reparse original SOFT independently of v1 JSON parser.
        block = blocks[gsm]
        assert re.search(r'^!Sample_title = (.+)$', block, re.M).group(1) == title
        chars = re.findall(r'^!Sample_characteristics_ch1 = (.+)$', block, re.M)
        fields = {a.split(':', 1)[0].strip().lower(): a.split(':', 1)[1].strip() for a in chars if ':' in a}
        assert fields['stage'] == sample['stage'] and fields['tissue'] == sample['tissue']
        rows.append(dict(dataset=r.dataset, library_id=r.library_id, official_accession=gsm, official_title=title,
                         official_stage=sample['stage'], official_tissue=sample['tissue'], author_sampling_phase='pre-treatment',
                         mapping_evidence='Unique exact GSM AND exact official title tokens in existing source filename, both consistent with original SOFT',
                         sampling_phase_evidence='explicit_series_overall_design_all32',
                         original_patient_identity_verified=r.patient_identity_verified, patient_identity_newly_verified=False,
                         clinical_outcome_adopted=False, original_authority_changed=False))
    table = pd.DataFrame(rows)
    assert len(table) == 32 and table.official_accession.is_unique and table.official_title.is_unique
    assert set(table.official_accession) == set(official)
    out = root / 'GSE305078_official_sampling_sidecar_20261016_v2'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__),
                official_cached_SOFT_sha256=sha(old / 'private_original_GEO_family.soft.gz'), original_library_registry_sha256=sha(source),
                prior_mapping_terminal_sha256=sha(old / 'private_terminal.json'),
                correction='Original source library_id contains both exact official GSM and title inside filename. V1 whole-string matching correctly held0; v2 requires both unique tokens and all32 bijective match, not fuzzy or prefix-based stage inference.',
                clinical_limit='Sampling phase and per-GSM author stage/tissue supported; no patient cross-dataset identity, response, therapy arm or clinical eligibility adopted. No authority mutation or new download.')
    (out / 'plan_and_official_source_provenance.json').write_text(json.dumps(plan, indent=2) + '\n')
    table.to_csv(out / 'private_exact_official_sampling_sidecar.tsv', sep='\t', index=False)
    proof = dict(status='PASS_ALL32_BIJECTIVE_DUAL_OFFICIAL_TOKENS_AND_INDEPENDENT_ORIGINAL_SOFT_STAGE_TISSUE_REPARSE_NO_OUTCOMES',
                 finished_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), exact_mapped_libraries=32,
                 original_authority_changed=False, patient_identity_or_outcome_adopted=False,
                 files={p.name: sha(p) for p in out.iterdir() if p.is_file()})
    (out / 'private_terminal.json').write_text(json.dumps(proof, indent=2) + '\n')
    print(json.dumps(proof))


if __name__ == '__main__':
    main()
