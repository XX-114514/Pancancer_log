#!/usr/bin/env python3
"""Versioned metadata overlay for GSE286827; original QC counts remain read-only.

Only resolves an official title when its patient token has one existing source
mapping across other GEX samples. Does not infer response or survival endpoints.
The TSV outputs are private server artifacts; publish only the summary JSON.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(); run = args.run.resolve()
    args.output_dir.mkdir(exist_ok=False)
    source = run/'manifests/GSE286827_sample_registry.tsv'
    proof = run/'manifests/accepted_GSE286827.json'
    official = run/'metadata/GSE286827/official_geo.json'
    frozen_hashes = {p.name: sha(p) for p in [source, proof, official]}
    reg = pd.read_csv(source, sep='\t', keep_default_na=False)
    geo = json.loads(official.read_text()); by_gsm = {s['gsm_id']: s for s in geo['samples']}
    pattern = re.compile(r'^(P\d+)_(pre|post)(\d*)_GEX$')
    matched = {}
    for _, row in reg[reg.patient_uid.ne('')].iterrows():
        m = pattern.fullmatch(by_gsm[row.library_id]['sample_title'])
        if m: matched.setdefault(m[1], set()).add(row.patient_uid)
    updates = {}
    for i, row in reg[reg.patient_uid.eq('')].iterrows():
        title = by_gsm[row.library_id]['sample_title']; m = pattern.fullmatch(title)
        if not m or len(matched.get(m[1], set())) != 1: continue
        patient = next(iter(matched[m[1]]))
        reg.loc[i, 'patient_uid'] = patient
        reg.loc[i, 'patient_identity_status'] = 'source_linked'
        reg.loc[i, 'timepoint'] = m[2]+m[3]
        reg.loc[i, 'patient_mapping_evidence'] = 'official_GSM_title_patient_token_unique_in_other_source_GEX_samples'
        reg.loc[i, 'timepoint_definition'] = 'official_title_code_no_interval_inferred'
        updates[row.library_id] = (patient, m[2]+m[3])
    if not updates: raise ValueError('No unambiguous official title mapping to repair')
    original = pd.read_csv(source, sep='\t', keep_default_na=False)
    for field in ['response', 'response_source', 'clinical_tier', 'treatment', 'tissue']:
        assert reg[field].equals(original[field]), 'No outcome/tier/treatment changes allowed'
    registry_output = args.output_dir/'GSE286827_sample_registry_patientfix_v2.tsv'
    reg.to_csv(registry_output, sep='\t', index=False)
    z = zarr.open_group(json.loads(proof.read_text())['output'], mode='r')
    obs = z['obs']; ids = np.asarray(read_elem(obs[obs.attrs['_index']])).astype(str)
    libraries = np.asarray(read_elem(obs['original_library_id'])).astype(str)
    affected = np.isin(libraries, list(updates))
    overlay = pd.DataFrame({'stable_cell_id': ids[affected], 'original_library_id': libraries[affected]})
    overlay['patient_uid'] = [updates[lib][0] for lib in overlay.original_library_id]
    overlay['patient_identity_status'] = 'source_linked'
    overlay['timepoint'] = [updates[lib][1] for lib in overlay.original_library_id]
    overlay['evidence_version'] = 'official_title_patientfix_v2'
    overlay_output = args.output_dir/'GSE286827_QC_cell_metadata_overlay_v2.tsv.gz'
    overlay.to_csv(overlay_output, sep='\t', index=False)
    reloaded = pd.read_csv(overlay_output, sep='\t', keep_default_na=False)
    assert reloaded.stable_cell_id.is_unique and np.array_equal(reloaded.stable_cell_id, ids[affected])
    assert len(ids) == 384628 and len(overlay) == 7463 and len(updates) == 1
    assert {p.name: sha(p) for p in [source, proof, official]} == frozen_hashes
    summary = {'schema': 'QC_metadata_overlay_acceptance_v1',
               'observed_utc': datetime.now(timezone.utc).isoformat(), 'dataset': 'GSE286827',
               'status': 'accepted_metadata_overlay_original_QC_preserved',
               'affected_QC_cells': len(overlay), 'libraries_repaired': len(updates),
               'remaining_unresolved_libraries': int(reg.patient_uid.eq('').sum()),
               'source_linked_patients': reg.loc[reg.patient_uid.ne(''), 'patient_uid'].nunique(),
               'original_QC_cells': len(ids), 'overlay_exact_stable_ID_subset': True,
               'outcome_or_tier_changed': False, 'original_counts_modified': False,
               'official_title_rule': 'GEX patient token must have exactly one existing official-source mapping',
               'timepoint_rule': 'Preserve literal official title code; no interval/event inferred',
               'original_evidence_sha256': frozen_hashes,
               'private_output_sha256': {p.name: sha(p) for p in [registry_output, overlay_output]},
               'private_output_directory': '${RUN_ROOT}/artifacts/QC_patient_mapping_repair_20261011_v2',
               'downstream_contract': 'Join overlay by stable_cell_id; retain original registry as historical v1.'}
    with (args.output_dir/'metadata_overlay_acceptance.json').open('x') as handle:
        json.dump(summary, handle, indent=2); handle.write('\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__': main()
