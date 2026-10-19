"""All cached source arms, exact versus coarse/marker-only endpoints."""
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_external_lineage_inputs_v2 import load, cell_arms
from NSCLC_external_endpoint_text_rules_v1 import CHANNELS, masks, fixture


def therapy_text(value):
    s = str(value).casefold()
    return bool(re.search(r'anti[ _-]?(?:pd[ _-]?1|pd[ _-]?l1|ctla[ _-]?4)|\bici\b|immunotherap|immune[ _-]?checkpoint|pembrolizumab|nivolumab|atezolizumab|ipilimumab|durvalumab|car[ _-]?t|tcr[ _-]?t|vaccine', s))


def summary(frame, columns):
    rows = []
    for key, group in frame.groupby(columns, sort=True, dropna=False):
        values = key if isinstance(key, tuple) else (key,)
        row = dict(zip(columns, values))
        reviewed = group.reviewed_patient.astype(bool)
        available = group.both_endpoints_atleast20
        row.update(libraries=len(group), libraries_receiver_atleast20=int(group.receiver_cells.ge(20).sum()),
                   libraries_both_endpoints_atleast20=int(available.sum()), sender_cells=int(group.sender_cells.sum()),
                   receiver_cells=int(group.receiver_cells.sum()), reviewed_patient_blocks=int(group.loc[reviewed, 'patient_block'].nunique()),
                   reviewed_patient_blocks_both_endpoints=int(group.loc[reviewed & available, 'patient_block'].nunique()),
                   unknown_identity_eligible_libraries=int((~reviewed & available).sum()),
                   any_response_field_in_current_qualification_table=False,
                   source_previously_discovered=bool(str(row['dataset']) == 'GSE243013'))
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 8 * 1024**2)
    out = root / 'NSCLC_external_endpoint_metadata_qualification_20261019_v1'
    assert not out.exists(); out.mkdir(mode=0o700)
    data = load(); lib = data['lib'].copy()
    proof = root / 'complete_external_lineage_opportunities_20261016_v2/private_independent_opportunity_audit.json'
    assert sha(proof) == '2d67723d6174a0301ea7157d038ef08f3aa715e695dd4414e7deafb602313fa7'
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    save('plan_before_all_source_endpoint_metadata.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), rules_sha256=sha(Path(__file__).with_name('NSCLC_external_endpoint_text_rules_v1.py')), resource_guard=guard, hashes=data['hashes'], channels=CHANNELS, scope='All61 source arms/1690 libraries, no outcome/expression/P filtering. Explicit gamma-delta, marker-only and broader myeloid are distinct; suffix grouping only counts, not annotation writes. Current metadata has no response field. Therapy keyword is not endpoint or clinical adoption.', all_individual_records_private=True, no_expression_download=True))
    save('private_rule_fixture.json', fixture())
    counts = np.zeros((len(lib), len(CHANNELS), 2), np.int64)
    diagnostics = {name: np.zeros(len(lib), np.int64) for name in ['missing_author_label', 'mixed_gamma_label', 'explicit_gamma_parent_discordant', 'marker_gamma_parent_discordant']}
    records = []; source_proofs = []
    for chunk in cell_arms(data):
        pairs, flags, classes, vcan, macro = masks(chunk, data['catalogs'][chunk['arm']])
        li = chunk['library_index']
        for channel, pair in enumerate(pairs):
            for endpoint, mask in enumerate(pair):
                counts[:, channel, endpoint] += np.bincount(li[mask], minlength=len(lib))
        for name, mask in flags.items(): diagnostics[name] += np.bincount(li[mask], minlength=len(lib))
        for i, label in enumerate(chunk['author'].categories):
            records.append(dict(arm=chunk['arm'], author_label=str(label), gamma_text_class=int(classes[i]), VCAN_named=bool(vcan[i]), unambiguous_macrophage_text=bool(macro[i])))
        source_proofs.append(dict(arm=chunk['arm'], author_values_sha256=chunk['author_values_sha256'], cells=chunk['cell_count']))
        save('private_progress.json', dict(stage='ALL_SOURCE_AUTHOR_METADATA_COUNTING', completed_arm=chunk['arm'], observed_UTC=datetime.now(timezone.utc).isoformat()))
    lib['therapy_keyword_flag'] = lib.treatment_arm.map(therapy_text)
    for forbidden in ['response', 'response_binary', 'clinical_response', 'outcome', 'response_label']:
        assert forbidden not in lib
    frames = []
    for channel, name in enumerate(CHANNELS):
        frame = lib.copy()
        frame['channel'] = name; frame['sender_cells'] = counts[:, channel, 0]; frame['receiver_cells'] = counts[:, channel, 1]
        frame['both_endpoints_atleast20'] = (counts[:, channel] >= 20).all(1)
        frames.append(frame)
    all_rows = pd.concat(frames, ignore_index=True)
    all_rows.to_csv(out / 'private_all_library_endpoint_counts.tsv.gz', sep='\t', index=False, compression='gzip')
    pd.DataFrame(records).to_csv(out / 'private_all_author_label_rule_classifications.tsv.gz', sep='\t', index=False, compression='gzip')
    pd.DataFrame(dict(arm=lib.arm, dataset=lib.dataset, library_id=lib.library_id, **diagnostics)).to_csv(out / 'private_all_library_label_diagnostics.tsv.gz', sep='\t', index=False, compression='gzip')
    bysource = summary(all_rows, ['arm', 'dataset', 'channel'])
    stratified = summary(all_rows, ['arm', 'dataset', 'channel', 'cancer_label', 'phase_class', 'material_class', 'therapy_keyword_flag'])
    bysource.to_csv(out / 'private_complete_source_endpoint_summary.tsv', sep='\t', index=False)
    stratified.to_csv(out / 'private_complete_source_stratified_endpoint_summary.tsv', sep='\t', index=False)
    save('private_source_metadata_read_proofs.json', source_proofs)
    files = {p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()}
    assert sum(d['bytes'] for d in files.values()) < 4 * 1024**2
    receipt = dict(status='COMPLETE_ALL_SOURCE_ENDPOINT_METADATA_QUALIFICATION_INDEPENDENT_AUDIT_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), files=files, source_metadata_hashes=data['hashes'], source_arms=len(source_proofs), channels=CHANNELS, libraries=len(lib), source_library_arm_rows=int(lib.groupby(['arm', 'dataset']).ngroups), total_metadata_cells=sum(d['cells'] for d in source_proofs), expression_matrix_read=False, clinical_response_adopted=False, P_values_computed=False, scientific_outputs_public=False)
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
