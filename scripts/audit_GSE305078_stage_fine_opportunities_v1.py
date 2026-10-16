"""Replay all stage/fine-pair counts with a pivot, no expression/effects."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'GSE305078_official_sampling_sidecar_20261016_v2'
    terminal = json.loads((out / 'private_terminal.json').read_text())
    assert terminal['status'].startswith('PASS_ALL32_BIJECTIVE_')
    for name, digest in terminal['files'].items():
        assert sha(out / name) == digest
    table_file = out / 'private_all_stage_fine_pair_opportunities.tsv'
    original = json.loads((out / 'private_stage_coverage_audit.json').read_text())
    assert sha(table_file) == original['output_sha256']
    side = pd.read_csv(out / 'private_exact_official_sampling_sidecar.tsv', sep='\t', dtype={'library_id': str})
    assert len(side) == 32 and side.library_id.is_unique and side.official_accession.is_unique
    assert side.author_sampling_phase.eq('pre-treatment').all()
    assert not side.patient_identity_newly_verified.any() and not side.clinical_outcome_adopted.any()
    units = pd.read_csv(root / 'existing_immune_fine_native_measurement_scope_20261016_v1/private_all_frozen_T_myeloid_node_units.tsv', sep='\t', dtype={'library_id': str})
    units = units[units.dataset.eq('GSE305078')]
    counts = units.pivot(index='library_id', columns='node_id', values='n_cells').fillna(0).reindex(side.library_id, fill_value=0)
    pairs = pd.read_csv(root / 'external_fine_pair_case_opportunities_20261016_v1/private_all_source_fine_pair_coverage.tsv', sep='\t')
    pairs = pairs[pairs.dataset.eq('GSE305078')]
    observed = pd.read_csv(table_file, sep='\t')
    keys = {(stage, r.extension_family, r.sender_node, r.receiver_node) for stage in set(side.official_stage) for r in pairs.itertuples()}
    assert len(observed) == len(keys) == 468
    assert keys == set(zip(observed.stage, observed.extension_family, observed.sender_node, observed.receiver_node))
    for row in observed.itertuples():
        ids = side.library_id[side.official_stage.eq(row.stage)]
        a = counts.reindex(columns=[row.sender_node, row.receiver_node], fill_value=0).loc[ids]
        expected = int(a.ge(20).all(axis=1).sum())
        assert row.both20_library_units == expected and row.all_stage_libraries == len(ids)
        assert not row.verified_patient_identity and row.no_effects
    result = dict(status='PASS_ALL468_STAGE_FINE_PAIR_OPPORTUNITIES_BY_INDEPENDENT_PIVOT_AND_ALL32_OFFICIAL_SIDECAR_FLAGS',
                  finished_UTC=datetime.now(timezone.utc).isoformat(), audit_source_sha256=sha(__file__),
                  official_sampling_terminal_sha256=sha(out / 'private_terminal.json'), coverage_proof_sha256=sha(out / 'private_stage_coverage_audit.json'),
                  source_unit_file_sha256=sha(root / 'existing_immune_fine_native_measurement_scope_20261016_v1/private_all_frozen_T_myeloid_node_units.tsv'),
                  stage_rows=len(observed), no_expression_or_effects_or_patient_promotion=True)
    p = out / 'private_independent_stage_coverage_audit.json'
    assert not p.exists()
    p.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
