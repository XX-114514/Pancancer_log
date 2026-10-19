"""Independent all-reference symbol, exclusion and coverage reconstruction."""
import json
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 512 * 1024)
    out = root / 'NSCLC92_receiver_reference_scope_audit_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    src = root / 'NSCLC92_receiver_function_reference_scope_20261019_v1'
    terminal = json.loads((src / 'private_terminal.json').read_text())
    assert sha(src / 'private_terminal.json') == 'fd4d775a3f4094f0359e020b72d4e85c2e546d13df50b9d0edf937cd3d909f2b'
    for name, item in terminal['files'].items():
        assert sha(src / name) == item['sha256']
    reference = Path(__file__).parent / 'hallmark50_public_reference_v3'
    manifest = reference / 'reference_provenance.json'
    plan = json.loads((src / 'plan_before_receiver_reference_gene_qualification.json').read_text())
    assert sha(manifest) == plan['public_reference_manifest_sha256']
    definitions = json.loads((src / 'private_all50_LR_disjoint_reference_definitions.json').read_text())
    recorded = {d['reference']: d for d in definitions}
    summary = pd.read_csv(src / 'private_all50_reference_gene_qualification.tsv', sep='\t').set_index('reference')
    gene_file = root / 'GSE243013_genes.csv.gz'
    assert sha(gene_file) == 'ac325c9ce6bfaf2270eb931e7d947bf31d1fb151b14cb166ed9f6c8f4b34a104'
    genes = pd.read_csv(gene_file, keep_default_na=False).geneSymbol.tolist()
    positions = defaultdict(list)
    for position, symbol in enumerate(genes):
        positions[str(symbol).strip().upper()].append(position)
    tensor = root / 'NSCLC_communication_tensor_family_20261019_v1'
    tt = json.loads((tensor / 'private_terminal.json').read_text())
    axis_path = tensor / 'private_tensor_axes.npz'
    assert sha(axis_path) == tt['files'][axis_path.name]['sha256']
    blocked = set(np.load(axis_path)['native_gene_columns'].astype(int).tolist())
    excluded = {str(genes[c]).strip().upper() for c in blocked}
    all_records = json.loads(manifest.read_text())['records']
    assert len(all_records) == len(recorded) == len(summary) == 50
    union = set()
    for record in all_records:
        name = record['name']
        path = reference / (name + '.json')
        assert sha(path) == record['sha256']
        raw = json.loads(path.read_text())[name]['geneSymbols']
        full = sorted(set(str(g).strip().upper() for g in raw))
        removed, remaining, unique, missing, ambiguous = [], [], {}, [], []
        for gene in full:
            candidates = positions.get(gene, [])
            if gene in excluded or (len(candidates) == 1 and candidates[0] in blocked):
                removed.append(gene)
                continue
            remaining.append(gene)
            if len(candidates) == 1:
                unique[gene] = candidates[0]
            elif not candidates:
                missing.append(gene)
            else:
                ambiguous.append(gene)
        d = recorded[name]
        for key, value in [('full_reference_genes', full), ('excluded_all_LR_symbols', removed),
                           ('remaining_genes', remaining), ('unique_native_columns', unique),
                           ('missing_source_symbols', missing), ('ambiguous_source_symbols', ambiguous)]:
            assert d[key] == value, (name, key)
        assert not blocked.intersection(unique.values())
        eligible = len(unique) >= 20 and len(unique) >= .8 * len(remaining)
        strong = eligible and len(remaining) >= .8 * len(full)
        assert d['complement_measurement_eligible'] == eligible
        assert d['strong_reference_representation'] == strong
        s = summary.loc[name]
        assert [s.original_genes, s.excluded_LR_genes, s.remaining_genes, s.uniquely_measured_genes] == [len(full), len(removed), len(remaining), len(unique)]
        np.testing.assert_allclose([s.retained_fraction, s.remaining_gene_measurement_fraction], [len(remaining) / len(full), len(unique) / len(remaining) if remaining else 0], atol=1e-15, rtol=1e-15)
        assert bool(s.complement_measurement_eligible) == eligible
        assert bool(s.strong_reference_representation) == strong
        union.update(unique.values())
    old = root / 'NSCLC_full_native_gene_unit_audit_20261016_v1/private_terminal.json'
    assert sha(old) == 'bb20365ee83758a7df8478ceadede2596fe749928dab5765e3bc0866d72dd8c1'
    proof = dict(status='ACCEPTED_ALL50_REFERENCE_NATIVE_MAPPING_EXCLUSIONS_AND_COVERAGE',
                 completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                 producer_terminal_sha256=sha(src / 'private_terminal.json'), reference_slots=50,
                 all_mapped_positions_checked=True, all_LR_exclusion_checked=True,
                 old_native_raw_audit_sha256=sha(old), old_raw_replay_units=64,
                 old_raw_replay_all_native_genes=True, all_native_units_raw_replayed=False,
                 resource_guard=guard, new_expression_scores=False, new_functional_P_values=False,
                 scientific_outputs_public=False)
    (out / 'private_terminal.json').write_text(json.dumps(proof, indent=2) + '\n')
    print(json.dumps(dict(status=proof['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
