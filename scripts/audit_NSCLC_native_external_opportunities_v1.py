"""Independently recheck every external gene/context opportunity, no effects."""
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha, VIEWS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    resource(root, 1024**2)
    out = root / 'NSCLC_native_interface_external_opportunities_20261018_v1'
    terminal = out / 'private_terminal.json'
    manifest = json.loads(terminal.read_text())
    for name, item in manifest['files'].items():
        assert sha(out / name) == item['sha256']
    plan = json.loads((out / 'plan_before_external_opportunity_read.json').read_text())
    assertions = [
        ('NSCLC_continuous_prior_native_association_inputs_20261016_v1',
         'private_independent_association_input_audit.json', 'native_input_audit_sha256'),
        ('NSCLC_external_programme_pair_inputs_20261016_v1',
         'private_independent_pair_input_audit.json', 'pair_input_audit_sha256'),
        ('NSCLC_native_MFA_external_RNA_universe_calibration_20261016_v1',
         'private_independent_native_RNA_universe_audit.json', 'RNA_universe_audit_sha256')]
    for directory, name, key in assertions:
        folder = root / directory
        assert sha(folder / name) == plan[key]
        assert json.loads((folder / name).read_text())['producer_terminal_sha256'] == sha(folder / 'private_terminal.json')
        original = json.loads((folder / 'private_terminal.json').read_text())
        for f, item in original['files'].items():
            assert sha(folder / f) == item['sha256']
    native = root / assertions[0][0]
    pairs = root / assertions[1][0]
    cal = root / assertions[2][0]
    active = pd.read_csv(native / 'private_all4800_pre_effect_eligibility.tsv.gz', sep='\t', keep_default_na=False)
    active = active[active.test_index.ge(0)].set_index('record')
    assert len(active) == 1338 and active.index.is_unique
    masks = pd.read_csv(cal / 'private_native_universe_masks.tsv.gz', sep='\t', index_col=0, keep_default_na=False)
    mapping = pd.read_csv(cal / 'private_all61_source_native_universe_map.tsv', sep='\t').set_index(['arm', 'dataset'])
    observed = pd.read_csv(out / 'private_all81618_interface_source_gene_opportunities.tsv.gz', sep='\t')
    assert len(observed) == 81618 and not observed.duplicated(['record', 'arm', 'dataset']).any()
    expected_genes = {}
    for (arm, dataset), m in mapping.iterrows():
        mask = masks['mask_%03d' % m.mask_index].to_dict()
        for record, a in active.iterrows():
            receptors = a.qualifying_reference_receptors.split('|')
            complete = sum(all(mask.get(g, False) for g in r.split('_')) for r in receptors)
            ligand = bool(mask.get(a.ligand, False))
            expected_genes[(record, arm, dataset)] = (ligand, len(receptors), complete, bool(complete), ligand and bool(complete), dataset == 'GSE243013')
    for row in observed.itertuples():
        got = (row.ligand_measured, row.original_qualifying_receptors,
               row.complete_measured_receptors, row.any_complete_receptor_measured,
               row.gene_opportunity, row.in_model_source)
        assert got == expected_genes[(row.record, row.arm, row.dataset)]
    unit = pd.read_csv(pairs / 'private_all_pair_analysis_units.tsv.gz', sep='\t',
        usecols=['slot', 'reviewed_patient', 'block_key'], keep_default_na=False)
    assert unit.reviewed_patient.dtype == bool
    all_counts = Counter(zip(unit.slot, unit.block_key))
    reviewed = unit[unit.reviewed_patient]
    reviewed_counts = Counter(zip(reviewed.slot, reviewed.block_key))
    totals = [Counter(), Counter()]
    for dest, counts in zip(totals, [all_counts, reviewed_counts]):
        for (slot, block), n in counts.items():
            if n >= 2:
                dest[slot] += n
    original = pd.read_csv(pairs / 'private_all_frozen_pair_slots_and_coverage.tsv', sep='\t', keep_default_na=False).set_index('slot')
    contexts = pd.read_csv(out / 'private_all173940_interface_context_opportunities.tsv.gz', sep='\t', keep_default_na=False)
    assert len(contexts) == 173940 and not contexts.duplicated(['record', 'context_id']).any()
    assert contexts.groupby('record').size().eq(130).all()
    contexts = contexts.merge(active[['sender_view', 'receiver_view']], left_on='record', right_index=True, validate='many_to_one')
    for c in original.columns:
        assert contexts[c].tolist() == original.loc[contexts.slot, c].tolist(), c
    for row in contexts.itertuples():
        assert sorted([row.sender_view, row.receiver_view], key=VIEWS.index) == [row.view_A, row.view_B]
        flags = expected_genes[(row.record, row.arm, row.dataset)]
        assert row.gene_opportunity == flags[4]
        for label, counts in zip(['all_units', 'reviewed_units'], totals):
            n = counts[row.slot]
            assert getattr(row, label + '_nonsingleton') == n
            assert getattr(row, label + '_coverage_opportunity') == (flags[4] and n >= 20 and not flags[5])
    report = dict(status='PASS_ALL81618_GENE_AND173940_CONTEXT_OPPORTUNITIES_NO_EFFECTS',
        finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        producer_terminal_sha256=sha(terminal), source_rows=len(observed), context_rows=len(contexts),
        no_new_counts_scores_P_or_outcomes_loaded=True,
        full_external_test_readiness_or_direct_communication=False)
    proof = out / 'private_independent_external_opportunity_audit.json'
    assert not proof.exists()
    proof.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(status=report['status'], finished_UTC=report['finished_UTC'], audit_sha256=sha(proof))))


if __name__ == '__main__':
    main()
