"""Outcome/effect-free external gene and paired-unit opportunities.

No count matrix, programme values, native association results or clinical
outcomes are loaded. Private endpoint lists and coverage tables stay server-side.
"""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha, VIEWS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def accepted(root, directory, proof_name, expected):
    folder = root / directory
    proof = folder / proof_name
    assert sha(proof) == expected
    terminal = folder / 'private_terminal.json'
    assert json.loads(proof.read_text())['producer_terminal_sha256'] == sha(terminal)
    manifest = json.loads(terminal.read_text())
    for name, item in manifest['files'].items():
        assert sha(folder / name) == item['sha256']
    return folder, proof


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 64 * 1024**2)
    native, native_proof = accepted(root,
        'NSCLC_continuous_prior_native_association_inputs_20261016_v1',
        'private_independent_association_input_audit.json',
        '18917c273d24aff4ed899abc95a5b6debffc893501f783e5d118eb14aa878500')
    pairs, pair_proof = accepted(root,
        'NSCLC_external_programme_pair_inputs_20261016_v1',
        'private_independent_pair_input_audit.json',
        'f5c46a83402c3258e4594aecd69135f690a9e80293b54b7e882eba49367e9ee0')
    cal, cal_proof = accepted(root,
        'NSCLC_native_MFA_external_RNA_universe_calibration_20261016_v1',
        'private_independent_native_RNA_universe_audit.json',
        'f9822656c7be74fa7953b04c8807b5da718bc2084bf3a946e7943504ed31d5b3')
    scope = pd.read_csv(native / 'private_all4800_pre_effect_eligibility.tsv.gz',
                        sep='\t', keep_default_na=False)
    active = scope[scope.test_index.ge(0)].copy()
    assert len(scope) == 4800 and len(active) == 1338
    assert active.test_index.tolist() == list(range(1338))
    masks = pd.read_csv(cal / 'private_native_universe_masks.tsv.gz',
                        sep='\t', index_col=0, keep_default_na=False)
    mapping = pd.read_csv(cal / 'private_all61_source_native_universe_map.tsv',
                          sep='\t', keep_default_na=False)
    assert len(mapping) == 61 and masks.isin([0, 1]).all().all()
    masks = masks.astype(bool)
    slots = pd.read_csv(pairs / 'private_all_frozen_pair_slots_and_coverage.tsv',
                        sep='\t', keep_default_na=False)
    # Read only metadata. Existing effect-valued columns are never loaded.
    units = pd.read_csv(pairs / 'private_all_pair_analysis_units.tsv.gz', sep='\t',
                        usecols=['slot', 'reviewed_patient', 'block_key'],
                        keep_default_na=False)
    assert units.reviewed_patient.dtype == bool
    opportunity = []
    for slot, g in units.groupby('slot', sort=True):
        row = {'slot': slot}
        for label, sub in [('all_units', g), ('reviewed_units', g[g.reviewed_patient])]:
            counts = sub.block_key.value_counts()
            row[label + '_nonsingleton'] = int(counts[counts.ge(2)].sum())
        opportunity.append(row)
    opportunities = pd.DataFrame(opportunity)
    assert slots.slot.is_unique and len(slots) == 7800
    assert set(opportunities.slot) <= set(slots.slot)
    slots = slots.merge(opportunities, on='slot', how='left', validate='one_to_one')
    # Zero here means no existing paired units, not missing RNA or absent biology.
    for c in ['all_units_nonsingleton', 'reviewed_units_nonsingleton']:
        slots[c] = slots[c].fillna(0).astype(int)
    order = {v: i for i, v in enumerate(VIEWS)}
    active['view_A'] = [min((a, b), key=order.get) for a, b in
                        zip(active.sender_view, active.receiver_view)]
    active['view_B'] = [max((a, b), key=order.get) for a, b in
                        zip(active.sender_view, active.receiver_view)]
    out = root / 'NSCLC_native_interface_external_opportunities_20261018_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__), native_input_audit_sha256=sha(native_proof),
        pair_input_audit_sha256=sha(pair_proof), RNA_universe_audit_sha256=sha(cal_proof),
        budget=budget, native_planned_positions=4800, native_effect_blind_eligible=1338,
        scope='Every pre-effect eligible native interface x all61 source-arms and all130 existing contexts. No effect-based candidate/source selection.',
        receptor_rule='At least one original native RNA-qualified reference receptor has every literal component measured. Gene measurement is not same-cell coexpression or protein assembly.',
        opportunity_rule='Existing frozen full-programme paired libraries and patient/library aggregation; remove singleton categorical blocks separately within all-unit and reviewed subsets. >=20 is only a coverage opportunity, not design identifiability or a new test.',
        limitations='Existing full-score universe gates do not certify the LR-masked external response. External pools are frozen subtype extensions, not exact native author views. Same-cell receptor RNA, masked-score calibration, complete model and study overlap remain future gates.',
        no_new_counts_scores_associations_outcomes_downloads_or_identity_changes=True)
    (out / 'plan_before_external_opportunity_read.json').write_text(json.dumps(plan, indent=2) + '\n')
    source_rows = []
    for m in mapping.itertuples():
        measured = set(masks.index[masks['mask_%03d' % m.mask_index]])
        for a in active.itertuples():
            receptors = a.qualifying_reference_receptors.split('|')
            complete = [r for r in receptors if set(r.split('_')) <= measured]
            source_rows.append(dict(record=a.record, arm=m.arm, dataset=m.dataset,
                ligand_measured=a.ligand in measured,
                original_qualifying_receptors=len(receptors),
                complete_measured_receptors=len(complete),
                any_complete_receptor_measured=bool(complete),
                gene_opportunity=(a.ligand in measured and bool(complete)),
                in_model_source=m.dataset == 'GSE243013'))
    sources = pd.DataFrame(source_rows)
    assert len(sources) == 1338 * 61
    contexts = active[['record', 'factor', 'view_A', 'view_B']].merge(
        slots, on=['factor', 'view_A', 'view_B'], how='left', validate='many_to_many')
    assert contexts.groupby('record').size().eq(130).all()
    contexts = contexts.merge(sources, on=['record', 'arm', 'dataset', 'in_model_source'],
                              how='left', validate='many_to_one')
    assert contexts.gene_opportunity.notna().all()
    for label in ['all_units', 'reviewed_units']:
        contexts[label + '_coverage_opportunity'] = (
            contexts.gene_opportunity & contexts[label + '_nonsingleton'].ge(20)
            & ~contexts.in_model_source)
    sources.to_csv(out / 'private_all81618_interface_source_gene_opportunities.tsv.gz',
                   sep='\t', index=False, compression='gzip')
    contexts.to_csv(out / 'private_all173940_interface_context_opportunities.tsv.gz',
                    sep='\t', index=False, compression='gzip')
    terminal = dict(status='COMPLETE_EFFECT_BLIND_EXTERNAL_OPPORTUNITIES_AUDIT_PENDING_NOT_TEST_READINESS',
        finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        source_rows=len(sources), context_rows=len(contexts),
        files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size)
               for p in out.iterdir() if p.is_file()},
        effects_or_P_or_new_expression_read=False, direct_communication_validated=False)
    assert sum(x['bytes'] for x in terminal['files'].values()) < budget['output_cap_bytes']
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2) + '\n')
    print(json.dumps({k: terminal[k] for k in ['status', 'finished_UTC', 'source_rows', 'context_rows']}))
    print(json.dumps({'terminal_sha256': sha(out / 'private_terminal.json')}))


if __name__ == '__main__':
    main()
