"""Freeze all external family/phase/fine-node inputs before association."""
import hashlib
import itertools
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from external_core_rank_model_v2 import geometry


CONTROL_NAMES = ['sender_fraction', 'sender_logN', 'sender_mean_logRNA', 'sender_mean_loggenes',
                 'receiver_fraction', 'receiver_logN', 'receiver_mean_logRNA', 'receiver_mean_loggenes',
                 'receiver_own_ligand', 'max_qualified_receptor_joint_fraction']


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as stream:
        for block in iter(lambda: stream.read(8388608), b''):
            h.update(block)
    return h.hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    config_file = root / 'external_ligand_programme_core_followup_config_20261016_v1.json'
    config = json.loads(config_file.read_text())
    native = root / 'existing_immune_fine_native_signatures_20261016_v4'
    side_root = root / 'GSE305078_official_sampling_sidecar_20261016_v2'
    opportunity = root / 'external_fine_pair_case_opportunities_20261016_v1'
    for p in [native / 'private_independent_audit.json', side_root / 'private_independent_stage_coverage_audit.json', opportunity / 'private_independent_audit.json']:
        assert json.loads(p.read_text())['status'].startswith('PASS_')
    terminal = json.loads((native / 'private_terminal.json').read_text())
    for name, digest in terminal['files'].items():
        assert sha(native / name) == digest
    units = pd.read_csv(native / 'private_all_original_fine_units_and_quality.tsv', sep='\t', dtype={'library_id': str})
    units = units[units.dataset.eq(config['source'])]
    side = pd.read_csv(side_root / 'private_exact_official_sampling_sidecar.tsv', sep='\t', dtype={'library_id': str})
    assert len(side) == 32 and side.library_id.is_unique and side.author_sampling_phase.eq('pre-treatment').all()
    assert not side.patient_identity_newly_verified.any() and not side.clinical_outcome_adopted.any()
    families = pd.read_csv(opportunity / 'private_all_deduplicated_extension_families.tsv', sep='\t')
    pairs = pd.read_csv(opportunity / 'private_all_source_fine_pair_coverage.tsv', sep='\t')
    pairs = pairs[pairs.dataset.eq(config['source'])].sort_values(['extension_family', 'sender_node', 'receiver_node'])
    assert len(families) == 4 and len(pairs) == 234
    family_lookup = families.set_index('extension_family')
    pm = pd.read_csv(native / 'private_all_source_parent_programme_scope.tsv.gz', sep='\t')
    pm = pm[pm.dataset.eq(config['source'])]
    tm = pd.read_csv(native / 'private_all_source_parent_TF_scope.tsv.gz', sep='\t')
    tm = tm[tm.dataset.eq(config['source'])]
    axis = pd.read_csv(native / 'private_all_source_panel_gene_axis.tsv.gz', sep='\t')
    axis = axis[axis.dataset.eq(config['source'])].set_index('gene')
    score_archive = np.load(native / 'private_frozen_programme_and_TF_scores.npz', allow_pickle=False)
    P = score_archive['programme']; T = score_archive['TF']
    panel_archive = np.load(native / 'private_limited_native_panel.npz', allow_pickle=False)
    L = panel_archive['mean_logCP10K']; C = panel_archive['positive_cell_counts']
    out = root / 'GSE305078_external_ligand_programme_core_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700); (out / 'private_designs').mkdir(mode=0o700)
    inputs = dict(config=config_file, native_audit=native / 'private_independent_audit.json', native_terminal=native / 'private_terminal.json',
                  official_sidecar=side_root / 'private_exact_official_sampling_sidecar.tsv', official_audit=side_root / 'private_independent_stage_coverage_audit.json',
                  all_pairs=opportunity / 'private_all_source_fine_pair_coverage.tsv', families=opportunity / 'private_all_deduplicated_extension_families.tsv')
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), kernel_sha256=sha(Path(__file__).with_name('external_core_rank_model_v2.py')),
                input_sha256={name: sha(p) for name, p in inputs.items()}, complete_family=config['complete_family'], config_sha256=sha(config_file),
                primary_slots=config['expected_primary_slots'], secondary_slots=config['expected_secondary_slots'], all_control_names=CONTROL_NAMES,
                controls='Complete10 technical/composition/own-ligand/receiver-RNA variables and intercept; SVD column-space projection, preserve original columns and record effective rank and strict full-column status. Do not delete variables according to effects.',
                interpretation='Lineage extension, not exact author-state transfer or clinical prediction. Source original Myeloid/T denominators retained; original discovery macrophage denominator is not an identical population.',
                clinical='Original clinical8 unavailable; no silent reduced version of old18 model. This is explicitly new technical-core analysis with unreviewed library-fallback units, no response read.',
                missingness='Gene availability and frozen score gates explicit, measured zeros retained. No complex joint fraction inferred from marginal RNA. Every planned stage/node/readout slot retained.',
                no_association_P_BH_or_effect_computed=True)
    (out / 'plan_before_any_external_association.json').write_text(json.dumps(plan, indent=2) + '\n')
    rows = []; membership = []; hashes = []
    for stage, pair in itertools.product(config['stages'], pairs.itertuples()):
        family = family_lookup.loc[pair.extension_family]
        parents = family.original_parent_records.split('|'); parent = parents[0]
        eligible_libraries = set(side.loc[side.official_stage.eq(stage), 'library_id'])
        u = units[units.library_id.isin(eligible_libraries) & units.n_cells.ge(20)]
        ss = u[u.node_id.eq(pair.sender_node)].set_index('library_id')
        rr = u[u.node_id.eq(pair.receiver_node)].set_index('library_id')
        libraries = sorted(set(ss.index) & set(rr.index))
        si = ss.loc[libraries, 'unit_index'].to_numpy(int); ri = rr.loc[libraries, 'unit_index'].to_numpy(int)
        kinds = ['masked_programme']
        if family.ligand == config['secondary_ligand']:
            kinds.append('fixed_functional_followup')
        for kind in kinds:
            model_id = f'EXT{len(rows)+1:04d}'
            row = dict(model_id=model_id, stage=stage, extension_family=pair.extension_family, original_parent_records=family.original_parent_records,
                       context=family.context, program=family.program, ligand=family.ligand, original_qualified_receptors=family.qualified_receptors,
                       sender_node=pair.sender_node, receiver_node=pair.receiver_node, readout=kind,
                       regulator=config['secondary_regulator'] if kind != 'masked_programme' else '', library_units=len(libraries),
                       status='READY_FOR_INDEPENDENT_INPUT_AUDIT', original_columns=12, effective_rank=np.nan, nuisance_rank=np.nan, residual_df=np.nan,
                       strict_full_column_rank=False, constant_controls='', no_patient_identity_or_endpoint_promotion=True)
            for k, lib in enumerate(libraries):
                membership.append(dict(model_id=model_id, library_id=lib, sender_unit=int(si[k]), receiver_unit=int(ri[k]), analysis_unit='unreviewed_library_fallback'))
            try:
                if len(libraries) < config['min_units']:
                    raise ValueError('INSUFFICIENT_LIBRARY_UNITS')
                ligand_axis = axis.loc[family.ligand]
                if not ligand_axis.source_native_unambiguous:
                    raise ValueError('LIGAND_UNMEASURED_OR_AMBIGUOUS')
                column = int(ligand_axis.panel_column)
                X = L[si, column]; own = L[ri, column]
                assert np.array_equal(X > 0, C[si, column] > 0)
                row['ligand_positive_libraries'] = int((X > 0).sum())
                if row['ligand_positive_libraries'] < config['min_ligand_positive_libraries']:
                    raise ValueError('INSUFFICIENT_LIGAND_POSITIVE_LIBRARIES')
                fractions = []
                for receptor in family.qualified_receptors.split('|'):
                    if '_' in receptor:
                        raise ValueError('COMPLEX_JOINT_RNA_NOT_AVAILABLE_FROM_MARGINAL_PANEL')
                    if receptor not in axis.index or not axis.loc[receptor, 'source_native_unambiguous']:
                        raise ValueError('ORIGINAL_RECEPTOR_UNMEASURED_OR_AMBIGUOUS')
                    rc = int(axis.loc[receptor, 'panel_column'])
                    fractions.append(C[ri, rc] / rr.loc[libraries, 'n_cells'].to_numpy())
                receptor_pool = [float((f * rr.loc[libraries, 'n_cells'].to_numpy()).sum() / rr.loc[libraries, 'n_cells'].sum()) for f in fractions]
                row['maximum_original_receptor_pooled_RNA_fraction'] = max(receptor_pool)
                if max(receptor_pool) < .01:
                    raise ValueError('NO_ORIGINAL_RECEPTOR_ONE_PERCENT_POOLED_RNA')
                roles = {}
                for role in ['sender', 'receiver_masked']:
                    candidates = pm[pm.record.isin(parents) & pm.role.eq(role)]
                    assert len(candidates) == len(parents) and candidates.score_column.nunique() == 1
                    if not candidates.status.str.startswith('MEASURED_').all():
                        raise ValueError('FROZEN_PROGRAMME_COVERAGE_FAILED')
                    roles[role] = int(candidates.score_column.iloc[0])
                sender_programme = P[si, roles['sender']]
                same_regulator = np.full(len(libraries), np.nan)
                if kind == 'masked_programme':
                    Y = P[ri, roles['receiver_masked']]
                else:
                    tf = tm[tm.record.isin(parents) & tm.TF.eq(config['secondary_regulator'])]
                    assert len(tf) == len(parents) and tf.score_column.nunique() == 1
                    if not tf.status.str.startswith('MEASURED_').all() or (tf.target_retained_fraction < .95).any() or (tf.background_retained_fraction < .95).any():
                        raise ValueError('INSUFFICIENT_MATCHED_FUNCTIONAL_TARGET_OR_BACKGROUND')
                    tc = int(tf.score_column.iloc[0]); Y = T[ri, tc]; same_regulator = T[si, tc]
                cv = []
                for frame in [ss.loc[libraries], rr.loc[libraries]]:
                    cv.extend([frame.within_original_major_fraction.to_numpy(), np.log1p(frame.n_cells.to_numpy()), frame.mean_logRNA.to_numpy(), frame.mean_loggenes.to_numpy()])
                controls = np.column_stack([*cv, own, np.max(np.column_stack(fractions), axis=1)])
                assert controls.shape == (len(libraries), 10)
                if not np.isfinite(Y).all() or np.ptp(Y) <= 1e-12:
                    raise ValueError('MISSING_OR_CONSTANT_RECEIVER_READOUT')
                g = geometry(X, controls, config['min_units'], config['min_residual_df'])
                yy = rankdata(Y).astype(float); yy -= g['q0'] @ (g['q0'].T @ yy)
                if np.linalg.norm(yy) <= 1e-8:
                    raise ValueError('NO_RECEIVER_VARIATION_AFTER_CONTROLS')
                row.update(effective_rank=g['rank'], nuisance_rank=g['nuisance_rank'], residual_df=g['residual_df'], condition=g['condition'],
                           strict_full_column_rank=g['strict_full_column_rank'], constant_controls='|'.join(np.asarray(CONTROL_NAMES)[g['constant_controls']]),
                           maximum_leverage=float(g['h'].max()))
                path = out / 'private_designs' / (model_id + '.npz')
                np.savez_compressed(path, X=X, Y=Y, C_raw=controls, covariate_names=np.asarray(CONTROL_NAMES), library_ids=np.asarray(libraries),
                                    sender_units=si, receiver_units=ri, sender_programme=sender_programme, sender_same_regulator=same_regulator)
                hashes.append(dict(model_id=model_id, sha256=sha(path)))
            except ValueError as error:
                row['status'] = str(error)
            rows.append(row)
    table = pd.DataFrame(rows)
    assert len(table) == config['complete_family'] and table.model_id.is_unique
    assert table.readout.eq('masked_programme').sum() == config['expected_primary_slots']
    assert table.readout.eq('fixed_functional_followup').sum() == config['expected_secondary_slots']
    table.to_csv(out / 'private_all676_frozen_model_scope.tsv', sep='\t', index=False)
    pd.DataFrame(membership).to_csv(out / 'private_all_model_library_membership.tsv.gz', sep='\t', index=False, compression='gzip')
    pd.DataFrame(hashes, columns=['model_id', 'sha256']).to_csv(out / 'private_ready_design_SHA256.tsv', sep='\t', index=False)
    result = dict(status='COMPLETE_ALL676_EXTERNAL_MODEL_INPUTS_CONTROL_SPACE_AND_UNAVAILABLE_REASONS_NO_EFFECTS',
                  finished_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), kernel_sha256=plan['kernel_sha256'],
                  complete_family=len(table), ready_designs=len(hashes), status_counts=table.status.value_counts().to_dict(),
                  files={p.name: sha(p) for p in out.iterdir() if p.is_file()}, no_external_associations_P_BH_or_clinical_outcomes=True)
    (out / 'private_input_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
