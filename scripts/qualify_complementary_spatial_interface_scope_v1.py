"""Complete existing spatial molecular/coordinate scope, no RNA effects."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 256 * 1024)
    request = root / 'private_complementary_spatial_interface_scope_request_20261019_v1.json'
    assert sha(request) == 'f31d71f1efdc8450f346020733b0edea354a962eeec57bcc58328ba36bd890d0'
    fixed = json.loads(request.read_text())
    manifest = root / 'immune_same_name_family_source_review_v1/spatial_native_input_repair_v2/analysis_input_manifest.tsv'
    assert sha(manifest) == '6fda649c8b2bf4abe010ba8d53d50e495bae29ce8d30648675d9be9bf91dcf83'
    files = pd.read_csv(manifest, sep='\t', keep_default_na=False)
    assert len(files) == 14 and not files.sample_id.duplicated().any()
    state = Path(json.loads((root / 'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    axis_path = state / 'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json'
    assert sha(axis_path) == 'b4818461998e0130d0db8b7a3ce2b75605c2dd13caeec1c39acd85b3a3caa589'
    axes = {v['sample_id']: v for v in json.loads(axis_path.read_text())}
    assert set(axes) == set(files.sample_id)
    old = root / 'immune_same_name_family_source_review_v1/spatial_contact_native_preflight_v3'
    accepted = json.loads((old / 'independent_audit.json').read_text())
    assert sha(old / 'inventory.json') == accepted['producer_inventory_sha256']
    old_inputs = {v['sample_id']: v for v in json.loads((old / 'inventory.json').read_text())['inputs']}
    repair = root / 'GSE299886_original_coordinate_sidecars_20261018_v1'
    assert sha(repair / 'private_terminal.json') == 'a1731ea5731ccc62fe8154cfb6b2d0b87d02cde12d235ed0461bad54a247e528'
    rt = json.loads((repair / 'private_terminal.json').read_text())
    restored = {v['sample_id']: v for v in json.loads((repair / 'private_all_source_coordinate_reconstruction_receipt.json').read_text())}
    out = root / 'complementary_spatial_interface_scope_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    save('plan_before_spatial_interface_input_reads.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), private_request_sha256=sha(request), manifest_sha256=sha(manifest), panel_axes_sha256=sha(axis_path), source_inputs=14, target_positions=len(fixed['targets']), resource_guard=guard, scope=fixed['scope'], expression_values_or_spatial_effects_read=False))
    records = []
    proofs = []
    for source in files.itertuples():
        path = Path(source.analysis_path)
        ax = axes[source.sample_id]
        current = sha(path)
        assert current == ax['input_sha256']
        lookup = ax['symbol_to_feature_index']
        symbols = set(ax['unambiguous_symbols'])
        assert symbols == set(lookup)
        with h5py.File(path, 'r') as h:
            obs = h['obs']
            idx = obs.attrs['_index']
            idx = idx.decode() if isinstance(idx, bytes) else str(idx)
            item = obs[idx]
            n = int(item.shape[0] if hasattr(item, 'shape') else item['codes'].shape[0])
            xy = np.asarray(read_elem(h['obsm']['spatial']))
            assert xy.shape == (n, 2)
            old_finite = int(np.isfinite(xy).all(axis=1).sum())
            fov_known = np.zeros(n, bool)
            if 'fov' in obs:
                fov = pd.Series(read_elem(obs['fov']))
                fov_known = (fov.notna() & ~fov.astype(str).str.strip().str.upper().isin(['', 'NA', 'NAN', 'NONE', 'NULL', '<NA>'])).to_numpy()
            authority = 'CURRENT_COORDINATE_SCHEMA_ONLY_RAW_GEOMETRY_NOT_ACCEPTED'
            if source.sample_id in old_inputs:
                assert current == old_inputs[source.sample_id]['sha256']
                authority = 'PREVIOUS_FIVE_SOURCE_NATIVE_COORDINATE_INPUT_PROOF'
            sidecar_sha = None
            if source.sample_id in restored:
                rr = restored[source.sample_id]
                sidecar = Path(rr['output'])
                assert sha(sidecar) == rt['files'][sidecar.name]['sha256']
                z = dict(np.load(sidecar, allow_pickle=False))
                assert current == rr['matrix_SHA256'] == str(z['existing_matrix_sha256'].item())
                np.testing.assert_array_equal(z['existing_row_index'], np.arange(n))
                xy = z['xy_source']
                fov_known = z['FOV_known']
                assert xy.shape == (n, 2) and np.isfinite(xy).all() and len(fov_known) == n
                authority = 'ORIGINAL_ARCHIVE_EXACT_ROW_REPAIRED_COORDINATE_SIDECAR'
                sidecar_sha = sha(sidecar)
            finite = int(np.isfinite(xy).all(axis=1).sum())
            proofs.append(dict(dataset_id=source.dataset_id, sample_id=source.sample_id, current_input_sha256=current, current_coordinate_finite_rows=old_finite, selected_coordinate_authority=authority, selected_coordinate_finite_rows=finite, total_rows=n, FOV_known_rows=int(fov_known.sum()), original_coordinate_sidecar_sha256=sidecar_sha, absolute_physical_unit_validated=False, patient_core_identity_not_upgraded=True))
            for target in fixed['targets']:
                ligand = target['original_ligand_genes']
                receptor = target['original_receptor_genes']
                strict = ligand + target['strict_ligand_assembly_addition']
                t = fixed['lineage_proxies']['T_like']
                my = fixed['lineage_proxies']['Myeloid_like']
                all_required = set(strict + receptor + t + my)
                columns = [lookup[v] for v in all_required if v in symbols]
                collisions = len(columns) - len(set(columns))
                complete = all_required.issubset(symbols) and collisions == 0
                records.append(dict(dataset_id=source.dataset_id, sample_id=source.sample_id, cancer=source.cancer, scope_id=target['axis_id'], native_cells=n, features=ax['features'], original_resource_endpoints_measurable=set(ligand + receptor).issubset(symbols), strict_ligand_assembly_measurable=set(strict).issubset(symbols), T_proxy_pair_measurable=set(t).issubset(symbols), Myeloid_proxy_pair_measurable=set(my).issubset(symbols), all_fixed_molecules_and_lineage_proxies_measurable=complete, physical_column_collisions=collisions, original_gamma_identity_validated=False, selected_coordinate_finite_rows=finite, selected_geometry_has_source_proof=authority != 'CURRENT_COORDINATE_SCHEMA_ONLY_RAW_GEOMETRY_NOT_ACCEPTED', FOV_known_rows=int(fov_known.sum()), missing_or_ambiguous_required_features=json.dumps(sorted(all_required - symbols)), independent_patient_mapping_reviewed=str(source.independent_patient_mapping_reviewed).lower() == 'true', actual_expression_or_neighbourhood_measured=False))
    assert len(records) == 14 * len(fixed['targets']) and len(proofs) == 14
    pd.DataFrame(records).to_csv(out / 'private_complete14_interface_molecular_and_geometry_scope.tsv', sep='\t', index=False)
    save('private_all14_coordinate_input_authority.json', proofs)
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file()) < 256 * 1024
    receipt = dict(status='COMPLETE_ALL14_COMPLEMENTARY_INTERFACE_MEASURABILITY_AND_COORDINATE_AUTHORITY_NO_SPATIAL_EFFECTS', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), private_request_sha256=sha(request), old_coordinate_parent_audit_sha256=sha(old / 'independent_audit.json'), restored_coordinate_terminal_sha256=sha(repair / 'private_terminal.json'), resource_guard=guard, full_current_input_hashes_checked=14, exact_gamma_identity_or_new_clinical_validation=False, scientific_outputs_public=False, files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
