"""Independent full descriptive review, role counts and prospective scope checks."""
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from measure_immune_RNA_MCP_fixed_lineage_diagnostic_v1 import PANELS


def bound(path, expected):
    assert sha(path / 'private_terminal.json') == expected
    terminal = json.loads((path / 'private_terminal.json').read_text())
    for name, value in terminal['files'].items():
        assert sha(path / name) == (value['sha256'] if isinstance(value, dict) else value)
    return terminal


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 8 * 1024**2)
    review = root / 'immune_RNA_MCP_function_and_atlas_scope_review_20261019_v1'
    bound(review, '4ff6fda25be9ae476692baadbfa3409e537ee9e4cb7be4337d75bf5f6b60e8ba')
    contract = root / 'immune_RNA_MCP_all_atlas_projection_contract_20261019_v1'
    bound(contract, '663c56f3aa8f01049b3d082b62ded4d4f51922cd5f714d56f5387a61ee950aca')
    opp = root / 'immune_RNA_MCP_all_atlas_opportunities_20261019_v1'
    bound(opp, '2a37f7af2f73ac5ceddcb065626df411ca6f016271fce7f0830a1bd658f83cd1')
    synthesis = json.loads((review / 'private_complete_function_scope_and_next_decision.json').read_text())
    functional = root / 'immune_RNA_MCP_functional_description_20261018_v1'
    comparisons = [('all600_functional_annotations', functional / 'private_all600_Hallmark_signed_annotations.tsv'),
                   ('all90_source_background_partitions', functional / 'private_all90_source_sentinel_RNA_component_descriptions.tsv'),
                   ('all54_source_pair_coverage', opp / 'private_all54_source_programme_pair_opportunities.tsv'),
                   ('all_context_coverage', opp / 'private_all_context_programme_pair_opportunities.tsv')]
    for key, path in comparisons:
        keep = key == 'all_context_coverage'
        expected = pd.read_csv(path, sep='\t', keep_default_na=not keep)
        pd.testing.assert_frame_equal(pd.DataFrame(synthesis[key]), expected, check_dtype=False)
    annotation = pd.DataFrame(synthesis['all600_functional_annotations'])
    assert pd.api.types.is_integer_dtype(annotation.signed_active_genes)
    pairs = pd.read_csv(opp / 'private_all2400_material_programme_pair_qualification.tsv.gz', sep='\t')
    scope = pd.read_csv(contract / 'private_all2400_frozen_projection_axis_scope.tsv', sep='\t')
    assert len(scope) == 2400 and not scope.duplicated(['sampleID', 'view', 'program']).any()
    expected_axes = set()
    for p in pairs.itertuples(index=False):
        if p.eligible:
            expected_axes.update((p.sampleID, view, p.program) for view in p.pair.split('|'))
    actual_axes = set(map(tuple, scope.loc[scope.selected_for_projection, ['sampleID', 'view', 'program']].to_numpy()))
    assert actual_axes == expected_axes
    assert scope.loc[scope.selected_for_projection, 'complete_gene_axis'].all()
    plan = json.loads((contract / 'plan_before_any_new_atlas_programme_expression.json').read_text())
    assert plan['no_new_scores_or_effects_or_expression_download']
    assert plan['prospective_assay_output_cap_bytes'] == 256 * 1024**2
    fixed = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    panels = json.loads(fixed.read_text())['panels']
    weights = pd.read_csv(root / 'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv', sep='\t')
    weights = weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    registry = pd.read_csv(root / 'immune_RNA_MCP_disjoint_receiver_preflight_20261018_v1/private_all_unique_directional_endpoint_slots.tsv', sep='\t')
    assert len(registry) == 394
    role_table = pd.read_csv(contract / 'private_frozen_all_programme_LR_and_diagnostic_gene_roles.tsv', sep='\t')
    actual_roles = {r.gene: set(r.roles.split('|')) for r in role_table.itertuples(index=False)}
    expected_roles = {}
    def add(genes, role):
        for gene in genes:
            expected_roles.setdefault(gene, set()).add(role)
    add(set(weights.gene), 'full_frozen_programme_gene')
    add({g for values in panels.values() for g in values}, 'original_fixed34_capture_gene')
    add({g for values in PANELS.values() for g in values}, 'original_fixed43_diagnostic_gene')
    for column, role in [('ligand', 'frozen394_ligand_component'), ('receptor', 'frozen394_receptor_component')]:
        add({g for value in registry[column].astype(str) for g in value.split('_')}, role)
    assert actual_roles == expected_roles
    target = 'CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea'
    eligible = set(pairs.loc[pairs.dataset.eq(target) & pairs.program.eq('MCP2') & pairs.pair.eq('R_FOXP3_capture|E_effector_capture') & pairs.eligible, 'sampleID'])
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'
    z = zarr.open_group(str(atlas), mode='r')
    obs = lambda field: np.asarray(read_elem(z['obs'][field])).astype(str)
    source, libraries, labels = obs('dataset'), obs('library_id'), obs('source_author_label')
    source_mask = source == target
    low = np.char.lower(labels[source_mask])
    dual_label_cells = int(np.sum((np.char.find(low, 'cd4') >= 0) & (np.char.find(low, 'cd8') >= 0)))
    assert dual_label_cells == 0
    old = root / 'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    cache = np.load(old / 'private_fixed34_gene_T_counts_and_depth.npz')
    raw, genes, rows = cache['raw_gene_counts'], cache['genes'].tolist(), cache['global_original_T_row_index']
    positive = lambda names: (raw[:, [genes.index(g) for g in names]] > 0).sum(axis=1)
    core = positive(panels['pan_T']) >= 2
    reg = (positive(panels['regulatory_associated']) >= 2) & (raw[:, genes.index('FOXP3')] > 0)
    eff = positive(['GZMB', 'PRF1', 'CTSW', 'NKG7']) >= 2
    take = (source[rows] == target) & np.isin(libraries[rows], [x.split('::', 1)[1] for x in eligible])
    def role(label):
        label = label.lower()
        for substrings, answer in [(['regulatory', 'treg'], 'author_regulatory_T_annotation'), (['cd4'], 'author_CD4_annotation'), (['cd8'], 'author_CD8_annotation')]:
            if any(s in label for s in substrings):
                return answer
        if 'natural killer' in label or label.startswith('nk'):
            return 'author_NK_annotation'
        if 't cell' in label or 't-cell' in label:
            return 'author_T_annotation_unspecified'
        return 'other_or_unresolved_author_annotation'
    expected_counts = {}
    for view, mask in [('parent_core_T', take & core), ('R_FOXP3_capture', take & core & reg & ~eff), ('E_effector_capture', take & core & eff & ~reg)]:
        expected_counts.update({(view, category): count for category, count in Counter(map(role, labels[rows[mask]])).items()})
    actual = pd.read_csv(review / 'private_CAR_product_author_role_capture_summary.tsv', sep='\t')
    assert {(r.view, r.author_category): r.cells for r in actual.itertuples(index=False)} == expected_counts
    for r in actual.itertuples(index=False):
        n = sum(c for (view, category), c in expected_counts.items() if view == r.view)
        assert r.view_cells == n and abs(r.fraction - r.cells / n) < 1e-14
    columns = ['library_id', 'patient_block', 'source_response', 'native_response_endpoint', 'clinical_endpoint_definition_status', 'clinical_validation_eligible']
    md = pd.DataFrame({field: obs(field)[source_mask] for field in columns})
    clinical = pd.read_csv(review / 'private_all59_CAR_metadata_concordance_no_score_read.tsv', sep='\t', keep_default_na=False).set_index('library_id')
    assert len(clinical) == 59
    for key, group in md.groupby('library_id'):
        row = clinical.loc[key]
        for field in ['patient_block', 'source_response', 'native_response_endpoint', 'clinical_endpoint_definition_status']:
            assert group[field].nunique(dropna=False) == 1
        response = group.source_response.iloc[0].strip().upper()
        assert row.source_response_class == (response if response in ['R', 'NR'] else 'OTHER_OR_UNREVIEWED')
        assert row.patient_block == group.patient_block.iloc[0]
        assert row.native_response_endpoint == group.native_response_endpoint.iloc[0]
        assert row.definition_status == group.clinical_endpoint_definition_status.iloc[0]
        assert row.all_cells_legacy_clinical_eligible == bool(group.clinical_validation_eligible.str.lower().eq('true').all())
        assert row.R_E_capture_eligible == (target + '::' + key in eligible)
    out = root / 'immune_RNA_MCP_atlas_review_and_contract_audit_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    receipt = dict(status='ACCEPTED_COMPLETE_FUNCTION_SCOPE_ROLE_METADATA_AND_PROSPECTIVE_CONTRACT_NO_EFFECTS', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), resource_guard=guard, review_terminal_sha256=sha(review / 'private_terminal.json'), contract_terminal_sha256=sha(contract / 'private_terminal.json'), full_descriptive_records_checked=1008, all_single_axis_positions=2400, selected_axes=len(actual_axes), all_gene_roles_checked=True, all_CAR_material_metadata_checked=59, CAR_dual_CD4_CD8_label_cells=dual_label_cells, no_new_scores_or_effects=True, new_clinical_endpoint_accepted=False)
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'), selected_axes=len(actual_axes))))


if __name__ == '__main__':
    main()
