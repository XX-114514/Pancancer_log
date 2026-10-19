"""Fixed public-model physical gates and complete context counts, no RNA."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
from external_functional_coverage_inputs_v1 import read_inputs
from NSCLC_external_lineage_inputs_v2 import load, STRATA
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def resolve(folder, name):
    master = pd.read_csv(folder / 'master_table.tsv', sep='\t', keep_default_na=False).set_index('name').signature
    model = pd.read_csv(folder / (name + '_scGate_Model.tsv'), sep='\t', keep_default_na=False)
    signatures = []
    for row in model.itertuples():
        value = row.signature if row.signature else master.loc[row.name]
        signatures.append({v.rstrip('+-') for v in value.split(';') if v})
    return set.union(*signatures), signatures


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); guard = resource(root, 1024**2)
    src = root / 'all61_fixed_gamma_scGate_native_opportunities_20261019_v1'
    t = json.loads((src / 'private_terminal.json').read_text())
    assert t['status'] == 'COMPLETE_ALL61_FIXED_MODEL_PHYSICAL_AXES_AND_ALL_MATERIAL_CONTEXTS_NO_RNA_OR_CLINICAL_ADOPTION'
    for name, d in t['files'].items(): assert sha(src / name) == d['sha256']
    plan = json.loads((src / 'plan_before_all_source_model_axis_qualification.json').read_text())
    runtime = root / 'gamma_scGate_official_runtime_qualification_20261019_v1'
    assert sha(runtime / 'private_terminal.json') == plan['runtime_sha256']
    modelroot = runtime / 'scGate_models-41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4/human'
    full, signatures = resolve(modelroot / 'HiTME', 'gdT')
    pbmc, _ = resolve(modelroot / 'PBMC', 'gdT')
    nk, _ = resolve(modelroot / 'HiTME', 'NK')
    own = set(pd.read_csv(root / 'GSE243013_genes.csv.gz', keep_default_na=False).geneSymbol.str.strip().str.upper())
    assert all(len(v & own) >= 1 and len(v & own) / len(v) >= .8 for v in signatures)
    required = dict(HiTME_full=full, PBMC_exact=pbmc, HiTME_native_subset=full & own, NK_HiTME_exact=nk)
    assert {k: sorted(v) for k, v in required.items()} == plan['models']
    data = load(); axes = read_inputs(root)
    assert data['hashes'] == plan['input_metadata_hashes'] and axes['input_sha256'] == plan['input_axis_hashes']
    tab = pd.read_csv(src / 'private_complete61_four_fixed_model_physical_opportunities.tsv', sep='\t', keep_default_na=False)
    assert len(tab) == 244 and not tab.duplicated(['arm', 'dataset', 'model']).any()
    for row in tab.itertuples():
        marker_map = axes['source_axes'][(row.arm, row.dataset)]['lookup']
        missing = []; ambiguous = []; columns = []
        for gene in required[row.model]:
            positions = marker_map.get(gene, set())
            if not positions: missing.append(gene)
            elif len(positions) != 1: ambiguous.append(gene)
            else: columns.extend(positions)
        unique = len(columns); collisions = unique - len(set(columns))
        assert row.fixed_features == len(required[row.model]) and row.uniquely_mapped_features == unique
        assert row.unmapped_features == len(missing) and row.ambiguous_features == len(ambiguous) and row.column_collisions == collisions
        assert row.all_fixed_features_uniquely_measurable == (not missing and not ambiguous and not collisions)
        assert json.loads(row.required_missing_or_ambiguous_features) == sorted(missing + ambiguous)
        assert not row.source_dependent_feature_deletion and not row.RNA_identity_or_clinical_accepted
    lib = data['lib'].copy()
    parents = data['counts'][data['counts'].major.isin(['T', 'T cell', 'T/NK lymphoid', 'Lymphoid cell', 'NK', 'T_NK_boundary'])]
    p = parents.groupby(['arm', 'dataset', 'library_id']).n_cells.sum().rename('frozen_lymphoid_parent_cells')
    lib = lib.merge(p.reset_index(), on=['arm', 'dataset', 'library_id'], how='left', validate='one_to_one')
    lib['frozen_lymphoid_parent_cells'] = lib.frozen_lymphoid_parent_cells.fillna(0).astype(int)
    coarse = ['arm', 'dataset', 'cancer_label', 'material_class', 'phase_class', 'capture_scope', 'ex_vivo']
    recorded = pd.read_csv(src / 'private_complete_source_material_capture_phase_contexts.tsv', sep='\t', keep_default_na=False).set_index(coarse)
    groups = lib.groupby(coarse, sort=True, dropna=False); assert len(groups) == len(recorded)
    for key, g in groups:
        row = recorded.loc[key]; reviewed = g.reviewed_patient.astype(bool)
        assert row.libraries == len(g) and row.canonical_native_cells == g.annotated_scope_cells.sum() and row.frozen_lymphoid_parent_cells == g.frozen_lymphoid_parent_cells.sum()
        assert row.reviewed_patient_blocks == g.loc[reviewed, 'patient_block'].nunique() and row.unreviewed_library_units == (~reviewed).sum()
        assert row.candidate_gamma_cells_not_yet_measured and not row.clinical_or_identity_adoption
    out = root / 'all61_fixed_gamma_scGate_native_opportunity_audit_20261019_v1'; assert not out.exists(); out.mkdir(mode=0o700)
    # Preserve exact original tissue/sample-class/treatment/timepoint contexts,
    # alongside the producer's explicitly coarser material/phase summaries.
    exact_fields = ['arm', 'dataset'] + STRATA
    rows = []
    for key, g in lib.groupby(exact_fields, sort=True, dropna=False):
        reviewed = g.reviewed_patient.astype(bool)
        rows.append(dict(zip(exact_fields, key), libraries=len(g), canonical_native_cells=int(g.annotated_scope_cells.sum()), frozen_lymphoid_parent_cells=int(g.frozen_lymphoid_parent_cells.sum()), reviewed_patient_blocks=int(g.loc[reviewed, 'patient_block'].nunique()), unreviewed_library_units=int((~reviewed).sum()), RNA_gamma_identity_not_measured=True, clinical_or_identity_adoption=False))
    exact = pd.DataFrame(rows); assert exact.libraries.sum() == len(lib) == 1690 and exact.canonical_native_cells.sum() == lib.annotated_scope_cells.sum()
    exact.to_csv(out / 'private_all_exact_tissue_sampleclass_treatment_phase_contexts.tsv', sep='\t', index=False)
    receipt = dict(status='ACCEPTED_ALL244_FIXED_PUBLIC_MODEL_NATIVE_GATES_COARSE_COUNTS_AND_COMPLETE_EXACT_SOURCE_CONTEXTS_NO_IDENTITY_ADOPTION', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(src / 'private_terminal.json'), fixed_model_gene_positions=244, library_arm_entries=1690, source_axes_previously_accepted=True, new_expression_read=False, clinical_or_identity_adoption=False, exact_contexts_sha256=sha(out / 'private_all_exact_tissue_sampleclass_treatment_phase_contexts.tsv'), resource_guard=guard, scientific_outputs_public=False)
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file()) < 1024**2
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
