"""Prospective count-only readiness for RNA-defined multicellular programme views.

No programme fitting, clinical test, cell identity promotion or LR selection.
"""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget = resource(root, 8 * 1024**2)
    source = root / 'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    audit = root / 'immune_CM_fixed_panel_same_cell_RNA_audit_20261018_v1/private_terminal.json'
    proof = json.loads(audit.read_text());assert proof['status'].startswith('ACCEPTED_ALL_ORIGINAL_T_CELLS_')
    assert sha(source / 'private_terminal.json') == proof['producer_terminal_sha256']
    for n, v in json.loads((source / 'private_terminal.json').read_text())['files'].items():
        assert sha(source / n) == v['sha256']
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1'
    mp = atlas / 'CM/T_whole_TME_primary_tumor/inputs/sample_metadata_private.tsv'
    meta = pd.read_csv(mp, sep='\t', keep_default_na=False)
    fixed_path = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    fixed = json.loads(fixed_path.read_text())
    units = pd.read_csv(source / 'private_original_T_unit_scope.tsv', sep='\t', keep_default_na=False)
    old_units_path = root / 'existing_immune_fine_native_measurement_scope_20261016_v1/private_all_frozen_T_myeloid_node_units.tsv'
    old_units = pd.read_csv(old_units_path, sep='\t', keep_default_na=False)
    z = np.load(source / 'private_fixed34_gene_T_counts_and_depth.npz', allow_pickle=False)
    raw = z['raw_gene_counts'];genes = z['genes'].tolist();unit = z['original_T_unit_index']
    library = pd.Index(meta.sampleID).get_indexer(units.dataset.astype(str) + '::' + units.library_id.astype(str))[unit]
    core = (raw[:, [genes.index(g) for g in fixed['panels']['pan_T']]] >= 1).sum(axis=1) >= 2
    gene_contract = json.loads((source / 'private_source_gene_measurement_contract.json').read_text())
    by_source = {r['dataset']: r for r in gene_contract}
    complete = np.array([all(by_source[s]['complete_panels'][k] for k in ['pan_T', 'regulatory_associated', 'cytotoxic']) for s in units.dataset])[unit]
    myeloid = old_units[old_units.major.eq('Myeloid')].copy()
    myeloid['sampleID'] = myeloid.dataset.astype(str) + '::' + myeloid.library_id.astype(str)
    my = myeloid.groupby('sampleID').n_cells.sum().reindex(meta.sampleID, fill_value=0).to_numpy()
    out = root / 'immune_CM_RNA_supported_multicell_views_preflight_20261018_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), budget=budget,
        raw_audit_sha256=sha(audit), original_CM_material_sha256=sha(mp), fixed_panel_sha256=sha(fixed_path), original_fine_unit_sha256=sha(old_units_path),
        scope='All118 original CM material libraries, both fixed count thresholds and pre separately. Primary threshold1 fixed before this readiness output; threshold2 sensitivity only.',
        views='Disjoint captured T patterns: FOXP3-supported R without>=2 effector genes; >=2 effector genes without FOXP3-supported R. All original Myeloid nodes form a third broad view. Both-positive and neither patterns recorded, not silently assigned to one view.',
        eligibility='Complete original panT/R/C genes; >=100 captured core-T and>=20 cells in each of the two disjoint T views and broad Myeloid. At least10 eligible libraries per source is source-fit feasibility only; no patient-independent N claim.',
        prospective_method='DIALOGUE or another joint expression-programme approach may be considered only after expression/design/source identifiability. If learning, exclude all34 grouping/clinical marker genes from programme features, preserve patient/source blocks and unknown timepoints. No programme is trained here.',
        limits='Markers define RNA capture views, not validated Treg/CD8 identities. A missing major view is uncaptured/unusable, not no biology. Existing CM/clinical/spatial families remain unchanged; no outcome inspected.')
    (out / 'plan_before_multicell_view_readiness.json').write_text(json.dumps(plan, indent=2) + '\n')
    rows = []
    for threshold in [1, 2]:
        r = ((raw[:, [genes.index(g) for g in fixed['panels']['regulatory_associated']]] >= threshold).sum(axis=1) >= 2) & (raw[:, genes.index('FOXP3')] >= threshold)
        e = (raw[:, [genes.index(g) for g in ['GZMB', 'PRF1', 'CTSW', 'NKG7']]] >= threshold).sum(axis=1) >= 2
        mask = (library >= 0) & core & complete
        state = r.astype(np.int8) + 2 * e.astype(np.int8)
        counts = np.bincount(library[mask] * 4 + state[mask], minlength=len(meta) * 4).reshape(len(meta), 4)
        for i, original in meta.iterrows():
            n = int(counts[i].sum());reg = int(counts[i, 1]);eff = int(counts[i, 2]);mye = int(my[i])
            reason = []
            if n < 100:reason.append('insufficient_complete_core_T_capture')
            if reg < 20:reason.append('insufficient_FOXP3_supported_R_only')
            if eff < 20:reason.append('insufficient_effector_supported_only')
            if mye < 20:reason.append('myeloid_not_captured' if mye == 0 else 'insufficient_myeloid_capture')
            rows.append(dict(sampleID=original.sampleID, dataset=original.dataset, cancer_type=original.cancer_type,
                timepoint=original.timepoint, material_class=original.material_class, threshold=threshold,
                original_CM_patient_verified=original.patient_verified, original_independent_unit=original.independent_unit,
                core_T=n, FOXP3_supported_R_only=reg, effector_only=eff, both_supported=int(counts[i, 3]),
                neither_supported=int(counts[i, 0]), original_myeloid_cells=mye, eligible=not reason, reason=';'.join(reason)))
    table = pd.DataFrame(rows);assert len(table) == len(meta) * 2
    table.to_csv(out / 'private_all_original_libraries_multicell_view_counts.tsv', sep='\t', index=False)
    summary = []
    for (source_name, threshold), group in table.groupby(['dataset', 'threshold'], sort=True):
        for stage, take in [('all_original_CM_material', group), ('explicit_pre_only', group[group.timepoint.eq('pre')])]:
            ok = take[take.eligible]
            reviewed = ok[ok.original_CM_patient_verified.astype(str).str.lower().isin(['true', '1'])]
            summary.append(dict(dataset=source_name, threshold=int(threshold), stage=stage, material_libraries=len(take),
                eligible_libraries=len(ok), source_fit_library_gate=len(ok) >= 10,
                eligible_reviewed_original_identity_units=int(reviewed.original_independent_unit.nunique()),
                independence='Original metadata review status only; repeated libraries may be same patient, no independently validated new cohort.',
                missing_reasons=take.loc[~take.eligible, 'reason'].value_counts().to_dict()))
    (out / 'private_source_level_multicell_readiness.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
    files = [p for p in out.iterdir() if p.is_file()]
    terminal = dict(status='COMPLETE_ALL118_LIBRARIES_BOTH_THRESHOLDS_MULTICELL_VIEW_COUNT_READINESS_NO_FIT',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in files}, no_programme_or_new_clinical_model_trained=True)
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2) + '\n')
    print(json.dumps(dict(status=terminal['status'], completed_UTC=terminal['completed_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
