"""Existing-panel feature feasibility only; no molecular names in stdout."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 512 * 1024)
    out = root / 'NSCLC92_existing_spatial_panel_qualification_20261019_v1'
    assert not out.exists(); out.mkdir(mode=0o700)
    state = Path(json.loads((root / 'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    source = state / 'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json'
    assert sha(source) == 'b4818461998e0130d0db8b7a3ce2b75605c2dd13caeec1c39acd85b3a3caa589'
    panels = json.loads(source.read_text()); assert len(panels) == 14
    targets = pd.read_csv(root / 'NSCLC92_receiver_function_inputs_20261019_v1/private_complete_old_BH_receiver_scope.tsv', sep='\t')
    vocab = json.loads((root / 'NSCLC_communication_tensor_family_20261019_v1/private_selected_LR_vocabulary.json').read_text())
    vocab = {v['canonical_resource_LR_id']: v for v in vocab}
    native = pd.read_csv(root / 'GSE243013_genes.csv.gz', keep_default_na=False).geneSymbol.str.strip().str.upper().tolist()
    references = json.loads((root / 'NSCLC92_receiver_function_reference_scope_20261019_v1/private_all50_LR_disjoint_reference_definitions.json').read_text())
    rows = []; detailed = []
    for number, panel in enumerate(panels, start=1):
        features = panel['features']
        assert isinstance(features, list) and all(isinstance(v, str) for v in features)
        feature_set = {v.strip().upper() for v in features}
        for target in targets.itertuples():
            lr = vocab[target.canonical_resource_LR_id]
            left = {native[c] for c in lr['native_ligand_columns']}
            right = {native[c] for c in lr['native_receptor_columns']}
            eligible = 0
            for d in references:
                genes = set(d['unique_native_columns'])
                overlap = genes.intersection(feature_set)
                good = len(overlap) >= 20 and len(overlap) >= .8 * len(genes)
                eligible += good
                detailed.append(dict(panel_index=number, scope_id=target.axis_id, reference=d['reference'], native_complement_genes=len(genes), panel_measured_genes=len(overlap), complement_panel_coverage=len(overlap) / len(genes) if genes else 0., complement_measurement_qualified=good))
            rows.append(dict(panel_index=number, scope_id=target.axis_id, measured_features=len(feature_set), all_resource_ligand_components_measured=left.issubset(feature_set), all_resource_receptor_components_measured=right.issubset(feature_set), both_resource_endpoints_measured=(left | right).issubset(feature_set), gamma_delta_lineage_marker_count=len({'TRDC', 'TRDV1', 'TRGC1', 'TRGC2'}.intersection(feature_set)), functional_complements_with_80percent_native_coverage=eligible, original_references=50, cell_level_resolution_and_lineage_not_qualified=True, actual_RNA_expression_not_tested=True))
    pd.DataFrame(rows).to_csv(out / 'private_all14_source_panel_feature_feasibility.tsv', sep='\t', index=False)
    pd.DataFrame(detailed).to_csv(out / 'private_all_panel_reference_feature_coverage.tsv', sep='\t', index=False)
    # Preserve source metadata only on the server, no dynamic identifiers in stdout.
    (out / 'private_source_metadata_pointer.json').write_text(json.dumps(dict(panel_axes_sha256=sha(source), source_register_sha256=sha(root / 'state_discovery_evidence_register_20261014_v1.json'), limits='Panel membership is measurability metadata only, not detected expression, cellular identity, cell contact or independent validation. No source/platform inferred from feature count. Encoded LR components do not establish full biological protein complex.'), indent=2) + '\n')
    receipt = dict(status='COMPLETE_EXISTING14_PANEL_FEATURE_FEASIBILITY_NO_SPATIAL_EFFECTS_OR_DOWNLOADS', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), panel_axes_sha256=sha(source), resource_guard=guard, no_dynamic_field_names_or_molecular_endpoints_exported=True, files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in receipt['files'].values()) < 512 * 1024
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
