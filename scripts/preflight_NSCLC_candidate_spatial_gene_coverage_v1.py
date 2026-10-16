"""All original candidate literal/literature groups on all existing spatial axes.

Only availability. No expression, proximity or programme transfer is tested.
"""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    state = Path(json.loads((root / 'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    axes_root = state / 'state_spatial_proxy_calibration_inputs_20261016_v1'
    inventory = json.loads((axes_root / 'inventory.json').read_text())
    axes_file = axes_root / 'fixed_all14_panel_gene_axes.json'
    assert sha(axes_file) == inventory['files'][axes_file.name]
    panels = json.loads(axes_file.read_text())
    assert len(panels) == 14
    native = root / 'NSCLC_candidate_lineage_cofactor_RNA_20261016_v1'
    assert json.loads((native / 'private_independent_audit.json').read_text())['status'].startswith('PASS_ALL_NATIVE_UNIT_MARGINALS_')
    evidence_root = root / 'NSCLC_automatic_ligand_response_review_20261016_v1'
    evidence_file = evidence_root / 'private_all_primary_BH_supported_rows.tsv'
    evidence_inventory = json.loads((evidence_root / 'private_terminal.json').read_text())
    assert sha(evidence_file) == evidence_inventory['files'][evidence_file.name]
    table = pd.read_csv(evidence_file, sep='\t')
    rule_file = root / 'NSCLC_automatic_ligand_molecular_review_20261016_v1/private_public_literature_component_rules.json'
    rules = json.loads(rule_file.read_text())['ligand_rules']
    rows = []
    for parent in table.itertuples():
        rule = rules[parent.ligand]
        for panel in panels:
            genes = set(panel['unambiguous_symbols'])
            assert sha_text(genes) == panel['gene_axis_sha256']
            records = [('ORIGINAL_LITERAL_REFERENCE', [parent.ligand], r.split('_')) for r in parent.qualified_receptors.split('|')]
            records += [('LITERATURE_DIAGNOSTIC_NOT_REPLACEMENT_HYPOTHESIS', s, r) for s in rule['sender_sets'] for r in rule['receiver_sets']]
            for basis, sender, receiver in records:
                rows.append(dict(record=parent.record, ligand=parent.ligand, dataset_id=panel['dataset_id'], sample_id=panel['sample_id'], gene_axis_sha256=panel['gene_axis_sha256'], input_sha256=panel['input_sha256'], basis=basis, sender_group='|'.join(sorted(sender)), receiver_group='|'.join(sorted(receiver)), sender_all_genes_measured=set(sender) <= genes, receiver_all_genes_measured=set(receiver) <= genes, both_all_genes_measured=(set(sender) | set(receiver)) <= genes, exact_author_subtypes_or_patients_validated=False, spatial_expression_or_proximity_tested=False))
    data = pd.DataFrame(rows).drop_duplicates(['record', 'sample_id', 'basis', 'sender_group', 'receiver_group'])
    out = root / 'NSCLC_candidate_spatial_gene_coverage_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    data.to_csv(out / 'private_all_candidate_spatial_group_coverage.tsv', sep='\t', index=False)
    summary = data.groupby(['ligand', 'basis', 'sender_group', 'receiver_group'], sort=True).agg(assays=('sample_id', 'nunique'), both_measured_assays=('both_all_genes_measured', 'sum')).reset_index()
    summary.to_csv(out / 'private_all_group_coverage_summary.tsv', sep='\t', index=False)
    result = dict(status='COMPLETE_ALL_ORIGINAL_PARENTS_ALL14_EXISTING_SPATIAL_AXES_LITERAL_AND_LITERATURE_GROUP_AVAILABILITY_ONLY', finished_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), axes_file_sha256=sha(axes_file), original_evidence_sha256=sha(evidence_file), molecular_rules_sha256=sha(rule_file), native_audit_sha256=sha(native / 'private_independent_audit.json'), original_parents=len(table), assays=len(panels), unique_assay_gene_axes=len({p['gene_axis_sha256'] for p in panels}), coverage_rows=len(data), literal_fully_measured_rows=int(data.loc[data.basis.eq('ORIGINAL_LITERAL_REFERENCE'), 'both_all_genes_measured'].sum()), files={p.name: sha(p) for p in out.iterdir() if p.is_file()}, new_download_or_expression_or_spatial_tests=False, missing_genes_not_negative=True, exact_cell_identity_or_independent_validation=False)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


def sha_text(genes):
    return hashlib.sha256(('\n'.join(sorted(genes)) + '\n').encode()).hexdigest()


if __name__ == '__main__':
    main()
