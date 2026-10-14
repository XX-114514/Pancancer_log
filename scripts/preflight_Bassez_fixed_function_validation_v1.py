"""Outcome-blind feature/capture readiness for fixed public functional readouts.

Counts stay on the server. Positive counts establish a conservative measured
gene lower bound, never prove that all-zero genes were unmeasured. No new model,
clinical statistic, old MCP transfer or fine-lineage mapping is performed.
"""
import hashlib, json, os, re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import anndata as ad
import h5py, numpy as np, pandas as pd, zarr
from anndata._io.specs import read_elem


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    project = Path(os.environ['PANCANCER_PROJECT_ROOT'])
    atlas = project / 'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    contract_path = atlas / 'config/source_contract_v2_private.json'
    contract = json.loads(contract_path.read_text())
    rec = next(r for r in contract['records'] if r['dataset'] == 'ICB_Bassez')
    assert rec['role'] == 'supplemental_retained_export'
    out = root / 'Bassez_fixed_function_readiness_20261014_v1'
    assert not out.exists()
    out.mkdir()
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(),
        question='Can independent fixed public IFN readouts be measured in existing Bassez myeloid/T cells?',
        reference_scope='All50 Hallmarks audited; prospective IFN-alpha and IFN-gamma readouts fixed before expression review.',
        minimum_reference_observed_fraction=.8, minimum_observed_genes=10,
        baseline_minimum_cells_per_author_lineage=20,
        source_measurement='Positive counts across all retained source cells give a conservative lower bound; global column existence is only an upper bound.',
        endpoint='Author E/NE T-cell clonal expansion; not pathological response, no Combined_outcome relabeling.',
        old95_weight_gate_unchanged=True, no_model_fitted=True, no_response_effect_test=True,
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_contract_sha256=hashlib.sha256(contract_path.read_bytes()).hexdigest())
    (out / 'plan_before_expression_review.json').write_text(json.dumps(plan, indent=2) + '\n')
    selected_positions = np.load(rec['row_positions_ref'])
    with h5py.File(rec['path']) as handle:
        obs = read_elem(handle['obs']).iloc[selected_positions]
        var = read_elem(handle['raw/var'])
        counts = ad.io.sparse_dataset(handle['raw/X'])
        observed = np.zeros(len(var), dtype=bool)
        valid_counts = True
        for start in range(0, len(selected_positions), 10000):
            block = counts[selected_positions[start:start+10000], :].tocsr()
            assert np.isfinite(block.data).all() and (block.data >= 0).all()
            assert np.allclose(block.data, np.rint(block.data), atol=1e-8)
            block.eliminate_zeros()
            observed[np.unique(block.indices)] = True
        featuremask_fields = [k for k in ['measured_in_source', 'feature_presence', 'source_gene_mask'] if k in handle.get('varm', {})]
    parent_rec = next(r for r in contract['records'] if r['dataset'] == 'GSE286827')
    parent_var = read_elem(zarr.open_group(parent_rec['path'], mode='r')[parent_rec['var_key']])
    id_pattern = re.compile(r'ENSG[0-9]{11}(?:\.[0-9]+)?')
    id_columns = defaultdict(set)
    lookup = defaultdict(set)
    for i, (name, row) in enumerate(var.iterrows()):
        if id_pattern.fullmatch(str(name)):
            id_columns[str(name).split('.')[0]].add(i)
        for value in [str(name)] + [str(row[c]) for c in ['gene_symbol', 'feature_name', 'gene_name_original'] if c in var]:
            if value not in ['', 'nan', 'None']:
                lookup[value].add(i)
    bridge = defaultdict(set)
    for name, row in parent_var.iterrows():
        if not id_pattern.fullmatch(str(name)) or row.get('feature_identity_status') != 'stable_Ensembl_ID':
            continue
        for alias in {str(row.get('gene_symbol', ''))} | set(str(row.get('source_gene_aliases', '')).split('|')):
            if alias:
                bridge[alias].add(str(name).split('.')[0])
    refroot = root / 'CRC_fine_pair_DIALOGUE_code_20261014_v1/hallmark50_public_reference_v3'
    provenance = json.loads((refroot / 'reference_provenance.json').read_text())
    assert len(provenance['records']) == 50
    coverage = []
    columns_private = []
    for reference in provenance['records']:
        path = refroot / (reference['name'] + '.json')
        assert hashlib.sha256(path.read_bytes()).hexdigest() == reference['sha256']
        genes = json.loads(path.read_text())[reference['name']]['geneSymbols']
        matched = []
        for gene in genes:
            possibilities = set(lookup.get(gene, set()))
            if not possibilities and len(bridge.get(gene, set())) == 1:
                possibilities = id_columns.get(next(iter(bridge[gene])), set())
            matched.append(next(iter(possibilities)) if len(possibilities) == 1 else -1)
        duplicates = {c for c in matched if c >= 0 and matched.count(c) > 1}
        matched = [c if c not in duplicates else -1 for c in matched]
        global_measured = sum(c >= 0 for c in matched)
        lower_bound = sum(c >= 0 and observed[c] for c in matched)
        coverage.append(dict(reference=reference['name'], original_genes=len(genes),
            global_unique_columns=global_measured, source_positive_confirmed_genes=lower_bound,
            source_positive_confirmed_fraction=lower_bound / len(genes),
            fixed_partial_readout_ready=lower_bound / len(genes) >= .8 and lower_bound >= 10,
            full_reference_measurement_known=lower_bound == len(genes)))
        for gene, column in zip(genes, matched):
            columns_private.append(dict(reference=reference['name'], gene=gene, physical_column=column,
                source_positive_confirmed=bool(column >= 0 and observed[column])))
    rows = []
    for (patient, phase), group in obs.groupby(['PMID_donor_id', 'pre_post'], observed=True):
        categories = group.author_cell_type_update.astype(str).value_counts()
        row = dict(patient_token=str(patient), phase=str(phase),
            author_T_cells=int(categories.get('T_cell', 0)),
            author_myeloid_cells=int(categories.get('Myeloid_cell', 0)))
        for field in ['outcome', 'Treatment.or.Mode.of.Action', 'Cancer_type_update', 'Primary_or_met']:
            values = sorted(set(group[field].astype(str)))
            row[field] = '|'.join(values)
            row[field + '_conflict'] = len(values) != 1
        rows.append(row)
    units = pd.DataFrame(rows)
    baseline = units[units.phase.eq('Pre')]
    eligible = baseline[(baseline.author_T_cells >= 20) & (baseline.author_myeloid_cells >= 20)]
    clean = eligible[~eligible[[c for c in eligible if c.endswith('_conflict')]].any(axis=1)]
    labelled = clean[clean.outcome.isin(['E', 'NE'])]
    pd.DataFrame(coverage).to_csv(out / 'private_all50_measurement_coverage.tsv', sep='\t', index=False)
    pd.DataFrame(columns_private).to_csv(out / 'private_reference_columns.tsv', sep='\t', index=False)
    units.to_csv(out / 'private_patient_phase_capture.tsv', sep='\t', index=False)
    np.save(out / 'private_source_positive_confirmed_columns.npy', np.flatnonzero(observed))
    report = dict(status='BASSEZ_FIXED_PUBLIC_FUNCTION_READINESS_AUDITED_NO_EFFECT_TEST',
        finished_UTC=datetime.now(timezone.utc).isoformat(), source_cells=len(obs),
        native_global_features=len(var), source_positive_confirmed_features=int(observed.sum()),
        per_source_measurement_mask_fields=featuremask_fields,
        author_baseline_tokens=len(baseline), both_lineages_atleast20=len(eligible),
        conflict_free_E_NE_tokens=len(labelled),
        treatment_by_expansion=pd.crosstab(labelled['Treatment.or.Mode.of.Action'], labelled.outcome).to_dict(),
        all50_fixed_partial_readouts_ready=sum(r['fixed_partial_readout_ready'] for r in coverage),
        prospective_IFN_readouts=[r for r in coverage if 'INTERFERON' in r['reference']],
        accepted_old_fine_node_mapping=False, clinical_response_validation=False,
        old_MCP_weights_changed=False, new_expression_download=False)
    (out / 'private_terminal.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
