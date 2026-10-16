"""Freeze three questions, then score fixed public partial IFN readouts.

All native baseline author T/myeloid cells are used. Patient data, source
columns and resulting means remain server-private. No effects are fitted.
"""
import hashlib, json, os, re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import anndata as ad
import h5py, numpy as np, pandas as pd, zarr
from anndata._io.specs import read_elem
from scipy.stats import rankdata


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    project = Path(os.environ['PANCANCER_PROJECT_ROOT'])
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    prior = root / 'Bassez_fixed_function_readiness_20261014_v1'
    terminal = json.loads((root / 'Bassez_fixed_function_readiness_20261014_v2/private_terminal.json').read_text())
    assert terminal['status'] == 'BASSEZ_FIXED_PUBLIC_FUNCTION_READINESS_AUDITED_NO_EFFECT_TEST'
    for record in terminal['prior_artifact_manifest']:
        assert digest(root / record['path']) == record['sha256']
    out = root / 'Bassez_IFN_validation_preparation_20261016_v1'
    assert not out.exists()
    out.mkdir()
    references = ['HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE']
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(),
        producer_sha256=digest(Path(__file__)),
        primary_questions=['myeloid_IFNa_vs_T_clonal_expansion', 'myeloid_IFNg_vs_T_clonal_expansion', 'T_IFNG_vs_myeloid_IFNg'],
        primary_BH_family=3, references=references,
        cell_selection='All native baseline author T_cell/Myeloid_cell cells in frozen eligible patient keys; no outcome-based subset or resampling cells.',
        feature_contract='Unique physical columns with positive-count source evidence; >=80% remaining reference genes and >=10. Partial functional sensitivity only.',
        interface_genes=['IFNG', 'IFNGR1', 'IFNGR2'],
        score='UCell scalar, mean ranks of positive counts, maxrank1500; interface columns excluded from signatures AND background.',
        expansion='Original E/NE labels only; E denotes T-clone expansion, not pathological response.',
        minimum_patients=20, minimum_residual_df=10,
        clinical_nuisance=['Treatment.or.Mode.of.Action', 'Cancer_type_update', 'My_logRNA', 'My_loggenes', 'My_logcells'],
        coupling_nuisance=['Treatment.or.Mode.of.Action', 'Cancer_type_update', 'T_logRNA', 'T_loggenes', 'My_logRNA', 'My_loggenes', 'T_logcells', 'My_logcells', 'My_IFNG'],
        statistical_contract='Rank-continuous readouts; E binary. 20000 reduced-model residual permutations and 2000 E/NE-stratified patient bootstraps, rerank each resample. No automatic deletion of required covariates.',
        categorical_encoding='Treatment and cancer labels one-hot with lexicographic reference; a constant categorical factor contributes no contrasts by definition, all levels reported.',
        Wilcoxon_sensitivity='Two fixed myeloid readouts, rank-biserial effect, 20000 E/NE label permutations within treatment, 2000 E/NE-stratified bootstrap; BH2 separately.',
        missing='Unavailable interface or failed reference measurement remains unavailable; planned primary slot retained internally P=1.',
        shared_genes='Do not interpret IFN-alpha and gamma as independent mechanisms; report their patient correlation descriptively.',
        scope='Author unsplit lineages; not accepted old172 mapping, C1QC identity, causal communication or full old MCP transfer.',
        no_new_expression_download=True)
    (out / 'plan_before_scores_or_effects.json').write_text(json.dumps(plan, indent=2) + '\n')
    units = pd.read_csv(prior / 'private_patient_phase_capture.tsv', sep='\t', keep_default_na=False)
    baseline = units[units.phase.eq('Pre')].copy()
    conflicts = baseline[[c for c in baseline if c.endswith('_conflict')]].apply(lambda s: s.astype(str).str.lower().eq('true')).any(axis=1)
    eligible = baseline[(baseline.author_T_cells >= 20) & (baseline.author_myeloid_cells >= 20) & ~conflicts & baseline.outcome.isin(['E', 'NE'])].set_index('patient_token').sort_index()
    assert len(eligible) == terminal['conflict_free_E_NE_tokens'] and eligible.index.is_unique
    contract_path = project / 'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json'
    contract = json.loads(contract_path.read_text())
    prior_plan = json.loads((prior / 'plan_before_expression_review.json').read_text())
    assert digest(contract_path) == prior_plan['source_contract_sha256']
    rec = next(r for r in contract['records'] if r['dataset'] == 'ICB_Bassez')
    assert rec['role'] == 'supplemental_retained_export'
    with h5py.File(rec['path']) as handle:
        obs = read_elem(handle['obs'])
        var = read_elem(handle['raw/var'])
        positions = np.load(rec['row_positions_ref'])
        md = obs.iloc[positions].copy()
        md['source_position'] = positions
        md['patient_key'] = md.PMID_donor_id.astype(str)
        md = md[md.patient_key.isin(eligible.index) & md.pre_post.astype(str).eq('Pre')]
        known = defaultdict(set)
        byid = defaultdict(set)
        for i, (name, row) in enumerate(var.iterrows()):
            if re.fullmatch(r'ENSG[0-9]{11}(?:\.[0-9]+)?', str(name)):
                byid[str(name).split('.')[0]].add(i)
            for value in [str(name)] + [str(row[c]) for c in ['gene_symbol', 'feature_name', 'gene_name_original'] if c in var]:
                if value not in ['', 'nan', 'None']:
                    known[value].add(i)
        parent_rec = next(r for r in contract['records'] if r['dataset'] == 'GSE286827')
        parent_var = read_elem(zarr.open_group(parent_rec['path'], mode='r')[parent_rec['var_key']])
        bridges = defaultdict(set)
        for name, row in parent_var.iterrows():
            if not re.fullmatch(r'ENSG[0-9]{11}(?:\.[0-9]+)?', str(name)) or row.get('feature_identity_status') != 'stable_Ensembl_ID':
                continue
            for alias in {str(row.get('gene_symbol', ''))} | set(str(row.get('source_gene_aliases', '')).split('|')):
                if alias:
                    bridges[alias].add(str(name).split('.')[0])
        observed = set(np.load(prior / 'private_source_positive_confirmed_columns.npy').tolist())
        interface = {}
        for gene in plan['interface_genes']:
            candidates = set(known.get(gene, set()))
            if not candidates and len(bridges.get(gene, set())) == 1:
                candidates = byid.get(next(iter(bridges[gene])), set())
            interface[gene] = next(iter(candidates)) if len(candidates) == 1 and next(iter(candidates)) in observed else None
        present = [v for v in interface.values() if v is not None]
        assert len(present) == len(set(present)), 'Interface physical-column collision'
        excluded = set(present)
        background = np.array(sorted(observed - excluded), int)
        physical_to_rank = {c: i for i, c in enumerate(background)}
        refcols = pd.read_csv(prior / 'private_reference_columns.tsv', sep='\t')
        signatures = {}; coverage = []
        for reference in references:
            frame = refcols[refcols.reference.eq(reference)]
            remaining = frame[~frame.gene.isin(plan['interface_genes']) & ~frame.physical_column.isin(excluded)]
            valid = remaining[remaining.source_positive_confirmed.astype(str).str.lower().eq('true') & remaining.physical_column.ge(0)]
            columns = sorted(set(valid.physical_column.astype(int)))
            fraction = len(columns) / len(remaining)
            ready = fraction >= .8 and len(columns) >= 10
            coverage.append(dict(reference=reference, original_genes=len(frame), after_interface_exclusion=len(remaining), confirmed_columns=len(columns), fraction=fraction, ready=ready))
            if ready:
                signatures[reference] = np.array([physical_to_rank[c] for c in columns], int)
        (out / 'private_score_measurement_contract.json').write_text(json.dumps(dict(interface=interface, coverage=coverage, background_columns=len(background)), indent=2) + '\n')
        counts = ad.io.sparse_dataset(handle['raw/X'])
        rows = []
        for patient in eligible.index:
            entry = eligible.loc[patient]
            row = dict(patient_token=patient, expanded=int(entry.outcome == 'E'), treatment=entry['Treatment.or.Mode.of.Action'], cancer_type=entry.Cancer_type_update)
            for lineage, prefix, expected in [('T_cell', 'T', entry.author_T_cells), ('Myeloid_cell', 'My', entry.author_myeloid_cells)]:
                take = md[md.patient_key.eq(patient) & md.author_cell_type_update.astype(str).eq(lineage)]
                assert len(take) == expected and len(take) >= 20
                sums = defaultdict(float)
                for start in range(0, len(take), 2000):
                    block = counts[np.sort(take.source_position.to_numpy(int)[start:start+2000]), :].tocsr().astype(float)
                    block.eliminate_zeros()
                    assert np.isfinite(block.data).all() and (block.data > 0).all()
                    total = np.asarray(block.sum(1)).ravel()
                    detected = np.diff(block.indptr)
                    assert (total > 0).all()
                    sums['logRNA'] += np.log1p(total).sum()
                    sums['loggenes'] += np.log1p(detected).sum()
                    for gene, column in interface.items():
                        if column is not None:
                            value = block[:, column].toarray().ravel()
                            sums[gene] += np.log1p(value * 10000 / total).sum()
                            sums[gene + '_positive'] += np.count_nonzero(value)
                    if all(interface.get(g) is not None for g in ['IFNGR1', 'IFNGR2']):
                        sums['joint_receptor_positive'] += np.count_nonzero((block[:, interface['IFNGR1']].toarray().ravel() > 0) & (block[:, interface['IFNGR2']].toarray().ravel() > 0))
                    sub = block[:, background].tocsr()
                    for i in range(sub.shape[0]):
                        lo, hi = sub.indptr[i:i+2]
                        ranks = np.full(len(background), 1500., dtype=float)
                        ranks[sub.indices[lo:hi]] = np.minimum(rankdata(-sub.data[lo:hi], method='average'), 1500)
                        for name, columns in signatures.items():
                            size = len(columns)
                            score = 1 - (ranks[columns].sum() - size*(size+1)/2) / (size*1500-size*(size+1)/2)
                            assert -1e-10 <= score <= 1+1e-10
                            sums[name] += score
                        if i == 0:
                            # Independent scalar first-cell check using count comparisons.
                            dense = sub[i].toarray().ravel(); positive = dense[dense > 0]
                            for name, columns in signatures.items():
                                values = []
                                for c in columns:
                                    v = dense[c]
                                    rank = (np.count_nonzero(positive > v) + (np.count_nonzero(positive == v)+1)/2) if v > 0 else 1500
                                    values.append(min(rank, 1500))
                                assert np.allclose(values, ranks[columns], atol=1e-12)
                row[prefix + '_cells'] = len(take)
                row[prefix + '_logcells'] = float(np.log1p(len(take)))
                for field, value in sums.items():
                    row[prefix + '_' + field] = float(value / len(take))
            rows.append(row)
    scores = pd.DataFrame(rows)
    assert len(scores) == len(eligible) and scores.patient_token.is_unique
    scores.to_csv(out / 'private_patient_scores_and_covariates.tsv', sep='\t', index=False)
    baseline.assign(analysis_eligible=baseline.patient_token.isin(eligible.index)).to_csv(out / 'private_all_baseline_inclusion.tsv', sep='\t', index=False)
    report = dict(status='BASSEZ_FIXED_IFN_SCORES_PREPARED_NO_EFFECTS', finished_UTC=datetime.now(timezone.utc).isoformat(),
        patient_blocks=len(scores), expansion_counts=scores.expanded.value_counts().to_dict(),
        all_native_baseline_lineage_cells=int(scores.T_cells.sum()+scores.My_cells.sum()),
        treatment_levels=sorted(scores.treatment.unique()), cancer_type_levels=sorted(scores.cancer_type.unique()),
        measurement=coverage, interface_source_positive_confirmed={g: c is not None for g, c in interface.items()},
        myeloid_joint_receptor_fraction_range=[float(scores.My_joint_receptor_positive.min()), float(scores.My_joint_receptor_positive.max())] if 'My_joint_receptor_positive' in scores else None,
        patient_table_sha256=digest(out / 'private_patient_scores_and_covariates.tsv'),
        new_expression_download=False, primary_effects_computed=False)
    (out / 'private_terminal.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
