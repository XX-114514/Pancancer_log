"""Independently reconstruct every supported official and heldout gene score."""
import gzip
import itertools
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def read_table(path):
    return pd.read_csv(path, sep='\t', dtype={'cell': str, 'sample': str})


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 8 * 1024**2)
    producer = root / 'immune_RNA_MCP_frozen_projection_20261018_v1'
    terminal = producer / 'private_terminal.json'
    record = json.loads(terminal.read_text())
    assert record['exit_code'] == 0
    for name, digest in record['files'].items():
        assert sha(producer / name) == digest
    contract_path = root / 'immune_RNA_MCP_postfit_contract_20261018_v1/plan_before_learned_outputs_or_projection.json'
    assert sha(contract_path) == '1172d887cffb7c3b0c28a21583154fbe3718a1f13d79edb85500343ead5e2b09'
    contract = json.loads(contract_path.read_text())
    inputs = root / 'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs'
    model = root / 'immune_RNA_MCP_official_models_20261018_v1'
    plan = json.loads((producer / 'plan_before_official_output_read.json').read_text())
    assert sha(model / 'private_terminal.json') == plan['original_model_terminal_sha256']
    for name, digest in plan['original_model_files_sha256'].items():
        assert sha(model / name) == digest
    output = root / 'immune_RNA_MCP_frozen_projection_audit_20261018_v1'
    assert not output.exists()
    output.mkdir(mode=0o700)
    (output / 'plan_before_complete_projection_replay.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        producer_terminal_sha256=sha(terminal), contract_sha256=sha(contract_path), budget=budget,
        scope='All33 planned slots. All supported raw gene values, discovery means/SD, signed coefficients, independently solved nuisance regression, every official saved/reconstructed/heldout cell score and every case mean. Model RDS/source SHA preservation; no new expression/model learning or scientific P.'), indent=2) + '\n')
    p = producer / 'private_R_projection'
    slots = read_table(p / 'private_all_planned_projection_slots.tsv')
    expected = []
    for fit in contract['planned_fits']:
        fit_id = 'K%d_seed%d_%s' % (fit['k'], fit['seed'], fit['source_mode'])
        expected.extend((fit_id, v, 'MCP%d' % k) for v, k in itertools.product(
            contract['view_to_official_name'], range(1, fit['k'] + 1)))
    keys = ['fit_id', 'view', 'program']
    assert len(slots) == 33 and not slots.duplicated(keys).any()
    assert set(map(tuple, slots[keys].to_numpy())) == set(expected)
    assert all(slots.official_view == slots.view.map(contract['view_to_official_name']))
    supported = slots.status.eq('OFFICIAL_PARITY_AND_FROZEN_SOURCE_PROJECTION_PASS')
    checks = []; maximum_error = 0.; total_cells = 0; total_cases = 0
    if supported.any():
        weights = read_table(p / 'private_frozen_gene_weights.tsv')
        cells = read_table(p / 'private_all_cell_score_reconstruction.tsv.gz')
        cases = read_table(p / 'private_all_case_mean_scores.tsv')
        assert not weights.duplicated(keys + ['gene']).any()
        assert not cells.duplicated(keys + ['partition', 'cell']).any()
        assert not cases.duplicated(keys + ['partition', 'sample']).any()
        supported_keys = set(map(tuple, slots.loc[supported, keys].to_numpy()))
        for data in [weights, cells, cases]:
            assert set(map(tuple, data[keys].to_numpy())) == supported_keys
        for view in contract['view_to_official_name']:
            matrices = {}; metadata = {}; genes = None
            for part in ['discovery', 'TNBC_holdout']:
                prefix = inputs / (view + '_' + part)
                md = read_table(Path(str(prefix) + '_metadata.tsv'))
                g = Path(str(prefix) + '_genes.tsv').read_text().splitlines()
                assert genes is None or genes == g
                genes = g
                with gzip.open(str(prefix) + '_tpm.mtx.gz', 'rb') as handle:
                    matrix = mmread(handle).tocsr()
                assert matrix.shape == (len(genes), len(md)) and not md.cell.duplicated().any()
                assert not md['sample'].isna().any() and np.isfinite(matrix.data).all()
                matrices[part] = matrix; metadata[part] = md
            train = metadata['discovery']; held = metadata['TNBC_holdout']
            assert train['sample'].nunique() == 38 and held['sample'].nunique() == 28
            assert train.partition.eq('reviewed_case_discovery').all()
            assert held.partition.eq('TNBC_source_holdout_sample_identity_pending').all()
            assert held.dataset.eq('GSE305078').all() and held.sourceCC.isna().all()
            assert not set(train.cell) & set(held.cell) and not set(train['sample']) & set(held['sample'])
            for row in slots[supported & slots.view.eq(view)].itertuples(index=False):
                keep = weights.fit_id.eq(row.fit_id) & weights.view.eq(view) & weights.program.eq(row.program)
                w = weights[keep].copy()
                assert len(w) == row.active_genes and len(w) > 0
                ix = pd.Index(genes).get_indexer(w.gene)
                assert (ix >= 0).all()
                raw_train = matrices['discovery'][ix].toarray()
                mu = raw_train.mean(axis=1); sd = raw_train.std(axis=1, ddof=1)
                assert np.all(sd > 0) and np.isfinite(w.weight).all()
                assert np.allclose(mu, w.discovery_mean, rtol=1e-10, atol=1e-10)
                assert np.allclose(sd, w.discovery_sd, rtol=1e-10, atol=1e-10)
                gene_weight = w.weight.to_numpy() / sd
                raw = np.asarray(gene_weight @ matrices['discovery'][ix]).ravel() - gene_weight @ mu
                raw_held = np.asarray(gene_weight @ matrices['TNBC_holdout'][ix]).ravel() - gene_weight @ mu
                adjusted = row.fit_id.endswith('source_adjusted')
                controls = ['cellQ', 'geneQ'] + (['sourceCC'] if adjusted else []) + ['fraction', 'logN']
                D = np.column_stack([np.ones(len(train)), train[controls].to_numpy(float)])
                assert np.linalg.matrix_rank(D) == D.shape[1]
                coefficients = np.linalg.lstsq(D, raw, rcond=None)[0]
                predicted = {'discovery': raw - D @ coefficients}
                common = [c for c in controls if c != 'sourceCC']
                ci = [0] + [controls.index(c) + 1 for c in common]
                Dh = np.column_stack([np.ones(len(held)), held[common].to_numpy(float)])
                predicted['TNBC_holdout'] = raw_held - Dh @ coefficients[ci]
                for part, values in predicted.items():
                    md = metadata[part]
                    s = cells[cells.fit_id.eq(row.fit_id) & cells.view.eq(view) &
                        cells.program.eq(row.program) & cells.partition.eq(part)].set_index('cell')
                    assert len(s) == len(md) and set(s.index) == set(md.cell)
                    s = s.loc[md.cell]
                    assert s['sample'].tolist() == md['sample'].tolist() and s.dataset.tolist() == md.dataset.tolist()
                    err = np.max(np.abs(values - s.score.to_numpy()))
                    if part == 'discovery':
                        err = max(err, np.max(np.abs(values - s.official_saved_score.to_numpy())))
                        assert s.source_offset_unknown.eq(False).all()
                    else:
                        assert s.official_saved_score.isna().all()
                        assert s.source_offset_unknown.eq(adjusted).all()
                    assert np.isfinite(values).all() and np.isfinite(err) and err < 1e-6
                    maximum_error = max(maximum_error, float(err)); total_cells += len(s)
                    expected_means = pd.DataFrame(dict(sample=md['sample'].to_numpy(),
                        dataset=md.dataset.to_numpy(), score=values)).groupby(['sample', 'dataset']).score.mean()
                    c = cases[cases.fit_id.eq(row.fit_id) & cases.view.eq(view) &
                        cases.program.eq(row.program) & cases.partition.eq(part)].set_index(['sample', 'dataset'])
                    assert len(c) == len(expected_means) and set(c.index) == set(expected_means.index)
                    assert np.allclose(c.loc[expected_means.index, 'score'], expected_means, rtol=1e-9, atol=1e-9)
                    total_cases += len(c)
                checks.append(dict(fit_id=row.fit_id, view=view, program=row.program,
                    positive_gene_count=len(w), all_official_and_projected_cell_scores_replayed=True))
        assert total_cells == len(cells) and total_cases == len(cases)
    result = dict(status='ALL_PLANNED_SLOTS_AND_ALL_SUPPORTED_PROJECTION_VALUES_ACCEPTED',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        producer_terminal_sha256=sha(terminal), planned_slots=33, supported_slots=int(supported.sum()),
        unavailable_slots=int((~supported).sum()), reconstructed_cell_score_records=total_cells,
        independently_reaggregated_case_score_records=total_cases, max_score_error=maximum_error,
        all_model_files_unchanged=True, seed_stability_and_source_association_inference='PENDING',
        independent_biological_or_clinical_validation=False)
    (output / 'private_complete_projection_acceptance.json').write_text(json.dumps(
        dict(result, all_supported_slot_checks=checks), indent=2) + '\n')
    (output / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'])


if __name__ == '__main__':
    main()
