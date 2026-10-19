"""Two predeclared complementary questions; human arrays remain in RAM."""
import hashlib
import json
import os
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
from immune_RNA_MCP_new_source_math_v1 import (
    design, independent_design, estimate, evaluate, wild_values, independent_wild,
    bootstrap_indices, bootstrap, leave_one, seed, adjust)
from HNSCC_author_state_clinical_kernel_v2 import fit


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_sha(value):
    a = np.ascontiguousarray(value)
    return hashlib.sha256(str(a.shape).encode() + str(a.dtype).encode() + a.tobytes()).hexdigest()


def close(a, b, tol=1e-8):
    assert np.allclose(a, b, rtol=0, atol=tol, equal_nan=True)


def synthetic_fixture():
    rng = np.random.default_rng(2026102007)
    C = rng.normal(size=(60, 3))
    names = np.array(['q0', 'q1', 'q2'])
    for kind in ['clinical', 'coordination']:
        x = np.tile([0., 1.], 30) if kind == 'clinical' else rng.normal(size=60)
        sign = 1 if kind == 'clinical' else -1
        y = sign * 1.5 * x + C @ np.array([.3, -.2, .1]) + rng.normal(size=60) * .2
        m, st = design(x, y, C, names, kind)
        im, ist = independent_design(x, y, C, names, kind)
        assert st == ist == 'ESTIMABLE'
        g, e = estimate(m); ig, iv = independent_wild(im, 999, 2026102007)
        v = wild_values(g, 999, 2026102007)
        close([g['beta'], g['t'], e], [ig['beta'], ig['t'], ig['effect']])
        close(v, iv)
        assert e * sign > 0
        bad, why = design(np.zeros(60), y, C, names, kind)
        assert bad is None and why == 'CONSTANT_EXPOSURE_OR_READOUT'
    return dict(positive_and_negative_and_constant_exposure=True, independent_QR_HC3_full999_replay=True)


def measure(root, config, save):
    project = root.parents[2]
    path = project / 'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json'
    assert sha(path) == '8d7f54a385389d0ecbc825155df727ed065cae2ed1c3f3ed584c1a7896647e71'
    contract = json.loads(path.read_text())
    rec = next(r for r in contract['records'] if r['dataset'] == 'ICB_Bassez')
    assert rec['role'] == 'supplemental_retained_export'
    ready = root / 'Bassez_fixed_function_readiness_20261014_v1'
    oldpath = root / 'Bassez_IFN_validation_preparation_20261016_v1/private_patient_scores_and_covariates.tsv'
    oldterm = json.loads((oldpath.parent / 'private_terminal.json').read_text())
    assert sha(oldpath) == oldterm['patient_table_sha256']
    old = pd.read_csv(oldpath, sep='\t', dtype={'patient_token': str}).set_index('patient_token').sort_index()
    assert old.index.is_unique and len(old) >= 20
    genes = config['sender_components'] + [config['receiver_component']]
    assert len(config['sender_components']) == 2 and len(set(genes)) == 3
    positive = set(np.load(ready / 'private_source_positive_confirmed_columns.npy').tolist())
    # Stable-ID bridge and source evidence are the previously accepted measurement contract.
    parent = next(r for r in contract['records'] if r['dataset'] == 'GSE286827')
    pv = read_elem(zarr.open_group(parent['path'], mode='r')[parent['var_key']])
    aliases = defaultdict(set)
    for name, r in pv.iterrows():
        if re.fullmatch(r'ENSG[0-9]{11}(?:\.[0-9]+)?', str(name)) and r.get('feature_identity_status') == 'stable_Ensembl_ID':
            for alias in {str(r.get('gene_symbol', ''))} | set(str(r.get('source_gene_aliases', '')).split('|')):
                if alias:
                    aliases[alias].add(str(name).split('.')[0])
    before = (Path(rec['path']).stat().st_size, Path(rec['path']).stat().st_mtime_ns)
    table = []
    manual_cells = 0
    total_cells = 0
    with h5py.File(rec['path'], 'r') as h:
        var = read_elem(h['raw/var'])
        lookup = defaultdict(set)
        ids = defaultdict(set)
        for i, (name, row) in enumerate(var.iterrows()):
            if re.fullmatch(r'ENSG[0-9]{11}(?:\.[0-9]+)?', str(name)):
                ids[str(name).split('.')[0]].add(i)
            for field in [str(name)] + [str(row[k]) for k in ['gene_symbol', 'feature_name', 'gene_name_original'] if k in var]:
                if field not in ['', 'nan', 'None']:
                    lookup[field].add(i)
        columns = []
        for gene in genes:
            possible = set(lookup.get(gene, set()))
            if not possible and len(aliases.get(gene, set())) == 1:
                possible = ids.get(next(iter(aliases[gene])), set())
            assert len(possible) == 1
            col = next(iter(possible))
            assert col in positive
            columns.append(col)
        assert len(set(columns)) == 3
        source_rows = np.load(rec['row_positions_ref'], allow_pickle=False)
        md = pd.DataFrame({k: np.asarray(read_elem(h['obs'][k]))[source_rows].astype(str)
                           for k in ['PMID_donor_id', 'pre_post', 'author_cell_type_update']})
        md['row'] = source_rows
        md = md[md.pre_post.eq('Pre') & md.PMID_donor_id.isin(old.index)].copy()
        totals = md.groupby('PMID_donor_id').size()
        raw = ad.io.sparse_dataset(h['raw/X'])
        for number, (token, prior) in enumerate(old.iterrows()):
            row = dict(expanded=int(prior.expanded), treatment=str(prior.treatment), cancer_type=str(prior.cancer_type))
            for prefix, label in [('T', 'T_cell'), ('My', 'Myeloid_cell')]:
                take = np.sort(md.loc[md.PMID_donor_id.eq(token) & md.author_cell_type_update.eq(label), 'row'].to_numpy(int))
                assert len(take) == int(prior[prefix + '_cells']) and len(take) >= 20
                sums = np.zeros(9)
                for lo in range(0, len(take), 1000):
                    selected = take[lo:lo + 1000]
                    block = raw[selected, :].tocsr().astype(float)
                    block.eliminate_zeros()
                    assert np.isfinite(block.data).all() and (block.data > 0).all()
                    assert np.array_equal(block.data, np.rint(block.data))
                    RNA = np.asarray(block.sum(1)).ravel()
                    detected = np.diff(block.indptr)
                    values = block[:, columns].toarray()
                    assert (RNA > 0).all()
                    value = np.log1p(values * 10000 / RNA[:, None])
                    # Independent log-difference normalization of every selected value.
                    alt = np.log(RNA[:, None] + values * 10000) - np.log(RNA[:, None])
                    close(value, alt, 1e-12)
                    sums += np.r_[np.log1p(RNA).sum(), np.log1p(detected).sum(), value.sum(0),
                                  (values > 0).sum(0), ((values[:, 0] > 0) & (values[:, 1] > 0)).sum()]
                    if lo == 0:
                        for offset, physical in enumerate(selected[:2]):
                            a, b = np.asarray(h['raw/X/indptr'][physical:physical + 2], dtype=int)
                            ix = np.asarray(h['raw/X/indices'][a:b], dtype=int)
                            data = np.asarray(h['raw/X/data'][a:b], dtype=float)
                            assert len(set(ix.tolist())) == len(ix)
                            close(RNA[offset], data.sum(), 0)
                            assert detected[offset] == np.count_nonzero(data)
                            for k, column in enumerate(columns):
                                close(values[offset, k], data[ix == column].sum(), 0)
                            manual_cells += 1
                avg = sums / len(take)
                close(avg[:2], [prior[prefix + '_logRNA'], prior[prefix + '_loggenes']], 1e-10)
                row[prefix + '_logRNA'], row[prefix + '_loggenes'] = avg[:2]
                row[prefix + '_logcells'] = float(np.log1p(len(take)))
                row[prefix + '_fraction'] = float(len(take) / totals.loc[token])
                row[prefix + '_cells'] = len(take)
                for k in range(3):
                    row[prefix + '_gene' + str(k)] = avg[2 + k]
                    row[prefix + '_positive' + str(k)] = avg[5 + k]
                row[prefix + '_joint_sender_positive'] = avg[8]
                total_cells += len(take)
            row['sender'] = min(row['My_gene0'], row['My_gene1'])
            row['receiver'] = row['T_gene2']
            row['interface'] = float(np.log1p(row['sender'] * row['receiver']))
            table.append(row)
            save('private_progress.json', dict(stage='NATIVE_FIXED_INTERFACE_MEASUREMENT', completed_author_units=number + 1,
                 total_author_units=len(old), observed_UTC=datetime.now(timezone.utc).isoformat()))
    assert before == (Path(rec['path']).stat().st_size, Path(rec['path']).stat().st_mtime_ns)
    frame = pd.DataFrame(table)
    proof = dict(source_contract_sha256=sha(path), old_accepted_unit_table_sha256=sha(oldpath),
                 unit_order_sha256=array_sha(old.index.to_numpy(str)), source_selection_sha256=array_sha(source_rows),
                 selected_source_stat_unchanged=True, unique_physical_interface_columns=len(columns),
                 all_selected_cells=total_cells, manually_reconstructed_full_CSR_rows=manual_cells,
                 current_native_quality_matches_prior_all_units=True,
                 per_unit_measurement_RAM_array_sha256=array_sha(frame.select_dtypes('number').to_numpy(float)),
                 per_unit_or_cell_RNA_scores_not_persisted=True)
    return frame, proof


def main():
    config = json.load(sys.stdin)
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    out.mkdir(mode=0o700)
    def save(name, value):
        (out / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    deps = {}
    for name in ['immune_RNA_MCP_new_source_math_v1.py', 'NSCLC_immune_RNA_MCP_cohort_math_v1.py',
                 'NSCLC_immune_RNA_MCP_cohort_audit_math_v1.py', 'HNSCC_author_state_clinical_kernel_v2.py']:
        p = Path(os.environ['PANCANCER_DEPENDENCY_ROOT']) / name
        assert sha(p) == config['dependency_sha256'][name]
        deps[name] = sha(p)
    save('private_plan_before_any_current_source_effect.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(),
         source_sha256=sha(__file__), dependency_sha256=deps, private_interface_definition=config,
         primary_family=2, null_scenarios=3, simulations_per_scenario=256, null_wild=999,
         main_wild=19999, bootstrap=2000, required_valid_bootstrap=1900, minimum_units=20,
         minimum_binary_groups=5, minimum_residual_df=10, null_rejection_rate_maximum=.10,
         clinical_interpretation='Author T-clone expansion; not clinical response or gamma identity',
         no_new_external_clinical_transfer_or_RNA_download=True, all_human_arrays_RAM_only=True))
    save('technical_synthetic_fixture.json', synthetic_fixture())
    frame, proof = measure(root, config, save)
    save('private_native_measurement_technical_proof.json', proof)
    numeric = [p + '_' + q for p in ['T', 'My'] for q in ['logRNA', 'loggenes', 'logcells', 'fraction']]
    C = frame[numeric].to_numpy(float)
    names = list(numeric)
    for field in ['treatment', 'cancer_type']:
        for value in sorted(frame[field].unique())[1:]:
            C = np.c_[C, frame[field].eq(value).to_numpy(float)]
            names.append('agent::' + field + '::' + str(value))
    names = np.asarray(names)
    continuous = np.array([not s.startswith('agent::') for s in names])
    clinicalC = C.copy()
    clinicalC[:, continuous] = rankdata(clinicalC[:, continuous], axis=0)
    axis = frame.interface.to_numpy(float)
    axis = (axis - axis.mean()) / axis.std(ddof=0) if np.ptp(axis) else axis
    specifications = [('BIC01', 'clinical', frame.expanded.to_numpy(float), axis, clinicalC),
                      ('BIC02', 'coordination', frame.sender.to_numpy(float), frame.receiver.to_numpy(float), C)]
    models = []
    for slot, kind, x, y, controls in specifications:
        m, st = design(x, y, controls, names, kind)
        im, ist = independent_design(x, y, controls, names, kind)
        assert st == ist
        if min((frame.expanded == 0).sum(), (frame.expanded == 1).sum()) < 5:
            m = im = None
            st = 'FEWER_THAN5_ORIGINAL_EXPANSION_GROUP_UNITS'
        models.append(dict(slot=slot, kind=kind, x=x, y=y, C=controls, m=m, im=im, status=st))
    nulls = []
    for entry in models:
        if entry['m'] is None:
            continue
        m, im = entry['m'], entry['im']
        D = m['D']; n = len(D)
        mu = D[:, :-1] @ np.linspace(-.2, .2, D.shape[1] - 1)
        quality = C[:, 0]
        quality = (quality - quality.mean()) / quality.std(ddof=0) if quality.std(ddof=0) > 0 else np.zeros(n)
        passed = True
        for scenario in ['normal', 't3_unit_variance', 'RNA_quality_heteroskedastic_normal']:
            rng = np.random.default_rng(seed(entry['slot'], 'null_errors_' + scenario))
            rejected = 0; maximum = 0.
            for simulation in range(256):
                err = rng.standard_t(3, n) / np.sqrt(3) if scenario == 't3_unit_variance' else rng.normal(size=n)
                if scenario == 'RNA_quality_heteroskedastic_normal':
                    err *= np.exp(.6 * quality)
                yy = mu + err
                g = fit(yy, D)
                random_seed = seed(entry['slot'], scenario + '|' + str(simulation))
                values = wild_values(g, 999, random_seed)
                ig, iv = independent_wild(im, 999, random_seed, y=yy)
                close([g['beta'], g['t']], [ig['beta'], ig['t']])
                close(values, iv)
                maximum = max(maximum, float(np.max(np.abs(values - iv))))
                count = int((np.abs(values) >= abs(g['t']) - 1e-12).sum())
                assert count == int((np.abs(iv) >= abs(ig['t']) - 1e-12).sum())
                rejected += (count + 1) / 1000 < .05
            rate = rejected / 256
            passed &= rate <= .10
            nulls.append(dict(slot=entry['slot'], scenario=scenario, simulations=256,
                              rejections=rejected, rejection_rate=rate, maximum_independent_error=maximum))
            save('private_progress.json', dict(stage='ALL_PREDECLARED_WORKING_NULLS', completed_scenarios=len(nulls),
                 observed_UTC=datetime.now(timezone.utc).isoformat()))
        entry['null_passed'] = bool(passed)
    save('private_all_working_null_diagnostics.json', nulls)
    results = []
    for e in models:
        m, im = e['m'], e['im']
        record = dict(slot=e['slot'], source='ICB_Bassez', kind=e['kind'], author_units=len(frame),
                      expansion1_units=int(frame.expanded.sum()), expansion0_units=int((frame.expanded == 0).sum()),
                      status=e['status'], effect=None, P=None, CI_low=None, CI_high=None,
                      reviewed_independent_patient_identity=False, same_fine_endpoint_validation=False)
        if m is not None and not e['null_passed']:
            record['status'] = 'PREDECLARED_NULL_DIAGNOSTIC_FAILED_NO_ADJUSTED_INFERENCE'
        elif m is not None:
            g, effect = estimate(m)
            ig, iv = independent_wild(im, 19999, seed(e['slot'], 'main_wild'))
            v = wild_values(g, 19999, seed(e['slot'], 'main_wild'))
            close([g['beta'], g['t'], effect], [ig['beta'], ig['t'], ig['effect']])
            close(v, iv)
            count = int((np.abs(v) >= abs(g['t']) - 1e-12).sum())
            assert count == int((np.abs(iv) >= abs(ig['t']) - 1e-12).sum())
            indices = bootstrap_indices(frame.expanded.to_numpy(int), e['slot'], 2000)
            b, reasons = bootstrap(e['x'], e['y'], e['C'], names, e['kind'], indices, m['nuisance_rank'])
            ib, ir = bootstrap(e['x'], e['y'], e['C'], names, e['kind'], indices, im['nuisance_rank'], independent=True)
            loo, lr = leave_one(e['x'], e['y'], e['C'], names, e['kind'], m['nuisance_rank'])
            iloo, ilr = leave_one(e['x'], e['y'], e['C'], names, e['kind'], im['nuisance_rank'], independent=True)
            close(b, ib); close(loo, iloo)
            assert np.array_equal(reasons, ir) and np.array_equal(lr, ilr)
            valid = np.isfinite(b); good = np.isfinite(loo)
            ci = np.quantile(b[valid], [.025, .975]) if valid.sum() >= 1900 else [None, None]
            record.update(status='ESTIMABLE_EXPLORATORY_SOURCE_COMPLEMENT', effect=effect, P=(count + 1) / 20000,
                          HC3_t=g['t'], CI_low=ci[0], CI_high=ci[1], valid_bootstrap=int(valid.sum()),
                          residual_df=m['df'], nuisance_rank=m['nuisance_rank'], condition=m['condition'],
                          max_leverage=m['max_leverage'], LOO_valid=int(good.sum()),
                          LOO_min=float(loo[good].min()) if good.any() else None,
                          LOO_max=float(loo[good].max()) if good.any() else None,
                          LOO_all_same_sign=bool(good.all() and (loo * effect > 0).all()),
                          all19999_independent_statistics_max_error=float(np.max(np.abs(v - iv))),
                          bootstrap_failure_counts=pd.Series(reasons[~valid]).value_counts().to_dict(),
                          main_random_stream_sha256=array_sha(v), bootstrap_stream_sha256=array_sha(b))
        points = []
        for mode, take in [('unadjusted', np.zeros(len(names), bool)), ('clinical_only', ~continuous),
                           ('RNA_quality_only', np.array(['logRNA' in s or 'loggenes' in s for s in names])),
                           ('capture_only', np.array(['logcells' in s or 'fraction' in s for s in names]))]:
            pm, ps = design(e['x'], e['y'], e['C'][:, take], names[take], e['kind'])
            ipm, ips = independent_design(e['x'], e['y'], e['C'][:, take], names[take], e['kind'])
            assert ps == ips
            point = None
            if pm is not None:
                _, point = estimate(pm); close(point, evaluate(ipm)['effect'])
            points.append(dict(mode=mode, status=ps, effect=point, no_new_P=True))
        record['prespecified_point_sensitivities'] = points
        results.append(record)
    p = np.array([r['P'] if r['P'] is not None else 1 for r in results])
    bh = adjust(p); by = adjust(p, by=True)
    close(bh, multipletests(p, method='fdr_bh')[1], 1e-12)
    close(by, multipletests(p, method='fdr_by')[1], 1e-12)
    for k, r in enumerate(results):
        r['BH2'] = float(bh[k]) if r['P'] is not None else None
        r['BY2'] = float(by[k]) if r['P'] is not None else None
    save('private_complete2_source_association_results.json', results)
    coverage = dict(source='ICB_Bassez', author_units=len(frame), original_counts_source_selection_preserved=True,
                    manually_reconstructed_cells=proof['manually_reconstructed_full_CSR_rows'],
                    total_selected_cells=proof['all_selected_cells'],
                    E_units=int(frame.expanded.sum()), NE_units=int((frame.expanded == 0).sum()),
                    sender_joint_positive_fraction_quantiles=frame.My_joint_sender_positive.quantile([0, .25, .5, .75, 1]).tolist(),
                    receiver_positive_fraction_quantiles=frame.T_positive2.quantile([0, .25, .5, .75, 1]).tolist(),
                    receiver_units_with_zero_mean=int(frame.receiver.eq(0).sum()),
                    sender_units_with_zero_mean=int(frame.sender.eq(0).sum()),
                    physical_patient_identity_and_gamma_specificity_not_validated=True)
    save('private_source_only_capture_summary.json', coverage)
    files = {p.name: dict(bytes=p.stat().st_size, sha256=sha(p)) for p in out.iterdir() if p.is_file()}
    assert sum(v['bytes'] for v in files.values()) < 4 * 1024 ** 2
    terminal = dict(status='COMPLETE_TWO_FIXED_COMPLEMENTARY_QUESTIONS_AND_ALL_NUMERIC_REPLAYS_PENDING_SCIENTIFIC_REVIEW',
                    completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                    files=files, registered_questions=2, human_expression_score_graph_arrays_persisted=False,
                    new_external_clinical_payload_or_RNA_download=False,
                    audit_scope='Same native selection/I-O; alternate normalization/all model QR/random replay; bounded manual CSR rows, not independent full raw archive')
    save('private_terminal.json', terminal)
    print(json.dumps(dict(status=terminal['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
