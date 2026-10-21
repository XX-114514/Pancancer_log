"""Frozen outcome-blind source-state family v2 with sequential bounded RAM replay.

Reuses the accepted reader, molecule-positive recipes and entire source grid.
Only source aggregates, technical proofs and simulation summaries may be saved.
"""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback

import numpy as np

import observed_T_marker_evidence_v1 as evidence
import qualify_GSE288199_external_state_measurement_v2 as qualification
import GSE288199_native_stream_v4 as source
import measure_GSE288199_native_state_capture_v3 as capture
import GSE288199_state_pair_math_v1 as math
import GSE288199_state_pair_R_audit_v1 as independent

R_STATES = ['R_TNFRSF9', 'R_LAG3', 'R_cycling']
E_STATES = ['E_dysfunction', 'E_cycling', 'E_memory']
MAIN_QC = 'human_RNA_QC500_250_MT20'
MAIN_PARENT = 'dual_regulatory_evidence'
MAIN_SORT = 'CD3_CD45_enriched'
CONTROL_NAMES = ['R_mean_log1p_human_UMI', 'R_mean_log1p_human_genes', 'R_log1p_N', 'R_fraction_of_T',
                 'E_mean_log1p_human_UMI', 'E_mean_log1p_human_genes', 'E_log1p_N', 'E_fraction_of_T']


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    return value


def quantile_ci(values):
    values = np.asarray(values)
    good = np.isfinite(values)
    return (np.quantile(values[good], [.025, .975]).tolist() if good.sum() >= 1900 else None), int(good.sum())


def measure_one(contract):
    (counts, totals, detected, mito, proof), full_sha = source.read_full(contract)
    rows = []
    for qi, qc in enumerate(capture.QC):
        eligible = np.ones(len(counts), bool) if qi == 0 else ((totals[:, 0] >= 500) & (detected >= 250) &
                    (mito <= .2*totals[:, 0]) & (totals[:, 0] > 0))
        depth = np.full(len(counts), -1, np.int8)
        if eligible.any():
            cuts = np.unique(np.quantile(totals[eligible, 0], [1/3, 2/3]))
            depth[eligible] = np.searchsorted(cuts, totals[eligible, 0], side='left')
        for parent in capture.PARENTS:
            flags = evidence.vector_evidence(counts[:, :len(evidence.GENES)], eligible, parent)
            np.testing.assert_array_equal(flags, evidence.scalar_evidence(counts[:, :len(evidence.GENES)], eligible, parent))
            np.testing.assert_array_equal(flags[:, 1:5].sum(axis=1), eligible.astype(int))
            other = flags[:, 0] & ~flags[:, 1:4].any(axis=1)
            for di, name in enumerate(capture.DEPTH):
                use = eligible.copy() if di == 0 else eligible & (depth == di-1)
                observed = flags[use].sum(axis=0).astype(np.int64)
                marker = {}
                for pname in ['regulatory_only', 'CD8_only']:
                    take = use & flags[:, evidence.METRICS.index(pname)]
                    marker[pname] = (counts[take] > 0).sum(axis=0).astype(int).tolist()
                quality = {key: float(np.mean(v[use])) if use.any() else None for key, v in
                    [('human_UMI', totals[:, 0]), ('other_GEX_UMI', totals[:, 1]),
                     ('antibody_UMI', totals[:, 2]), ('human_genes', detected)]}
                pq = {}
                for label, parent_mask in [('R', flags[:, 1]), ('E', flags[:, 2]), ('O', other)]:
                    take = use & parent_mask
                    n = int(take.sum())
                    pq[label] = [float(np.mean(np.log1p(totals[take, 0]))) if n else None,
                                 float(np.mean(np.log1p(detected[take]))) if n else None,
                                 float(np.log1p(n)), n / int(observed[0]) if observed[0] else None]
                ox = counts[use & other]
                pos = {g: ox[:, evidence.GENES.index(g)] > 0 for g in
                       ['TNFRSF9', 'LAG3', 'MKI67', 'TOP2A', 'UBE2C']}
                other_states = [int(pos['TNFRSF9'].sum()), int(pos['LAG3'].sum()),
                                int((pos['MKI67'] & (pos['TOP2A'] | pos['UBE2C'])).sum())]
                rows.append(dict(case=contract['case'], phase=contract['phase'], sort=contract['sort'], QC=qc,
                    parent=parent, depth=name, n=int(use.sum()), metrics=observed.tolist(),
                    parent_marker_positive=marker, quality=quality, parent_quality=pq,
                    other_T_cells=int((use & other).sum()), other_T_states=other_states))
    return dict(rows=rows, full_SHA=full_sha, entries=proof['entries'])


def inputs(all_rows, phase, qc=MAIN_QC, parent=MAIN_PARENT, depth='all_capture', other=False):
    selected = [r for r in all_rows if r['phase'] == phase and r['sort'] == MAIN_SORT and
                r['QC'] == qc and r['parent'] == parent and r['depth'] == depth]
    assert len({r['case'] for r in selected}) == len(selected)
    eligible = [r for r in selected if r['metrics'][1] >= 20 and r['metrics'][2] >= 20 and
                (not other or r['other_T_cells'] >= 20)]
    return selected, eligible


def arrays(rows, rstate, estate, other=False):
    n = len(rows)
    j, k = evidence.METRICS.index(rstate), evidence.METRICS.index(estate)
    x = np.asarray([r['metrics'][j]/r['metrics'][1] for r in rows])
    y = np.asarray([r['metrics'][k]/r['metrics'][2] for r in rows])
    c = np.asarray([r['parent_quality']['R']+r['parent_quality']['E']+
                    (r['parent_quality']['O'] if other else []) for r in rows], float).reshape(n, 12 if other else 8)
    if other:
        o = np.asarray([r['other_T_states'][R_STATES.index(rstate)]/r['other_T_cells'] for r in rows])
        return x, y, c, o
    return x, y, c


def block(label, x, y, c, f, indices=None, boots=None, signs=None, wild=None):
    return clean(dict(label=label, x=x, y=y, C=c,
        expected=[f[a] for a in ['beta', 'se', 't', 'rho']] if f is not None else None,
        indices=indices, bootstrap=boots, signs=signs, wild=wild))


def main():
    started = time.monotonic()
    task = Path(os.environ['PANCANCER_STATE_PAIR_TASK']).resolve()
    frozen = json.loads((task/'private_frozen_implementation_v1.json').read_text())
    code = Path(__file__).parent
    for name, digest in frozen['code_SHA256'].items():
        assert sha(code/name) == digest
    out = task/'production'
    assert not out.exists()
    out.mkdir(mode=0o700)
    def save(name, value):
        (out/name).write_text(json.dumps(clean(value), indent=2, allow_nan=False)+'\n')
    prior = Path(os.environ['PANCANCER_STATE_PAIR_PRIOR_CAPTURE'])
    assert sha(prior/'private_terminal.json') == frozen['prior_terminal_SHA256']
    assert sha(prior/'private_all96_native_state_capture_opportunities.json') == frozen['prior_all96_SHA256']
    assert sha(prior/'private_complete_stream_conservation.json') == frozen['prior_conservation_SHA256']
    old_grid = json.loads((prior/'private_all96_native_state_capture_opportunities.json').read_text())
    old_proof = json.loads((prior/'private_complete_stream_conservation.json').read_text())
    save('plan_before_new_state_pair_input.json', frozen)
    raw = qualification.fetch(qualification.BASE+'soft/GSE288199_family.soft.gz', 1024**2)
    assert qualification.digest(raw) == qualification.SOFT_SHA
    samples, series = qualification.soft_rows(raw)
    gex = {}
    for row in samples.values():
        match = re.fullmatch(r'(HN[0-9]+) (Pre|Post) Gene expression data', row['title'])
        if match:
            key = (match[1], match[2].lower())
            assert key not in gex
            typ = row['traits']['cell type']
            assert typ in ['CD45+ CD3+ T cells', 'CD45+ cells']
            gex[key] = MAIN_SORT if typ == 'CD45+ CD3+ T cells' else 'CD45_enriched'
    items = []
    for url in sorted(series['supplementary_file']):
        match = re.search(r'/GSE288199_(HN[0-9]+)_(pre|post)_T_matrix.mtx.gz$', url)
        if match:
            key = match.groups()
            items.append(dict(case=key[0], phase=key[1], sort=gex[key],
                              url=qualification.BASE+'suppl/'+url.rsplit('/', 1)[-1]))
    assert len(items) == len(gex) == 60
    with ThreadPoolExecutor(max_workers=3) as pool:
        contracts = []
        for k, (contract, _) in enumerate(pool.map(lambda a: capture.retry(capture.make_contract, a), items), 1):
            contracts.append(contract)
            save('private_progress.json', dict(stage='PINNED_HEADER_CENSUS', completed=k, total=len(items)))
    assert capture.signature(contracts) == old_proof['source_header_manifest_SHA256']
    assert sum(c['bytes'] for c in contracts) <= 16*1024**3
    assert max(c['bytes'] for c in contracts) <= 4*1024**3
    assert sys.platform == 'darwin'
    vm = subprocess.run(['/usr/bin/vm_stat'], capture_output=True, text=True, timeout=10)
    assert vm.returncode == 0
    page = int(re.search(r'page size of ([0-9]+) bytes', vm.stdout)[1])
    mem = {m[1]: int(m[2]) for line in vm.stdout.splitlines()
           if (m := re.fullmatch(r'([^:]+):\s+([0-9]+)\.', line.strip()))}
    available = sum(mem[k] for k in ['Pages free', 'Pages inactive', 'Pages speculative'])*page
    required = max(c['dims'][1] for c in contracts)*(len(capture.GENES)*8+8*12+30)+512*1024**2
    assert available > required + 3*1024**3 and shutil.disk_usage(task).free > 128*1024**2
    save('private_input_budget_gate.json', dict(accepted=True, header_signature=capture.signature(contracts),
        RAM_estimate=required, available_RAM=available, compressed_source_bytes=sum(c['bytes'] for c in contracts),
        individual_payloads_or_arrays_saved=False))
    rows = []; results = []
    with ThreadPoolExecutor(max_workers=1) as pool:
        futures = [pool.submit(capture.retry, measure_one, c) for c in contracts]
        for k, future in enumerate(as_completed(futures), 1):
            result, _ = future.result(); rows.extend(result['rows']); results.append(result)
            save('private_progress.json', dict(stage='FULL_NATIVE_RNA_REMEASUREMENT', completed=k, total=len(items)))
    assert sum(r['entries'] for r in results) == old_proof['observed_entries']
    digest = qualification.digest(json.dumps(sorted(r['full_SHA'] for r in results)).encode())
    assert digest == old_proof['all_compressed_payload_digests_SHA256']
    assert capture.aggregate(rows) == old_grid
    save('private_complete_source_parity.json', dict(status='PASS_ALL_ORIGINAL_SOURCE_OPPORTUNITIES_AND_FULL_STREAM_DIGESTS',
         original_terminal_SHA256=frozen['prior_terminal_SHA256'], old_source_grid_exactly_reproduced=True,
         full_compressed_payload_multiset_SHA256=digest, all_chunks_R_Matrix_and_cell_rules_replayed=True,
         no_raw_cell_or_coded_unit_arrays_written=True))
    # Construct the entire fixed design family before inspecting pair effects.
    calibration = []; calibrated = set()
    for phase, rstate in product(capture.PHASES, R_STATES):
        _, use = inputs(rows, phase)
        x, _, c = arrays(use, rstate, E_STATES[0])
        m, why = math.design(x, c)
        label = phase+'|'+rstate
        if m is None:
            calibration.append(dict(design=label, status=why, scenarios=[]))
        else:
            scenarios = math.calibration(m, label)
            calibration.append(dict(design=label, status='CALIBRATED_WORKING_MODELS', scenarios=scenarios))
            calibrated.add(label)
        save('private_progress.json', dict(stage='FIXED_DESIGN_WORKING_NULL_CALIBRATION', completed=len(calibration), total=6))
    calibration_pass = all(t['passed_finite_alarm'] for d in calibration for t in d['scenarios'])
    save('private_working_null_calibration.json', dict(passed=calibration_pass, designs=calibration,
        meaning='Finite fixed-design linear-rank-score alarm only, not proof of unit independence or all-null validity'))
    if not calibration_pass:
        save('private_terminal.json', dict(status='CONTROLLED_STOP_WORKING_NULL_CALIBRATION_FAILED_NO_REAL_PAIR_EFFECTS',
             completed_UTC=datetime.now(timezone.utc).isoformat(), source_parity_accepted=True,
             original_CM_labels_or_closed_inference_changed=False))
        return
    main_rows = []; reference_rows = []; sensitivity_rows = []; audit_proofs = []
    def audit_one(b):
        proof = independent.audit([b])
        audit_proofs.extend(proof['blocks'] if isinstance(proof['blocks'], list) else [proof['blocks']])
    for phase, rstate, estate in product(capture.PHASES, R_STATES, E_STATES):
        label = phase+'|'+rstate+'|'+estate
        selected, use = inputs(rows, phase)
        x, y, c = arrays(use, rstate, estate)
        m, f, why = math.association(x, y, c)
        result = dict(source='GSE288199', phase=phase, R_state=rstate, E_state=estate, status=why,
            source_units=len(selected), eligible_source_units=len(use), reviewed_independent_patients=False,
            control_names=CONTROL_NAMES, rho=None, HC3_t=None, P=None, CI=None)
        rng = np.random.default_rng(math.seed('bootstrap|'+phase))
        takes = rng.integers(0, len(use), size=(2000, len(use))) if len(use) else None
        if f is not None:
            assert phase+'|'+rstate in calibrated
            rng = np.random.default_rng(math.seed('wild|'+phase))
            signs = np.where(rng.random((19999, len(use))) < .5, -1., 1.)
            wild = math.wild_statistics(m, f, signs)
            boots, reasons = math.bootstrap(x, y, c, takes)
            ci, valid = quantile_ci(boots)
            loo = []
            for j in range(len(use)):
                keep = np.arange(len(use)) != j
                _, ff, _ = math.association(x[keep], y[keep], c[keep])
                loo.append(ff['rho'] if ff is not None else np.nan)
            finite = np.asarray(loo)[np.isfinite(loo)]
            result.update(rho=f['rho'], HC3_t=f['t'], P=float((1+np.sum(np.abs(wild)>=abs(f['t'])-1e-12))/20000),
                CI=ci, valid_bootstrap=valid, bootstrap_reasons=dict(Counter(reasons)), residual_df=m['df'],
                nuisance_rank=m['rank'], max_leverage=float(m['h'].max()), leave_one_valid=len(finite),
                leave_one_min=float(finite.min()) if len(finite) else None,
                leave_one_max=float(finite.max()) if len(finite) else None,
                leave_one_all_same_direction=bool(len(finite)==len(use) and np.all(np.sign(finite)==np.sign(f['rho']))),
                unadjusted_spearman=float(np.corrcoef(math.ranks(x),math.ranks(y))[0,1]))
            audit_one(block(label, x, y, c, f, takes, boots, signs, wild))
        else:
            audit_one(block(label, x, y, c, None))
        main_rows.append(result)
        # Same records and same complete R/E/other-T nuisance space; descriptive contrasts only.
        _, matched = inputs(rows, phase, other=True)
        xx, yy, cc, oo = arrays(matched, rstate, estate, other=True)
        contrast = dict(source='GSE288199', phase=phase, R_state=rstate, E_state=estate,
                        eligible_source_units=len(matched), status='UNESTIMABLE', R_rho=None, other_T_rho=None,
                        difference=None, pointwise_CI=None, new_P_or_FDR=False)
        mr, fr, wr = math.association(xx, yy, cc)
        mo, fo, wo = math.association(oo, yy, cc)
        contrast.update(R_status=wr, other_T_status=wo)
        if fr is not None and fo is not None:
            rng = np.random.default_rng(math.seed('matched_bootstrap|'+phase))
            ix = rng.integers(0, len(matched), size=(2000, len(matched)))
            br, _ = math.bootstrap(xx, yy, cc, ix); bo, _ = math.bootstrap(oo, yy, cc, ix)
            ci, nv = quantile_ci(br-bo)
            contrast.update(status='DESCRIPTIVE_MATCHED_OTHER_T_COMPARISON', R_rho=fr['rho'], other_T_rho=fo['rho'],
                difference=fr['rho']-fo['rho'], pointwise_CI=ci, valid_paired_bootstraps=nv)
            audit_one(block(label+'|matched_R',xx,yy,cc,fr,ix,br))
            audit_one(block(label+'|matched_other_T',oo,yy,cc,fo,ix,bo))
        else:
            audit_one(block(label+'|matched_R',xx,yy,cc,fr))
            audit_one(block(label+'|matched_other_T',oo,yy,cc,fo))
        reference_rows.append(contrast)
        settings = ([('alternate_parent', MAIN_QC, 'costate_regulatory_sensitivity', 'all_capture'),
                     ('all_released_columns', capture.QC[0], MAIN_PARENT, 'all_capture')] +
                    [(depth, MAIN_QC, MAIN_PARENT, depth) for depth in capture.DEPTH[1:]])
        for setting, qc, parent, depth in settings:
            _, scope = inputs(rows, phase, qc=qc, parent=parent, depth=depth)
            sx, sy, sc = arrays(scope, rstate, estate)
            _, sf, sw = math.association(sx, sy, sc)
            chosen = {v['case'] for v in scope}
            original_same = [v for v in use if v['case'] in chosen]
            bx, by, bc = arrays(original_same, rstate, estate)
            _, bf, bw = math.association(bx, by, bc)
            shared = {v['case'] for v in original_same}
            modified_same = [v for v in scope if v['case'] in shared]
            assert {v['case'] for v in original_same} == {v['case'] for v in modified_same}
            sxm, sym, scm = arrays(modified_same, rstate, estate)
            _, mf, mw = math.association(sxm, sym, scm)
            sensitivity_rows.append(dict(source='GSE288199', phase=phase, R_state=rstate, E_state=estate,
                setting=setting, source_units=len(scope), status=sw, rho=sf['rho'] if sf else None,
                shared_source_units=len(original_same), modified_measurement_same_units_status=mw,
                modified_measurement_same_units_rho=mf['rho'] if mf else None,
                original_measurement_same_units_status=bw, original_measurement_same_units_rho=bf['rho'] if bf else None,
                new_P_or_CI=False))
        save('private_progress.json', dict(stage='COMPLETE_FIXED18_AND_MATCHED_REFERENCES_AND_SENSITIVITIES', completed=len(main_rows), total=18))
    assert len(main_rows) == len(reference_rows) == 18 and len(sensitivity_rows) == 90
    p = [r['P'] if r['P'] is not None else 1 for r in main_rows]
    bh, by = math.adjustments(p), math.adjustments(p, by=True)
    for r, q, qq in zip(main_rows, bh, by):
        r['BH18'] = float(q) if r['P'] is not None else None
        r['BY18'] = float(qq) if r['P'] is not None else None
    assert len(audit_proofs) == 54
    save('private_progress.json', dict(stage='INDEPENDENT_R_ALL_DESIGNS_AND_COMPLETE_RANDOM_STREAM_REPLAY', submitted_blocks=len(audit_proofs)))
    accepted = dict(status='PASS_ALL54_SUBMITTED_SOURCE_DESIGNS_AND_FULL_RANDOM_STREAM_QR_REPLAY',blocks=audit_proofs,
                    RAM_mode='one statistical block per R invocation; no combined individual payload retained')
    save('private_complete_independent_R_QR_audit.json', accepted)
    save('private_complete18_state_pair_results.json', main_rows)
    save('private_complete18_matched_other_T_contrasts.json', reference_rows)
    save('private_complete90_measurement_sensitivities.json', sensitivity_rows)
    save('private_terminal.json', dict(status='COMPLETE_FIXED18_SOURCE_STATE_ASSOCIATIONS_AND_ALL_REFERENCES_AND_SENSITIVITIES_R_AUDITED',
        completed_UTC=datetime.now(timezone.utc).isoformat(), elapsed_seconds=time.monotonic()-started,
        independent_patients_treatment_or_outcomes_not_qualified=True,
        not_direct_communication_or_original_CM_projection_or_clinical_validation=True,
        original_CM_labels_and_closed_families_unchanged=True,
        all_individual_raw_score_quality_and_resampling_arrays_RAM_only=True,
        files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()}))


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        out = Path(os.environ['PANCANCER_STATE_PAIR_TASK'])/'production'
        if out.is_dir():
            (out/'private_failed_terminal.json').write_text(json.dumps(dict(status='FAILED_NOT_SCIENTIFICALLY_ACCEPTED',
                observed_UTC=datetime.now(timezone.utc).isoformat(), error_class=type(error).__name__,
                frames=[dict(file=Path(f.filename).name,line=f.lineno,function=f.name) for f in traceback.extract_tb(error.__traceback__)],
                no_individual_data_in_error_output=True),indent=2)+'\n')
        raise SystemExit(1)
