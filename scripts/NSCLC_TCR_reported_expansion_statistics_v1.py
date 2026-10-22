"""Finite local-RAM source associations of reported expansion within RNA types.

No sequence reconstruction, patient identity export or outcome-chosen labels.
These postoperative author material units are not reviewed independent patients.
"""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import subprocess
import numpy as np
import GSE288199_state_pair_math_v1 as math
import NSCLC_observed_T_state_statistics_v1 as established

NAMESPACE = 'NSCLC_reported_TCR_expansion_v1'
KINDS = ['coordination', 'conditional_R', 'conditional_E']


def select(payload, r, e, minimum=20, ids=None, cancer=None):
    assert len(payload['rows']) == 243
    assert len({v['case_index'] for v in payload['rows']}) == 243
    base = [v for v in payload['rows'] if v['clinical_eligible']
            and min(v['TCR_groups'][name]['TCR'] for name in [r, e, 'other_CD4_POOL']) >= minimum
            and (ids is None or v['case_index'] in ids)
            and (cancer is None or v['cancer_type'] == cancer)]
    counts = Counter(established.agent(v['drug']) for v in base)
    rows = sorted([v for v in base if counts[established.agent(v['drug'])] >= 10], key=lambda v: v['case_index'])
    return rows, dict(eligible_before_agent10=len(base), retained_after_agent10=len(rows),
        removed_for_agent_group_below10=len(base)-len(rows), response1=sum(v['response'] == 1 for v in rows),
        response0=sum(v['response'] == 0 for v in rows), reviewed_independent_patients=0)


def fraction(rows, group, ge2=False):
    key = 'expanded_ge2' if ge2 else 'expanded'
    return np.asarray([v['TCR_groups'][group][key]/v['TCR_groups'][group]['TCR'] for v in rows])


def controls(rows, r, e, mode='full', ge2=False):
    if not rows:
        return np.empty((0, 0))
    assert mode in ['unadjusted', 'clinical_only', 'capture', 'full']
    if mode == 'unadjusted':
        return np.empty((len(rows), 0))
    names = sorted({established.agent(v['drug']) for v in rows})
    c = np.c_[[v['clinical'] for v in rows],
              [[established.agent(v['drug']) == name for name in names] for v in rows]]
    if mode == 'clinical_only':
        return c
    quality = []
    for v in rows:
        values = []
        for group in [r, e, 'other_CD4_POOL']:
            g = v['TCR_groups'][group]
            assert 0 < g['TCR'] <= g['RNA']
            values += [np.log1p(g['RNA']), np.log1p(g['TCR']), g['TCR']/g['RNA'],
                       g['called_mean_logRNA'], g['called_mean_loggenes']]
        values += [v['TCR_groups'][r]['RNA']/v['TCR_groups']['Treg_POOL']['RNA'],
                   v['TCR_groups'][e]['RNA']/v['TCR_groups']['CD8_POOL']['RNA']]
        quality.append(values)
    c = np.c_[c, quality]
    return np.c_[c, fraction(rows, 'other_CD4_POOL', ge2)] if mode == 'full' else c


def vectors(kind, rows, r, e, mode='full', ge2=False):
    c = controls(rows, r, e, mode, ge2)
    rv, ev = fraction(rows, r, ge2), fraction(rows, e, ge2)
    response = np.asarray([v['response'] for v in rows], float)
    if kind == 'coordination':
        return rv, ev, c
    assert kind in KINDS
    return (response, rv, np.c_[c, ev]) if kind == 'conditional_R' else (response, ev, np.c_[c, rv])


def complete(kind, rows, r, e, indices):
    x, y, c = vectors(kind, rows, r, e)
    m, f, why = established.association(rows, x, y, c)
    label = '|'.join([NAMESPACE, kind, r, e])
    result = dict(kind=kind, R_group=r, E_group=e, n=len(rows), status=why)
    if f is None:
        audit = established.replay(label, x, y, c, None) if established.qualification(rows) == 'ELIGIBLE_SOURCE_CASE_DESIGN' else None
        return result, audit, None
    result.update({key: f[key] for key in ['rho', 'beta', 'se', 't']})
    signs = np.where(np.random.default_rng(math.seed(label+'|wild')).random((19999, len(rows))) < .5, -1, 1).astype(np.int8)
    wild = math.wild_statistics(m, f, signs)
    boot, reasons = math.bootstrap(x, y, c, indices)
    valid = np.isfinite(boot)
    result.update(raw_P=float((1+np.sum(np.abs(wild) >= abs(f['t'])-1e-12))/20000),
        CI95=np.quantile(boot[valid], [.025, .975]).tolist() if valid.sum() >= 1900 else None,
        bootstrap_valid=int(valid.sum()), bootstrap_status_counts=dict(Counter(reasons)),
        nuisance_rank=m['rank'], residual_df=m['df'])
    audit = established.replay(label, x, y, c, f, signs, wild, indices, boot)
    trials = math.calibration(m, label)
    for trial in trials:
        trial.pop('simulation_pvalues')
    effects, missing = [], Counter()
    for j in range(len(rows)):
        keep = np.arange(len(rows)) != j
        _, ff, why = established.association([v for i, v in enumerate(rows) if i != j], x[keep], y[keep], c[keep])
        if ff is None:
            missing[why] += 1
        else:
            effects.append(ff['rho'])
    result['leave_one_author_material_out'] = dict(estimable=len(effects),
        same_direction=sum(np.sign(v) == np.sign(f['rho']) for v in effects),
        rho_range=[min(effects), max(effects)] if effects else None, missing_status_counts=dict(missing))
    return result, audit, dict(label=label, scenarios=trials)


def families(results, pairs):
    assert len(results) == 3*len(pairs)
    coord = [v for v in results if v['kind'] == 'coordination']
    conditional = [v for v in results if v['kind'] != 'coordination']
    primary, regulatory = [], []
    lookup = {(v['R_group'], v['E_group'], v['kind']): v for v in results}
    for r, e in pairs:
        three = [lookup[(r, e, k)] for k in KINDS]
        primary.append(dict(R_group=r, E_group=e, n=three[0]['n'],
            status='ESTIMABLE' if all(v['status'] == 'ESTIMABLE' for v in three) else 'FIXED_COMPONENT_NOT_ESTIMABLE',
            joint_IUT_P=max(v.get('raw_P', 1.) for v in three),
            all_CI_qualified=all(v.get('CI95') is not None for v in three),
            interpretation='Coordination and both conditional pathology associations; not predictive gain or interaction'))
        regulatory.append(dict(R_group=r, E_group=e, n=three[0]['n'],
            regulatory_IUT_P=max(v.get('raw_P', 1.) for v in three[:2]),
            interpretation='Coordination and regulatory conditional pathology association only'))
    tables = [('coordination40', coord, 'raw_P'), ('conditional80', conditional, 'raw_P'),
              ('primary40', primary, 'joint_IUT_P'), ('regulatory40', regulatory, 'regulatory_IUT_P')]
    for name, table, key in tables:
        values = [v.get(key, 1.) for v in table]
        for v, bh, by in zip(table, math.adjustments(values), math.adjustments(values, by=True)):
            v[name+'_BH'] = float(bh)
            v[name+'_BY'] = float(by)
    return tables


def run(payload, task, plan):
    assert payload['semantics_qualified']
    code = Path(__file__).parent
    for name, digest in plan['code_SHA256'].items():
        assert hashlib.sha256((code/name).read_bytes()).hexdigest() == digest
    os.environ['PANCANCER_R_BLOCK_REFERENCE'] = plan['Rscript_path']
    out = task/'local_statistics_v1'
    out.mkdir(exist_ok=False)
    save = lambda name, v: (out/name).write_text(json.dumps(established.clean(v), indent=2, allow_nan=False)+'\n')
    started = datetime.now(timezone.utc).isoformat()
    pairs = plan['all40_pairs']
    assert len(pairs) == 40 and len(set(map(tuple, pairs))) == 40
    measurement = json.loads((task/'private_complete_local_RNA_bridge_v1.json').read_text())
    reference = {(v['R_group'], v['E_group']): v for v in measurement['all40_clinical_capture_positions']}
    coverage, jobs = [], []
    for r, e in pairs:
        rows, cov = select(payload, r, e)
        assert len(rows) == reference[r, e]['qualified_author_materials']
        assert cov['response1'] == reference[r, e]['response1'] and cov['response0'] == reference[r, e]['response0']
        coverage.append(dict(R_group=r, E_group=e, **cov))
        indices = established.stratified_indices([v['response'] for v in rows], 2000, NAMESPACE+'|'+r+'|'+e) if established.qualification(rows) == 'ELIGIBLE_SOURCE_CASE_DESIGN' else None
        jobs.extend((k, rows, r, e, indices) for k in KINDS)
    save('private_all40_coverage.json', coverage)
    results, audits, calibrations = [], [], []
    with ThreadPoolExecutor(max_workers=2) as executor:
        for result, audit, cal in executor.map(lambda a: complete(*a), jobs):
            results.append(result)
            if audit is not None:
                audits.append(audit)
            if cal is not None:
                calibrations.append(cal)
            save('private_checkpoint.json', dict(completed=len(results), expected=120,
                stage='COMPLETE_FIXED_COMPONENT_AND_FULL_R_REPLAY', observed_UTC=datetime.now(timezone.utc).isoformat()))
    tables = families(results, pairs)
    pv = {name: [v.get(key, 1.) for v in table] for name, table, key in tables}
    rcode = "suppressPackageStartupMessages(library(jsonlite));a<-fromJSON(file('stdin'));cat(toJSON(lapply(a,function(x)list(BH=p.adjust(x,'BH'),BY=p.adjust(x,'BY'))),digits=16))"
    check = subprocess.run([plan['Rscript_path'], '--vanilla', '-e', rcode], input=json.dumps(pv), text=True, capture_output=True, check=True, timeout=60)
    adjusted = json.loads(check.stdout)
    for name, values in pv.items():
        for field, by in [('BH', False), ('BY', True)]:
            np.testing.assert_allclose(adjusted[name][field], math.adjustments(values, by), atol=1e-13, rtol=0)
    for name, table, key in tables:
        save('private_'+name+'.json', table)
    points = []
    for r, e in pairs:
        main, _ = select(payload, r, e)
        ids = {v['case_index'] for v in main}
        settings = [('unadjusted', 'unadjusted', False, {}, None),
                    ('clinical_only', 'clinical_only', False, {}, None),
                    ('capture_no_other_CD4_expansion', 'capture', False, {}, None),
                    ('full', 'full', False, {}, None), ('reported_size_ge2', 'full', True, {}, None),
                    ('TCR_min10', 'full', False, dict(minimum=10), None),
                    ('TCR_min30', 'full', False, dict(minimum=30), None),
                    ('LUAD_same_main', 'full', False, dict(ids=ids), 'LUAD'),
                    ('LUSC_same_main', 'full', False, dict(ids=ids), 'LUSC')]
        for setting, mode, ge2, kw, cancer in settings:
            rows, _ = select(payload, r, e, cancer=cancer, **kw)
            for kind in KINDS:
                x, y, c = vectors(kind, rows, r, e, mode, ge2)
                _, f, why = established.association(rows, x, y, c)
                points.append(dict(setting=setting, kind=kind, R_group=r, E_group=e, n=len(rows), status=why, rho=f['rho'] if f else None))
    assert len(points) == 1080
    save('private_all120_components.json', results)
    save('private_all1080_point_diagnostics.json', points)
    save('private_all_real_random_stream_R_audits.json', audits)
    save('private_all_fixed_working_null_diagnostics.json', calibrations)
    accepted = bool(calibrations) and all(s['passed_finite_alarm'] for c in calibrations for s in c['scenarios'])
    save('private_terminal.json', dict(status='COMPLETE_LOCAL_TCR_SOURCE_FAMILY_WORKING_NULL_ACCEPTED' if accepted else 'COMPLETE_SOURCE_FAMILY_NULL_UNACCEPTED_OR_NO_FITS',
        started_UTC=started, completed_UTC=datetime.now(timezone.utc).isoformat(),
        all120_slots_and40_pairs_retained=True, all_fitted_actual19999wild2000bootstrap_R_replayed=True,
        all_complete_families_R_adjusted=True, individual_and_resampling_arrays_local_RAM_only=True,
        source_postoperative_exploratory_not_independent_validation_or_prediction=True,
        no_sequence_reconstruction_direct_communication_or_original_CM_validation=True))
    print(json.dumps(dict(status='LOCAL_TCR_SOURCE_FAMILY_COMPLETE', working_null_accepted=accepted)), flush=True)
