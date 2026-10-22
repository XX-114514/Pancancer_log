"""Synthetic numerical, denominator and coverage checks; no human data."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from datetime import datetime, timezone
import argparse
import copy
import hashlib
import json
import os
import numpy as np
import NSCLC_TCR_reported_expansion_statistics_v1 as family
import NSCLC_observed_T_state_statistics_v1 as established


def fixture():
    rng = np.random.default_rng(486531)
    rs = ['Treg_POOL', 'synthetic_R1', 'synthetic_R2', 'synthetic_R3']
    es = ['synthetic_E'+str(j) for j in range(10)]
    rows = []
    for i in range(243):
        response = float(i % 2)
        latent = rng.normal()+.7*response
        groups = {}
        for name in rs+es+['other_CD4_POOL', 'CD8_POOL']:
            n = int(rng.integers(140, 260))
            rn = n+int(rng.integers(20, 140))
            shift = .35*latent if name in ['Treg_POOL', 'synthetic_R1'] else -.3*latent if name == 'synthetic_E0' else rng.normal(0, .2)
            p = 1/(1+np.exp(1.5-shift))
            k = int(rng.binomial(n, p))
            groups[name] = dict(RNA=rn, TCR=n, expanded=k, expanded_ge2=k+int(rng.binomial(n-k, .1)),
                called_mean_logRNA=float(rng.uniform(6, 8)), called_mean_loggenes=float(rng.uniform(4, 6)))
        groups['Treg_POOL']['RNA'] = sum(groups[n]['RNA'] for n in rs[1:])
        groups['CD8_POOL']['RNA'] = sum(groups[n]['RNA'] for n in es)
        groups['synthetic_R3']['TCR'] = 5
        groups['synthetic_R3']['expanded'] = 1
        groups['synthetic_R3']['expanded_ge2'] = 2
        rows.append(dict(case_index=i, clinical_eligible=i<160, response=response,
            clinical=rng.normal(size=8).tolist(), drug='synthetic_agent(C3)', cancer_type='synthetic_cancer', TCR_groups=groups))
    return dict(rows=rows, semantics_qualified=True), [(r, e) for r in rs for e in es]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--Rscript', required=True)
    a = p.parse_args()
    assert not a.output.exists()
    os.environ['PANCANCER_R_BLOCK_REFERENCE'] = a.Rscript
    payload, pairs = fixture()
    assembled = []
    for r, e in pairs:
        rows, cov = family.select(payload, r, e)
        assert len(rows) == (0 if r == 'synthetic_R3' else 160)
        for kind in family.KINDS:
            x, y, c = family.vectors(kind, rows, r, e)
            _, f, why = established.association(rows, x, y, c)
            assert (f is None) == (r == 'synthetic_R3')
            assembled.append(dict(kind=kind, R_group=r, E_group=e, n=len(rows), status=why, raw_P=.5))
    tables = family.families(assembled, pairs)
    assert [len(v) for _, v, _ in tables] == [40, 80, 40, 40]
    missing = copy.deepcopy(payload)
    for row in missing['rows']:
        if row['case_index'] >= 30:
            row['TCR_groups']['synthetic_R1']['TCR'] = 5
    rr, _ = family.select(missing, 'synthetic_R1', 'synthetic_E0')
    assert established.qualification(rr) == 'FEWER_THAN50_AUTHOR_CASE_UNITS'
    rows, _ = family.select(payload, 'synthetic_R1', 'synthetic_E0')
    indices = established.stratified_indices([v['response'] for v in rows], 2000, family.NAMESPACE+'|pure_fixture')
    with ThreadPoolExecutor(max_workers=2) as executor:
        full = list(executor.map(lambda kind: family.complete(kind, rows, 'synthetic_R1', 'synthetic_E0', indices), family.KINDS))
    assert all(v['status'] == 'ESTIMABLE' and audit['full_wild_statistics_checked'] == 19999 and audit['full_bootstrap_statistics_checked'] == 2000 for v, audit, cal in full)
    assert all(s['passed_finite_alarm'] for v, audit, cal in full for s in cal['scenarios'])
    assert full[0][0]['rho'] < 0
    assert full[1][0]['rho'] > 0
    a.output.write_text(json.dumps(dict(status='PASS_ALL120_ASSEMBLY_THREE_COMPLETE_RANDOM_R_BRANCHES_AND_UNESTIMABLE',
        completed_UTC=datetime.now(timezone.utc).isoformat(), no_human_data=True,
        actual_real_fitted_designs_all_full_R_required=True,
        source_code_SHA256=hashlib.sha256(Path(family.__file__).read_bytes()).hexdigest(),
        audits=[audit for v, audit, cal in full], working_nulls=[cal for v, audit, cal in full]), indent=2)+'\n')
    print('PASS_SYNTHETIC_LOCAL_TCR_FAMILY', flush=True)


if __name__ == '__main__':
    main()
