"""Synthetic numerical and source-unit scope checks; no human data."""
import json
import numpy as np
import GSE288199_state_pair_math_v1 as kernel
import GSE288199_state_pair_R_audit_v1 as independent
import run_GSE288199_state_pairs_v1 as worker


def main():
    rng = np.random.default_rng(2026102201)
    np.testing.assert_allclose(kernel.ranks([3,1,1,2]), [4,1.5,1.5,3])
    c = rng.normal(size=(80, 8))
    x = rng.normal(size=80)
    blocks = []
    for label, sign in [('positive', 1), ('negative', -1)]:
        y = sign*x + .3*c[:,0] + rng.normal(size=80)*.2
        m, f, why = kernel.association(x, y, c)
        assert why == 'ESTIMABLE' and np.sign(f['rho']) == sign
        signs = np.where(rng.random((511,80)) < .5, -1., 1.)
        wild = kernel.wild_statistics(m, f, signs)
        takes = rng.integers(0,80,size=(128,80))
        boots, reasons = kernel.bootstrap(x,y,c,takes)
        assert np.isfinite(boots).all()
        blocks.append(worker.block(label,x,y,c,f,takes,boots,signs,wild))
    for label, xx, yy, cc in [('constant_exposure', np.zeros(80), x, c),
        ('constant_readout', x, np.zeros(80), c),
        ('too_few_units',x[:19],x[:19],c[:19]),
        ('controlled_exposure',c[:,0],x,c)]:
        _, f, _ = kernel.association(xx,yy,cc)
        assert f is None
        blocks.append(worker.block(label,xx,yy,cc,None))
    # Exact nuisance collinearity preserves the full span, rather than dropping a scientific control.
    cc = np.c_[c, c[:,0]]
    _, f1, _ = kernel.association(x,y,c)
    _, f2, why = kernel.association(x,y,cc)
    assert why == 'ESTIMABLE'
    np.testing.assert_allclose(f1['rho'],f2['rho'],atol=1e-12)
    blocks.append(worker.block('equivalent_nuisance_space',x,y,cc,f2))
    np.testing.assert_allclose(kernel.adjustments([.01,.04,1]),[.03,.06,1])
    assert worker.quantile_ci(np.r_[np.arange(1899),np.full(101,np.nan)])[0] is None
    assert worker.quantile_ci(np.arange(2000))[1] == 2000
    # Reconstruct the original capture summary from a synthetic complete read without network.
    counts = rng.integers(0,3,size=(120,len(worker.capture.GENES)))
    totals = np.c_[rng.integers(500,3000,size=120),np.zeros((120,2))]
    detected = rng.integers(250,600,size=120); mito = np.zeros(120)
    contract = dict(case='synthetic_only',phase='pre',sort=worker.MAIN_SORT)
    original = worker.source.read_full
    worker.source.read_full = lambda _: ((counts,totals,detected,mito,dict(entries=1)), 'synthetic')
    try:
        out = worker.measure_one(contract)
    finally:
        worker.source.read_full = original
    assert len(out['rows']) == 16
    for row in out['rows']:
        assert sum(row['metrics'][1:5]) == row['n']
        assert row['other_T_cells'] <= row['metrics'][0]
        assert all(v <= row['other_T_cells'] for v in row['other_T_states'])
    proof = independent.audit(blocks)
    assert proof['status'].startswith('PASS_ALL_')
    print(json.dumps(dict(status='PASS_SYNTHETIC_RANK_SPACE_HC3_BOOTSTRAP_R_QR_AND_MEASUREMENT_SCOPE',
        checks=['ties','positive','negative','constant_exposure','constant_readout','coverage','confounded_exposure',
                'equivalent_nuisance_space','BH','CI_valid_draw_gate','synthetic_capture_conservation',
                'independent_R_all_fixture_statistics_and_resampling']),indent=2))


if __name__ == '__main__':
    main()
