"""Verify unchanged v2 numerics and single-block R acceptance; synthetic only."""
import numpy as np
import test_GSE288199_state_pair_math_v1 as fixtures
import run_GSE288199_state_pairs_v2 as worker
import GSE288199_state_pair_math_v1 as kernel
import GSE288199_state_pair_R_audit_v1 as independent


def main():
    fixtures.worker = worker
    fixtures.main()
    rng = np.random.default_rng(2026102202)
    x = rng.normal(size=60); c = rng.normal(size=(60, 8)); y = x + rng.normal(size=60)
    _, f, why = kernel.association(x, y, c)
    assert why == 'ESTIMABLE'
    answer = independent.audit([worker.block('single_block_synthetic',x,y,c,f)])
    rows = answer['blocks'] if isinstance(answer['blocks'], list) else [answer['blocks']]
    assert len(rows) == 1 and rows[0]['estimable']
    print('PASS_SINGLE_BLOCK_R_RUNTIME_AND_UNCHANGED_NUMERICAL_FIXTURES')


if __name__ == '__main__':
    main()
