"""Synthetic full-stream validation, independent pairwise tie ranks and replay."""
import hashlib
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.stats import beta
from GSE176021_disjoint_function_math_v1 import observed,permutation_numerators,probability,bootstrap_indices,bootstrap,response_quality_point,adjust
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def ranks(x):
    x = np.asarray(x)
    return 2*(x[:,None]>x[None,:]).sum(axis=1)+(x[:,None]==x[None,:]).sum(axis=1)-len(x)


def independent(a,b):
    aa,bb = ranks(a),ranks(b)
    num = int(sum(int(x)*int(y) for x,y in zip(aa,bb)))
    den = float(np.sqrt(sum(aa*aa)*sum(bb*bb)))
    return num/den if den else np.nan,num,den


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root,2*1024**2)
    out = root/'GSE176021_disjoint_function_math_check_20261019_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    original = np.arange(10,dtype=float)
    fixtures = [('positive',original,original**2),('negative_ties',np.array([0,0,1,1,2,2,3,3,4,4]),np.array([5,4,4,3,3,2,2,1,1,0])),('nonlinear',original,np.array([2,8,3,5,0,1,9,4,6,7])),('zero_variation',original,np.zeros(10)),('two_levels',np.array([0]*5+[1]*5),np.array([1,0,1,0,1,0,1,0,1,0]))]
    indexes = bootstrap_indices()
    idxseed = int.from_bytes(hashlib.sha256(b'2026101908|all10_shared_case_pairs').digest()[:8],'little')
    assert np.array_equal(indexes,np.random.default_rng(idxseed).integers(0,10,size=(2000,10)))
    permutation_positions = 0
    bootstrap_positions = 0
    for label,x,y in fixtures:
        e,n,d = observed(x,y)
        ie,inn,idd = independent(x,y)
        np.testing.assert_allclose([e,n,d],[ie,inn,idd],rtol=0,atol=1e-12,equal_nan=True)
        bs = bootstrap(x,y,indexes)
        other = np.array([independent(x[z],y[z])[0] for z in indexes])
        np.testing.assert_allclose(bs,other,atol=1e-12,rtol=0,equal_nan=True)
        bootstrap_positions += len(bs)
        if d:
            draws = 99999 if label in ['positive','negative_ties'] else 19999
            actual = permutation_numerators(x,y,label,draws)
            s = int.from_bytes(hashlib.sha256(('2026101908|'+label+'|full_case_permutation').encode()).digest()[:8],'little')
            rng = np.random.default_rng(s)
            a,b = ranks(x),ranks(y)
            expect = []
            for lo in range(0,draws,717):
                u = rng.random((min(717,draws-lo),10))
                order = np.argsort(u,axis=1)
                expect.extend(np.sum(b[order]*a[None,:],axis=1).tolist())
            assert np.array_equal(actual,expect)
            count = sum(abs(v)>=abs(inn) for v in expect)
            upper = 1. if count == draws else float(beta.ppf(1-.01/400,count+1,draws-count))
            np.testing.assert_allclose(probability(actual,n),[(count+1)/(draws+1),min(1.,(1+draws*upper)/(draws+1)),count],rtol=0,atol=1e-12)
            permutation_positions += draws
    for quality in [np.c_[original,original],np.c_[original,original**2],np.c_[original%3,(original*3)%7]]:
        x,y = fixtures[2][1:]
        got,status = response_quality_point(x,y,quality)
        D = np.column_stack([np.ones(10)]+[ranks(quality[:,j]) for j in range(2)])
        if np.linalg.matrix_rank(D) < 3:
            assert np.isnan(got) and status.startswith('COLLINEAR')
        else:
            Q,_ = np.linalg.qr(D)
            a,b = ranks(x).astype(float),ranks(y).astype(float)
            a -= Q@(Q.T@a);b -= Q@(Q.T@b)
            if min(np.linalg.norm(a),np.linalg.norm(b)) < 1e-10:
                assert np.isnan(got) and status.startswith('ZERO')
            else:
                np.testing.assert_allclose(got,np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b)),atol=1e-12,rtol=0)
    p = np.r_[np.array([1e-5,.0002,.003,.05,1.]),np.ones(395)]
    for by in [False,True]:
        expected = np.array([min(1.,min(400*value/(p<=value).sum() for value in sorted(set(p)) if value>=target)) for target in p])
        if by:
            # Apply the harmonic factor before clipping, independently of BH clipping.
            expected = np.array([min(1.,min(400*sum(1/k for k in range(1,401))*value/(p<=value).sum() for value in sorted(set(p)) if value>=target)) for target in p])
        np.testing.assert_allclose(adjust(p,by),expected,atol=1e-12,rtol=0)
    receipt = {'status':'ACCEPTED_SYNTHETIC_COMPACT_RANK_PERMUTATION_BOOTSTRAP_MC400_QUALITY_AND_BH_BY_MATH','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'math_module_sha256':sha(Path(__file__).parent/'GSE176021_disjoint_function_math_v1.py'),'resource_guard':guard,'synthetic_permutation_positions':permutation_positions,'synthetic_bootstrap_positions':bootstrap_positions,'actual_biological_tests':False,'complete_future_family_denominator':400}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__ == '__main__':main()
