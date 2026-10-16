"""All input rows, nuisance spans and80k bootstrap design draws independently."""
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.linalg import qr
import statsmodels.api as sm

from NSCLC_native_MFA_interpretation_inputs_v1 import load, sha, VIEWS
from NSCLC_native_MFA_clinical_inputs_v1 import read, matched_quota, COMMON
from HNSCC_author_state_clinical_kernel_v2 import seed,fit


def main():
    data=load();root,pre,plan,scope,pairs,arrays=read();matched_quota(root)
    proof=pre/'private_independent_clinical_input_audit.json';assert not proof.exists()
    old=root/'NSCLC_native_MFA_clinical_coverage_20261016_v1/private_all208_clinical_eligibility.tsv'
    c=pd.read_csv(old,sep='\t').set_index('sampleID').loc[data['cases'].sampleID]
    mroot=root/'NSCLC_native_MFA_interpretation_20261016_v1';mt=json.loads((mroot/'private_terminal.json').read_text())
    assert sha(mroot/'private_terminal.json')=='91f42e9d7422f0cf878bf15ee22fa1f42ab0f695674db38de073508cea8d9c13'
    assert sha(mroot/'private_all_mixture_models_and_scores.npz')==mt['files']['private_all_mixture_models_and_scores.npz']['sha256']
    with np.load(mroot/'private_all_mixture_models_and_scores.npz',allow_pickle=False) as z:mix={k:z[k] for k in z.files if k.startswith('mixture_residual_')}
    checked=0;draws=0;max_projection_error=0.;max_synthetic_error=0.
    for (part,view),a in arrays.items():
        vi=VIEWS.index(view);original=np.ones(len(c),bool) if part=='all_frozen_four_view_cases' else data['cases'].partition.eq(part).to_numpy()
        take=original&c.clinical_eligible.to_numpy(bool);p=data['prepared'][vi];scores=data['model']['own_view_scores_'+str(vi)]
        assert np.array_equal(a['cases'],data['cases'].sampleID.to_numpy(str)[take])
        assert np.array_equal(a['strata'],data['cases'].holdout_block.to_numpy(str)[take])
        assert np.array_equal(a['X'],c.loc[take,'response'].to_numpy())
        sd=scores[data['train']].std(axis=0,ddof=1);mean=scores[data['train']].mean(axis=0)
        assert np.allclose(a['Y'],(scores[take]-mean)/sd,atol=1e-12,rtol=1e-12)
        assert np.allclose(a['Y_mixture'],mix['mixture_residual_'+str(vi)][take]/sd,atol=1e-12,rtol=1e-12)
        cov=np.column_stack([c.loc[take,COMMON].to_numpy(float),p['quality'][take]])
        raw=np.column_stack([np.ones(take.sum()),(cov-cov.mean(axis=0))/cov.std(axis=0)])
        assert np.allclose(raw,a['Zraw'],atol=1e-12,rtol=1e-12)
        q,r,pivot=qr(raw,mode='economic',pivoting=True)
        s=np.linalg.svd(raw,compute_uv=False);rank=int((s>s[0]*1e-12).sum());q=q[:,:rank]
        expected_projection=q@q.T;actual_projection=a['D'][:,:-1]@a['D'][:,:-1].T
        error=float(np.max(abs(expected_projection-actual_projection)));assert error<1e-10
        max_projection_error=max(max_projection_error,error)
        # Synthetic vector only, not an observed clinical outcome fit.
        y=np.sin(np.arange(len(raw))*.37)+np.cos(np.arange(len(raw))*.19)
        reference=sm.OLS(y,np.column_stack([raw,a['X']])).fit(cov_type='HC3')
        actual=fit(y,a['D']);error=float(max(abs(actual['beta']-reference.params[-1]),abs(actual['t']-reference.tvalues[-1])))
        assert error<1e-9;max_synthetic_error=max(max_synthetic_error,error)
        checked+=a['Y'].size
        if part=='internal_holdout':
            for row in scope[scope.view==view].itertuples():
                assert row.design_sha256==sha(pre/('private_'+part+'_'+view+'.npz'))
                rng=np.random.default_rng(seed(row.record,'bootstrap'));indices=np.empty((2000,len(raw)),int);offset=0
                for level in sorted(set(a['strata'])):
                    ids=np.where(a['strata']==level)[0];indices[:,offset:offset+len(ids)]=ids[rng.integers(0,len(ids),(2000,len(ids)))];offset+=len(ids)
                valid=0
                for start in range(0,2000,128):
                    sv=np.linalg.svd(a['D'][indices[start:start+128]],compute_uv=False)
                    rank1=(sv>sv[:,[0]]*1e-12).sum(axis=1);rank2=(sv>sv[:,[0]]*np.finfo(float).eps*max(a['D'].shape)).sum(axis=1)
                    good=(rank1==a['D'].shape[1])&(rank1==rank2)&(sv[:,0]/np.maximum(sv[:,-1],np.finfo(float).tiny)<1e8);valid+=int(good.sum())
                assert valid==row.bootstrap_design_valid;draws+=2000
    assert checked==16160 and draws==80000
    lookup=scope.set_index('record')
    for row in pairs.itertuples():
        a=lookup.loc[row.record_A];b=lookup.loc[row.record_B]
        assert a.factor==b.factor==row.factor and a.view==row.view_A and b.view==row.view_B
        assert np.array_equal(arrays[('internal_holdout',a.view)]['cases'],arrays[('internal_holdout',b.view)]['cases'])
    report=dict(status='PASS_ALL40_NATIVE_CLINICAL_AND60_PAIR_INPUTS_FULL_CONTROL_SPANS_AND80000_BOOTSTRAP_DESIGN_DRAWS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(pre/'private_terminal.json'),standardized_score_values_checked=checked,
        bootstrap_design_draws_checked=draws,maximum_projection_error=max_projection_error,maximum_synthetic_HC3_error=max_synthetic_error,
        actual_clinical_effects_computed=False,source_patients_reclassified=False)
    proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
