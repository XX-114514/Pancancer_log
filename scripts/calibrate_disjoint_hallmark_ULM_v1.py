"""Synthetic official-source and independent two-group OLS checks before RNA scores."""
import ast,hashlib,json,os,warnings
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,statsmodels.api as sm
from disjoint_hallmark_ULM_kernel_v1 import score


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);p=r/'fullgene_receiver_target_code_v1/collectri_reference_v1/official_decoupler192_ULM_source.txt'
    assert sha(p)=='17954dbf35dbb10f358a330d271588867fc54be11552997726992353157664a1'
    ns={'np':np};f=[n for n in ast.parse(p.read_text()).body if isinstance(n,ast.FunctionDef) and n.name in {'mat_cov','mat_cor','t_val'}]
    assert len(f)==3;exec(compile(ast.Module(body=f,type_ignores=[]),str(p),'exec'),ns)
    out=r/'NSCLC_disjoint_hallmark_ULM_calibration_20261016_v1';assert not out.exists();out.mkdir(mode=0o700);rng=np.random.default_rng(2026101614)
    maximum=0.;checked=0;ols_checks=0
    for g in [80,501,3000]:
        for trial in range(8):
            x=rng.lognormal(0,1,(17,g));x[rng.random(x.shape)<.5]=0
            targets=[np.sort(rng.choice(g,size=k,replace=False)) for k in [20,25,min(70,g-1)]]
            got=score(x,targets);w=np.zeros((g,3))
            for j,ix in enumerate(targets):w[ix,j]=1
            expected=ns['t_val'](ns['mat_cor'](w,x.T),g-2)
            assert np.allclose(got['ULM_score'],expected,rtol=1e-9,atol=1e-9)
            maximum=max(maximum,float(np.max(np.abs(got['ULM_score']-expected))));checked+=expected.size
            for i in range(17):
                for j,ix in enumerate(targets):
                    fit=sm.OLS(x[i],np.column_stack([np.ones(g),w[:,j]])).fit()
                    assert np.isclose(fit.tvalues[-1],got['ULM_score'][i,j],rtol=1e-8,atol=1e-8)
                    rest=np.ones(g,bool);rest[ix]=False
                    assert np.isclose(x[i,ix].mean(),got['target_mean'][i,j]) and np.isclose(x[i,rest].mean(),got['nontarget_mean'][i,j])
                    assert np.count_nonzero(x[i,ix])==got['expressed_target_genes'][i,j];ols_checks+=1
    x=np.ones((3,80));x[1,:20]=2;x[2,:20]=0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore');z=score(x,[np.arange(20)])
    assert np.isnan(z['ULM_score'][0,0]) and z['ULM_score'][1,0]>0 and z['ULM_score'][2,0]<0
    assert z['expressed_target_genes'][2,0]==0 and z['target_mean'][2,0]==0
    # Excluding learned physical columns is independent of their expression values.
    full=rng.normal(size=(7,110));bg=np.arange(10,110);a=score(full[:,bg],[np.arange(20)])
    full[:,:10]=1e12;z=score(full[:,bg],[np.arange(20)]);assert np.array_equal(a['ULM_score'],z['ULM_score'])
    result=dict(status='PASS_SYNTHETIC_SPARSE_ULM_OFFICIAL192_ALL_SCALAR_OLS_MEANS_DETECTION_ZERO_AND_EXCLUDED_COLUMN_CHECKS',finished_UTC=datetime.now(timezone.utc).isoformat(),
                code_sha256=sha(__file__),kernel_sha256=sha(Path(__file__).with_name('disjoint_hallmark_ULM_kernel_v1.py')),official_sha256=sha(p),official_scores_checked=checked,scalar_OLS_checks=ols_checks,maximum_official_difference=maximum,real_expression_used=False)
    (out/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
