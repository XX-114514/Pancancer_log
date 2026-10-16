"""Freeze complete inference and diagnose small-sample fixed-design null size."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from HNSCC_function_clinical_kernel_v1 import fit,wild_count,bootstrap,seed


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);pre=r/'HNSCC_author_state_clinical_preflight_20261016_v1';rc=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/preflight_HNSCC_author_state_clinical_v1/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (pre/'private_terminal.json').exists():
        p=Path('/proc')/str(rc['pid'])/'stat';assert p.exists();a=p.read_text().rsplit(')',1)[1].split();assert a[0]!='Z' and a[19]==str(rc['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(5)
    frozen=json.loads((pre/'plan_before_author_state_clinical_effects.json').read_text());assert frozen['kernel_sha256']==sha(Path(__file__).with_name('HNSCC_function_clinical_kernel_v1.py'))
    models=pd.read_csv(pre/'private_unique_numeric_designs.tsv',sep='\t');assert len(models)>=3;out=r/'HNSCC_author_state_clinical_calibration_20261016_v1';assert not out.exists();out.mkdir(mode=0o700);kernel=Path(__file__).with_name('HNSCC_function_clinical_kernel_v1.py')
    smallest=models.sort_values(['n','model_id']).model_id.tolist();leverage=models.sort_values(['max_leverage','model_id'],ascending=[False,True]).model_id.tolist();middle=models.sort_values(['max_leverage','model_id']).model_id.tolist();middle=middle[len(middle)//2:]+middle[:len(middle)//2];chosen=[]
    for candidates in [smallest,leverage,middle]:chosen.append(next(x for x in candidates if x not in chosen))
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),kernel_sha256=sha(kernel),preflight_terminal_sha256=sha(pre/'private_terminal.json'),unique_design_scope_sha256=sha(pre/'private_unique_numeric_designs.tsv'),selected_designs=chosen,selection='Smallest n, highest leverage, median leverage; unique tie-broken numeric SHA; no observed clinical coefficients.',formal_inference=dict(primary_family=221,secondary_single_family=26,supplementary_nested_pair=20,supplementary_nested_single=9,wild_draws=199999,bootstrap_draws=2000,bootstrap='Whole patient resampling within fixed reviewed treatment x HPV strata; fixed original rank scale, no outcome selection; invalid full-rank/condition draws retained missing, CI only>=1900 valid.',clinical_test='Single endpoint HC3 t; restricted null HC2-adjusted residual Rademacher wild bootstrap, approximate working-model inference, not randomized treatment or exact small-sample proof.',primary_pair='max of A/B clinical P on identical patients; shared/marginal association, no counterpart state fraction in own model; whole221 primary BH/BY and MCupper BH; auxiliary26, nested20 and nested9 separately corrected.',MC='Plus-one P; simultaneous .99 Clopper-Pearson upper bound across unique numerical endpoint designs; corresponding family correction. No precision escalation after results.',interval='2000 conditional-stratum patient bootstrap percentile coefficient interval, all draws kept with invalid flags. LOO influence descriptive; no new P.',robust_candidate='Both endpoints CI exclude0, all LOO valid/same sign and primary BH/BY/MCupper as separate recorded gates. Not direct LR or independent validation.'),diagnostics=dict(simulations_each=1024,wild_draws_each=1999,error_models=['normal','t3_unit_variance','heteroskedastic_normal'],null_mean='Fixed analysis-rank-scale nuisance linear mean, not re-ranking generated Y; evaluates working-model bootstrap behaviour only.',gate='Every one of9 null rejection rates at alpha .05 <=.075; diagnostic tolerance, no exact-size guarantee; any fail prevents real effects.',power='Additional strong target coefficient1.0 normal-error scenario, descriptive, not a gate.',seed=2026101631),new_observed_clinical_effects=False)
    (out/'plan_before_any_real_clinical_effect.json').write_text(json.dumps(plan,indent=2)+'\n');rng=np.random.default_rng(plan['diagnostics']['seed']);rows=[]
    for number,identifier in enumerate(chosen):
        p=pre/'private_designs'/(identifier+'.npz');assert sha(p)==models.set_index('model_id').loc[identifier,'design_sha256'];archive=np.load(p,allow_pickle=False);d=archive['D'];strata=archive['strata'];n=len(d);mean=d[:,:-1]@np.linspace(-.2,.2,d.shape[1]-1)
        for distribution in plan['diagnostics']['error_models']+['strong_target_normal']:
            rejected=0
            for k in range(1024):
                noise=rng.standard_t(3,n)/np.sqrt(3) if distribution=='t3_unit_variance' else rng.normal(size=n)
                if distribution=='heteroskedastic_normal':noise*=np.exp(.6*d[:,1])
                y=mean+noise+(d[:,-1] if distribution=='strong_target_normal' else 0.);g=fit(y,d,min_df=10);pv=(1+wild_count(g,1999,2026101631+number*100000+(['normal','t3_unit_variance','heteroskedastic_normal','strong_target_normal'].index(distribution))*2000+k))/2000;rejected+=pv<.05
                if (k+1)%128==0:(out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),design_number=number+1,error_model=distribution,simulations=k+1))+'\n')
            rows.append(dict(model_id=identifier,n=n,error_model=distribution,simulations=1024,rejections=int(rejected),rejection_rate=rejected/1024))
        # Meaningful coefficient parity and fixed-stratum sampling fixture.
        import statsmodels.api as sm
        test=mean+.5*d[:,-1]+rng.normal(size=n);g=fit(test,d);m=sm.OLS(test,d).fit(cov_type='HC3');assert np.allclose([g['beta'],g['t']],[m.params[-1],m.tvalues[-1]],atol=1e-9,rtol=1e-9)
        b=bootstrap(test,d,strata,16,seed(identifier,'calibration_boot'));assert np.isfinite(b).sum()>=12
    pd.DataFrame(rows).to_csv(out/'private_all12_null_and_power_diagnostics.tsv',sep='\t',index=False);passed=all(x['rejection_rate']<=.075 for x in rows if x['error_model']!='strong_target_normal')
    terminal=dict(status='PASS_SMALL_N_COMPLETE_NUISANCE_WORKING_MODEL_NULL_DIAGNOSTICS' if passed else 'FAIL_SMALL_N_NULL_DIAGNOSTICS_DO_NOT_RUN_REAL_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),kernel_sha256=sha(kernel),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_observed_clinical_effects=False,exact_small_sample_or_causal_guarantee=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
