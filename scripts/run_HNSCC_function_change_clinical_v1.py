"""Complete fixed clinical families; atomically publish each private numeric model."""
import hashlib,json,os,subprocess,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import beta as beta_dist
from statsmodels.stats.multitest import multipletests
from HNSCC_function_clinical_kernel_v1 import fit,wild_count,bootstrap,loo,seed


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def atomic(path,value):
    assert not path.exists();tmp=path.with_suffix(path.suffix+'.tmp');assert not tmp.exists();tmp.write_text(json.dumps(value,indent=2)+'\n');os.replace(tmp,path)


def wait(r,directory,proof,job):
    p=r/directory/proof;rc=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2'/job/'launch_receipt.json').read_text());deadline=time.monotonic()+28800
    while not p.exists():
        stat=Path('/proc')/str(rc['pid'])/'stat';assert stat.exists();a=stat.read_text().rsplit(')',1)[1].split();assert a[0]!='Z' and a[19]==str(rc['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(5)
    assert json.loads(p.read_text())['status'].startswith('PASS_');return p


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);pre=r/'HNSCC_function_change_clinical_preflight_20261016_v1';cal=r/'HNSCC_function_clinical_calibration_20261016_v1';input_proof=wait(r,pre.name,'private_independent_clinical_input_audit.json','audit_HNSCC_function_change_clinical_inputs_v1');calibration=wait(r,cal.name,'private_terminal.json','calibrate_HNSCC_function_clinical_null_v1');contract=json.loads((cal/'plan_before_any_real_clinical_effect.json').read_text());parameters=contract['formal_inference'];kernel=Path(__file__).with_name('HNSCC_function_clinical_kernel_v1.py');assert sha(kernel)==contract['kernel_sha256']==json.loads(calibration.read_text())['kernel_sha256']
    terminal=json.loads((pre/'private_terminal.json').read_text());assert sha(pre/'private_terminal.json')==json.loads(input_proof.read_text())['producer_terminal_sha256']
    for name,h in terminal['files'].items():assert sha(pre/name)==h
    model=pd.read_csv(pre/'private_unique_numeric_designs.tsv',sep='\t').sort_values('model_id').reset_index(drop=True);scope=pd.read_csv(pre/'private_all8000_clinical_scope.tsv',sep='\t',keep_default_na=False);assert len(scope)==8000
    out=r/'HNSCC_function_change_clinical_results_20261016_v1';assert not out.exists();fs=subprocess.check_output(['df','-P',str(r.resolve())],text=True).splitlines()[-1].split()[0];quota=subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True);qr=[s.split() for s in quota.splitlines() if s.split() and s.split()[0]==fs];assert len(qr)==1;free=(int(qr[0][3])-int(qr[0][1].rstrip('*')))*1024;assert free>500_000_000_000+1_073_741_824+268_435_456
    out.mkdir(mode=0o700);(out/'private_models').mkdir(mode=0o700);plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),kernel_sha256=sha(kernel),formal_contract_sha256=sha(cal/'plan_before_any_real_clinical_effect.json'),input_audit_sha256=sha(input_proof),calibration_sha256=sha(calibration),preflight_terminal_sha256=sha(pre/'private_terminal.json'),parameters=parameters,unique_numeric_designs=len(model),primary_family=4800,secondary_family=3200,quota_free_bytes=free,output_cap_bytes=268_435_456,no_old_family_or_identity_change=True,clinical_or_direct_communication_validation_accepted=False)
    atomic(out/'plan_before_observed_clinical_effects.json',plan);all_models=[];size=0
    for number,row in enumerate(model.itertuples()):
        identifier=row.model_id;p=pre/'private_designs'/(identifier+'.npz');assert sha(p)==row.design_sha256;archive=np.load(p,allow_pickle=False);d=archive['D'];y=archive['Y'];strata=archive['strata'];B=int(parameters['wild_draws']);draws=int(parameters['bootstrap_draws']);record=dict(model_id=identifier,design_sha256=row.design_sha256,n=int(len(y)),residual_df=int(len(y)-d.shape[1]),status='UNESTIMABLE_AFTER_OUTCOME',beta=np.nan,HC3_t=np.nan,HC3_se=np.nan,P=np.nan,MC_upper_P=np.nan,wild_exceedances=-1,valid_bootstraps=0,CI_low=np.nan,CI_high=np.nan,LOO_valid=0,LOO_same_sign=False,LOO_min=np.nan,LOO_max=np.nan,incremental_partial_R2=np.nan)
        boots=np.full(draws,np.nan);leave=np.full(len(y),np.nan)
        try:g=fit(y,d,min_df=10)
        except ValueError as error:record['status']='UNESTIMABLE_'+str(error)
        else:
            exceed=wild_count(g,B,seed(identifier,'wild'));pv=(exceed+1)/(B+1);upper=1. if exceed==B else float(beta_dist.ppf(1-.01/len(model),exceed+1,B-exceed));upper=min(1.,(1+B*upper)/(B+1));boots=bootstrap(y,d,strata,draws,seed(identifier,'bootstrap'));leave=loo(y,d);valid=np.isfinite(boots);lo,hi=np.quantile(boots[valid],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan];z=d[:,:-1];reduced=y-z@np.linalg.lstsq(z,y,rcond=None)[0];error=y-d@(np.linalg.pinv(d,rcond=1e-12)@y)
            record.update(status='ESTIMATED_EXPLORATORY_FUNCTION_CHANGE_CLINICAL_ASSOCIATION',beta=g['beta'],HC3_t=g['t'],HC3_se=float(np.sqrt(g['variance'])),P=pv,MC_upper_P=upper,wild_exceedances=exceed,valid_bootstraps=int(valid.sum()),CI_low=float(lo),CI_high=float(hi),LOO_valid=int(np.isfinite(leave).sum()),LOO_same_sign=bool(np.isfinite(leave).all() and (np.sign(leave)==np.sign(g['beta'])).all()),LOO_min=float(np.nanmin(leave)) if np.isfinite(leave).any() else np.nan,LOO_max=float(np.nanmax(leave)) if np.isfinite(leave).any() else np.nan,incremental_partial_R2=float(1-(error@error)/(reduced@reduced)))
        resample=out/'private_models'/(identifier+'.npz');np.savez_compressed(resample,bootstrap_coefficients=boots,LOO_coefficients=leave);record['resamples_sha256']=sha(resample);record['completed_UTC']=datetime.now(timezone.utc).isoformat();atomic(out/'private_models'/(identifier+'.json'),record);all_models.append(record);size+=resample.stat().st_size+(out/'private_models'/(identifier+'.json')).stat().st_size;assert size<plan['output_cap_bytes']
        (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),unique_models_finished=number+1,unique_models_planned=len(model),last_model_id=identifier,private_model_bytes=size))+'\n')
    models=pd.DataFrame(all_models).set_index('model_id');models.to_csv(out/'private_all_unique_clinical_model_results.tsv',sep='\t',na_rep='NA');rows=[]
    for original in scope.to_dict('records'):
        result={**original,'association_status':'INPUT_UNAVAILABLE','clinical_P':np.nan,'MC_upper_clinical_P':np.nan}
        if original['status']=='READY_NO_OBSERVED_CLINICAL_EFFECT':
            endpoints=['A'] if original['family']=='single3200' else ['A','B'];ready=True
            for key in endpoints:
                item=models.loc[original[key+'_model']];ready &= str(item.status).startswith('ESTIMATED_')
                for name in ['beta','HC3_t','HC3_se','P','MC_upper_P','CI_low','CI_high','valid_bootstraps','LOO_valid','LOO_same_sign','LOO_min','LOO_max','incremental_partial_R2']:result[key+'_'+name]=item[name]
            if ready:
                result.update(association_status='ESTIMATED_SHARED_MARGINAL_CLINICAL_ASSOCIATION' if len(endpoints)==2 else 'ESTIMATED_SINGLE_POOL_CLINICAL_ASSOCIATION',clinical_P=max(result[key+'_P'] for key in endpoints),MC_upper_clinical_P=max(result[key+'_MC_upper_P'] for key in endpoints))
        rows.append(result)
    results=pd.DataFrame(rows);assert len(results)==8000
    for family,index in results.groupby('family').groups.items():
        assert len(index)==(4800 if family=='pair4800' else 3200);valid=results.loc[index,'association_status'].str.startswith('ESTIMATED_')
        for name,method,column in [('BH_q','fdr_bh','clinical_P'),('BY_q','fdr_by','clinical_P'),('MC_upper_BH_q','fdr_bh','MC_upper_clinical_P')]:
            results.loc[index,name]=multipletests(results.loc[index,column].fillna(1).to_numpy(float),method=method)[1];results.loc[np.array(index)[~valid.to_numpy()],name]=np.nan
    results.to_csv(out/'private_all8000_longitudinal_function_clinical_results.tsv',sep='\t',index=False,na_rep='NA');assert sum(p.stat().st_size for p in out.rglob('*') if p.is_file())<plan['output_cap_bytes']
    terminal=dict(status='COMPLETE_ALL_UNIQUE_CLINICAL_MODELS_AND4800_SHARED3200_SINGLE_FULL_FAMILIES_INDEPENDENT_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),unique_numeric_models=len(model),logical_slots=len(results),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_direct_communication_or_independent_validation_or_prediction=False)
    atomic(out/'private_terminal.json',terminal);print(json.dumps(terminal))


if __name__=='__main__':main()
