"""Complete matched-case conditional clinical family after input/calibration PASS."""
import hashlib,json,os,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.stats import beta as beta_dist
from statsmodels.stats.multitest import multipletests
from paired_clinical_increment_kernel_v1 import fit,wild_count,resample,loo


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def wait_proof(r,path,job):
    p=r/path;j=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2'/job/'launch_receipt.json').read_text());deadline=time.monotonic()+7200
    while not p.exists():
        stat=Path('/proc')/str(j['pid'])/'stat';assert stat.exists();s=stat.read_text().rsplit(')',1)[1].split()
        assert s[19]==str(j['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    assert json.loads(p.read_text())['status'].startswith('PASS_');return p


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);pre=r/'NSCLC_paired_programme_clinical_increment_preflight_20261016_v1'
    input_proof=wait_proof(r,pre.name+'/private_independent_input_audit.json','audit_NSCLC_paired_clinical_increment_inputs_v1')
    calibration=wait_proof(r,'paired_clinical_increment_calibration_20261016_v2/private_terminal.json','calibrate_paired_clinical_increment_v2')
    oldplan=json.loads((pre/'plan_before_any_paired_clinical_increment_effect.json').read_text());kernel=Path(__file__).with_name('paired_clinical_increment_kernel_v1.py')
    assert sha(kernel)==oldplan['kernel_sha256']==json.loads(calibration.read_text())['kernel_sha256']
    terminal=json.loads((pre/'private_terminal.json').read_text());assert sha(pre/'private_terminal.json')==json.loads(input_proof.read_text())['producer_terminal_sha256']
    for n,h in terminal['files'].items():assert sha(pre/n)==h
    scope=pd.read_csv(pre/'private_all572_scope.tsv',sep='\t');hashes=pd.read_csv(pre/'private_ready_design_SHA256.tsv',sep='\t').set_index('record').sha256
    out=r/'NSCLC_paired_programme_clinical_increment_results_20261016_v1';assert not out.exists();out.mkdir(mode=0o700);(out/'private_resamples').mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),kernel_sha256=sha(kernel),preflight_plan_sha256=sha(pre/'plan_before_any_paired_clinical_increment_effect.json'),
              input_audit_sha256=sha(input_proof),calibration_sha256=sha(calibration),family=572,wild_draws=oldplan['wild_draws'],bootstrap_draws=2000,seed=oldplan['seed'],
              primary='IUT max of BOTH conditional HC3 restricted-HC2 wild P, one complete572 BH; fixed full standardized design and original binary response.',
              no_direct_communication_or_causality_or_prediction_claim=True)
    (out/'plan_before_true_paired_clinical_effects.json').write_text(json.dumps(plan,indent=2)+'\n');rows=[];resample_hashes=[]
    for i,row in enumerate(scope.itertuples()):
        z=row._asdict();z.update(association_status='INPUT_UNAVAILABLE',pair_IUT_P=np.nan,pair_MC_upper_IUT_P=np.nan);rows.append(z)
        if row.status!='READY_NO_CLINICAL_EFFECT_COMPUTED':continue
        p=pre/'private_designs'/(row.record+'.npz');assert sha(p)==hashes.loc[row.record];a=np.load(p,allow_pickle=False);y=a['Y'];d=a['D']
        g=fit(y,d,min_df=10);seed=plan['seed']+i*10;B=plan['wild_draws']
        for key,target,offset in [('A',-2,0),('B',-1,1)]:
            exceed=wild_count(g,target,B,seed+offset);pv=(1+exceed)/(B+1)
            upper=1. if exceed==B else float(beta_dist.ppf(1-.01/1144,exceed+1,B-exceed))
            upper=min(1.,(1+B*upper)/(B+1))
            z.update({key+'_beta':float(g['beta'][target]),key+'_HC3_t':float(g['t'][target]),key+'_HC3_se':float(np.sqrt(g['variance'][target])),
                      key+'_wild_exceedances':exceed,key+'_P':pv,key+'_MC_upper_P':upper})
        z['pair_IUT_P']=max(z['A_P'],z['B_P']);z['pair_MC_upper_IUT_P']=max(z['A_MC_upper_P'],z['B_MC_upper_P'])
        boot=resample(y,d,2000,seed+2);leave=loo(y,d);p=out/'private_resamples'/(row.record+'.npz');np.savez_compressed(p,bootstrap_coefficients=boot,LOO_coefficients=leave)
        resample_hashes.append(dict(record=row.record,sha256=sha(p)))
        for j,key in enumerate(['A','B']):
            v=boot[:,j];valid=np.isfinite(v);z[key+'_valid_bootstraps']=int(valid.sum())
            lo,hi=np.quantile(v[valid],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan]
            z[key+'_CI_low']=lo;z[key+'_CI_high']=hi;u=leave[:,j]
            z[key+'_LOO_valid']=int(np.isfinite(u).sum());z[key+'_LOO_same_sign']=bool(np.isfinite(u).all() and (np.sign(u)==np.sign(g['beta'][-2+j])).all())
            z[key+'_LOO_min']=float(np.nanmin(u));z[key+'_LOO_max']=float(np.nanmax(u))
            single=np.delete(d,-1 if j==0 else -2,axis=1);bf=np.linalg.lstsq(single,y,rcond=None)[0]
            z[key+'_matched_single_end_beta']=float(bf[-1])
            reduced=np.delete(d,-2+j,axis=1);er=y-reduced@np.linalg.lstsq(reduced,y,rcond=None)[0]
            z[key+'_incremental_partial_R2']=float(1-(g['residual']@g['residual'])/(er@er))
        mask=np.column_stack([a['masked_A'],a['masked_B']]);mask_ok=np.isfinite(mask).all() and (mask.std(0)>1e-12).all()
        z['nonshared_conditional_status']='NONSHARED_READOUT_UNAVAILABLE'
        if mask_ok:
            dm=np.column_stack([d[:,:-2],(mask-mask.mean(0))/mask.std(0)])
            try:
                mg=fit(y,dm,check_hc3=False);z['nonshared_A_beta']=float(mg['beta'][-2]);z['nonshared_B_beta']=float(mg['beta'][-1]);z['nonshared_conditional_status']='DESCRIPTIVE_ONLY_NO_NEW_P'
            except ValueError as error:z['nonshared_conditional_status']='UNESTIMABLE_'+str(error)
        z['association_status']='ESTIMATED_PAIRED_CONDITIONAL_CLINICAL_EXPLORATORY'
        checkpoint=dict(observed_UTC=datetime.now(timezone.utc).isoformat(),models_completed=len(resample_hashes),last_record=row.record,total_planned=572)
        (out/'private_checkpoint.json').write_text(json.dumps(checkpoint)+'\n')
    t=pd.DataFrame(rows);assert len(t)==572 and t.record.is_unique
    valid=t.association_status.eq('ESTIMATED_PAIRED_CONDITIONAL_CLINICAL_EXPLORATORY')
    for name,method,pvals in [('BH_q_pair572','fdr_bh',t.pair_IUT_P.fillna(1)),('BY_q_pair572','fdr_by',t.pair_IUT_P.fillna(1)),('MC_upper_BH_q_pair572','fdr_bh',t.pair_MC_upper_IUT_P.fillna(1))]:
        t[name]=multipletests(pvals,method=method)[1];t.loc[~valid,name]=np.nan
    t.to_csv(out/'private_all572_paired_clinical_results.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(resample_hashes).to_csv(out/'private_resample_SHA256.tsv',sep='\t',index=False)
    terminal=dict(status='COMPLETE_ALL572_PAIRED_CLINICAL_INFERENCE_INDEPENDENT_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),planned_models=572,estimated_models=int(valid.sum()),
                  code_sha256=sha(__file__),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},independent_patient_or_direct_communication_or_baseline_prediction=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
