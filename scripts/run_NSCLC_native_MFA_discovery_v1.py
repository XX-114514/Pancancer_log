"""Official primary fit and two fixed stability fits; clinical outcomes never read."""
import contextlib,importlib.metadata,itertools,json,os,subprocess,sys,time
from datetime import datetime,timezone
from pathlib import Path
import h5py,numpy as np,pandas as pd
from scipy.optimize import linear_sum_assignment
from NSCLC_native_MFA_math_v1 import sha,ridge_score


def wait(root,directory,filename,job):
    p=root/directory/filename;receipt=json.loads((root/'HNSCC_receiver_followup_job_control_20261016_v2'/job/'launch_receipt.json').read_text());deadline=time.monotonic()+7200
    while not p.exists():
        stat=Path('/proc')/str(receipt['pid'])/'stat';assert stat.exists();s=stat.read_text().rsplit(')',1)[1].split();assert s[0]!='Z' and s[19]==str(receipt['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(3)
    assert json.loads(p.read_text())['status'].startswith('PASS_');return p


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);source=root/'NSCLC_native_MFA_inputs_20261016_v1';input_proof=wait(root,source.name,'private_independent_MFA_input_audit.json','audit_NSCLC_native_MFA_inputs_v1');runtime_proof=wait(root,'NSCLC_MFA_runtime_20261016_v1','private_terminal.json','prepare_NSCLC_MFA_runtime_v1');terminal=json.loads((source/'private_terminal.json').read_text());assert sha(source/'private_terminal.json')==json.loads(input_proof.read_text())['producer_terminal_sha256']
    for name,h in terminal['files'].items():assert sha(source/name)==h
    contract=json.loads((source/'plan_before_native_expression_and_MFA_fit.json').read_text());options=contract['model_before_any_fit'];assert options['factors']==10 and options['primary_seed']==20261017 and options['stability_seeds']==[20261018,20261019]
    sys.path.insert(0,str(root/'NSCLC_MFA_runtime_20261016_v1/library'));assert importlib.metadata.version('mofapy2')=='0.7.2';from mofapy2.run.entry_point import entry_point
    views=contract['primary_views'];arrays={}
    for v in views:
        a=np.load(source/('private_'+v+'_prepared.npz'),allow_pickle=False);arrays[v]={k:a[k] for k in a.files}
    train=arrays[views[0]]['train'];cases=arrays[views[0]]['cases'].astype(str);assert train.sum()==145
    for v in views:assert np.array_equal(arrays[v]['cases'],cases) and np.array_equal(arrays[v]['train'],train)
    out=root/'NSCLC_native_MFA_discovery_20261016_v1';assert not out.exists()
    fs=subprocess.check_output(['df','-P',str(root.resolve())],text=True).splitlines()[-1].split()[0];qt=subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True);qr=[s.split() for s in qt.splitlines() if s.split() and s.split()[0]==fs];assert len(qr)==1;free=(int(qr[0][3])-int(qr[0][1].rstrip('*')))*1024;assert free>500_000_000_000+1_073_741_824+268_435_456;out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),projection_helper_sha256=sha(Path(__file__).with_name('NSCLC_native_MFA_math_v1.py')),input_audit_sha256=sha(input_proof),runtime_proof_sha256=sha(runtime_proof),frozen_input_and_model_contract_sha256=sha(source/'plan_before_native_expression_and_MFA_fit.json'),options=options,scope='Primary20261017 and stability20261018/19 all fixed; four views145 training cases only. No holdout or clinical choice, all10 factor slots preserved.',projection='Own-view frozen-weight precision-weighted ridge, prior penalty1; not exact variational posterior, no cross-view information used in a view score.',quota_free_bytes=free,output_cap_bytes=268435456,no_direct_LR_or_clinical_or_patient_independence_claim=True)
    (out/'plan_before_any_native_MFA_fit.json').write_text(json.dumps(plan,indent=2)+'\n')
    # A meaningful numerical fixture compares independent augmented least squares.
    rng=np.random.default_rng(2026101741)
    for _ in range(4):
        w=rng.normal(size=(120,10));t=rng.lognormal(size=120);x=rng.normal(size=(63,120));left=np.vstack([np.sqrt(t)[:,None]*w,np.eye(10)]);right=np.vstack([(x*np.sqrt(t)).T,np.zeros((10,63))]);expected=np.linalg.lstsq(left,right,rcond=None)[0].T;assert np.allclose(expected,ridge_score(x,w,t),rtol=1e-10,atol=1e-11)
    records=[];posterior={}
    for seed in [options['primary_seed']]+options['stability_seeds']:
        directory=out/('seed_'+str(seed));directory.mkdir(mode=0o700)
        with (directory/'private_official_fit.log').open('x') as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
            ent=entry_point();ent.set_data_options(scale_views=False,scale_groups=False,center_groups=False,use_float32=False);ent.set_data_matrix([[arrays[v]['X'][train].copy()] for v in views],likelihoods=['gaussian']*4,views_names=views,groups_names=['NSCLC_discovery'],samples_names=[cases[train].tolist()],features_names=[arrays[v]['selected_genes'].astype(str).tolist() for v in views]);ent.set_model_options(factors=10,spikeslab_factors=False,spikeslab_weights=False,ard_factors=False,ard_weights=True);ent.set_train_options(iter=1000,convergence_mode='medium',dropR2=None,quiet=True,seed=seed,gpu_mode=False);ent.build();ent.run();e=ent.model.getExpectations();z=e['Z']['E'].copy();assert z.shape==(145,10) and np.isfinite(z).all();weights=[e['W'][m]['E'].copy() for m in range(4)];precision=[]
            for m,w in enumerate(weights):
                assert w.shape==(arrays[views[m]]['X'].shape[1],10) and np.isfinite(w).all();t=np.asarray(e['Tau'][m]['E'],float)
                if t.ndim==2:assert t.shape==(145,w.shape[0]) and np.allclose(t,t[:1],rtol=1e-10,atol=1e-12);t=t[0].copy()
                assert t.shape==(w.shape[0],) and np.isfinite(t).all() and (t>0).all();precision.append(t)
            ent.save(str(directory/'private_official_model.hdf5'),save_data=False,save_parameters=True,expectations=['W','Z','Tau'])
        # Canonical sign depends only on discovery weights, never on outcomes.
        stack=np.vstack(weights);maximum=np.argmax(abs(stack),axis=0);orientation=np.sign(stack[maximum,np.arange(10)]);orientation[orientation==0]=1;z*=orientation
        weights=[w*orientation for w in weights];save=dict(discovery_factor_scores=z,orientation_sign=orientation,cases=cases,train=train,views=np.array(views));scores=[]
        for m,view in enumerate(views):
            score=ridge_score(arrays[view]['X'],weights[m],precision[m]);scores.append(score);save['weights_'+str(m)]=weights[m];save['precision_'+str(m)]=precision[m];save['own_view_scores_'+str(m)]=score
        np.savez_compressed(directory/'private_frozen_parameters_and_own_view_scores.npz',**save)
        with h5py.File(directory/'private_official_model.hdf5') as f:elbo=f['training_stats/elbo'][()];times=f['training_stats/time'][()]
        finite=np.flatnonzero(np.isfinite(elbo));assert len(finite)>1;last=int(finite[-1]);record=dict(seed=seed,role='primary' if seed==options['primary_seed'] else 'stability_only',discovery_cases=145,factors=10,elbo_observations=len(finite),last_recorded_iteration=last,maximum_iteration_reached=last>=999,final_ELBO=float(elbo[last]),official_model_sha256=sha(directory/'private_official_model.hdf5'),parameter_archive_sha256=sha(directory/'private_frozen_parameters_and_own_view_scores.npz'),all_four_views_finite=True,holdout_used_for_fit=False)
        (directory/'private_fit_terminal.json').write_text(json.dumps(record,indent=2)+'\n');records.append(record);posterior[seed]=z;(out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),fixed_seeds_finished=len(records),fixed_seeds_planned=3))+'\n');assert sum(p.stat().st_size for p in out.rglob('*') if p.is_file())<plan['output_cap_bytes']
    primary=posterior[options['primary_seed']];stability=[]
    for seed in options['stability_seeds']:
        other=posterior[seed];corr=np.corrcoef(primary.T,other.T)[:10,10:];a,b=linear_sum_assignment(-abs(corr))
        for i,j in zip(a,b):stability.append(dict(primary_factor='NSCLC_MFA_20261016_F'+str(i+1).zfill(2),stability_seed=seed,matched_other_factor=int(j+1),discovery_score_correlation=float(corr[i,j]),absolute_correlation=float(abs(corr[i,j])),holdout_used_for_alignment=False))
    pd.DataFrame(records).to_csv(out/'private_all_three_fixed_seed_fit_status.tsv',sep='\t',index=False);pd.DataFrame(stability).to_csv(out/'private_all20_discovery_only_factor_stability_matches.tsv',sep='\t',index=False)
    result=dict(status='COMPLETE_THREE_FIXED_NATIVE_MOFA_FITS_AND_OWN_VIEW_PROJECTIONS_NO_HOLDOUT_TEST_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),files={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()},primary_factor_slots=10,stability_fits=2,no_clinical_or_direct_communication_result=True);(out/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
