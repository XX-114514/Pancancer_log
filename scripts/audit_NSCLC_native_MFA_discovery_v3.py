"""Every saved posterior, own-view QR projection, sample key and stability match."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import h5py,numpy as np,pandas as pd
from scipy.optimize import linear_sum_assignment


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as stream:
        for b in iter(lambda:stream.read(16777216),b''):h.update(b)
    return h.hexdigest()


def strings(dataset):return [x.decode() if isinstance(x,bytes) else str(x) for x in dataset[()]]


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);source=root/'NSCLC_native_MFA_discovery_20261016_v3';receipt=json.loads((root/'HNSCC_receiver_followup_job_control_20261016_v2/run_NSCLC_native_MFA_discovery_v3/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (source/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[0]!='Z' and s[19]==str(receipt['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(3)
    terminal=json.loads((source/'private_terminal.json').read_text());assert terminal['code_sha256']==receipt['source_sha256']
    for name,h in terminal['files'].items():assert sha(source/name)==h
    plan=json.loads((source/'plan_before_any_native_MFA_fit.json').read_text());inputs=root/'NSCLC_native_MFA_inputs_20261016_v1';contract=json.loads((inputs/'plan_before_native_expression_and_MFA_fit.json').read_text());views=contract['primary_views'];assert plan['frozen_input_and_model_contract_sha256']==sha(inputs/'plan_before_native_expression_and_MFA_fit.json')
    axis=pd.read_csv(inputs/'private_case_axis.tsv',sep='\t',dtype={'sampleID':str});train=axis.partition.eq('discovery').to_numpy();cases=axis.sampleID.tolist();data={v:np.load(inputs/('private_'+v+'_prepared.npz'),allow_pickle=False) for v in views};status=pd.read_csv(source/'private_all_three_fixed_seed_fit_status.tsv',sep='\t');assert status.seed.tolist()==[20261017,20261018,20261019];maximum=0.;all_z={};checks=[]
    for seed in status.seed:
        directory=source/('seed_'+str(seed));a=np.load(directory/'private_frozen_parameters_and_own_view_scores.npz',allow_pickle=False);assert a['cases'].astype(str).tolist()==cases and np.array_equal(a['train'],train) and a['views'].astype(str).tolist()==views;canonical=np.vstack([a['weights_'+str(m)] for m in range(4)]);z=a['discovery_factor_scores'];assert z.shape==(145,10);all_z[seed]=z
        with h5py.File(directory/'private_official_model.hdf5') as f:
            assert strings(f['samples/NSCLC_discovery'])==axis.loc[train,'sampleID'].tolist();saved=[]
            for view in views:
                assert strings(f['features/'+view])==[view+'::'+g for g in data[view]['selected_genes'].astype(str)];w=f['expectations/W/'+view][()];assert w.shape==(10,data[view]['X'].shape[1]);saved.append(w.T)
            stack=np.vstack(saved);cos=canonical.T@stack/(np.linalg.norm(canonical,axis=0)[:,None]*np.linalg.norm(stack,axis=0)[None,:]);ii,jj=linear_sum_assignment(-abs(cos));assert np.array_equal(ii,np.arange(10)) and (abs(cos[ii,jj])>1-1e-10).all();sign=np.sign(cos[ii,jj]);raw_z=f['expectations/Z/NSCLC_discovery'][()].T;assert np.allclose(raw_z[:,jj]*sign,z,rtol=1e-9,atol=1e-10)
            for m,view in enumerate(views):
                w=a['weights_'+str(m)];t=a['precision_'+str(m)];assert np.allclose(saved[m][:,jj]*sign,w,rtol=1e-9,atol=1e-10);noise=f['expectations/Tau/'+view][()]
                if noise.shape==(w.shape[0],145):noise=noise.T
                if noise.ndim==2:assert noise.shape==(145,w.shape[0]) and np.allclose(noise,noise[:1],rtol=1e-10,atol=1e-12);noise=noise[0]
                assert np.allclose(noise,t,rtol=1e-10,atol=1e-12) and (t>0).all();left=np.vstack([np.sqrt(t)[:,None]*w,np.eye(10)]);q,rr=np.linalg.qr(left,mode='reduced');right=np.vstack([(data[view]['X']*np.sqrt(t)).T,np.zeros((10,len(cases)))]);expected=np.linalg.solve(rr,q.T@right).T;scores=a['own_view_scores_'+str(m)];assert scores.shape==(208,10) and np.allclose(scores,expected,rtol=1e-9,atol=1e-10);maximum=max(maximum,float(np.max(abs(scores-expected))))
            if seed==20261017:
                with h5py.File(root/'NSCLC_native_MFA_discovery_20261016_v2/seed_20261017/private_official_model.hdf5') as old:
                    for path in ['expectations/W/'+v for v in views]+['expectations/Tau/'+v for v in views]+['expectations/Z/NSCLC_discovery']:
                        assert np.allclose(f[path][()],old[path][()],rtol=1e-10,atol=1e-11), 'Same-seed primary parameter parity failed'
            elbo=f['training_stats/elbo'][()];finite=np.flatnonzero(np.isfinite(elbo));row=status[status.seed.eq(seed)].iloc[0];assert row.elbo_observations==len(finite) and row.last_recorded_iteration==finite[-1] and np.isclose(row.final_ELBO,elbo[finite[-1]],rtol=1e-12) and not row.holdout_used_for_fit
        checks.append(dict(seed=int(seed),saved_W_Z_Tau_verified=True,all832_view_case10_factor_projections_replayed=True,discovery_only_sample_and_feature_order_exact=True))
    stability=pd.read_csv(source/'private_all20_discovery_only_factor_stability_matches.tsv',sep='\t');assert len(stability)==20
    for seed in [20261018,20261019]:
        x=all_z[20261017];y=all_z[seed];xx=x-x.mean(axis=0);yy=y-y.mean(axis=0);corr=xx.T@yy/np.sqrt(np.sum(xx**2,axis=0)[:,None]*np.sum(yy**2,axis=0)[None,:]);ii,jj=linear_sum_assignment(-abs(corr));g=stability[stability.stability_seed.eq(seed)].set_index('primary_factor')
        for i,j in zip(ii,jj):
            row=g.loc['NSCLC_MFA_20261016_F'+str(i+1).zfill(2)];assert row.matched_other_factor==j+1 and np.isclose(row.discovery_score_correlation,corr[i,j],atol=1e-12) and np.isclose(row.absolute_correlation,abs(corr[i,j]),atol=1e-12) and not row.holdout_used_for_alignment
    proof=dict(status='PASS_ALL_THREE_OFFICIAL_MODELS_SAVED_POSTERIORS_24960_OWN_VIEW_SCORES_AND20_TRAIN_ONLY_STABILITY_MATCHES',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),all_three_models_checked=checks,maximum_QR_projection_difference=maximum,same_seed_primary_parameter_parity_with_preserved_v2=True,primary_maximum_iteration_reached=bool(status[status.seed.eq(20261017)].maximum_iteration_reached.iloc[0]),independent_patient_or_direct_LR_validation=False);p=source/'private_independent_MFA_discovery_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
