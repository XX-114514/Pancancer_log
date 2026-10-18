"""Rebuild all eight current-RNA clinical inputs and every random value."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from audit_immune_RNA_MCP_stability_and_coordination_v1 import independent_fit, family_q, rho, equal
from audit_immune_RNA_MCP_directional_associations_v2 import boot_svd, independent_leave


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    source=root/'HNSCC_immune_RNA_MCP_clinical_20261019_v1'
    terminal=json.loads((source/'private_terminal.json').read_text())
    assert terminal['source_sha256']==sha(Path(__file__).with_name('test_HNSCC_immune_RNA_MCP_clinical_v1.py'))
    for name,item in terminal['files'].items():assert sha(source/name)==item['sha256']
    out=root/'HNSCC_immune_RNA_MCP_clinical_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_full_input_and_random_replay.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        producer_terminal_sha256=sha(source/'private_terminal.json'),
        scope='All8 frozen scores, clinical keys, nuisance values, full permutation/bootstrap streams,24 sensitivities,LOO and BH/BY; independent pandas ranks/OLS and batched SVD.'),indent=2)+'\n')
    base=pd.read_csv(root/'HNSCC_immune_RNA_MCP_clinical_metadata_20261019_v1/private_baseline_case_clinical_keys.tsv',sep='\t',keep_default_na=False).sort_values('patient_registry')
    assert len(base)==26 and base.patient_registry.is_unique
    decomposition=pd.read_csv(root/'immune_RNA_MCP_original_node_decomposition_20261018_v1/private_all396_case_view_decompositions.tsv',sep='\t')
    metadata=pd.read_csv(root/'immune_RNA_MCP_source_coordination_20261018_v1/private_all_case_view_controls.tsv',sep='\t').set_index(['sample','view'])
    main=pd.read_csv(source/'private_all8_clinical_statistics.tsv',sep='\t');assert main.slot.tolist()==list(range(8))
    sens=pd.read_csv(source/'private_all24_point_sensitivities.tsv',sep='\t');assert len(sens)==24
    saved_inputs=pd.read_csv(source/'private_all_patient_model_inputs.tsv.gz',sep='\t')
    scale_table=pd.read_csv(source/'private_all_fixed_discovery_scales.tsv',sep='\t')
    views=['R_FOXP3_capture','E_effector_capture','Myeloid']
    slots=[(p,[v]) for p in ['MCP1','MCP2'] for v in views]+[(p,views[1:]) for p in ['MCP1','MCP2']]
    perms=0;boots=0;input_rows=0;maximum_error=0.
    for slot,(program,selected_views) in enumerate(slots):
        C=[base.neo_type.eq('D+T').to_numpy(int),base['HPV/p16'].eq('positive').to_numpy(int)]
        names=['treatment_DplusT','HPV_positive'];values=[]
        for view in selected_views:
            d=decomposition[decomposition.program.eq(program)&decomposition.view.eq(view)].set_index('sample')
            training_samples=metadata[metadata.partition.eq('discovery')].reset_index()
            ids=training_samples.loc[training_samples.view.eq(view),'sample'];assert len(ids)==38
            train=d.loc[ids,'score'].to_numpy(float)
            mean=np.sum(train)/len(train);sd=np.sqrt(np.sum((train-mean)**2)/len(train))
            scale=scale_table[(scale_table.slot==slot)&scale_table.view.eq(view)].iloc[0]
            equal(scale['mean'],mean);equal(scale.sd,sd)
            value=(d.loc[base.fit_sample,'score'].to_numpy()-mean)/sd;values.append(value)
            for col in ['fraction','logN','cellQ','geneQ','composition_prediction']:
                names.append(view+'__'+col)
                if col=='composition_prediction':v=d.loc[base.fit_sample,col].to_numpy(float)
                else:v=np.array([metadata.loc[(s,view),col] for s in base.fit_sample],float)
                C.append(v)
        x=values[0] if len(values)==1 else np.where(values[0]<values[1],values[0],values[1])
        y=base.regression.to_numpy(float);C=np.column_stack(C)
        valid=np.isfinite(np.c_[x,y,C]).all(axis=1);x=x[valid];y=y[valid];C=C[valid];b=base[valid]
        saved=saved_inputs[saved_inputs.slot.eq(slot)]
        assert saved.patient_registry.tolist()==b.patient_registry.tolist()
        assert np.allclose(saved.score,x) and np.allclose(saved.regression,y) and np.allclose(saved[names],C)
        input_rows+=len(b)
        f,status=independent_fit(x,y,C);r=main[main.slot.eq(slot)].iloc[0]
        assert r.status==status and r.n_patients==len(b) and r.program==program and r.views=='|'.join(selected_views)
        equal(r.adjusted_rho,f['effect'] if f else np.nan)
        if f:
            z=np.load(source/('private_slot_%02d_complete_random_stream.npz'%slot),allow_pickle=False)
            assert z['patients'].tolist()==b.patient_registry.tolist() and z['control_names'].tolist()==names
            assert np.allclose(z['x'],x) and np.allclose(z['y'],y) and np.allclose(z['controls'],C)
            rng=np.random.default_rng(np.random.SeedSequence([2026101813,slot,0]))
            indexes=np.stack([rng.permutation(len(b)) for _ in range(19999)]).astype(np.int16)
            assert np.array_equal(indexes,z['permutation_indices'])
            effects=np.empty(19999)
            for start in range(0,19999,1024):
                Y=f['fit_y'][:,None]+f['ry'][indexes[start:start+1024]].T
                residual=Y-f['D']@np.linalg.lstsq(f['D'],Y,rcond=None)[0]
                effects[start:start+residual.shape[1]]=f['rx']@residual/(np.linalg.norm(f['rx'])*np.linalg.norm(residual,axis=0))
            error=float(np.max(np.abs(effects-z['permutation_effects'])));assert error<1e-9;maximum_error=max(maximum_error,error)
            equal(r.permutation_p,(1+np.count_nonzero(np.abs(effects)>=abs(f['effect'])-1e-12))/20000)
            rng=np.random.default_rng(np.random.SeedSequence([2026101813,slot,1]));takes=rng.integers(0,len(b),size=(2000,len(b)),dtype=np.int16)
            assert np.array_equal(takes,z['bootstrap_indices'])
            boot,reasons=boot_svd(x,y,C,takes)
            assert reasons.tolist()==z['bootstrap_reasons'].tolist()
            assert np.allclose(boot,z['bootstrap_effects'],rtol=0,atol=1e-9,equal_nan=True)
            good=np.isfinite(boot);assert int(good.sum())==r.valid_bootstrap
            ci=np.quantile(boot[good],[.025,.975]) if good.sum()>=1800 else [np.nan,np.nan]
            equal(r.lower,ci[0]);equal(r.upper,ci[1]);assert r.residual_df==f['residual_df']
            count,low,high=independent_leave(x,y,C)
            assert count==r.leave_one_estimable;equal(low,r.leave_one_min);equal(high,r.leave_one_max)
            perms+=19999;boots+=2000
        else:assert pd.isna(r.permutation_p) and pd.isna(r.lower) and pd.isna(r.upper)
        for s in sens[sens.slot.eq(slot)].itertuples(index=False):
            xx=x;yy=y;cc=C;n=len(b)
            if s.mode=='same_case_unadjusted':
                effect=rho(x,y);reason='DESCRIPTIVE_UNADJUSTED' if np.isfinite(effect) else 'FEWER_THAN3_COMPLETE_UNITS'
            else:
                if s.mode=='treatment_HPV_only':cc=C[:,:2]
                elif s.mode=='author_concordant_full_controls':
                    keep=b.author_binary_concordant.astype(str).str.lower().eq('true').to_numpy();xx=x[keep];yy=y[keep];cc=C[keep];n=int(keep.sum())
                else:raise AssertionError('UNDECLARED_SENSITIVITY')
                ff,reason=independent_fit(xx,yy,cc);effect=ff['effect'] if ff else np.nan
            assert n==s.n_patients and reason==s.status;equal(s.point,effect)
    p=main.permutation_p.fillna(1).to_numpy()
    for col,by in [('BH8',False),('BY8',True)]:
        expected=family_q(p,by)
        for i,r in enumerate(main.itertuples(index=False)):equal(getattr(r,col),expected[i] if np.isfinite(r.permutation_p) else np.nan)
    t=dict(status='ACCEPTED_ALL8_CURRENT_RNA_CLINICAL_INPUTS_AND_COMPLETE_RANDOM_STREAMS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),
        positions=8,all_point_sensitivities=24,input_rows_checked=input_rows,permutation_values_replayed=perms,
        bootstrap_values_replayed=boots,max_permutation_error=maximum_error,independent_clinical_validation=False)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n')
    print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
