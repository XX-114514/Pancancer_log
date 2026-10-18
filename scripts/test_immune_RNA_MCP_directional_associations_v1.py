"""Complete fixed 394-position directional RNA association family."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from test_immune_RNA_MCP_source_coordination_v1 import fit_rank_partial, resample, adjust_family


def point(x,y,C):
    f,status=fit_rank_partial(x,y,C)
    return (f['effect'] if f is not None else np.nan),status


def raw_rho(x,y):
    return float(spearmanr(x,y).statistic) if len(x)>=3 and min(np.ptp(x),np.ptp(y))>0 else np.nan


def leave_one(x,y,C):
    effects=[]
    for i in range(len(x)):
        take=np.arange(len(x))!=i;effect,_=point(x[take],y[take],C[take]);effects.append(effect)
    effects=np.asarray(effects)
    return dict(leave_one_estimable_units=int(np.isfinite(effects).sum()),
        leave_one_min=float(np.nanmin(effects)) if np.isfinite(effects).any() else np.nan,
        leave_one_max=float(np.nanmax(effects)) if np.isfinite(effects).any() else np.nan)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,512*1024**2)
    source=root/'immune_RNA_MCP_disjoint_receiver_preflight_20261018_v1'
    acceptance=root/'immune_RNA_MCP_disjoint_receiver_preflight_audit_20261018_v1/private_terminal.json'
    accepted=json.loads(acceptance.read_text());assert accepted['status']=='ACCEPTED_ALL_DISJOINT_READOUTS_NODE_BASELINES_ENDPOINTS_AND_COVERAGE'
    assert accepted['producer_terminal_sha256']==sha(source/'private_terminal.json')=='469df1863a6770168a7c67fd78078455b586f654ef7be3e350ede9a7f3592f2b'
    for name,item in json.loads((source/'private_terminal.json').read_text())['files'].items():assert sha(source/name)==item['sha256']
    registry=pd.read_csv(source/'private_all_unique_directional_endpoint_slots.tsv',sep='\t');assert len(registry)==394
    inputs=pd.read_csv(source/'private_all_slot_case_model_inputs.tsv.gz',sep='\t')
    coverage=pd.read_csv(source/'private_all_source_model_coverage.tsv',sep='\t').set_index(['slot','dataset'])
    controls=[c+'_'+s for s in ['S','R'] for c in ['fraction','logN','cellQ','geneQ']]+['sender_composition','receiver_composition','receptor_RNA']
    out=root/'immune_RNA_MCP_directional_associations_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        preflight_terminal_sha256=sha(source/'private_terminal.json'),independent_preflight_acceptance_sha256=sha(acceptance),
        scope='All394 unique explicit endpoint positions; fixed120 selection slots, no receptor/sender/candidate replacement. Original 2026-10-18 RNA MCP v1,K2 source-adjusted.',
        main='TNBC sample units with original receptor1pct/complete-design/fidelity qualification; ligand limiting-component RNA vs disjoint receiver programme, original11 required controls. Two-sided partial rank correlation.',
        random=dict(seed=2026101811,residual_permutations=19999,block_bootstrap=2000,minimum_valid_bootstrap_for_CI=1000),
        multiplicity='BH394 and BY394, all394 planned slots; unavailable p=1 only for denominator, displayed P/q remain missing.',
        sensitivities='Same-case original8,original8+two-node point effects, LR RNA product proxy,10pct receptor support,reviewed-patient subset,all-unit leave-one, all2 discovery source descriptions. No additional P/CI family.',
        limitations='Discovery descriptions not validation. No clinical endpoint, independent patient or causal receptor-binding claim. Alternative receptor positions sharing ligand/programme not independent biological evidence. Conditional fixed-weight/node-baseline uncertainty.',
        privacy='No learned candidates, coefficients, case rows or scientific effects public Git; effect review after complete independent random-stream acceptance.')
    (out/'plan_before_any_directional_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    rows=[];sens=[];discovery=[];modes=['original8_same_cases','original10_same_cases','LR_product_proxy','receptor10pct','reviewed_patient_only']
    expected_status={
        'FEWER_THAN20_COMPLETE_UNITS':'FEWER_THAN20_COMPLETE_UNITS',
        'NONFINITE_COMPLETE_DESIGN':'NONFINITE_DESIGN','CONSTANT_SCORE':'CONSTANT_PREDICTOR_OR_READOUT',
        'CONSTANT_REQUIRED_COVARIATE_NO_AUTOMATIC_DELETION':'CONSTANT_REQUIRED_COVARIATE',
        'COLLINEAR_REQUIRED_DESIGN_NO_AUTOMATIC_DELETION':'COLLINEAR_REQUIRED_DESIGN',
        'ILL_CONDITIONED_REQUIRED_DESIGN':'ILL_CONDITIONED_REQUIRED_DESIGN',
        'FEWER_THAN10_FULL_MODEL_RESIDUAL_DF':'FEWER_THAN10_FULL_MODEL_RESIDUAL_DF',
        'SCORE_EXPLAINED_BY_REQUIRED_COVARIATES':'SCORE_EXPLAINED_BY_REQUIRED_CONTROLS','ESTIMABLE':'ESTIMABLE'}
    for item in registry.itertuples(index=False):
        all_cases=inputs[inputs.slot.eq(item.slot)]
        hold=all_cases[all_cases.dataset.eq('GSE305078')];d=hold[hold.primary_eligible]
        x=d.ligand_RNA.to_numpy();y=d.receiver_score.to_numpy();C=d[controls].to_numpy(float)
        fitted,status=fit_rank_partial(x,y,C);prior=coverage.loc[(item.slot,'GSE305078')]
        assert expected_status[status]==prior.design_status_before_fidelity
        if not item.all3_source_fidelity_pass:fitted=None;status='DISJOINT_RECEIVER_NOT_FAITHFUL'
        row=dict(slot=item.slot,prototype_id=item.prototype_id,program=item.program,sender_view=item.sender_view,receiver_view=item.receiver_view,
            dataset='GSE305078',source_units=len(hold),complete_units=len(d),reviewed_patient_units=int(d.reviewed_patient.sum()),status=status,
            unadjusted_rho=raw_rho(x,y),adjusted_rho=np.nan,lower=np.nan,upper=np.nan,permutation_p=np.nan,valid_bootstrap=0,
            design_condition=np.nan,residual_df=np.nan,leave_one_estimable_units=0,leave_one_min=np.nan,leave_one_max=np.nan)
        if fitted is not None:
            result,draws=resample(x,y,C,fitted,int(item.slot),2026101811)
            row.update(result,adjusted_rho=fitted['effect'],design_condition=fitted['condition'],residual_df=fitted['residual_df'],**leave_one(x,y,C))
            np.savez_compressed(out/('private_slot_%04d_random_stream.npz'%item.slot),samples=d['sample'].to_numpy(str),
                x=x,y=y,controls=C,control_names=np.asarray(controls),**draws)
        rows.append(row)
        for mode in modes:
            s=d.copy();names=controls.copy();xx=x
            if mode=='original8_same_cases':names=controls[:8]
            elif mode=='original10_same_cases':names=controls[:10]
            elif mode=='LR_product_proxy':xx=np.log1p(s.ligand_RNA.to_numpy()*s.receptor_RNA.to_numpy())
            elif mode=='receptor10pct':s=s[s.receptor_support10pct];xx=s.ligand_RNA.to_numpy()
            elif mode=='reviewed_patient_only':s=s[s.reviewed_patient];xx=s.ligand_RNA.to_numpy()
            effect,reason=point(xx,s.receiver_score.to_numpy(),s[names].to_numpy(float))
            sens.append(dict(slot=item.slot,dataset='GSE305078',mode=mode,complete_units=len(s),point_effect=effect,status=reason,no_new_P_or_CI=True))
        for dataset,g in all_cases[~all_cases.dataset.eq('GSE305078')].groupby('dataset',sort=True):
            s=g[g.primary_eligible];xx=s.ligand_RNA.to_numpy();yy=s.receiver_score.to_numpy();cc=s[controls].to_numpy(float)
            effect,reason=point(xx,yy,cc)
            discovery.append(dict(slot=item.slot,dataset=dataset,source_units=len(g),complete_units=len(s),reviewed_patient_units=int(s.reviewed_patient.sum()),
                unadjusted_rho=raw_rho(xx,yy),adjusted_rho=effect,status=reason,**leave_one(xx,yy,cc),fitted_source_not_validation=True))
        if len(rows)%10==0 or len(rows)==394:
            (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),planned_slots=394,
                completed_slots=len(rows),completed_estimable_slots=sum(r['status']=='ESTIMABLE' for r in rows),partial_effects_review_allowed=False))+'\n')
    main=pd.DataFrame(rows);assert len(main)==394 and len(sens)==1970 and len(discovery)==788
    main['family_adjustment_p']=main.permutation_p.fillna(1)
    main['BH394']=adjust_family(main.family_adjustment_p);main['BY394']=adjust_family(main.family_adjustment_p,by=True)
    main.loc[main.permutation_p.isna(),['BH394','BY394']]=np.nan
    main.to_csv(out/'private_all394_main_directional_tests.tsv',sep='\t',index=False)
    pd.DataFrame(sens).to_csv(out/'private_all1970_point_sensitivities.tsv',sep='\t',index=False)
    pd.DataFrame(discovery).to_csv(out/'private_all788_discovery_source_descriptions.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()};assert sum(p['bytes'] for p in files.values())<512*1024**2
    (out/'private_terminal.json').write_text(json.dumps(dict(status='COMPLETE_ALL394_DIRECTIONAL_POSITIONS_AND_RANDOM_STREAMS_ACCEPTANCE_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,
        planned_main_slots=394,estimable_main_slots=int(main.status.eq('ESTIMABLE').sum()),sensitivity_positions=1970,discovery_descriptions=788,
        clinical_or_causal_validation=False),indent=2)+'\n')
    print('ALL394_DIRECTIONAL_TEST_POSITIONS_COMPLETE_RANDOM_ACCEPTANCE_PENDING')


if __name__=='__main__':main()
