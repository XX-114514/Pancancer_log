"""Eight frozen current-RNA exploratory clinical associations, all private."""
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


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,32*1024**2)
    qualification=root/'HNSCC_immune_RNA_MCP_clinical_metadata_20261019_v1'
    accepted=qualification/'private_independent_key_acceptance_20261019_v1.json'
    assert sha(accepted)=='d3b462a3017fe214b8c4c700cea87702a1b46f874fdd6a7a8bc3179c7b242751'
    metadata=root/'immune_RNA_MCP_source_coordination_20261018_v1/private_all_case_view_controls.tsv'
    decomposition=root/'immune_RNA_MCP_original_node_decomposition_20261018_v1'
    audit=root/'immune_RNA_MCP_original_node_decomposition_audit_20261018_v1/private_terminal.json'
    assert sha(audit)=='f109cecf33be6ad8a937bce8c50d567aa1426f1a1e1ea05a0c20335ce4e123d7'
    assert sha(metadata)=='4a13825e483a8256dd150a061513a0024874b9963e0e8caf4826557bc1849bf6'
    out=root/'HNSCC_immune_RNA_MCP_clinical_20261019_v2';assert not out.exists();out.mkdir(mode=0o700)
    views=['R_FOXP3_capture','E_effector_capture','Myeloid'];programs=['MCP1','MCP2']
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        key_acceptance_sha256=sha(accepted),node_acceptance_sha256=sha(audit),
        slots=[dict(program=p,views=[v]) for p in programs for v in views]+[dict(program=p,views=views[1:]) for p in programs],
        primary='All8 two-sided baseline pathology rank partial correlations; treatment/HPV plus fraction,logN,cellQ,geneQ,discovery-node prediction per contributing view.',
        joint_score='Minimum of fixed38-discovery standardized E and Myeloid case means. No clinical score, cutoff or direction tuning.',
        rank_design='>=20 unique baseline patients,full rank,condition<1e8,>=10 full residual df; no automatic covariate deletion.',
        seed=2026101813,permutations=19999,bootstrap=2000,valid_bootstrap_for_CI=1800,
        multiplicity='All8 BH/BY; unavailable P=1 internally, q displayed missing.',
        point_sensitivities=['same_case_unadjusted','treatment_HPV_only','author_concordant_full_controls'],
        limitations='Internal exploratory association: source used for RNA discovery; project previously saw outcomes. Clinical association is not independent prediction or direct cell communication. Conditional fixed-weight/node/scale intervals.',
        no_new_expression_download_or_gene_weight_fit=True)
    (out/'plan_before_any_current_RNA_clinical_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    term=json.loads((decomposition/'private_terminal.json').read_text())
    filename='private_all396_case_view_decompositions.tsv'
    assert sha(decomposition/filename)==term['files'][filename]
    scores=pd.read_csv(decomposition/filename,sep='\t')
    meta=pd.read_csv(metadata,sep='\t')
    joined=scores.merge(meta,on=['sample','dataset','view'],validate='many_to_one')
    assert len(joined)==396
    base=pd.read_csv(qualification/'private_baseline_case_clinical_keys.tsv',sep='\t',keep_default_na=False).sort_values('patient_registry')
    assert len(base)==26 and base.patient_registry.is_unique
    assert set(base.neo_type)=={'D','D+T'} and set(base['HPV/p16'])=={'negative','positive'}
    scales=[];rows=[];sens=[];modelrows=[]
    for slot,item in enumerate(plan['slots']):
        frame=base.copy();control_names=['treatment_DplusT','HPV_positive']
        frame['treatment_DplusT']=frame.neo_type.eq('D+T').astype(int)
        frame['HPV_positive']=frame['HPV/p16'].eq('positive').astype(int)
        z=[]
        for v in item['views']:
            one=joined[joined.program.eq(item['program'])&joined.view.eq(v)].set_index('sample')
            discovery=one[one.partition.eq('discovery')]
            assert len(discovery)==38
            mean=float(discovery.score.mean());sd=float(discovery.score.std(ddof=0));assert sd>0
            scales.append(dict(slot=slot,program=item['program'],view=v,mean=mean,sd=sd,training_cases=38))
            b=one.loc[frame.fit_sample]
            assert b.dataset.eq('GSE286827').all() and b.reviewed_patient.all()
            z.append((b.score.to_numpy()-mean)/sd)
            for col in ['fraction','logN','cellQ','geneQ','composition_prediction']:
                name=v+'__'+col;frame[name]=b[col].to_numpy();control_names.append(name)
        frame['score']=z[0] if len(z)==1 else np.minimum(z[0],z[1])
        complete=np.isfinite(frame[['score','regression']+control_names].to_numpy(float)).all(axis=1)
        d=frame[complete].copy();x=d.score.to_numpy();y=d.regression.to_numpy();C=d[control_names].to_numpy(float)
        fitted,status=fit_rank_partial(x,y,C)
        row=dict(slot=slot,program=item['program'],views='|'.join(item['views']),n_patients=len(d),status=status,
            adjusted_rho=np.nan,lower=np.nan,upper=np.nan,permutation_p=np.nan,valid_bootstrap=0,residual_df=np.nan,
            leave_one_estimable=0,leave_one_min=np.nan,leave_one_max=np.nan)
        if fitted is not None:
            result,draws=resample(x,y,C,fitted,slot,2026101813)
            row.update(result,adjusted_rho=fitted['effect'],residual_df=fitted['residual_df'])
            valid=np.isfinite(draws['bootstrap_effects'])
            ci=np.quantile(draws['bootstrap_effects'][valid],[.025,.975]) if valid.sum()>=1800 else [np.nan,np.nan]
            row.update(lower=float(ci[0]),upper=float(ci[1]))
            np.savez_compressed(out/('private_slot_%02d_complete_random_stream.npz'%slot),
                patients=d.patient_registry.to_numpy(str),x=x,y=y,controls=C,control_names=np.asarray(control_names),**draws)
            loo=[]
            for i in range(len(d)):
                keep=np.arange(len(d))!=i;f,_=fit_rank_partial(x[keep],y[keep],C[keep]);loo.append(f['effect'] if f else np.nan)
            loo=np.asarray(loo);row['leave_one_estimable']=int(np.isfinite(loo).sum())
            if np.isfinite(loo).any():row.update(leave_one_min=float(np.nanmin(loo)),leave_one_max=float(np.nanmax(loo)))
        rows.append(row)
        for mode in plan['point_sensitivities']:
            sub=d.copy();effect=np.nan;reason='FEWER_THAN3_COMPLETE_UNITS'
            if mode=='same_case_unadjusted':
                if len(sub)>=3 and min(np.ptp(sub.score),np.ptp(sub.regression))>0:
                    effect=float(spearmanr(sub.score,sub.regression).statistic);reason='DESCRIPTIVE_UNADJUSTED'
            else:
                cols=control_names[:2] if mode=='treatment_HPV_only' else control_names
                if mode=='author_concordant_full_controls':sub=sub[sub.author_binary_concordant.astype(str).str.lower().eq('true')]
                f,reason=fit_rank_partial(sub.score.to_numpy(),sub.regression.to_numpy(),sub[cols].to_numpy(float))
                if f:effect=f['effect']
            sens.append(dict(slot=slot,mode=mode,n_patients=len(sub),point=effect,status=reason,no_new_P_or_CI=True))
        modelrows.append(d.assign(slot=slot))
    result=pd.DataFrame(rows);assert len(result)==8
    internal=result.permutation_p.fillna(1).to_numpy();result['BH8']=adjust_family(internal);result['BY8']=adjust_family(internal,True)
    result.loc[result.permutation_p.isna(),['BH8','BY8']]=np.nan
    result.to_csv(out/'private_all8_clinical_statistics.tsv',sep='\t',index=False)
    pd.DataFrame(sens).to_csv(out/'private_all24_point_sensitivities.tsv',sep='\t',index=False)
    pd.DataFrame(scales).to_csv(out/'private_all_fixed_discovery_scales.tsv',sep='\t',index=False)
    pd.concat(modelrows,ignore_index=True).to_csv(out/'private_all_patient_model_inputs.tsv.gz',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
    assert sum(v['bytes'] for v in files.values())<32*1024**2
    t=dict(status='COMPLETE_ALL8_INTERNAL_CURRENT_RNA_CLINICAL_ASSOCIATIONS_INDEPENDENT_ACCEPTANCE_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n')
    print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
