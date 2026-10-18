"""Independent complete-case, receptor support and nuisance-span eligibility."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from NSCLC_native_MFA_clinical_inputs_v1 import COMMON
from NSCLC_external_exact_rank_guards_v1 import rational_rank,any_exactly_fitted
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def inspect_design(a,b,quality,blocks):
    if len(a)<50:raise ValueError('INSUFFICIENT_INFORMATIVE_UNITS')
    if not (np.isfinite(a).all() and np.isfinite(b).all()):raise ValueError('NONFINITE_PROGRAMME')
    if not np.isfinite(quality).all():raise ValueError('NONFINITE_QUALITY')
    qr=rankdata(quality,axis=0);sd=qr.std(axis=0);standard=(qr-qr.mean(axis=0))/np.where(sd>0,sd,1)
    dummies=[(blocks==s).astype(int) for s in sorted(set(blocks))]
    z=np.column_stack([np.ones(len(a)),*dummies,standard]);y=rankdata(np.column_stack([a,b]),axis=0)
    coef,_,rank,singular=np.linalg.lstsq(z,y,rcond=1e-10)
    if len(a)-rank<10:raise ValueError('QUALITY_AND_BLOCK_DESIGN_INSUFFICIENT_RESIDUAL_DF')
    residual=y-z@coef
    if np.linalg.norm(residual,axis=0).min()<1e-6:
        zi=np.column_stack([np.ones(len(a)),*dummies,2*qr]).astype(np.int64)
        if rational_rank(zi)==rank and any_exactly_fitted(zi,(2*y).astype(np.int64)):
            raise ValueError('INDIVIDUAL_PROGRAMME_RESIDUAL_CONSTANT')
    if (np.linalg.norm(residual,axis=0)<1e-10).any():raise ValueError('INDIVIDUAL_PROGRAMME_RESIDUAL_CONSTANT')
    machine=int((singular>singular[0]*np.finfo(float).eps*max(z.shape)).sum())
    if rank!=machine or singular[rank-1]/singular[0]<1e-8:raise ValueError('UNSTABLE_FULL_NUISANCE_SPAN')
    return int(rank),len(a)-int(rank)


def main():
    data=load();root=data['root'];resource(root,1048576)
    out=root/'NSCLC_continuous_prior_native_association_inputs_20261016_v1'
    terminal=json.loads((out/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(out/name)==item['sha256']
    plan=json.loads((out/'plan_before_any_native_interface_association.json').read_text())
    assert (plan['permutations'],plan['bootstrap'],plan['min_valid_bootstrap'])==(199999,2000,1900)
    frozen=json.loads((out/'private_frozen_family_before_effects.json').read_text())
    src=root/'NSCLC_continuous_prior_native_interfaces_20261016_v2'
    proof=src/'private_independent_native_interface_audit.json';assert sha(proof)==frozen['input_audit_sha256']
    accepted=json.loads(proof.read_text());assert accepted['producer_terminal_sha256']==sha(src/'private_terminal.json')
    with np.load(src/'private_all_native_interface_measurements.npz',allow_pickle=False) as saved:z={k:saved[k] for k in saved.files}
    control_path=root/'NSCLC_native_MFA_clinical_coverage_20261016_v1/private_all208_clinical_eligibility.tsv'
    assert sha(control_path)==frozen['clinical_covariate_source_sha256']
    c=pd.read_csv(control_path,sep='\t',usecols=['sampleID']+COMMON).set_index('sampleID')
    clinical=c.loc[data['cases'].sampleID,COMMON].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    with np.load(out/'private_all208_covariates_only.npz',allow_pickle=False) as saved:
        assert np.array_equal(clinical,saved['clinical'],equal_nan=True) and saved['clinical_names'].tolist()==COMMON
        assert np.array_equal(saved['cases'],z['cases']) and np.array_equal(saved['train'],z['train'])
    interfaces=pd.read_csv(src/'private_all4800_native_interface_readiness.tsv.gz',sep='\t',keep_default_na=False)
    scope=pd.read_csv(out/'private_all4800_pre_effect_eligibility.tsv.gz',sep='\t',keep_default_na=False).set_index('interface_index')
    observed_edges=pd.read_csv(out/'private_all_fixed_receptor_case_coverage.tsv.gz',sep='\t',keep_default_na=False)
    observed_parts=pd.read_csv(out/'private_all36_covariate_partition_coverage.tsv',sep='\t')
    views={v:i for i,v in enumerate(VIEWS)};receptors={r:i for i,r in enumerate(z['receptors'])};parts={};expected_parts=[]
    for sender in VIEWS:
        for receiver in VIEWS:
            if sender==receiver:continue
            q=np.column_stack([clinical,data['prepared'][views[sender]]['quality'],data['prepared'][views[receiver]]['quality']])
            complete=np.isfinite(q).all(axis=1);block=data['cases'].holdout_block.to_numpy(str)
            for part,mask in [('internal_holdout',~data['train']),('discovery',data['train']),('all_frozen_cases',np.ones(len(block),bool))]:
                candidate=np.flatnonzero(complete&mask);selected=[]
                for b in sorted(set(block[candidate])):
                    ids=candidate[block[candidate]==b]
                    if len(ids)>=2:selected.extend(ids.tolist())
                take=np.array(sorted(selected),int);parts[(sender,receiver,part)]=(take,q[take],block[take])
                expected_parts.append(dict(sender_view=sender,receiver_view=receiver,partition=part,informative_cases=len(take),category_blocks=len(set(block[take]))))
    pd.testing.assert_frame_equal(pd.DataFrame(expected_parts),observed_parts,check_dtype=False)
    all_edges=[];eligible=0
    for interface in interfaces.itertuples():
        actual=scope.loc[interface.interface_index];sv=views[interface.sender_view];rv=views[interface.receiver_view]
        hi,hq,hb=parts[(interface.sender_view,interface.receiver_view,'internal_holdout')];di,dq,db=parts[(interface.sender_view,interface.receiver_view,'discovery')]
        assert actual.record=='NATIVE_PRIOR_INTERFACE|%04d'%interface.interface_index
        for key in ['slot_index','sender_view','receiver_view','factor','direction','prospective_position','ligand','ligand_axis_index']:assert actual[key]==getattr(interface,key)
        assert actual.active_cases==len(hi) and actual.discovery_cases==len(di)
        reasons=[];compatible=[];rank=df=0
        if interface.status!='RNA_INPUT_READY_EFFECT_MODEL_NOT_RUN':reasons.append('INPUT:'+interface.status)
        else:
            j=interface.ligand_axis_index
            if min((z['ligand_positive'][sv,hi,j]>0).sum(),(z['ligand_positive'][sv,di,j]>0).sum())<10:reasons.append('SENDER_POSITIVE_IN_FEWER_THAN10_CASES_IN_EITHER_PARTITION')
            for receptor in interface.receiver_RNA_compatible_receptors.split('|'):
                k=receptors[receptor];hn=int((z['receptor_positive'][rv,hi,k]>0).sum());dn=int((z['receptor_positive'][rv,di,k]>0).sum());ok=min(hn,dn)>=10
                all_edges.append(dict(interface_index=interface.interface_index,receptor=receptor,holdout_joint_positive_cases=hn,discovery_joint_positive_cases=dn,both_partition10case_gate=ok))
                if ok:compatible.append(receptor)
            if not compatible:reasons.append('NO_FIXED_RECEPTOR_WITH10_CASES_IN_BOTH_PARTITIONS')
            for label,ix,q,b in [('PRIMARY',hi,hq,hb),('DISCOVERY',di,dq,db)]:
                try:
                    rr,dd=inspect_design(z['ligand_mean'][sv,ix,j],interface.direction*z['response'][interface.slot_index,ix],q,b)
                    if label=='PRIMARY':rank,df=rr,dd
                except ValueError as error:reasons.append(label+':'+str(error))
        expected=';'.join(reasons) if reasons else 'ELIGIBLE_FOR_FIXED_PRIMARY_TEST'
        assert actual.status==expected and actual.qualifying_reference_receptors=='|'.join(compatible)
        assert actual.nuisance_rank==rank and actual.residual_df==df
        if not reasons:assert actual.test_index==eligible;eligible+=1
        else:assert actual.test_index==-1
    pd.testing.assert_frame_equal(pd.DataFrame(all_edges),observed_edges,check_dtype=False)
    assert eligible==frozen['eligible_primary_records'] and len(scope)==frozen['all_planned_positions']==4800
    assert frozen['no_response_clinical_eligibility_or_pathology_columns_loaded'] and not frozen['associations_or_P_calculated']
    result=dict(status='PASS_ALL4800_OUTCOME_FREE_COVARIATE_CASES_FIXED_RECEPTOR_COVERAGE_OLS_NUISANCE_AND_FROZEN_FAMILY',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),eligible_primary_records=eligible,all_planned_positions=4800,clinical_outcomes_or_association_effects_loaded=False)
    path=out/'private_independent_association_input_audit.json';assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],finished_UTC=result['finished_UTC'],audit_sha256=sha(path))))


if __name__=='__main__':main()
