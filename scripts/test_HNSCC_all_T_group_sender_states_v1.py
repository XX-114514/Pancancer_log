"""All predeclared T groups, common frozen receiver and global families."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from test_HNSCC_frozen_T_sender_states_v1 import family,bh


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    inputs=root/'HNSCC_all_T_group_sender_states_20261014_v1'
    out=root/'HNSCC_all_T_group_sender_associations_20261014_v1';assert not out.exists();out.mkdir()
    contract=pd.read_csv(inputs/'private_complete_group_candidate_contract.tsv',sep='\t')
    group_names=sorted(contract.T_candidate_group.unique());assert len(group_names)==3
    (out/'execution_contract_before_complete_family_effects.json').write_text(json.dumps(dict(
      created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      group_names=group_names,scope_plan_sha256=hashlib.sha256((root/'HNSCC_all_T_subtype_sender_capture_20261014_v1/full_T_scope_before_sender_effect_interpretation.json').read_bytes()).hexdigest(),
      BH='all T groups *admissible unique signatures in one family per model, including missing/insufficient internally P1',
      all_candidate_family='all T groups *all unique signatures, separately diagnostic',
      original_CD8='reuse all original192 coefficient/permutation/bootstrap results; only add complete-family multiplicity',
      new_group_seed='20261014 + stable SHA of candidate group and model',
      receiver_selected_after_previous_clinical_analysis=True,not_confirmatory_clinical_or_causal=True),indent=2)+'\n')
    old=pd.read_csv(root/'HNSCC_frozen_T_sender_associations_20261014_v1/all_private_sender_state_associations.tsv',sep='\t')
    columns=['T_logRNA','T_loggenes','Mac_logRNA','Mac_loggenes','T_fraction','Mac_fraction','T_logcells','Mac_logcells']
    all_rows=[]
    for group in group_names:
        candidates=contract[contract.T_candidate_group.eq(group)]
        if group=='CD8_Effector':
            for row in old.to_dict('records'):row['T_candidate_group']=group;all_rows.append(row)
            continue
        directory=inputs/group
        md=pd.read_csv(directory/'private_patient_receiver_covariates.tsv',sep='\t',index_col=0)
        raw=pd.read_csv(directory/'private_patient_T_state_means.tsv',sep='\t',index_col=0).loc[md.index]
        orth=pd.read_csv(directory/'private_patient_T_orthogonal_state_means.tsv',sep='\t',index_col=0).loc[md.index]
        models=[('primary_RNA_capture',raw,md[columns].to_numpy(float)),('sensitivity_unadjusted',raw,np.empty((len(md),0))),
          ('sensitivity_demographic',raw,md[columns+['D_plus_T','HPV_positive','male','age']].to_numpy(float)),
          ('sensitivity_sender_gene_overlap_removed',orth,md[columns].to_numpy(float))]
        y=md.receiver_programme.to_numpy()
        for name,frame,cov in models:
            group_col='orthogonal_group' if name.endswith('overlap_removed') else 'observed_signature_group'
            representatives=candidates.sort_values('programme_id').drop_duplicates(group_col)
            existing=representatives[representatives.programme_id.isin(frame.columns)]
            random_seed=20261014+int(hashlib.sha256((group+'|'+name).encode()).hexdigest()[:8],16)
            fit=family(frame[existing.programme_id].to_numpy(),y,cov,random_seed)
            by_group={g:r for g,r in zip(existing[group_col],fit)}
            for signature,meta in candidates.groupby(group_col,sort=True):
                flags=bool(meta.foreign_RNA_sensitive.any() or meta.basic_expression_review_flag.any())
                eligible=bool(meta.primary_admissible.any());assert not (eligible and flags)
                row=dict(T_candidate_group=group,model=name,signature_group=signature,programme_ids='|'.join(sorted(meta.programme_id)),
                  primary_admissible=eligible,n_patients=len(md),nuisance_columns=cov.shape[1],
                  foreign_RNA_flag=bool(meta.foreign_RNA_sensitive.any()),basic_expression_flag=bool(meta.basic_expression_review_flag.any()))
                row.update(by_group.get(signature,dict(status='signature_not_measurable',beta=np.nan,t=np.nan,p=np.nan,lower=np.nan,upper=np.nan,valid_bootstrap=0)))
                all_rows.append(row)
    result=pd.DataFrame(all_rows);assert len(result)==576 and not result.duplicated(['T_candidate_group','model','signature_group']).any()
    result['BH_all_T_groups_admissible_family']=np.nan;result['BH_all_T_groups_all_candidate_family']=np.nan
    summaries=[]
    for model,index in result.groupby('model').groups.items():
        f=result.loc[index];assert len(f)==144
        result.loc[index,'BH_all_T_groups_all_candidate_family']=np.where(f.p.notna(),bh(f.p.fillna(1)),np.nan)
        eligible=f.primary_admissible.to_numpy(bool);keep=f.index[eligible]
        q=bh(result.loc[keep,'p'].fillna(1));result.loc[keep,'BH_all_T_groups_admissible_family']=np.where(result.loc[keep,'p'].notna(),q,np.nan)
        summaries.append(dict(model=model,all_unique_tests=len(f),admissible_unique_tests=len(keep),estimated_tests=int(f.p.notna().sum())))
    result.to_csv(out/'all576_private_T_group_associations.tsv',sep='\t',index=False,na_rep='NA')
    terminal=dict(status='FULL_ALL_T_GROUP_FROZEN_STATE_ASSOCIATION_FAMILIES_COMPLETE_PRIVATE',finished_UTC=datetime.now(timezone.utc).isoformat(),
      groups=group_names,families=summaries,original_CD8_results_preserved=True,clinical_P_values_added=False,formal_states_or_direct_interaction_verified=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
