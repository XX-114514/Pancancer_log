"""Complete all already frozen T candidate groups; preserve original CD8 run."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import numpy as np,pandas as pd
from diagnose_frozen_MCP_two_end_functions_v1 import hnscc_inputs,score_functions,seed


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    scope=root/'HNSCC_all_T_subtype_sender_capture_20261014_v1'
    scope_plan=json.loads((scope/'full_T_scope_before_sender_effect_interpretation.json').read_text())
    assert scope_plan['keep_CD8_screen_as_original_nested_round']
    dictionary_path=root/'state_programme_dictionary_20261014_v2/frozen_programme_dictionary.json'
    programmes=[r for r in json.loads(dictionary_path.read_text())['programmes'] if 'T cell' in r['allowed_majors']]
    assert len(programmes)==51
    receiver=pd.read_csv(root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1/private_patient_subtype_scores.tsv',sep='\t')
    receiver=receiver[receiver.subtype_candidate.eq('Macro_C1QC')].set_index('patient_id').sort_index()
    groups=[r['T_subtype_candidate'] for r in json.loads((scope/'aggregate_terminal.json').read_text())['capture']]
    out=root/'HNSCC_all_T_group_sender_states_20261014_v1';assert not out.exists();out.mkdir()
    (out/'plan_before_all_group_scores.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),
      producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      scope_plan_sha256=hashlib.sha256((scope/'full_T_scope_before_sender_effect_interpretation.json').read_bytes()).hexdigest(),
      groups=groups,candidates='all51 existing signatures, same flags and exact observed signature dedup',
      min_cells=20,max_cells=50,min_primary_patients=20,
      family='all observed T groups *all admissible unique physical signatures including unestimable internally P1',
      preserve_original_CD8_cell_draw_and_result=True,no_new_clinical_outcome_test=True,
      T_group_labels_are_candidate_annotations_not_accepted_states=True),indent=2)+'\n')
    g,rows,cells,var,eligible,measured=hnscc_inputs(root,project)
    lookup={}
    for i,(name,row) in enumerate(var.iterrows()):
        if not measured[i] or str(name).startswith('AMBIGUOUS::'):continue
        names={str(name),str(row.get('gene_symbol',''))}
        if row.get('feature_identity_status')=='stable_Ensembl_ID':names.update(str(row.get('source_gene_aliases','')).split('|'))
        for alias in names:
            if alias not in ['','nan','None']:lookup.setdefault(alias,set()).add(i)
    excluded=set(np.flatnonzero(~measured|var.index.astype(str).str.startswith('AMBIGUOUS::')))
    receiver_gene=pd.read_csv(root/'HNSCC_frozen_MCP_gene_axis_20261014_v2/Macnode18_private_frozen_gene_columns.tsv',sep='\t')
    receiver_columns=set(receiver_gene.loc[receiver_gene.native_column.ge(0),'native_column'].astype(int))
    references={r['programme_id']:r['primary_genes'] for r in programmes}
    old=root/'HNSCC_frozen_T_sender_states_20261014_v1'
    contract=pd.read_csv(root/'HNSCC_frozen_T_sender_associations_20261014_v1/private_signature_dedup_contract.tsv',sep='\t')
    major=cells[cells.timepoint.eq('pre')].groupby(['patient_id','major']).size()
    base_md=pd.read_csv(old/'private_patient_receiver_covariates.tsv',sep='\t',index_col=0)
    summaries=[];counts=[]
    for group in groups:
        dest=out/group;dest.mkdir()
        chosen=cells[cells.patient_id.isin(receiver.index)&cells.timepoint.eq('pre')&cells.major.eq('T')&cells.subtype.eq(group)&cells.evidence.eq('multigene_subtype_candidate')]
        number=chosen.groupby('patient_id').size();patients=sorted(number[number>=20].index)
        assert set(patients).issubset(receiver.index)
        if group=='CD8_Effector':
            assert patients==list(base_md.index)
            md=base_md.copy();primary=pd.read_csv(old/'private_patient_T_state_means.tsv',sep='\t',index_col=0).loc[patients]
            orth=pd.read_csv(old/'private_patient_T_orthogonal_state_means.tsv',sep='\t',index_col=0).loc[patients]
            counts.append(dict(group=group,patients=len(patients),source='exact cached previous draw'))
        else:
            rng=np.random.default_rng(seed('HNSCC_all_T_group_'+group));positions=[]
            for patient,frame in chosen[chosen.patient_id.isin(patients)].groupby('patient_id',sort=True):
                ix=cells.index.get_indexer(frame.index);positions.extend(rng.choice(ix,min(50,len(ix)),replace=False))
            positions=np.sort(np.asarray(positions,int));cm=cells.iloc[positions]
            raw=ad.io.sparse_dataset(g['layers/counts'])[rows[positions],:].tocsr()
            total=np.asarray(raw.sum(1)).ravel();detected=np.asarray((raw>0).sum(1)).ravel()
            assert (total>0).all() and np.isfinite(raw.data).all() and (raw.data>=0).all()
            values,names,cov=score_functions(raw,lookup,excluded,references)
            ov,onames,ocov=score_functions(raw,lookup,excluded|receiver_columns,references)
            p=pd.DataFrame(values,columns=names);p['patient_id']=cm.patient_id.to_numpy();primary=p.groupby('patient_id').mean().loc[patients]
            p=pd.DataFrame(ov,columns=onames);p['patient_id']=cm.patient_id.to_numpy();orth=p.groupby('patient_id').mean().loc[patients]
            q=pd.DataFrame(dict(patient_id=cm.patient_id.to_numpy(),T_logRNA=np.log1p(total),T_loggenes=np.log1p(detected))).groupby('patient_id').mean().loc[patients]
            md=base_md.loc[patients].copy();md[['T_logRNA','T_loggenes']]=q
            md['T_fraction']=[number.loc[p]/major.loc[p,'T'] for p in patients];md['T_logcells']=np.log1p(number.loc[patients])
            pd.DataFrame(cov).to_csv(dest/'private_signature_coverage.tsv',sep='\t',index=False)
            pd.DataFrame(ocov).to_csv(dest/'private_orthogonal_coverage.tsv',sep='\t',index=False)
            counts.append(dict(group=group,patients=len(patients),source='fixed subtype draw',selected_cells=len(cm),scalar_UCell_parity=True))
        assert md.notna().all().all() and set(primary.index)==set(md.index)==set(orth.index)
        primary.to_csv(dest/'private_patient_T_state_means.tsv',sep='\t');orth.to_csv(dest/'private_patient_T_orthogonal_state_means.tsv',sep='\t')
        md.to_csv(dest/'private_patient_receiver_covariates.tsv',sep='\t')
        for record in programmes:
            row=contract.set_index('programme_id').loc[record['programme_id']].to_dict()
            row.update(T_candidate_group=group,programme_id=record['programme_id'],measurable=record['programme_id'] in primary,orthogonal_measurable=record['programme_id'] in orth)
            summaries.append(row)
    pd.DataFrame(summaries).to_csv(out/'private_complete_group_candidate_contract.tsv',sep='\t',index=False)
    report=dict(status='COMPLETE_ALL_T_CANDIDATE_GROUP_FIXED_STATE_INPUTS_NO_NEW_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
      groups=counts,programme_ids_per_group=51,clinical_outcomes_tested=False,CD8_original_preserved=True)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
