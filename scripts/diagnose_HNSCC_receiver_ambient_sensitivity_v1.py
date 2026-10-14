"""Matched raw/adopted ambient-removal sensitivity of fixed receiver scores.
No refit, new clinical P or reinterpretation of non-adopted correction as done.
"""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from diagnose_frozen_MCP_two_end_functions_v1 import hnscc_inputs
from diagnose_HNSCC_clinical_composition_v1 import beta


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    out=root/'HNSCC_receiver_adopted_ambient_sensitivity_20261014_v1';assert not out.exists();out.mkdir()
    (out/'plan_before_correction_diagnostics.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),
      scope='all previously selected baseline macrophage subtype samples, exact fixed cells',
      models=['raw_original','removed_counts_with_raw_denominator_and_quality','corrected_counts_corrected_denominator_and_quality'],
      correction='only already adopted ambient_removed_counts, not all estimated contamination',
      weights_and_training_normalization='unchanged CRC weights/means/SD/quality coefficients, >=95% measured weight retained',
      no_new_clinical_P_or_BH=True,clinical_variants=['arm_HPV','arm_HPV_age_sex'],
      interpretation='measurement sensitivity; changing count correction changes input domain, not a refitted programme'),indent=2)+'\n')
    g,source_rows,cells,var,eligible,measured=hnscc_inputs(root,project)
    selected=pd.read_csv(root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1/private_selected_cells.tsv',sep='\t')
    local=cells.index.get_indexer(selected.cell);assert (local>=0).all()
    md=cells.iloc[local];raw=ad.io.sparse_dataset(g['layers/counts'])[source_rows[local],:].tocsr().astype(float)
    assert 'ambient_removed_counts' in g['layers']
    removed=ad.io.sparse_dataset(g['layers/ambient_removed_counts'])[source_rows[local],:].tocsr().astype(float)
    corrected=(raw-removed).tocsr();assert (corrected.data>=-1e-7).all()
    corrected.data=np.maximum(corrected.data,0);corrected.eliminate_zeros()
    total=np.asarray(raw.sum(1)).ravel();ctotal=np.asarray(corrected.sum(1)).ravel()
    assert (ctotal>0).all() and (removed.data>=0).all() and (ctotal<=total+1e-7).all()
    detected=np.asarray((raw>0).sum(1)).ravel();cdetected=np.asarray((corrected>0).sum(1)).ravel()
    w=pd.read_csv(root/'HNSCC_frozen_MCP_gene_axis_20261014_v2/Macnode18_private_frozen_gene_columns.tsv',sep='\t')
    valid=w.native_column.ge(0);indices=w.loc[valid,'native_column'].to_numpy(int)
    assert np.abs(w.loc[valid,'weight']).sum()/np.abs(w.weight).sum()>=.95
    crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1'
    q=np.load(crc/'private_inputs/Macnode18_frozen_projection.npz')
    qc=pd.read_csv(crc/'private_transfer_fixed_weights_v1/Macnode18_MCP1_quality.tsv',sep='\t').set_index('term').value
    clinical=pd.read_csv(root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_eligible_libraries_clinical.tsv',sep='\t')
    clinical=clinical[clinical.timepoint.eq('pre')].set_index('patient_id')
    author=pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',sep='\t')
    pre=author[author.time.eq('Pre')].set_index('patient_registry');clinical['age']=pre.Age.loc[clinical.index]
    clinical['male']=pre.Gender.str.strip().str.lower().map({'m':1,'male':1,'f':0,'female':0}).loc[clinical.index]
    old=pd.read_csv(root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1/private_patient_subtype_scores.tsv',sep='\t')
    scores=[]
    for model,matrix,denom,genes in [('raw_original',raw,total,detected),('removed_counts_with_raw_denominator_and_quality',corrected,total,detected),('corrected_counts_corrected_denominator_and_quality',corrected,ctotal,cdetected)]:
        values=np.log1p(matrix[:,indices].toarray()*10000/denom[:,None])
        z=(values-w.loc[valid,'mean'].to_numpy())/w.loc[valid,'sd'].to_numpy()
        quality=(np.column_stack((np.log1p(denom),np.log1p(genes)))-q['quality_mean'])/q['quality_sd']
        programme=z@w.loc[valid,'weight'].to_numpy()-qc['(Intercept)']-quality[:,0]*qc.cellQ-quality[:,1]*qc.geneQ
        assert np.isfinite(programme).all()
        frame=pd.DataFrame(dict(patient_id=md.patient_id.to_numpy(),subtype_candidate=selected.subtype_candidate.to_numpy(),score=programme,
          removed_RNA_fraction=(total-ctotal)/total,any_adopted_removal=(total-ctotal)>1e-7))
        means=frame.groupby(['patient_id','subtype_candidate']).mean().reset_index();means['model']=model;scores.extend(means.to_dict('records'))
    scores=pd.DataFrame(scores);summary=[]
    for subtype in sorted(scores.subtype_candidate.unique()):
        group=scores[scores.subtype_candidate.eq(subtype)]
        baseline=group[group.model.eq('raw_original')].set_index('patient_id').sort_index()
        expected=old[old.subtype_candidate.eq(subtype)].set_index('patient_id').loc[baseline.index]
        assert np.allclose(baseline.score,expected.score,rtol=1e-9,atol=1e-9)
        cm=clinical.loc[baseline.index];C=np.column_stack((np.ones(len(cm)),cm.treatment_x.eq('D+T').to_numpy(float),cm.HPV.eq('positive').to_numpy(float)))
        y=cm.regression.to_numpy(float)
        for model,frame in group.groupby('model'):
            f=frame.set_index('patient_id').loc[baseline.index];x=f.score.to_numpy()
            for variant,N,cov in [('arm_HPV',C,np.empty((len(cm),0))),('arm_HPV_age_sex',np.column_stack((C,cm.male)),cm.age.to_numpy(float)[:,None])]:
                point=beta(x,y,N,cov);loo=[]
                for i in range(len(x)):
                    keep=np.arange(len(x))!=i;value=beta(x[keep],y[keep],N[keep],cov[keep])
                    if np.isfinite(value):loo.append(value)
                summary.append(dict(subtype_candidate=subtype,model=model,variant=variant,n_patients=len(cm),beta=point,
                  raw_corrected_patient_score_rho=float(spearmanr(baseline.score,f.score).statistic),
                  mean_selected_cell_adopted_removal_fraction=float(f.any_adopted_removal.mean()),mean_removed_RNA_fraction=float(f.removed_RNA_fraction.mean()),
                  LOO_valid=len(loo),LOO_min=min(loo) if loo else np.nan,LOO_max=max(loo) if loo else np.nan,
                  status='descriptive_estimated' if np.isfinite(point) else 'rank_or_df_insufficient'))
    assert len(summary)==18
    scores.to_csv(out/'private_matched_patient_subtype_scores.tsv',sep='\t',index=False)
    pd.DataFrame(summary).to_csv(out/'all18_private_adopted_ambient_diagnostics.tsv',sep='\t',index=False,na_rep='NA')
    terminal=dict(status='MATCHED_ADOPTED_AMBIENT_RECEIVER_SENSITIVITY_COMPLETE_NO_NEW_P',finished_UTC=datetime.now(timezone.utc).isoformat(),
      diagnostics=summary,raw_score_replay=True,estimated_contamination_not_assumed_adopted=True)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
