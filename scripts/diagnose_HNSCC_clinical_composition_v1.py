"""Post-result descriptive composition and leave-one-patient diagnosis.
All variants retained; no new P/BH, no replacement of the primary model.
"""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata


def beta(x,y,clinical,cov):
    x=rankdata(x);y=rankdata(y)
    if min(x.std(),y.std())<1e-12:return np.nan
    x=(x-x.mean())/x.std();y=(y-y.mean())/y.std()
    extras=np.apply_along_axis(rankdata,0,cov) if cov.shape[1] else cov
    X=np.column_stack((clinical,extras,x))
    if len(y)-X.shape[1]<10 or np.linalg.matrix_rank(X)!=X.shape[1]:return np.nan
    return float(np.linalg.lstsq(X,y,rcond=None)[0][-1])


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    inputs=root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1'
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    out=root/'HNSCC_clinical_descriptive_composition_20261014_v1';assert not out.exists();out.mkdir()
    (out/'plan_before_diagnostics.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),
      scope='all3 supported baseline endpoints, all variants, leave1patient',
      variants=['clinical_only','clinical_plus_both_candidate_fractions_and_counts','clinical_plus_all_macro_subtype_fractions_except_lexical_reference','clinical_plus_joint_capture_and_subtypes','clinical_only_author_concordant_subset'],
      min_residual_df=10,no_new_P_or_BH=True,
      limitations='post-result descriptive diagnosis; composition adjustment can condition away a biological mediator; RNA quality still frozen training correction only'),indent=2)+'\n')
    md=pd.read_csv(inputs/'private_eligible_libraries_clinical.tsv',sep='\t')
    baseline=md[md.timepoint.eq('pre')].set_index('patient_id')
    scores=pd.read_csv(inputs/'private_patient_library_frozen_scores.tsv',sep='\t')
    s=scores[scores.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='score')
    z=scores[scores.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='CRC_training_z')
    s['equal_mean_CRC_training_z']=z.mean(axis=1);baseline=baseline.loc[s.index]
    chunks=[]
    for chunk in pd.read_csv(atlas/'tables/all_core_annotation_private.tsv.gz',sep='\t',usecols=['dataset','library_id','major','subcluster'],dtype=str,keep_default_na=False,chunksize=200000):
        take=chunk[chunk.dataset.eq('GSE286827')&chunk.library_id.isin(baseline.library_id)]
        if len(take):chunks.append(take)
    cells=pd.concat(chunks,ignore_index=True)
    provisional=pd.read_csv(atlas/'tables/subcluster_catalog_private.tsv',sep='\t')
    refined=pd.read_csv(atlas/'tables/subcluster_catalog_refined_private.tsv',sep='\t')
    bridge=provisional[['subcluster','node_id']].merge(refined[['node_id','subtype_candidate','subtype_evidence_status']],on='node_id',validate='one_to_one').set_index('subcluster')
    cells['subtype']=cells.subcluster.map(bridge.subtype_candidate)
    cells['candidate_status']=cells.subcluster.map(bridge.subtype_evidence_status)
    assert cells.candidate_status.notna().all()
    major=cells.groupby(['library_id','major']).size()
    macro=cells[cells.major.eq('Myeloid')&cells.subtype.str.startswith('Macro_')&cells.candidate_status.eq('multigene_subtype_candidate')]
    table=macro.groupby(['library_id','subtype']).size().unstack(fill_value=0).sort_index(axis=1)
    assert table.shape[1]>=2
    table=table.loc[baseline.library_id]
    assert np.array_equal(table.sum(axis=1).to_numpy(),baseline.Macrophage_candidate_cells.to_numpy())
    composition=table.div(table.sum(axis=1),axis=0).iloc[:,1:].to_numpy()
    capture=np.array([[row.CD8_candidate_cells/major.loc[row.library_id,'T'],row.Macrophage_candidate_cells/major.loc[row.library_id,'Myeloid'],np.log1p(row.CD8_candidate_cells),np.log1p(row.Macrophage_candidate_cells)] for row in baseline.itertuples()])
    clinical=np.column_stack((np.ones(len(baseline)),baseline.treatment_x.eq('D+T').to_numpy(float),baseline.HPV.eq('positive').to_numpy(float)))
    y=baseline.regression.to_numpy(float)
    variants=[('clinical_only',np.empty((len(y),0)),np.ones(len(y),bool)),
      ('clinical_plus_both_candidate_fractions_and_counts',capture,np.ones(len(y),bool)),
      ('clinical_plus_all_macro_subtype_fractions_except_lexical_reference',composition,np.ones(len(y),bool)),
      ('clinical_plus_joint_capture_and_subtypes',np.column_stack((capture,composition)),np.ones(len(y),bool)),
      ('clinical_only_author_concordant_subset',np.empty((len(y),0)),baseline.author_binary_concordant.to_numpy(bool))]
    records=[];leave=[]
    for endpoint in s.columns:
        for name,cov,mask in variants:
            x=s.loc[mask,endpoint].to_numpy();yy=y[mask];cc=clinical[mask];cv=cov[mask]
            effect=beta(x,yy,cc,cv);lo=[]
            for i in range(len(x)):
                keep=np.arange(len(x))!=i;value=beta(x[keep],yy[keep],cc[keep],cv[keep])
                leave.append(dict(endpoint=endpoint,variant=name,left_out_ordinal=i+1,beta=value))
                if np.isfinite(value):lo.append(value)
            records.append(dict(endpoint=endpoint,variant=name,n_patients=len(x),covariates=cc.shape[1]+cv.shape[1]-1,
              beta=effect,status='descriptive_estimated' if np.isfinite(effect) else 'rank_or_df_insufficient',
              LOO_valid=len(lo),LOO_min=min(lo) if lo else np.nan,LOO_max=max(lo) if lo else np.nan))
    result=pd.DataFrame(records);assert len(result)==15
    result.to_csv(out/'all15_private_descriptive_diagnostics.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(leave).to_csv(out/'private_leave_one_patient.tsv',sep='\t',index=False,na_rep='NA')
    (out/'private_terminal.json').write_text(json.dumps(dict(status='ALL15_CLINICAL_COMPOSITION_DIAGNOSTICS_COMPLETE_NO_NEW_P',
      finished_UTC=datetime.now(timezone.utc).isoformat(),macro_subtype_dimensions=table.shape[1]),indent=2)+'\n')
    print(result.to_json(orient='records'))


if __name__=='__main__':main()
