"""Independent full-point/LOO audit of adopted ambient-removal diagnostics."""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from audit_HNSCC_subtype_function_followup_v1 import effect,check_loo,same
from audit_HNSCC_fixed_MCP_clinical_v1 import ranks


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=root/'HNSCC_receiver_adopted_ambient_sensitivity_20261014_v1'
    scores=pd.read_csv(out/'private_matched_patient_subtype_scores.tsv',sep='\t')
    result=pd.read_csv(out/'all18_private_adopted_ambient_diagnostics.tsv',sep='\t');assert len(result)==18
    clinical=pd.read_csv(root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_eligible_libraries_clinical.tsv',sep='\t')
    clinical=clinical[clinical.timepoint.eq('pre')].set_index('patient_id')
    author=pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',sep='\t')
    pre=author[author.time.eq('Pre')].set_index('patient_registry')
    clinical['age']=pre.Age.loc[clinical.index]
    clinical['male']=pre.Gender.str.strip().str.lower().map({'m':1,'male':1,'f':0,'female':0}).loc[clinical.index]
    for row in result.itertuples():
        s=scores[scores.subtype_candidate.eq(row.subtype_candidate)&scores.model.eq(row.model)].set_index('patient_id').sort_index()
        raw=scores[scores.subtype_candidate.eq(row.subtype_candidate)&scores.model.eq('raw_original')].set_index('patient_id').loc[s.index]
        c=clinical.loc[s.index];N=np.column_stack((np.ones(len(c)),c.treatment_x.eq('D+T').to_numpy(float),c.HPV.eq('positive').to_numpy(float)))
        extra=c.age.to_numpy(float)[:,None] if row.variant.endswith('age_sex') else np.empty((len(c),0))
        if row.variant.endswith('age_sex'):N=np.column_stack((N,c.male.to_numpy(float)))
        check_loo(effect,s.score.to_numpy(),c.regression.to_numpy(float),[N,extra],row)
        same(np.corrcoef(ranks(raw.score.to_numpy()),ranks(s.score.to_numpy()))[0,1],row.raw_corrected_patient_score_rho)
        same(s.any_adopted_removal.mean(),row.mean_selected_cell_adopted_removal_fraction)
        same(s.removed_RNA_fraction.mean(),row.mean_removed_RNA_fraction)
    report=dict(status='PASS_ALL18_ADOPTED_AMBIENT_DESCRIPTIVE_EFFECTS_LOO_AND_SCORE_RANKS',finished_UTC=datetime.now(timezone.utc).isoformat(),
      no_new_clinical_P_or_BH=True,correction_extrapolated_to_nonadopted_cells=False)
    (out/'private_independent_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
