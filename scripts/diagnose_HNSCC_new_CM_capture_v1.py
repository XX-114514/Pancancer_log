"""All30 cohort capture and projection diagnostics; no new statistical tests."""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    inputs=root/'HNSCC_new_CM_baseline_preparation_20261016_v1'
    original=root/'HNSCC_new_CM_frozen_model_export_20261016_v2'
    out=root/'HNSCC_new_CM_capture_diagnostics_20261016_v1';assert not out.exists()
    scores=pd.read_csv(inputs/'private_patient_baseline_component_scores.tsv',sep='\t')
    projection=pd.read_csv(inputs/'private_patient_projection_diagnostics.tsv',sep='\t')
    catalogue=pd.read_csv(inputs/'private_all30_component_scope.tsv',sep='\t')
    rows=[];scope_rows=[]
    for scope in sorted(scores.scope.unique()):
        fq=pd.read_csv(atlas/'CM'/scope/'inputs/frequency_private.tsv',sep='\t',index_col=0)
        meta=pd.read_csv(atlas/'CM'/scope/'inputs/sample_metadata_private.tsv',sep='\t')
        selected=meta[meta.dataset.eq('GSE286827')&meta.timepoint.eq('pre')]
        capture=fq[selected.sampleID].gt(0).any(axis=1)
        W=pd.read_csv(original/scope/'W_private.tsv',sep='\t',index_col=0)
        frame=scores[scores.scope.eq(scope)].pivot(index='patient_id',columns='component',values='activity')
        rank_groups={}
        for component in frame:
            ranks=tuple(rankdata(frame[component]).tolist());rank_groups.setdefault(ranks,[]).append(component)
        for component in frame:
            x=frame[component].to_numpy();threshold=1e-10*max(1,float(np.max(abs(x))))
            others=next(names for ranks,names in rank_groups.items() if component in names)
            record=catalogue[catalogue.scope.eq(scope)&catalogue.component.eq(component)].iloc[0]
            row=dict(scope=scope,component=component,baseline_patients=len(x),positive_activity_patients=int((x>threshold).sum()),
                activity_IQR=float(np.quantile(x,.75)-np.quantile(x,.25)),activity_SD=float(np.std(x)),
                maximum_patient_share_of_activity=float(np.max(x)/np.sum(x)) if np.sum(x)>0 else np.nan,
                captured_raw_W_weight_fraction=float(W.loc[capture,component].sum()/W[component].sum()),
                raw_W_nodes_observed_in_any_baseline=int(capture.sum()),raw_W_nodes_total=len(W),
                identical_activity_rank_components='|'.join(sorted(others)),
                network_members=int(record.supported_network_members),weighted_major_composition=record.weighted_major_composition)
            rows.append(row)
        residual=projection[projection.scope.eq(scope)].relative_residual
        scope_rows.append(dict(scope=scope,patient_blocks=len(frame),components=len(frame.columns),
            distinct_within_cohort_activity_rank_patterns=len(rank_groups),
            relative_projection_residual_min=float(residual.min()),median=float(residual.median()),max=float(residual.max())))
    assert len(rows)==30
    out.mkdir();pd.DataFrame(rows).to_csv(out/'all30_private_capture_and_activity_diagnostics.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(scope_rows).to_csv(out/'private_scope_projection_diagnostics.tsv',sep='\t',index=False)
    report=dict(status='ALL30_COHORT_CAPTURE_AND_ACTIVITY_DIAGNOSTICS_COMPLETE_NO_NEW_P',finished_UTC=datetime.now(timezone.utc).isoformat(),
        scopes=scope_rows,components_with_fewer_than5_positive_patients=sum(r['positive_activity_patients']<5 for r in rows),
        post_result_exclusions_or_FDR_refits=False,
        interpretation='Sparse/duplicate within-cohort activities are measurement diagnostics, not a new negative test or independent modules. Original all30 primary family unchanged.')
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    print(pd.DataFrame(rows)[['scope','component','positive_activity_patients','maximum_patient_share_of_activity','captured_raw_W_weight_fraction','identical_activity_rank_components']].to_json(orient='records'))


if __name__=='__main__':main()
