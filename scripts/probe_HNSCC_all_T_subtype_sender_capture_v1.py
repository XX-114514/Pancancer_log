"""Audit all existing T subtype candidates before viewing sender effects."""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    receiver=pd.read_csv(root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1/private_patient_subtype_scores.tsv',sep='\t')
    receiver=receiver[receiver.subtype_candidate.eq('Macro_C1QC')];patients=set(receiver.patient_id)
    chunks=[]
    for chunk in pd.read_csv(atlas/'tables/all_core_annotation_private.tsv.gz',sep='\t',usecols=['dataset','patient_id','timepoint','major','subcluster','library_id'],dtype=str,keep_default_na=False,chunksize=200000):
        take=chunk[chunk.dataset.eq('GSE286827')&chunk.patient_id.isin(patients)&chunk.timepoint.eq('pre')&chunk.major.eq('T')]
        if len(take):chunks.append(take)
    cells=pd.concat(chunks,ignore_index=True)
    provisional=pd.read_csv(atlas/'tables/subcluster_catalog_private.tsv',sep='\t')
    refined=pd.read_csv(atlas/'tables/subcluster_catalog_refined_private.tsv',sep='\t')
    bridge=provisional[['subcluster','node_id']].merge(refined[['node_id','subtype_candidate','subtype_evidence_status']],on='node_id',validate='one_to_one').set_index('subcluster')
    cells['candidate']=cells.subcluster.map(bridge.subtype_candidate);cells['evidence']=cells.subcluster.map(bridge.subtype_evidence_status)
    counts=cells[cells.evidence.eq('multigene_subtype_candidate')].groupby(['candidate','patient_id']).size()
    rows=[]
    for candidate in sorted(counts.index.get_level_values('candidate').unique()):
        g=counts.xs(candidate,level='candidate').reindex(sorted(patients),fill_value=0)
        rows.append(dict(T_subtype_candidate=candidate,receiver_baseline_patients=len(patients),patients_with_any_T_subtype=int(g.gt(0).sum()),
          patients_both_endpoints_atleast20=int(g.ge(20).sum()),primary20_patient_capture_gate=bool(g.ge(20).sum()>=20)))
    out=root/'HNSCC_all_T_subtype_sender_capture_20261014_v1';assert not out.exists();out.mkdir()
    counts.to_csv(out/'private_patient_T_subtype_counts.tsv',sep='\t')
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),scope='all existing multigene T candidate groups to same frozen C1QC receiver',
      candidate_programmes='all51 frozen T signatures; same stability/foreignRNA/basic flags, physical signature dedup and RNA/capture nuisance contract',
      family='all observed T groups times all admissible unique signatures, including insufficient capture as internally P1',
      no_current_sender_effects_viewed_for_expansion=True,clinical_target_selected_after_prior_clinical_results=True,
      keep_CD8_screen_as_original_nested_round=True,new_family_not_independent_clinical_confirmation=True)
    (out/'full_T_scope_before_sender_effect_interpretation.json').write_text(json.dumps(plan,indent=2)+'\n')
    report=dict(status='ALL_T_CANDIDATE_GROUP_CAPTURE_AUDITED_SCOPE_FROZEN_NO_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),capture=rows)
    (out/'aggregate_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
