"""Outcome-free exact reviewed-patient/pre-post RNA capture qualification."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    author=root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    capture=root/'immune_CM_RNA_supported_multicell_views_preflight_20261018_v1/private_all_original_libraries_multicell_view_counts.tsv'
    original=root/'immune_RNA_MCP_discovery_scope_20261018_v1/private_metadata_only_case_source_partition.tsv'
    assert sha(author)=='9942706870518c20c87a45fed44caff6bec22a1a512b6f194dc17d406df5cd22'
    assert sha(capture)=='b4d7af78d6346b55cd0d383c5d786664d0f5386b7a065912c04b916df10e5e94'
    assert sha(capture.parent/'private_independent_count_readiness_acceptance_20261018_v2.json')=='3f1261b9083f1cefcf0e072daecdef7c3db5d1fbf7f44b4ac1ff237efb126c15'
    scope=root/'immune_RNA_MCP_discovery_scope_20261018_v1';assert sha(scope/'private_terminal.json')=='0f4833a07bd5cdf5e8f6806d06805777dddf78adad60a6151eca991dfdda6196'
    assert sha(original)==json.loads((scope/'private_terminal.json').read_text())['files'][original.name]['sha256']
    out=root/'HNSCC_immune_RNA_MCP_paired_materials_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_any_paired_expression_or_effect.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        scope='All reviewed HNSCC author patients, exactpre/post library keys and same frozen rawRNA view counts; omit recurrence from pre/post comparison. No response or expression score fields read.',
        holdout='Exclude any pair whose post library was used in originalRNA discovery. Baseline may have contributed training; not independent patient or external gene-programme validation.',
        question='Feasibility of future paired within-patient programme dynamics, distinguishing stable person-associated RNA from time-associated coordinated changes. No new P family, score fitting or cell identity change.'),indent=2)+'\n')
    a=pd.read_csv(author,sep='\t',keep_default_na=False,usecols=['library_id','patient_registry','timepoint_registry'])
    assert len(a)==57 and a.library_id.is_unique
    a['sampleID']='GSE286827::'+a.library_id.astype(str)
    c=pd.read_csv(capture,sep='\t',keep_default_na=False);c=c[c.threshold.eq(1)].set_index('sampleID')
    d=pd.read_csv(original,sep='\t',keep_default_na=False);discovery=set(d.loc[d.partition.eq('reviewed_case_discovery'),'sampleID'])
    fields={'R_FOXP3_capture':'FOXP3_supported_R_only','E_effector_capture':'effector_only','Myeloid':'original_myeloid_cells'}
    rows=[]
    for patient,g in a.groupby('patient_registry',sort=True):
        pre=g[g.timepoint_registry.eq('pre')];post=g[g.timepoint_registry.eq('post')]
        assert len(pre)<=1 and len(post)<=1
        if len(pre)!=1 or len(post)!=1:continue
        prekey=pre.sampleID.iloc[0];postkey=post.sampleID.iloc[0];assert prekey in c.index and postkey in c.index
        row=dict(patient_registry=patient,pre_sampleID=prekey,post_sampleID=postkey,post_used_in_original_discovery=postkey in discovery,baseline_used_in_original_discovery=prekey in discovery)
        for view,field in fields.items():
            row[view+'_pre_cells']=int(c.loc[prekey,field]);row[view+'_post_cells']=int(c.loc[postkey,field])
            row[view+'_paired_capture_ge20']=min(row[view+'_pre_cells'],row[view+'_post_cells'])>=20
        rows.append(row)
    paired=pd.DataFrame(rows);assert paired.patient_registry.is_unique
    paired.to_csv(out/'private_all_exact_pre_post_patient_material_pairs.tsv',sep='\t',index=False)
    held=paired[~paired.post_used_in_original_discovery]
    summary=dict(dataset='GSE286827',reviewed_author_patients=int(a.patient_registry.nunique()),exact_pre_post_pairs=len(paired),post_training_excluded_pairs=int(paired.post_used_in_original_discovery.sum()),
        post_not_used_in_discovery_pairs=len(held),paired_view_coverage={v:int(held[v+'_paired_capture_ge20'].sum()) for v in fields},
        paired_pair_coverage={va+'|'+vb:int((held[va+'_paired_capture_ge20']&held[vb+'_paired_capture_ge20']).sum()) for va,vb in [('R_FOXP3_capture','E_effector_capture'),('R_FOXP3_capture','Myeloid'),('E_effector_capture','Myeloid')]},
        phase_comparison_not_causal_therapy_effect=True,programme_discovery_patient_overlap=True,no_new_response_or_score_effects=True)
    (out/'private_source_paired_material_qualification.json').write_text(json.dumps(summary,indent=2)+'\n')
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
    t=dict(status='COMPLETE_EXACT_REVIEWED_PATIENT_PRE_POST_CAPTURE_QUALIFICATION_NO_EFFECTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
