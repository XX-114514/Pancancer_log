#!/usr/bin/env python3
"""Export an audited, aggregate preparation index; never perform integration."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


CAR='CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea'


def build(audit_path, history_path):
    audit=json.loads(audit_path.read_text())
    if not audit['all_checks_pass']:
        raise ValueError('Independent count/identity audit has not passed')
    history=json.loads(history_path.read_text())
    context={r['dataset']:r for r in history['records']}
    candidates={}
    for row in audit['integration_candidate_index']:
        candidates.setdefault(row['dataset'],{})[row['library_token']]=row['n_cells']
    records=[]
    for dataset in audit['datasets']:
        name=dataset['dataset']; prior=context.get(name)
        branch='manufactured_CAR_T_product' if name==CAR else 'T_enriched_nuclei' if name=='GSE292046' else 'tissue_or_comparator_context_review'
        libraries=[]
        for library in dataset['libraries']:
            token=library['library_token']
            if not library['all_checks_pass']:
                raise ValueError('Unaudited shard cannot enter preparation manifest')
            strict=candidates.get(name,{}).get(token,0)
            selected='counts' if prior else 'counts_decontx' if strict else 'counts'
            gates=['baseline_study_alias_and_library_overlap_pending','user_integration_scope_pending']
            if prior:gates.append('author_correction_count_slot_lineage_pending')
            else:gates.append('author_processing_history_and_ambient_model_review_pending')
            if name==CAR:gates.append('author_DP_state_doublet_retention_sensitivity_pending')
            if not strict:gates.append('strict_technical_eligibility_not_established')
            libraries.append(dict(library_token=token,n_prepared_cells=library['n_output'],
                n_strict_technical_candidate_cells=strict,planned_analysis_count_layer=selected,
                prepared_shard_logical_ref='objects/'+name+'/'+token+'/prepared.zarr',
                strict_technical_mask_role='diagnostic_candidate_not_final_cell_selection',
                prepared_normalized_X_may_use_different_layer=True,
                planned_layer_is_final=False,normalization_must_be_recomputed_from_selected_layer=True,
                counts_semantics='immutable_published_upstream_counts_not_proof_of_uncorrected_droplet_counts',
                release_gates=gates,integration_allowed=False))
        records.append(dict(dataset=name,branch=branch,n_prepared_cells=dataset['n_output'],
            author_history_reviewed_from_primary=bool(prior),
            author_ambient_method=prior['author_ambient_method'] if prior else 'pending',
            subset_scope=prior.get('published_subset_scope') if prior else 'context_review_pending',libraries=libraries))
    return dict(schema='public_preintegration_preparation_manifest_v1',observed_utc=datetime.now(timezone.utc).isoformat(),
        audit_observed_utc=audit['observed_utc'],audit_sha256=hashlib.sha256(audit_path.read_bytes()).hexdigest(),
        author_history_sha256=hashlib.sha256(history_path.read_bytes()).hexdigest(),
        n_sources=len(records),n_libraries=sum(len(r['libraries']) for r in records),
        n_audited_prepared_cells=sum(r['n_prepared_cells'] for r in records),
        scope='completed_audited_sources_only_not_full_queue_acceptance',
        source_final_AC=False,final_accepted_new_cells=0,cross_study_integrated=False,records=records)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--audit',required=True);p.add_argument('--author-history',required=True)
    p.add_argument('--output',required=True);args=p.parse_args()
    report=build(Path(args.audit),Path(args.author_history))
    with open(args.output,'x') as f:json.dump(report,f,indent=2,allow_nan=False)
    print(json.dumps({k:v for k,v in report.items() if k!='records'}))
