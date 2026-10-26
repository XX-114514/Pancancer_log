"""All-archive unit census and explicit qualified/held/unknown analysis views.

Preserves frozen reviewed patient blocks, native units and source roles. New
clinical eligibility is never inferred from observed metadata flag names alone.
This reads metadata and accepted identity sidecars, not expression matrices.
"""
import argparse,collections,datetime,gc,hashlib,json,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
from annotation_representation_pilot_v1 import save
from atlas_analysis_view_eligibility_v1 import qualify_unit_views
from audit_full_source_local_graphs_v1 import index_at
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
from summarize_full_native_annotation_evidence_v1 import read_source

VIEWS=['all_archive','expression_reference','whole_tissue_composition','within_selected_lineage_composition','cell_product','treatment_clinical']


def main(root):
    began=time.monotonic();cfg=json.loads((root/'config/archive_unit_view_spec_v1.json').read_text())
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Unit view input binding changed')
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_atlas_analysis_view_eligibility_v1.py')],capture_output=True,text=True,timeout=120);(root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 4 tests' not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Actual source/unit view counterexamples failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':4,'no_skips':True})
    candidates,roles,extension=map(Path,[cfg['candidates'],cfg['source_roles'],cfg['design_extension']]);ci=index_at(candidates,'full_native_cell_identity_candidate_index_v1.json');reader=NativeArchiveCollection(cfg['archive_manifest'])
    rt=json.loads((roles/'provenance/terminal_v1.json').read_text());rp=roles/'outputs/source_roles_and_coverage_v1.json'
    et=json.loads((extension/'provenance/terminal_v1.json').read_text());ep=extension/'outputs/source_design_qualification_extension_v1.json'
    if sha256_file(rp)!=rt['output_sha256'] or sha256_file(ep)!=et['output_sha256'] or et['original_source_role_output_sha256']!=sha256_file(rp):raise ValueError('Accepted source review output changed')
    rr={r['source']:r for r in json.loads(rp.read_text())['sources']};er={r['source']:r for r in json.loads(ep.read_text())['sources']}
    if set(rr)!=set(ci) or set(er)!=set(ci) or set(ci)!=set(reader.archive_partitions):raise ValueError('All archive unit/source scope differs')
    meta=pd.read_csv(reader.legacy_config['sample_metadata_path'],sep='\t',dtype=str,keep_default_na=False)
    if not meta.sampleID.is_unique:raise ValueError('Frozen unit metadata contains duplicate sample keys')
    inherited=meta.set_index('sampleID').to_dict('index');results=[];summaries=[];overall=collections.Counter();totals={v:collections.Counter() for v in VIEWS}
    for ds in sorted(ci):
        storage_gate(root,cfg);cm,ca=read_source(candidates,ci[ds]);rows=reader.rows(ds);raw_units=reader.frozen_expression_units(ds).astype(str);units,code,number=np.unique(raw_units,return_inverse=True,return_counts=True);fit=ca['fit_pool_candidate'];nfit=np.bincount(code[fit],minlength=len(units));np.testing.assert_array_equal(ca['native_row'],rows)
        if len(rows)!=rr[ds]['n_archive_cells'] or len(units)!=rr[ds]['n_archive_expression_unit_keys']:raise ValueError('Frozen source/unit denominator changed')
        capture=rr[ds]['design_decision']['capture_class'];unit_records=[];blocks=[];view_code=np.zeros((len(units),len(VIEWS)),np.uint8)
        for j,unit in enumerate(units):
            old=inherited.get(str(unit)) if reader.is_legacy(ds) else None
            if old is not None and old['gse_id']!=ds:raise ValueError('Frozen unit/source join is not exact')
            usable=old is not None or str(unit).strip().lower() not in {'','nan','none','na','n/a','unknown','<na>','unreviewed','not_available'}
            reviewed=old is not None and old['patient_review_status']=='reviewed_identity' and old['unit_kind']=='reviewed_patient' and bool(old['analysis_unit'])
            block_kind='inherited_reviewed_patient' if reviewed else 'native_expression_unit_fallback' if usable else 'unresolved_unit_key_archive_bucket'
            block_key=old['analysis_unit'] if reviewed else str(unit);blocks.append(json.dumps([block_kind,block_key],ensure_ascii=False))
            # Accepted fields here have narrow technical meanings: exact local
            # native source/asset binding, measured count-domain evidence and
            # native expression-unit identity. They do not certify paper asset
            # versions, capture census, prior QC, phase or endpoints.
            review={'archive_membership':'accepted','source_asset_join':'accepted','native_measurement':'accepted','expression_unit_join':'accepted','capture_class':capture,'prior_QC_scope':'unknown','assay_scope':'unknown','annotation_release':'unknown','source_alias_scope':'unknown','material_and_anatomy':'unknown','capture_design':'unknown','unselected_tissue_sampling':'unknown','lineage_selection_scope':'unknown','within_lineage_capture_comparability':'unknown','product_design':'unknown','patient_identity':'accepted' if reviewed else 'unknown','patient_sample_library_join':'unknown','treatment_phase':'unknown','clinical_endpoint':'unknown','endpoint_time_origin':'unknown','repeated_sampling_plan':'unknown','clinical_cohort_scope':'unknown'}
            if not usable:review['expression_unit_join']='unknown'
            restriction=er[ds]['current_conservative_analysis_roles']['whole_tissue_composition']
            if restriction.startswith(('not_eligible','do_not_use')):review['unselected_tissue_sampling']='rejected'
            qualified=qualify_unit_views(review)
            for k,view in enumerate(VIEWS):
                status=qualified[view]['status'];view_code[j,k]={'unknown':0,'accepted':1,'rejected':2}[status];totals[view][status+'_units']+=1;totals[view][status+'_archive_rows']+=int(number[j])
            unit_records.append({'source':ds,'native_expression_unit':str(unit),'native_unit_key_sha256':hashlib.sha256(json.dumps([ds,str(unit)],ensure_ascii=False).encode()).hexdigest(),'native_unit_key_usable_or_frozen_reviewed':usable,'n_archive_rows':int(number[j]),'n_fit_pool_candidate_rows':int(nfit[j]),'inherited_frozen_unit_review':None if old is None else {f:old[f] for f in ['patient_review_status','analysis_unit','unit_kind','sorting_enrichment_status','ecology_eligibility','lineage_eligibility','eligible_min100_cells','sample_class','cancerType']},'inherited_reviewed_patient_block_preserved':reviewed,'analysis_block_kind':block_kind,'analysis_block_original_key':block_key,'technical_review':review,'technical_source_asset_scope':'existing_local_native_source_axis_and_hash_binding_not_new_public_paper_asset_or_patient_join','qualified_views':{v:qualified[v] for v in VIEWS},'current_source_roles':er[ds]['current_conservative_analysis_roles'],'native_flag_values_not_automatically_clinical_qualification':True,'source_alias_and_cross_source_patient_independence_not_accepted':True})
        block_names,unit_block=np.unique(np.asarray(blocks),return_inverse=True);block_dictionary=[dict(zip(['kind','original_key'],json.loads(b))) for b in block_names]
        arrays={'native_row':rows,'native_expression_unit_code':code.astype(np.int32),'analysis_block_code':unit_block[code].astype(np.int32),'fit_pool_candidate':fit,'unit_view_status_code':view_code,'expression_unit_dictionary':units}
        # Unit-level arrays have their own axis, distinct from per-cell arrays.
        output=root/'outputs'/(ds+'_archive_unit_views_v1.npz')
        with output.open('xb') as handle:np.savez_compressed(handle,**arrays)
        with np.load(output,allow_pickle=False) as f:
            for key,value in arrays.items():np.testing.assert_array_equal(f[key],value)
        mp=output.with_suffix('.json');save(mp,{'source':ds,'n_archive_rows':len(rows),'native_ids_sha256':cm['native_ids_sha256'],'native_unit_axis_sha256':text_axis_sha256(raw_units),'n_expression_units_not_patients':len(units),'unit_records':unit_records,'analysis_block_dictionary':block_dictionary,'view_order':VIEWS,'view_status_decoder':{'0':'unknown','1':'accepted','2':'rejected'},'unit_axis_fields':['unit_view_status_code','expression_unit_dictionary'],'all_other_array_fields_have_native_archive_row_axis':True,'science_view_cell_selection_additionally_requires_fit_pool_and_final_identity_gate':True,'original_source_role_record':rr[ds],'source_design_extension_record':er[ds],'new_patient_identity_clinical_eligibility_or_independence_assigned':False})
        reviewed_count=sum(r['inherited_reviewed_patient_block_preserved'] for r in unit_records);summary={'source':ds,'n_archive_rows':len(rows),'n_units':len(units),'n_fit_rows':int(fit.sum()),'n_inherited_reviewed_patient_units':reviewed_count,'n_inherited_reviewed_patient_blocks_not_independent_across_sources':sum(r['kind']=='inherited_reviewed_patient' for r in block_dictionary),'n_frozen_unit_metadata_matches':sum(r['inherited_frozen_unit_review'] is not None for r in unit_records)};summaries.append(summary);overall.update({k:v for k,v in summary.items() if k!='source'})
        results.append({'source':ds,'n_archive_rows':len(rows),'n_units':len(units),'arrays':output.name,'arrays_sha256':sha256_file(output),'metadata':mp.name,'metadata_sha256':sha256_file(mp)})
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_complete_sources':len(results),**summary});reader._fields.clear();del ca,raw_units,arrays,unit_records;gc.collect()
    if overall['n_archive_rows']!=sum(r['n_archive_rows'] for r in ci.values()):raise ValueError('Unit view changed whole archive denominator')
    ip=root/'outputs/archive_unit_view_index_v1.json';save(ip,results);summary=root/'outputs/archive_unit_view_coverage_v1.json';save(summary,{'source_summaries':summaries,'totals':dict(overall),'views':{k:dict(v) for k,v in totals.items()},'frozen_unit_metadata_rows':len(meta),'legacy_reviewed_status_units_in_frozen_table':int(meta.patient_review_status.eq('reviewed_identity').sum()),'existing_reviewed_blocks_preserved_not_new_patient_count':True,'unknown_scientific_gates_not_filled_by_technical_success':True})
    save(root/'provenance/terminal_v1.json',{'status':'ALL_ARCHIVE_UNIT_CENSUS_INHERITED_REVIEWED_BLOCKS_AND_EXPLICIT_VIEW_GATES_MATERIALIZED_INDEPENDENT_AUDIT_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(results),'n_archive_rows':overall['n_archive_rows'],'n_expression_units_not_patients':overall['n_units'],'index_sha256':sha256_file(ip),'summary_sha256':sha256_file(summary),'new_clinical_view_eligibility_assigned':False,'native_expression_read':False,'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
