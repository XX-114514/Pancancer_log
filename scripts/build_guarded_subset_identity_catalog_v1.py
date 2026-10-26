"""Materialize audited exact-group IDs, display markers and archive crosswalks.

The catalog preserves frozen guarded identity proposals, uncertainty and states.
It does not promote structural subgroup IDs to new biological cell types or
change source, clinical, malignancy or joint-representation qualification.
"""
import argparse,datetime,gc,json,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
from annotation_representation_pilot_v1 import save
from atlas_identity_release_keys_v1 import membership_key,assign_display_numbers,select_display_marker,display_label,version_crosswalk
from audit_full_source_local_graphs_v1 import index_at
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256,row_axis_sha256
from summarize_full_native_annotation_evidence_v1 import read_source


def main(root):
    began=time.monotonic();cfg=json.loads((root/'config/guarded_identity_catalog_spec_v1.json').read_text())
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Identity catalog bound input changed')
    marker,auditor,candidates=map(Path,[cfg['markers'],cfg['marker_auditor'],cfg['candidates']])
    for stage in [marker,auditor]:
        s=json.loads((stage/'provenance/supervisor_terminal_v1.json').read_text())
        if s['actual_child_exit_code'] or s['wall_timeout']:raise ValueError('Actual full marker production and audit required')
    audit=json.loads((auditor/'provenance/terminal_v1.json').read_text());mt=json.loads((marker/'provenance/terminal_v1.json').read_text())
    if audit['scope']!='full_completed_archive' or audit['producer_terminal_sha256']!=sha256_file(marker/'provenance/terminal_v1.json'):raise ValueError('Exact full marker audit provenance differs')
    ap=auditor/'outputs/guarded_subset_native_marker_audit_v1.json'
    if sha256_file(ap)!=audit['output_sha256']:raise ValueError('Full marker audit output changed')
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_atlas_identity_release_keys_v1.py')],capture_output=True,text=True,timeout=120);(root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 5 tests' not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Actual key/marker/crosswalk tests failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':5,'no_skips':True})
    mi=index_at(marker,'full_source_local_marker_evidence_index_v1.json');ci=index_at(candidates,'full_native_cell_identity_candidate_index_v1.json');reader=NativeArchiveCollection(cfg['archive_manifest'])
    if set(mi)!=set(ci) or set(ci)!=set(reader.archive_partitions):raise ValueError('Catalog all-source scope differs')
    global_catalog=[];source_catalog_records={}
    # First freeze exact group keys from source-namespaced native IDs. Neither
    # marker names nor candidate type strings enter these stable group keys.
    for ds in sorted(mi):
        storage_gate(root,cfg);record=mi[ds];p=marker/'outputs'/record['metadata']
        if sha256_file(p)!=record['metadata_sha256']:raise ValueError('Accepted own-marker metadata changed')
        meta=json.loads(p.read_text());cm,ca=read_source(candidates,ci[ds]);_,groups=read_source(marker,meta['source_graph_record']);ids=reader.native_ids(ds);rows=reader.rows(ds)
        if len(set(ids))!=len(ids):raise ValueError('Native IDs require reviewed library namespaces before catalog keys')
        np.testing.assert_array_equal(rows,ca['native_row']);np.testing.assert_array_equal(rows,groups['native_row']);np.testing.assert_array_equal(ca['fit_pool_candidate'],groups['fit_pool_candidate'])
        if text_axis_sha256(ids)!=meta['native_ids_sha256'] or meta['n_archive_rows']!=len(rows):raise ValueError('Catalog archive identity changed')
        support={(r['source_local_group_code'],r['gene']):r for r in meta['marker_unit_direction_summaries']};catalog=[]
        for desc,group in zip(meta['guarded_subsets'],meta['groups']):
            code=desc['guarded_subset_code'];positions=np.flatnonzero(groups['source_local_group_code']==code)
            if code!=group['source_local_group_code'] or len(positions)!=desc['n_subset_fit_rows'] or row_axis_sha256(rows[positions])!=desc['native_subset_rows_sha256']:raise ValueError('Audited exact catalog membership differs')
            # The top bound marker alone supplies the display token; all20
            # candidates remain in the original own-expression evidence file.
            key=membership_key(ds,ids[positions].astype(str));own=[{'gene':m['gene'],'marker_evidence_type':'own_subset_expression_specificity','membership_key':key,'review_status':'accepted','evidence_record_sha256':record['metadata_sha256']} for m in group['own_group_markers'][:1]]
            selected=select_display_marker(own,key,allow_parent=False)
            top=support.get((code,selected['marker']));units=top['n_comparable_expression_units'] if top else 0;positive=top['n_units_positive_mean_ratio'] if top else 0
            status='own_marker_multiple_expression_unit_direction_supported' if units>=2 and positive>units/2 else 'own_marker_expression_support_without_multiple_unit_majority'
            if len(positions)<30:status='small_group_own_expression_descriptive_only'
            if selected['marker'] is None:status='no_accepted_own_marker'
            catalog.append({'source':ds,'source_subset_code':code,'stable_group_sha256':key,'parent_source_local_group_code':desc['parent_source_local_group_code'],'stable_parent_group_sha256':desc['stable_parent_group_sha256'],'native_subset_rows_sha256':desc['native_subset_rows_sha256'],'n_subset_fit_rows':len(positions),'consensus_celltype_candidate':desc['guarded_celltype_candidate'],**selected,'group_expression_evidence_status':status,'n_comparable_expression_units_not_patients':units,'n_units_positive_mean_ratio':positive,'membership_identical_to_parent':desc['membership_identical_to_parent'],'marker_DEG_status':'no_qualified_exact_membership_inferential_DEG_input_supplied','new_biological_celltype_not_inferred_from_marker_or_number':True})
        if len(catalog)!=len(meta['groups']):raise ValueError('Catalog subset denominator differs')
        global_catalog.extend(catalog);source_catalog_records[ds]=catalog;reader._fields.clear();del ca,groups,ids,meta;gc.collect()
    keys=[r['stable_group_sha256'] for r in global_catalog]
    prior=json.loads(Path(cfg['previous_display_registry']).read_text()) if cfg.get('previous_display_registry') else {}
    numbers=assign_display_numbers(keys,prior)
    for record in global_catalog:
        record['global_subcluster_number']=numbers['current'][record['stable_group_sha256']];record['global_subcluster_id']='atlas_group_'+str(record['global_subcluster_number']);record['subcluster']=display_label(record['consensus_celltype_candidate'],record['marker'],record['global_subcluster_number']);record['label_version']=cfg['label_version']
    catalog_file=root/'outputs/guarded_subset_identity_catalog_v1.json';save(catalog_file,{'schema':'guarded_subset_identity_catalog_v1','groups':global_catalog,'fields_are_native_guarded_proposals_with_audited_structural_keys':True,'confidence_is_evidence_provenance_not_probability':True,'new_authoritative_identity':False,'global_IDs_are_exact_groups_not_automatically_cross_source_biological_types':True})
    registry=root/'outputs/global_display_number_registry_v1.json';save(registry,numbers['registry']);index=[]
    # Archive rows without a fit group keep their direct guarded proposals and
    # an explicit absent group number. No per-cell synthetic fine group is made.
    for ds in sorted(mi):
        storage_gate(root,cfg);cm,ca=read_source(candidates,ci[ds]);meta=json.loads((marker/'outputs'/mi[ds]['metadata']).read_text());_,groups=read_source(marker,meta['source_graph_record']);rows=reader.rows(ds);ids=reader.native_ids(ds);code=groups['source_local_group_code'];fit=groups['fit_pool_candidate'];numbers_for_source=np.asarray([r['global_subcluster_number'] for r in source_catalog_records[ds]],np.int32);global_number=np.full(len(rows),-1,np.int32);global_number[fit]=numbers_for_source[code[fit]]
        values={**ca,'global_subcluster_number':global_number,'guarded_source_subset_code':code}
        output=root/'outputs'/(ds+'_guarded_identity_catalog_rows_v1.npz')
        with output.open('xb') as handle:np.savez_compressed(handle,**values)
        with np.load(output,allow_pickle=False) as f:
            if set(f.files)!=set(values):raise ValueError('Catalog cell serialization keys differ')
            for name,value in values.items():np.testing.assert_array_equal(f[name],value)
        label_by_code=np.asarray([r['subcluster'] for r in source_catalog_records[ds]],object);display=np.full(len(rows),'no_fitted_subcluster',object);display[fit]=label_by_code[code[fit]];status=np.where(fit,'audited_native_guarded_candidate','archive_excluded_direct_candidate_no_fitted_subcluster');old=reader.frozen_labels(ds)
        before=np.asarray(old['fine']).astype(str)
        missing=np.isin(np.char.lower(before),['','nan','none','<na>'])
        before=np.where(missing,'old_label_missing',before)
        crosswalk=version_crosswalk(ds,ids.astype(str),before,display,status)
        mp=output.with_suffix('.json');save(mp,{'source':ds,'n_archive_rows':len(rows),'n_fit_group_rows':int(fit.sum()),'n_archive_excluded_rows_retained':int((~fit).sum()),'native_ids_sha256':text_axis_sha256(ids),'arrays':output.name,'arrays_sha256':sha256_file(output),'label_version':cfg['label_version'],'candidate_dictionaries':cm['dictionaries'],'parent_local_subcluster_key':'source + parent_source_local_group_code; parent stable SHA in bound catalog','global_number_dictionary':catalog_file.name,'marker_origin_and_evidence_type_in_catalog':True,'old_to_new_display_crosswalk':crosswalk,'frozen_author_labels_and_RNA_state_vectors_retained_in_bound_native_inputs':True,'clinical_and_composition_views_not_newly_qualified':True,'new_authoritative_identity':False})
        index.append({'source':ds,'n_archive_rows':len(rows),'arrays':output.name,'arrays_sha256':sha256_file(output),'metadata':mp.name,'metadata_sha256':sha256_file(mp)})
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'n_complete_sources':len(index),'stage':'all_archive_candidate_rows_group_registry_and_crosswalk_written'});reader._fields.clear();del ca,groups,values,ids,old,display,crosswalk;gc.collect()
    if sum(r['n_archive_rows'] for r in index)!=mt['n_archive_rows'] or len(global_catalog)!=mt['n_guarded_subsets']:raise ValueError('Whole archive catalog scope differs')
    ip=root/'outputs/full_guarded_identity_catalog_index_v1.json';save(ip,index)
    save(root/'provenance/terminal_v1.json',{'status':'ALL_ARCHIVE_GUARDED_IDENTITY_GROUP_KEYS_MARKERS_AND_CROSSWALKS_MATERIALIZED_INDEPENDENT_AUDIT_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(index),'n_archive_rows':sum(r['n_archive_rows'] for r in index),'n_exact_guarded_groups':len(global_catalog),'index_sha256':sha256_file(ip),'catalog_sha256':sha256_file(catalog_file),'display_registry_sha256':sha256_file(registry),'marker_full_audit_terminal_sha256':sha256_file(auditor/'provenance/terminal_v1.json'),'new_authoritative_identity':False,'new_clinical_composition_or_malignancy_qualification':False,'expression_matrix_copied':False,'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
