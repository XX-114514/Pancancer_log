"""Recompute own native expression markers for every guarded local subset.

Subsets are exact parent-group x guarded-celltype memberships. They remain
candidate partitions, with state/identity and archive/model eligibility separate.
The accepted native-expression arithmetic is reused; no parent marker is copied.
"""
import argparse,datetime,gc,json,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from extract_full_source_local_marker_evidence_v1 import source_evidence,storage_gate,save_arrays
from full_context_atlas_v1 import sha256_file
from guarded_native_marker_partition_v1 import partition_guarded_subsets
from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256,row_axis_sha256
from summarize_full_native_annotation_evidence_v1 import read_source


def main(root):
    began=time.monotonic();cfg=json.loads((root/'config/native_marker_evidence_spec_v1.json').read_text())
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Guarded marker bound input changed')
    tests=[(['test_guarded_native_marker_partition_v1.py'],3),(['test_native_marker_unit_support_v1.py'],3),(['test_native_annotation_review_v1.py','NativeReviewTests.test_native_total_zero_RNA_and_group_denominators','NativeReviewTests.test_marker_fallback_is_own_group_expression_not_DEG'],2)]
    for arguments,n in tests:
        r=subprocess.run([sys.executable,'-B',str(root/'tests'/arguments[0]),*arguments[1:]],capture_output=True,text=True,timeout=120)
        (root/'provenance'/(arguments[0]+'.stderr.txt')).write_text(r.stderr)
        if r.returncode or ('Ran %d tests'%n) not in r.stderr or 'skipped' in r.stderr.lower():raise ValueError('Actual own-subset marker numerical tests failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':8,'no_skips':True})
    candidates,graphs,reference=map(Path,[cfg['candidates'],cfg['graphs'],cfg['reference']])
    ci=index_at(candidates,'full_native_cell_identity_candidate_index_v1.json');gi=index_at(graphs,'full_native_source_local_graph_index_v1.json');ri=index_at(reference,'full_archive_reference_index_v1.json')
    reader=NativeArchiveCollection(cfg['archive_manifest']);expected=set(reader.archive_partitions)
    if set(ci)!=expected or set(gi)!=expected or set(ri)!=expected:raise ValueError('Guarded marker all-source scope differs')
    results=[]
    for ds in sorted(expected,key=lambda s:(ci[s]['n_archive_rows'],s)):
        storage_gate(root,cfg);cm,ca=read_source(candidates,ci[ds]);gm,ga=read_source(graphs,gi[ds]);rows=reader.rows(ds)
        np.testing.assert_array_equal(rows,ca['native_row']);np.testing.assert_array_equal(rows,ga['native_row'])
        for name in ['fit_pool_candidate','fine_graph_fitted']:
            np.testing.assert_array_equal(ca[name],ga[name])
        np.testing.assert_array_equal(ca['parent_source_local_group_code'],ga['source_local_group_code'])
        if cm['native_ids_sha256']!=text_axis_sha256(reader.native_ids(ds)) or gm['native_ids_sha256']!=cm['native_ids_sha256']:raise ValueError('Guarded marker native identity differs')
        names=cm['dictionaries']['consensus_celltype'];types=ca['consensus_celltype_code'];fit=ga['fit_pool_candidate']
        if np.any(types>=len(names)):raise ValueError('Guarded subset type dictionary capacity differs')
        code,pairs=partition_guarded_subsets(ga['source_local_group_code'],types,fit)
        if not len(pairs):raise ValueError('No fitted native subset to characterize')
        descriptors=[];parents={r['source_local_group_code']:r for r in cm['parent_group_guarded_counts']}
        for group,(parent,t) in enumerate(pairs):
            member=code==group;original=ga['source_local_group_code']==parent
            descriptors.append({'guarded_subset_code':group,'parent_source_local_group_code':int(parent),'guarded_celltype_candidate':names[t],'n_subset_fit_rows':int(member.sum()),'n_parent_fit_rows':int(original.sum()),'membership_identical_to_parent':bool(np.array_equal(member,original)),'stable_parent_group_sha256':parents[parent]['stable_parent_group_sha256'],'native_subset_rows_sha256':row_axis_sha256(rows[member]),'marker_origin':'recomputed_from_exact_own_native_subset','candidate_not_final_identity':True})
        arrays={'native_row':rows,'source_local_group_code':code,'fine_graph_fitted':ga['fine_graph_fitted'],'fit_pool_candidate':fit,'provisional_family_code':ga['provisional_family_code']}
        ap=root/'outputs'/(ds+'_guarded_marker_input_axes_v1.npz');save_arrays(ap,**arrays)
        mp=ap.with_suffix('.json');save(mp,{'source':ds,'native_ids_sha256':cm['native_ids_sha256'],'n_archive_rows':len(rows),'n_source_local_groups':len(pairs),'partition_kind':'exact_frozen_parent_by_guarded_celltype_candidate','family_comparator_origin':'frozen_parent_graph_provisional_family_not_final_consensus','guarded_subsets':descriptors,'candidate_source_record':ci[ds],'original_graph_record':gi[ds],'authoritative_identity_changed':False})
        record={'source':ds,'n_archive_rows':len(rows),'arrays':ap.name,'arrays_sha256':sha256_file(ap),'metadata':mp.name,'metadata_sha256':sha256_file(mp)}
        result=source_evidence(root,cfg,reader,root,record,reference,ri[ds])
        base_output=root/'outputs'/result['metadata'];metadata=json.loads(base_output.read_text())
        metadata.update(partition_kind='exact_frozen_parent_by_guarded_celltype_candidate',guarded_subsets=descriptors,all_subset_markers_recomputed_from_exact_membership=True,parent_marker_copied_into_subset=False,provisional_family_comparator_origin='frozen_parent_graph_not_new_consensus',new_final_identity=False,independent_biological_validation=False)
        output=root/'outputs'/(ds+'_guarded_subset_native_marker_evidence_v1.json')
        save(output,metadata);result['metadata']=output.name;result['metadata_sha256']=sha256_file(output);result['guarded_marker_axis_record']=record
        results.append(result)
        save(root/'checkpoints'/(ds+'_complete.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'stage':'all_guarded_local_subsets_own_native_expression_and_unit_support_complete','n_guarded_subsets':len(pairs),'n_complete_sources':len(results)})
        reader._fields.clear();del ca,ga,arrays,metadata,descriptors,code,pairs;gc.collect()
    output=root/'outputs/full_source_local_marker_evidence_index_v1.json';save(output,results)
    save(root/'provenance/terminal_v1.json',{'status':'ALL_GUARDED_LOCAL_SUBSET_OWN_NATIVE_MARKERS_COMPUTED_INDEPENDENT_AUDIT_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(results),'n_archive_rows':sum(r['n_archive_rows'] for r in results),'n_guarded_subsets':sum(r['n_local_groups'] for r in results),'index_sha256':sha256_file(output),'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'inferential_DEG_computed':False,'parent_markers_borrowed':False,'new_final_identity':False})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True);root=Path(parser.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
