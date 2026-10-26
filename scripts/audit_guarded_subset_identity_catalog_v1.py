"""Independent all-row/group checks for the native guarded display catalog."""
import argparse,collections,datetime,gc,hashlib,json,resource,time
from pathlib import Path
import numpy as np
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
from summarize_full_native_annotation_evidence_v1 import read_source


def main(root):
    began=time.monotonic();cfg=json.loads((root/'config/catalog_audit_spec_v1.json').read_text());producer=Path(cfg['producer'])
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Catalog audit bound input changed')
    index=index_at(producer,'full_guarded_identity_catalog_index_v1.json');terminal=json.loads((producer/'provenance/terminal_v1.json').read_text());spec=json.loads((producer/'config/guarded_identity_catalog_spec_v1.json').read_text())
    for b in spec['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Catalog upstream input changed')
    cp=producer/'outputs/guarded_subset_identity_catalog_v1.json';rp=producer/'outputs/global_display_number_registry_v1.json'
    if sha256_file(cp)!=terminal['catalog_sha256'] or sha256_file(rp)!=terminal['display_registry_sha256']:raise ValueError('Catalog/registry changed')
    catalog=json.loads(cp.read_text());registry=json.loads(rp.read_text());groups=catalog['groups'];by_number={r['global_subcluster_number']:r for r in groups}
    if len(by_number)!=len(groups) or len({r['stable_group_sha256'] for r in groups})!=len(groups):raise ValueError('Duplicate group number or stable key')
    if any(type(v) is not int or v<1 for v in registry.values()) or len(set(registry.values()))!=len(registry):raise ValueError('Invalid append-only number registry')
    if catalog['new_authoritative_identity'] or not catalog['confidence_is_evidence_provenance_not_probability']:raise ValueError('Catalog evidence promoted')
    candidates=Path(spec['candidates']);marker=Path(spec['markers']);ci=index_at(candidates,'full_native_cell_identity_candidate_index_v1.json');mi=index_at(marker,'full_source_local_marker_evidence_index_v1.json');reader=NativeArchiveCollection(spec['archive_manifest'])
    if set(index)!=set(ci) or set(ci)!=set(mi) or set(mi)!=set(reader.archive_partitions):raise ValueError('Full catalog source scope differs')
    reports=[];verified_keys=set()
    for ds in sorted(index):
        storage_gate(root,cfg);record=index[ds];m,a=read_source(producer,record);cm,ca=read_source(candidates,ci[ds]);mp=marker/'outputs'/mi[ds]['metadata']
        if sha256_file(mp)!=mi[ds]['metadata_sha256']:raise ValueError('Audited own-marker metadata changed')
        native=json.loads(mp.read_text());_,partition=read_source(marker,native['source_graph_record']);rows=reader.rows(ds);ids=reader.native_ids(ds)
        np.testing.assert_array_equal(a['native_row'],rows)
        if m['native_ids_sha256']!=text_axis_sha256(ids) or m['new_authoritative_identity'] or m['clinical_and_composition_views_not_newly_qualified'] is not True:raise ValueError('Catalog source identity or qualification differs')
        if set(a)!=set(ca)|{'global_subcluster_number','guarded_source_subset_code'}:raise ValueError('Catalog altered retained candidate schema')
        for key,value in ca.items():np.testing.assert_array_equal(a[key],value)
        code=partition['source_local_group_code'];fit=partition['fit_pool_candidate'];np.testing.assert_array_equal(a['guarded_source_subset_code'],code)
        if not np.all(a['global_subcluster_number'][~fit]==-1):raise ValueError('Excluded row acquired fitted subcluster')
        if m['n_fit_group_rows']!=int(fit.sum()) or m['n_archive_excluded_rows_retained']!=int((~fit).sum()):raise ValueError('Catalog denominator changed')
        support={(r['source_local_group_code'],r['gene']):r for r in native['marker_unit_direction_summaries']}
        source_groups={r['source_subset_code']:r for r in groups if r['source']==ds}
        if set(source_groups)!=set(range(len(native['groups']))):raise ValueError('Source exact-group coverage differs')
        number_to_label={}
        for d,g in zip(native['guarded_subsets'],native['groups']):
            local=d['guarded_subset_code'];positions=np.flatnonzero(code==local);r=source_groups[local];h=hashlib.sha256()
            for token in ['atlas_exact_membership_v1',ds,*sorted(str(v) for v in ids[positions])]:
                data=token.encode('utf-8');h.update(len(data).to_bytes(8,'big'));h.update(data)
            key=h.hexdigest();verified_keys.add(key);number=r['global_subcluster_number']
            if r['stable_group_sha256']!=key or registry.get(key)!=number or r['global_subcluster_id']!='atlas_group_'+str(number):raise ValueError('Independent stable membership/number replay differs')
            np.testing.assert_array_equal(a['global_subcluster_number'][positions],np.full(len(positions),number,np.int32))
            if r['n_subset_fit_rows']!=len(positions) or r['consensus_celltype_candidate']!=d['guarded_celltype_candidate'] or r['stable_parent_group_sha256']!=d['stable_parent_group_sha256']:raise ValueError('Group type, parent or denominator changed')
            if r['native_subset_rows_sha256']!=hashlib.sha256(np.asarray(rows[positions],dtype='<i8').tobytes()).hexdigest() or r['membership_identical_to_parent']!=d['membership_identical_to_parent']:raise ValueError('Native group row provenance differs')
            expected=g['own_group_markers'][0]['gene'] if g['own_group_markers'] else None
            if r['marker']!=expected or (expected is not None and (r['marker_origin']!=key or r['marker_evidence_type']!='own_subset_expression_specificity' or r['evidence_record_sha256']!=mi[ds]['metadata_sha256'])) or r['parent_marker_borrowed']:raise ValueError('Own marker origin/ranking differs')
            label=d['guarded_celltype_candidate']+'_'+(expected or 'marker_unresolved')+'_'+str(number)
            if r['subcluster']!=label or r['label_version']!=spec['label_version']:raise ValueError('Single-marker display/version replay differs')
            evidence=support.get((local,expected));nu=evidence['n_comparable_expression_units'] if evidence else 0;positive=evidence['n_units_positive_mean_ratio'] if evidence else 0
            status='own_marker_multiple_expression_unit_direction_supported' if nu>=2 and positive>nu/2 else 'own_marker_expression_support_without_multiple_unit_majority'
            if len(positions)<30:status='small_group_own_expression_descriptive_only'
            if expected is None:status='no_accepted_own_marker'
            if r['n_comparable_expression_units_not_patients']!=nu or r['n_units_positive_mean_ratio']!=positive or r['group_expression_evidence_status']!=status or r['marker_DEG_status']!='no_qualified_exact_membership_inferential_DEG_input_supplied':raise ValueError('Marker unit/small-group evidence status differs')
            number_to_label[number]=label
        old=np.asarray(reader.frozen_labels(ds)['fine']).astype(str);old=np.where(np.isin(np.char.lower(old),['','nan','none','<na>']),'old_label_missing',old)
        new=np.asarray([number_to_label[int(n)] if n>=0 else 'no_fitted_subcluster' for n in a['global_subcluster_number']],str);status=np.where(fit,'audited_native_guarded_candidate','archive_excluded_direct_candidate_no_fitted_subcluster');counts=collections.Counter(zip(old,new,status));before_to_after=collections.defaultdict(set);after_to_before=collections.defaultdict(set)
        for before,after in zip(old,new):before_to_after[before].add(after);after_to_before[after].add(before)
        observed={(r['old_label'],r['new_label'],r['new_identity_status']):r for r in m['old_to_new_display_crosswalk']}
        if set(observed)!=set(counts) or len(observed)!=len(m['old_to_new_display_crosswalk']):raise ValueError('Crosswalk drops or duplicates archive row categories')
        for key,count in counts.items():
            r=observed[key];split=len(before_to_after[key[0]])>1;merge=len(after_to_before[key[1]])>1;relation='many_to_many' if split and merge else 'split' if split else 'merge' if merge else 'one_to_one'
            if r['source']!=ds or r['n_archive_records']!=count or r['old_label_splits']!=split or r['new_label_merges']!=merge or r['relation']!=relation or r['mapping_is_identity_truth']:raise ValueError('Independent complete archive crosswalk differs')
        reports.append({'source':ds,'n_archive_rows':len(rows),'n_exact_groups':len(source_groups),'n_crosswalk_records':len(observed),'all_candidate_values_exactly_preserved':True,'all_membership_keys_independently_replayed':True})
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_complete_sources':len(reports),**reports[-1]});reader._fields.clear();del a,ca,partition,ids,old,new;gc.collect()
    if verified_keys!={r['stable_group_sha256'] for r in groups} or len(groups)!=terminal['n_exact_guarded_groups'] or sum(r['n_archive_rows'] for r in reports)!=terminal['n_archive_rows']:raise ValueError('Full catalog denominator/key scope differs')
    output=root/'outputs/guarded_identity_catalog_audit_v1.json';save(output,{'sources':reports,'new_authoritative_identity':False,'clinical_and_biological_transfer_not_promoted':True})
    save(root/'provenance/terminal_v1.json',{'status':'ALL_ARCHIVE_GUARDED_GROUP_KEYS_MARKER_DISPLAYS_AND_CROSSWALKS_INDEPENDENTLY_AUDITED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(reports),'n_archive_rows':sum(r['n_archive_rows'] for r in reports),'n_exact_groups':len(groups),'output_sha256':sha256_file(output),'producer_terminal_sha256':sha256_file(producer/'provenance/terminal_v1.json'),'new_authoritative_identity':False,'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
