"""Bind numeric acceptance, exercise the native-backed reader, export summaries."""
import argparse,collections,csv,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def table(p,rows):
    with Path(p).open('x',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def main(root):
    start=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);A=Path(cfg['auditor']);J=Path(cfg['joined'])
    for stage in [P,A]:assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    manifest_name=cfg.get('producer_manifest_name','TME_scoped_working_annotation_manifest_v1.json')
    if Path(manifest_name).name!=manifest_name:raise ValueError('Manifest name must be a filename within producer outputs')
    m=P/'outputs'/manifest_name;at=A/'provenance/terminal_v1.json';t=load(at)
    assert sha(m)==t['producer_manifest_SHA256'];manifest=load(m);groups=load(manifest['group_catalog']['path'])['groups'];summary=load(manifest['summary']['path'])
    acceptance={'schema':'accepted_TME_scoped_working_annotation_v1','accepted_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'producer_manifest':{'path':str(m),'SHA256':sha(m)},'independent_audit_terminal':{'path':str(at),'SHA256':sha(at)},
        'acceptance_scope':'complete_numeric_working_annotation_not_biological_accuracy_or_final_taxonomy',
        'final_module_input_ready':False,'new_CM_started':False,'all_source_qualified_records_and_unresolved_retained':True,
        'source_counts_are_not_independent_patients':True}
    ap=root/'outputs/accepted_TME_scoped_working_annotation_v1.json';save(ap,acceptance)
    sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(root/'scripts'))
    from TME_scoped_annotation_view_v1 import TMEScopedAnnotationView
    reader=TMEScopedAnnotationView(ap);roles={};source_rows=[];major=collections.Counter();fine=collections.Counter();review=collections.Counter()
    for source in reader.sources:
        s=reader.source(source);a=s['arrays'];meta=s['metadata'];nr=len(a['native_row'])
        for role,mask in [('new_pDC_context',a['annotation_review_code']==1),('new_TCR_context',a['annotation_review_code']==6),('competing_or_incomplete_context',np.isin(a['annotation_review_code'],[2,3,4])),('without_fitted_group',~a['fit_pool_candidate_not_final_QC'])]:
            if role not in roles and mask.any():roles[role]=(source,np.array([int(np.flatnonzero(mask)[0])],np.int64))
        if 'source_filtered_boundary' not in roles and len(a['source_native_position'])>1:
            gaps=np.flatnonzero(np.diff(a['source_native_position'])>1)
            if len(gaps):roles['source_filtered_boundary']=(source,np.array([int(gaps[0]),int(gaps[0]+1),int(gaps[0])],np.int64))
        major.update(meta['major_counts']);fine.update(meta['celltype_counts']);review.update(meta['annotation_review_counts'])
        source_rows.append({'source':source,'n_source_qualified_records':nr,'n_fit_candidates_not_final_QC':meta['n_fit_candidate_not_final_QC'],
            'n_without_fitted_group':meta['n_without_fitted_group'],'n_major_unresolved':meta['major_counts'].get('unresolved',0),
            'n_new_pDC_related_context_supported':meta['annotation_review_counts'].get('1',0),'n_structural_groups_not_celltypes':meta['n_source_annotation_groups']})
    tested=[]
    required={'new_pDC_context','competing_or_incomplete_context','without_fitted_group','source_filtered_boundary'}
    assert required.issubset(roles) and set(roles).issubset(required|{'new_TCR_context'})
    if summary['totals'].get('n_TCR_coarse_updates',0):assert 'new_TCR_context' in roles
    for role,(source,p) in roles.items():
        s=reader.source(source);a=s['arrays'];result=reader.metadata(source,p);columns=result['columns']
        np.testing.assert_array_equal(columns['source_native_position'],a['source_native_position'][p]);np.testing.assert_array_equal(columns['native_row'],a['native_row'][p]);assert (columns['source_qualification_code']==1).all()
        if role=='new_pDC_context':
            assert (columns['consensus_celltype']=='pDC_related').all() and (columns['consensus_major']=='Myeloid').all()
            assert (columns['coarse_evidence']=='additional_native_context_pDC_related_working_proposal').all()
            assert (columns['fine_evidence']=='additional_native_context_not_orthogonal_identity_validation').all()
            assert 'baseline__coarse_evidence' in columns and 'baseline__fine_evidence' in columns
        if role=='competing_or_incomplete_context':assert (columns['consensus_major']=='unresolved').all()
        if role=='new_TCR_context':
            assert (columns['consensus_major']=='T').all() and (columns['consensus_celltype']=='T').all()
            assert (columns['coarse_evidence']=='additional_TCR_core_reference_and_local_context_supported').all()
            assert (columns['fine_evidence']=='coarse_T_RNA_context_only_no_new_fine_identity').all()
        if role=='without_fitted_group':assert (columns['global_subcluster_number']==-1).all() and all(v is None for v in columns['subcluster'])
        genes=['CD3D','LILRA4','__UNMEASURED_SENTINEL__'];actual=reader.read_counts(source,p,genes);expected=reader.atlas.read_counts(source,a['source_native_position'][p],genes)
        assert (actual['counts']!=expected['counts']).nnz==0
        for k in ['native_row','native_total_counts','measured']:np.testing.assert_array_equal(actual[k],expected[k])
        assert not actual['measured'][-1] and actual['missing_features_are_not_biological_zeros']
        tested.append({'role':role,'n_query_rows':len(p),'metadata_axis_and_native_counts_and_missingness_matched':True})
    try:reader.source('__not_in_TME_scope__')
    except KeyError:pass
    else:raise AssertionError('Out-of-scope source was accepted')
    source=reader.sources[0]
    for bad in [[-1],[len(reader.source(source)['arrays']['native_row'])]]:
        try:reader.metadata(source,bad)
        except IndexError:pass
        else:raise AssertionError('Out-of-range scoped row accepted')
    totals=summary['totals'];assert sum(major.values())==sum(fine.values())==totals['n_source_qualified_rows']
    major_rows=[{'major':k,'n_source_qualified_records':v,'fraction_of_source_qualified_archive':v/totals['n_source_qualified_rows']} for k,v in major.most_common()]
    fine_rows=[{'celltype':k,'n_source_qualified_records':v,'fraction_of_source_qualified_archive':v/totals['n_source_qualified_rows'],'identity_level':'unresolved' if k=='unresolved' else 'working_RNA_identity_or_coarse_fallback_not_accuracy_certified'} for k,v in fine.most_common()]
    table(root/'outputs/TME_current_major_summary_v1.tsv',major_rows);table(root/'outputs/TME_current_celltype_summary_v1.tsv',fine_rows);table(root/'outputs/TME_current_source_annotation_coverage_v1.tsv',source_rows)
    group_rows=[{k:g.get(k) for k in ['source','global_subcluster_number','consensus_celltype','marker','subcluster','n_scope_group_rows','marker_evidence_type','parent_marker_borrowed','borrow_reason','stable_group_sha256','marker_origin','parent_archive_group_sha256','parent_archive_display_number','marker_native_mean_CP10K','marker_detection_fraction','marker_comparator_scope','marker_comparator_cell_count','marker_DEG_status','marker_small_group_status']} for g in groups]
    table(root/'outputs/TME_current_subcluster_catalog_v1.tsv',group_rows)
    if 'annotation_membership_crosswalk' in manifest:
        cross=manifest['annotation_membership_crosswalk'];assert sha(cross['path'])==cross['SHA256']
        table(root/'outputs/TME_annotation_membership_crosswalk_v3.tsv',load(cross['path']))
    if 'marker_display_crosswalk' in manifest:
        cross=manifest['marker_display_crosswalk'];assert sha(cross['path'])==cross['SHA256']
        table(root/'outputs/TME_marker_display_crosswalk_v2.tsv',load(cross['path']))
        for g in groups:
            assert not g['parent_marker_borrowed'] and g['marker_origin']==g['stable_group_sha256']
        # Exercise actual metadata with the new catalog, including at least
        # one changed marker, rather than accepting only the manifest hash.
        changed=next(g for g in groups if g.get('previous__marker',g['marker'])!=g['marker'])
        src=reader.source(changed['source']);pos=np.flatnonzero(src['arrays']['annotation_group_code']==changed['source_annotation_group_code'])[:2]
        columns=reader.metadata(changed['source'],pos)['columns']
        assert len(pos)>0 and all(v==changed['marker'] for v in columns['marker'])
        assert all(v==changed['subcluster'] for v in columns['subcluster'])
        assert all(v==manifest['label_version'] for v in columns['label_version'])
        assert all(v==changed['stable_group_sha256'] for v in columns['marker_origin'])
        tested.append({'role':'changed_marker_current_catalog_readback','n_query_rows':len(pos),'metadata_marker_display_and_origin_matched':True})
    marker_counts=collections.Counter(g['marker_evidence_type'] for g in groups);marker_cell_counts=collections.Counter()
    for g in groups:marker_cell_counts[g['marker_evidence_type']]+=g['n_scope_group_rows']
    report={'status':'ACCEPTED_WORKING_ANNOTATION_READER_NATIVE_COUNT_MAPPING_AND_PRIVATE_SUMMARIES_DELIVERED',
        'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':totals,'major_counts':dict(major),'celltype_counts':dict(fine),
        'annotation_review_counts':dict(review),'marker_group_counts':dict(marker_counts),'marker_cell_counts':dict(marker_cell_counts),
        'bounded_reader_cases':tested,'rejected_out_of_scope_source_and_invalid_positions':True,
        'all_count_readbacks_equal_original_source_counts':True,'native_counts_read_only_for_bounded_interface_cases':True,
        'acceptance_contract_SHA256':sha(ap),'producer_manifest_SHA256':sha(m),'independent_audit_terminal_SHA256':sha(at),
        'new_CM_started':False,'final_module_input_ready':False,'biological_annotation_accuracy_not_estimated':True,
        'wall_seconds':time.monotonic()-start}
    rp=root/'outputs/TME_annotation_delivery_report_v1.json';save(rp,report)
    files={p.name:sha(p) for p in (root/'outputs').glob('*') if p.is_file()}
    save(root/'provenance/terminal_v1.json',{'status':report['status'],'completed_UTC':report['completed_UTC'],'totals':totals,
        'delivery_report_SHA256':sha(rp),'files':files,'new_CM_started':False,'final_module_input_ready':False,'wall_seconds':report['wall_seconds']})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except BaseException:
        import traceback
        save(root/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
