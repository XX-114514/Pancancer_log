"""Materialize a conservative working annotation on the accepted TME scope.

Native counts, the archived annotation and source decisions stay immutable.
Structural subgroups, biological identity, states and scientific eligibility
remain separate; this working release is not the final module input.
"""
import argparse,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def hist(x):
    a,b=np.unique(np.asarray(x).astype(str),return_counts=True);return dict(zip(a.tolist(),map(int,b)))

REVIEW=['baseline_guarded_identity_retained','additional_native_pDC_context_supported',
        'pDC_competing_specific_context_unresolved','pDC_required_context_unmeasured_unresolved',
        'pDC_additional_context_not_met_unresolved','baseline_pDC_conservative_related_alias']

def main(root):
    start=time.monotonic();cfg=load(root/'config/spec_v1.json');J=Path(cfg['joined']);P=Path(cfg['pDC_proposals']);A=Path(cfg['pDC_proposal_audit'])
    for b in cfg['bindings'].values():assert sha(b['path'])==b['SHA256']
    for stage in [P,A]:assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    sys.path.insert(0,str(J/'scripts'))
    sys.path.insert(0,str(root/'scripts'))
    from joined_native_atlas_delivery_v1 import JoinedNativeAtlas
    from atlas_identity_release_keys_v1 import membership_key,assign_display_numbers,display_label
    from review_source_qualified_identity_context_v1 import source_status
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    atlas=JoinedNativeAtlas(J/'outputs/joined_native_atlas_delivery_v1.json');view=load(cfg['source_view']);cp=Path(view['work_units']['path']);assert sha(cp)==view['work_units']['sha256'];core=load(cp)['units']
    eligible_sources=sorted({u['source'] for u in core if u['source_qualification']=='来源资格通过'})
    pidx={r['source']:r for r in load(P/'outputs/pDC_context_proposal_index_v1.json')}
    context=Path(load(P/'config/spec_v1.json')['context']);cidx={r['source']:r for r in load(context/'outputs/disputed_pDC_context_index_v1.json')}
    catalog=[];source_records=[];totals=collections.Counter();previous_registry={r['stable_group_sha256']:r['global_subcluster_number'] for r in atlas.groups.values()}
    for ds in eligible_sources:
        storage_gate(root,cfg);s=atlas.source(ds);ca,ea=s['ca'],s['ea'];n=len(s['rows']);gate=source_status(s,core,view);positions=np.flatnonzero(gate==1);nr=len(positions)
        assert nr==sum(u['n_archive_rows'] for u in core if u['source']==ds and u['source_qualification']=='来源资格通过')
        olddict=s['cm']['candidate_dictionaries'];major_dict=list(olddict['consensus_major']);fine_dict=list(olddict['consensus_celltype'])
        for name in ['Myeloid','unresolved']:
            if name not in major_dict:major_dict.append(name)
        for name in ['pDC_related','unresolved']:
            if name not in fine_dict:fine_dict.append(name)
        oldmajor=ca['consensus_major_code'][positions];oldfine=ca['consensus_celltype_code'][positions]
        major=oldmajor.copy();fine=oldfine.copy();review=np.zeros(nr,np.uint8)
        baseline_fine=np.asarray(olddict['consensus_celltype'],object)[oldfine]
        alias=baseline_fine=='pDC';fine[alias]=fine_dict.index('pDC_related');review[alias]=5
        raw_context=None;proposal=None;proposed_input_positions=np.array([],np.int64);new_supported=np.zeros(nr,bool)
        if ds in pidx:
            pi=pidx[ds];assert sha(pi['path'])==pi['SHA256'];pm=load(pi['path']);pa=Path(pi['path']).parent/pm['arrays'];assert sha(pa)==pm['arrays_SHA256']
            with np.load(pa,allow_pickle=False) as f:proposal={k:f[k] for k in f.files}
            loc=np.searchsorted(positions,proposal['source_native_position']);inside=loc<len(positions);inside[inside]&=positions[loc[inside]]==proposal['source_native_position'][inside]
            psel=np.flatnonzero(inside);target=loc[inside];assert np.array_equal(proposal['source_qualification_code'][psel],np.ones(len(psel),np.int8))
            assert np.all(np.asarray(olddict['consensus_major'],object)[oldmajor[target]]=='unresolved')
            np.testing.assert_array_equal(proposal['native_row'][psel],s['rows'][positions[target]])
            disposition=proposal['disposition_code'][psel];review[target]=np.where(disposition==0,1,disposition+1)
            supported=target[disposition==0];new_supported[supported]=True;major[supported]=major_dict.index('Myeloid');fine[supported]=fine_dict.index('pDC_related')
            if len(supported):
                ci=cidx[ds];assert sha(ci['path'])==ci['SHA256'];cm=load(ci['path']);cf=Path(ci['path']).parent/cm['arrays'];assert sha(cf)==cm['arrays_SHA256']
                with np.load(cf,allow_pickle=False) as f:
                    cases=f['stratum_code']==0;raw_context={'positions':f['source_native_position'][cases],'values':f['native_gene_values'][cases],'genes':cm['genes'],'measured':f['measured']}
                np.testing.assert_array_equal(raw_context['positions'],proposal['source_native_position'])
        current_fine=np.asarray(fine_dict,object)[fine];current_major=np.asarray(major_dict,object)[major];fit=ca['fit_pool_candidate'][positions].astype(bool)
        original_group=ca['global_subcluster_number'][positions];group_code=np.full(nr,-1,np.int32);source_groups=[]
        # Each reviewed identity split and source-scope intersection gets its
        # own exact membership key. The old parent number is never reused for
        # a changed membership, and numbers do not establish biological types.
        for parent_number in np.unique(original_group[fit]):
            assert parent_number>=1;parent=atlas.groups[int(parent_number)];members=np.flatnonzero(fit&(original_group==parent_number))
            for label in sorted(set(current_fine[members])):
                take=members[current_fine[members]==label];native_positions=positions[take];key=membership_key(ds,s['ids'][native_positions].astype(str).tolist())
                old_key=parent['stable_group_sha256'];identical=key==old_key
                record={'source':ds,'source_annotation_group_code':len(source_groups),'stable_group_sha256':key,'parent_archive_group_sha256':old_key,
                    'parent_archive_display_number':int(parent_number),'parent_source_local_group_code':parent['parent_source_local_group_code'],
                    'n_scope_group_rows':len(take),'consensus_celltype':label,'major_distribution':hist(current_major[take]),
                    'native_rows_SHA256':hashlib.sha256(np.asarray(s['rows'][native_positions],'<i8').tobytes()).hexdigest(),
                    'membership_identical_to_archive_parent':identical,'new_biological_type_not_inferred_from_group_number':True,
                    'marker_DEG_status':'no_qualified_exact_membership_inferential_DEG_supplied'}
                if new_supported[take].any():
                    assert new_supported[take].all() and raw_context is not None
                    hit=np.searchsorted(raw_context['positions'],native_positions);np.testing.assert_array_equal(raw_context['positions'][hit],native_positions)
                    j=raw_context['genes'].index('LILRA4');assert raw_context['measured'][j];v=raw_context['values'][hit,j];assert (v>=2).all()
                    record.update(marker='LILRA4',marker_origin=key,marker_evidence_type='own_subset_native_expression_fixed_identity_context',parent_marker_borrowed=False,
                        marker_selection_rule='fixed identity-context marker; exact subgroup native expression measured; not genome-wide DEG ranking',
                        marker_native_mean=float(v.mean()),marker_detection_fraction=float((v>0).mean()),marker_native_minimum=float(v.min()),
                        borrow_reason=None)
                elif identical:
                    record.update({k:parent.get(k) for k in ['marker','marker_origin','marker_evidence_type','parent_marker_borrowed','borrow_reason']})
                    record['marker_selection_rule']='accepted exact-membership archive own-marker evidence reused without changing membership'
                else:
                    record.update(marker=parent['marker'],marker_origin=old_key,marker_evidence_type='parent_group_borrowed_explicit',parent_marker_borrowed=True,
                        borrow_reason='source-qualified intersection differs from the measured full archive parent; own expression/DEG not recomputed on this intersection',
                        marker_selection_rule='explicit borrowed parent marker; not relabeled as exact-subset expression')
                group_code[take]=len(source_groups);source_groups.append(record);catalog.append(record)
        assert np.all(group_code[fit]>=0) and np.all(group_code[~fit]==-1)
        author=atlas.native.author_labels(ds);author_values=np.asarray(author['values'][positions]).astype(str);author_dict,author_code=np.unique(author_values,return_inverse=True)
        arrays={'source_native_position':positions,'native_row':s['rows'][positions],
            'baseline_consensus_major_code':oldmajor,'baseline_consensus_celltype_code':oldfine,
            'consensus_major_code':major,'consensus_celltype_code':fine,'annotation_review_code':review,
            'baseline_coarse_evidence_code':ca['coarse_evidence_code'][positions],'baseline_fine_evidence_code':ca['fine_evidence_code'][positions],
            'fit_pool_candidate_not_final_QC':fit,'source_qualification_code':gate[positions],
            'parent_source_local_group_code':ca['parent_source_local_group_code'][positions],
            'baseline_global_subcluster_number':original_group,'annotation_group_code':group_code,
            'native_family_bitmask':ca['native_family_bitmask'][positions],'reference_family_bitmask':ca['reference_family_bitmask'][positions],
            'author_label_code':author_code.astype(np.int32),'RNA_state_native_logCP10K':ea['state_native_logCP10K'][positions],
            'native_total_RNA':ea['native_total_RNA'][positions]}
        op=root/'outputs'/(ds+'_TME_scoped_annotation_v1.npz')
        with op.open('xb') as f:np.savez_compressed(f,**arrays)
        with np.load(op,allow_pickle=False) as f:
            for k,v in arrays.items():np.testing.assert_array_equal(f[k],v)
        metadata={'source':ds,'arrays':op.name,'arrays_SHA256':sha(op),'n_source_qualified_rows':nr,'n_archive_source_rows':n,
            'n_fit_candidate_not_final_QC':int(fit.sum()),'n_without_fitted_group':int((~fit).sum()),'native_ids_source_SHA256':s['em']['native_ids_sha256'],
            'baseline_candidate_dictionaries':olddict,'candidate_dictionaries':{'consensus_major':major_dict,'consensus_celltype':fine_dict,'annotation_review':REVIEW},
            'author_label_dictionary':author_dict.tolist(),'author_label_status':author['status'],'state_order':s['em']['state_order'],
            'state_measurement_coverage':s['em']['state_measurement_coverage'],'n_source_annotation_groups':len(source_groups),
            'major_counts':hist(current_major),'celltype_counts':hist(current_fine),'annotation_review_counts':hist(review),
            'old_to_new_celltype_counts':[{'baseline':a,'current':b,'n_cells':int(c)} for (a,b),c in sorted(collections.Counter(zip(baseline_fine,current_fine)).items())],
            'source_view_SHA256':cfg['bindings']['source_view']['SHA256'],'baseline_candidate_metadata_SHA256':atlas.ci[ds]['metadata_sha256'],
            'pDC_proposal_metadata_SHA256':pidx[ds]['SHA256'] if ds in pidx else None,
            'label_version':'TME_source_qualified_working_annotation_20261027_v1','identity_status':'conservative_working_annotation_not_final_authoritative_taxonomy',
            'confidence_scale':'categorical_evidence_not_calibrated_probability','no_new_QC_malignancy_patient_or_module_eligibility':True}
        source_records.append(metadata);totals.update(n_sources=1,n_source_qualified_rows=nr,n_fit_candidates_not_final_QC=int(fit.sum()),
            n_without_fitted_group=int((~fit).sum()),n_context_pDC_updates=int(new_supported.sum()),n_existing_pDC_related_aliases=int(alias.sum()),n_exact_scoped_groups=len(source_groups))
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(totals)})
        del s,ca,ea,arrays,raw_context,author_values;atlas._cache.clear();atlas.native._fields.clear();gc.collect()
    numbers=assign_display_numbers([r['stable_group_sha256'] for r in catalog],previous_registry)
    for r in catalog:
        number=numbers['current'][r['stable_group_sha256']];r['global_subcluster_number']=number;r['global_subcluster_id']='TME_annotation_group_'+str(number);r['subcluster']=display_label(r['consensus_celltype'],r['marker'],number)
    group_path=root/'outputs/TME_scoped_annotation_group_catalog_v1.json';save(group_path,{'groups':catalog,'label_numbers_are_exact_groups_not_shared_biological_types':True})
    registry=root/'outputs/TME_annotation_display_number_registry_v1.json';save(registry,numbers['registry']);index=[]
    for m in source_records:
        m['group_catalog']={'path':str(group_path),'SHA256':sha(group_path)};mp=root/'outputs'/(m['source']+'_TME_scoped_annotation_v1.json');save(mp,m)
        index.append({'source':m['source'],'metadata':mp.name,'metadata_SHA256':sha(mp),'arrays':m['arrays'],'arrays_SHA256':m['arrays_SHA256'],'n_source_qualified_rows':m['n_source_qualified_rows']})
    assert totals['n_source_qualified_rows']==sum(u['n_archive_rows'] for u in core if u['source_qualification']=='来源资格通过')
    ip=root/'outputs/TME_scoped_annotation_index_v1.json';save(ip,index)
    sp=root/'outputs/TME_scoped_annotation_summary_v1.json';save(sp,{'totals':dict(totals),'sources':[{k:m[k] for k in ['source','n_source_qualified_rows','major_counts','celltype_counts','annotation_review_counts','n_without_fitted_group']} for m in source_records],
        'scope_is_source_qualified_not_final_QC':True,'all_unresolved_cells_retained':True,'native_counts_copied':False,'new_CM_started':False})
    manifest={'schema':'TME_scoped_working_annotation_v1','joined_manifest':str(J/'outputs/joined_native_atlas_delivery_v1.json'),
        'source_view':cfg['source_view'],'bindings':cfg['bindings'],'index':{'path':str(ip),'SHA256':sha(ip)},'group_catalog':{'path':str(group_path),'SHA256':sha(group_path)},
        'summary':{'path':str(sp),'SHA256':sha(sp)},'label_version':'TME_source_qualified_working_annotation_20261027_v1',
        'release_status':'WORKING_ANNOTATION_INDEPENDENT_AUDIT_REQUIRED_NOT_MODULE_FIT_READY','source_records':totals['n_source_qualified_rows']}
    mp=root/'outputs/TME_scoped_working_annotation_manifest_v1.json';save(mp,manifest)
    save(root/'provenance/terminal_v1.json',{'status':'ALL_SOURCE_QUALIFIED_TME_WORKING_ANNOTATION_MATERIALIZED_INDEPENDENT_AUDIT_REQUIRED',
        'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(totals),'manifest_SHA256':sha(mp),'index_SHA256':sha(ip),'catalog_SHA256':sha(group_path),
        'summary_SHA256':sha(sp),'native_count_matrix_read':False,'native_counts_labels_old_CM_unchanged':True,
        'output_bytes':sum(p.stat().st_size for p in root.rglob('*') if p.is_file()),'wall_seconds':time.monotonic()-start,
        'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
