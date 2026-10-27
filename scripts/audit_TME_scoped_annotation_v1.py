"""Independently check every qualified row and exact scoped annotation group."""
import argparse,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def hist(x):
    return dict(sorted(collections.Counter(map(str,x)).items()))
def exact_key(source,ids):
    assert len(ids)==len(set(ids));h=hashlib.sha256()
    for text in ['atlas_exact_membership_v1',source]+sorted(map(str,ids)):
        data=text.encode();h.update(len(data).to_bytes(8,'big'));h.update(data)
    return h.hexdigest()

def main(root):
    start=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    pc=load(P/'config/spec_v1.json');J=Path(pc['joined']);sys.path.insert(0,str(J/'scripts'))
    from joined_native_atlas_delivery_v1 import JoinedNativeAtlas
    atlas=JoinedNativeAtlas(J/'outputs/joined_native_atlas_delivery_v1.json');view=load(pc['source_view']);corep=Path(view['work_units']['path']);assert sha(corep)==view['work_units']['sha256'];core=load(corep)['units'];by={u['qualification_unit_key_sha256']:u for u in core}
    older=load(view['base_v4_manifest']['path']);origin=older['bindings']['author_bridge_registry'];assert sha(origin['path'])==origin['sha256'];original_source=load(origin['path'])['source']
    overrides={original_source:{'path':older['override_arrays_path'],'sha256':older['override_arrays_sha256']}};overrides.update(view['new_source_overrides'])
    manifest_path=P/'outputs/TME_scoped_working_annotation_manifest_v1.json';terminal=load(P/'provenance/terminal_v1.json');assert sha(manifest_path)==terminal['manifest_SHA256'];manifest=load(manifest_path)
    for b in manifest['bindings'].values():assert sha(b['path'])==b['SHA256']
    ip=Path(manifest['index']['path']);gp=Path(manifest['group_catalog']['path']);assert sha(ip)==manifest['index']['SHA256'] and sha(gp)==manifest['group_catalog']['SHA256'];index=load(ip);groups=load(gp)['groups']
    old_registry={g['stable_group_sha256']:g['global_subcluster_number'] for g in atlas.groups.values()};newkeys=sorted({g['stable_group_sha256'] for g in groups}-set(old_registry));expected_numbers=dict(old_registry)
    for i,key in enumerate(newkeys,max(old_registry.values())+1):expected_numbers[key]=i
    assert len(groups)==len({g['stable_group_sha256'] for g in groups})==len({g['global_subcluster_number'] for g in groups})
    pp=Path(pc['pDC_proposals']);pidx={r['source']:r for r in load(pp/'outputs/pDC_context_proposal_index_v1.json')}
    context=Path(load(pp/'config/spec_v1.json')['context']);cidx={r['source']:r for r in load(context/'outputs/disputed_pDC_context_index_v1.json')}
    totals=collections.Counter();checked_keys=set()
    assert {r['source'] for r in index}=={u['source'] for u in core if u['source_qualification']=='来源资格通过'}
    for record in index:
        ds=record['source'];s=atlas.source(ds);ca,ea=s['ca'],s['ea'];n=len(s['rows'])
        if ds in overrides:
            b=overrides[ds];assert sha(b['path'])==b['sha256']
            with np.load(b['path'],allow_pickle=False) as f:np.testing.assert_array_equal(f['native_row'],s['rows']);eligible=f['source_qualification_code']==1
        else:
            eligible=np.zeros(n,bool)
            for j,u in enumerate(s['um']['unit_records']):
                if by[u['native_unit_key_sha256']]['source_qualification']=='来源资格通过':eligible|=s['ua']['native_expression_unit_code']==j
        expected_positions=np.flatnonzero(eligible);mp=P/'outputs'/record['metadata'];ap=P/'outputs'/record['arrays'];assert sha(mp)==record['metadata_SHA256'] and sha(ap)==record['arrays_SHA256'];m=load(mp)
        with np.load(ap,allow_pickle=False) as f:a={k:f[k] for k in f.files}
        np.testing.assert_array_equal(a['source_native_position'],expected_positions);pos=a['source_native_position'];np.testing.assert_array_equal(a['native_row'],s['rows'][pos]);assert (a['source_qualification_code']==1).all()
        for name,key in [('baseline_consensus_major_code','consensus_major_code'),('baseline_consensus_celltype_code','consensus_celltype_code'),('baseline_coarse_evidence_code','coarse_evidence_code'),('baseline_fine_evidence_code','fine_evidence_code'),('parent_source_local_group_code','parent_source_local_group_code'),('baseline_global_subcluster_number','global_subcluster_number'),('native_family_bitmask','native_family_bitmask'),('reference_family_bitmask','reference_family_bitmask')]:np.testing.assert_array_equal(a[name],ca[key][pos])
        fit=ca['fit_pool_candidate'][pos].astype(bool);np.testing.assert_array_equal(a['fit_pool_candidate_not_final_QC'],fit)
        np.testing.assert_array_equal(a['RNA_state_native_logCP10K'],ea['state_native_logCP10K'][pos]);np.testing.assert_array_equal(a['native_total_RNA'],ea['native_total_RNA'][pos])
        author=atlas.native.author_labels(ds);np.testing.assert_array_equal(np.asarray(m['author_label_dictionary'],object)[a['author_label_code']],np.asarray(author['values'][pos]).astype(str))
        md=s['cm']['candidate_dictionaries'];oldmajor=np.asarray(md['consensus_major'],object)[ca['consensus_major_code'][pos]];oldfine=np.asarray(md['consensus_celltype'],object)[ca['consensus_celltype_code'][pos]]
        major=oldmajor.copy();fine=oldfine.copy();review=np.zeros(len(pos),np.uint8);alias=oldfine=='pDC';fine[alias]='pDC_related';review[alias]=5;changed=np.zeros(len(pos),bool)
        raw=None
        if ds in pidx:
            pm=load(pidx[ds]['path']);fp=Path(pidx[ds]['path']).parent/pm['arrays'];assert sha(fp)==pm['arrays_SHA256']
            with np.load(fp,allow_pickle=False) as f:
                proposal_positions=f['source_native_position'];status=f['disposition_code'];gate=f['source_qualification_code']
                valid=gate==1;target=np.searchsorted(pos,proposal_positions[valid]);np.testing.assert_array_equal(pos[target],proposal_positions[valid]);assert (oldmajor[target]=='unresolved').all()
                for code in [0,1,2,3]:
                    z=target[status[valid]==code];review[z]=1 if code==0 else code+1
                    if code==0:major[z]='Myeloid';fine[z]='pDC_related';changed[z]=True
            if changed.any():
                cm=load(cidx[ds]['path']);cf=Path(cidx[ds]['path']).parent/cm['arrays'];assert sha(cf)==cm['arrays_SHA256']
                with np.load(cf,allow_pickle=False) as f:
                    rows=f['stratum_code']==0;raw=(f['source_native_position'][rows],f['native_gene_values'][rows,cm['genes'].index('LILRA4')])
        actual_major=np.asarray(m['candidate_dictionaries']['consensus_major'],object)[a['consensus_major_code']];actual_fine=np.asarray(m['candidate_dictionaries']['consensus_celltype'],object)[a['consensus_celltype_code']]
        np.testing.assert_array_equal(actual_major,major);np.testing.assert_array_equal(actual_fine,fine);np.testing.assert_array_equal(a['annotation_review_code'],review)
        source_groups=sorted([g for g in groups if g['source']==ds],key=lambda g:g['source_annotation_group_code']);assert len(source_groups)==m['n_source_annotation_groups'];assert (a['annotation_group_code'][~fit]==-1).all() and (a['annotation_group_code'][fit]>=0).all()
        seen=np.zeros(len(pos),bool)
        for group in source_groups:
            code=group['source_annotation_group_code'];mask=a['annotation_group_code']==code;members=np.flatnonzero(mask);assert len(members)==group['n_scope_group_rows'] and (fit[members]).all() and not seen[members].any();seen[members]=True
            assert set(fine[members])=={group['consensus_celltype']};parent=atlas.groups[group['parent_archive_display_number']];assert (a['baseline_global_subcluster_number'][members]==group['parent_archive_display_number']).all()
            key=exact_key(ds,s['ids'][pos[members]].astype(str).tolist());assert key==group['stable_group_sha256'];assert expected_numbers[key]==group['global_subcluster_number'];checked_keys.add(key)
            assert hashlib.sha256(np.asarray(s['rows'][pos[members]],'<i8').tobytes()).hexdigest()==group['native_rows_SHA256'];assert hist(major[members])==group['major_distribution']
            identical=key==parent['stable_group_sha256'];assert identical==group['membership_identical_to_archive_parent']
            if changed[members].any():
                assert changed[members].all() and raw is not None;hits=np.searchsorted(raw[0],pos[members]);np.testing.assert_array_equal(raw[0][hits],pos[members]);values=raw[1][hits]
                assert (values>=2).all() and group['marker']=='LILRA4' and group['marker_origin']==key and not group['parent_marker_borrowed']
                assert group['marker_evidence_type']=='own_subset_native_expression_fixed_identity_context';assert group['marker_native_mean']==float(values.mean()) and group['marker_detection_fraction']==float((values>0).mean()) and group['marker_native_minimum']==float(values.min())
            elif identical:
                for k in ['marker','marker_origin','marker_evidence_type','parent_marker_borrowed','borrow_reason']:assert group[k]==parent.get(k)
            else:
                assert group['marker']==parent['marker'] and group['marker_origin']==parent['stable_group_sha256'] and group['parent_marker_borrowed'] and group['marker_evidence_type']=='parent_group_borrowed_explicit' and group['borrow_reason']
            expected=group['consensus_celltype']+'_'+(group['marker'] or 'marker_unresolved')+'_'+str(group['global_subcluster_number']);assert expected==group['subcluster']
        np.testing.assert_array_equal(seen,fit);assert hist(major)==m['major_counts'] and hist(fine)==m['celltype_counts'] and hist(review)==m['annotation_review_counts']
        assert sum(r['n_cells'] for r in m['old_to_new_celltype_counts'])==len(pos)
        totals.update(n_sources=1,n_source_qualified_rows=len(pos),n_fit_candidates_not_final_QC=int(fit.sum()),n_without_fitted_group=int((~fit).sum()),n_context_pDC_updates=int(changed.sum()),n_existing_pDC_related_aliases=int(alias.sum()),n_exact_scoped_groups=len(source_groups))
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(totals)})
        del s,ca,ea,a,raw;atlas._cache.clear();atlas.native._fields.clear();gc.collect()
    assert checked_keys=={g['stable_group_sha256'] for g in groups} and dict(totals)==terminal['totals']
    out={'status':'COMPLETE_TME_SCOPE_ALL_ANNOTATION_ROWS_NATIVE_STATES_EXACT_GROUPS_MARKER_ORIGINS_AND_PARENT_CROSSWALKS_INDEPENDENTLY_ACCEPTED',
         'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(totals),'producer_manifest_SHA256':sha(manifest_path),
         'producer_terminal_SHA256':sha(P/'provenance/terminal_v1.json'),'producer_annotation_or_group_key_helpers_imported':False,
         'all_source_qualified_rows_and_unresolved_retained':True,'numerical_acceptance_not_identity_accuracy_or_final_module_qualification':True,
         'native_count_matrix_read':False,'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024}
    op=root/'outputs/TME_scoped_annotation_audit_v1.json';save(op,out);save(root/'provenance/terminal_v1.json',dict(out,output_SHA256=sha(op)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except BaseException:
        import traceback
        save(root/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
