"""Full working-view delta audit using accepted proposal and subset evidence."""
import argparse,collections,datetime,hashlib,json,re,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')

def main(root):
    start=time.monotonic();U=Path(load(root/'config/spec_v1.json')['producer']);assert load(U/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    uc=load(U/'config/spec_v1.json');S=Path(uc['native_subsets']);SA=Path(uc['native_subset_audit']);sc=load(S/'config/spec_v1.json');P=Path(sc['producer']);pc=load(P/'config/spec_v1.json');J=Path(pc['joined'])
    for s in [S,SA,P,Path(sc['auditor'])]:assert load(s/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    assert load(SA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(S/'provenance/terminal_v1.json')
    mp=U/'outputs/TME_TCR_working_annotation_manifest_v3.json';assert sha(mp)==load(U/'provenance/terminal_v1.json')['manifest_SHA256'];m=load(mp)
    for b in m['bindings'].values():assert sha(b['path'])==b['SHA256']
    old=load(m['bindings']['previous_working_annotation']['path'])
    for c in [m,old]:
        for k in ['index','group_catalog','summary']:assert sha(c[k]['path'])==c[k]['SHA256']
    oi={r['source']:r for r in load(old['index']['path'])};ni={r['source']:r for r in load(m['index']['path'])};assert set(oi)==set(ni);oldbase=Path(old['index']['path']).parent;newbase=Path(m['index']['path']).parent
    oldgroups=load(old['group_catalog']['path'])['groups'];groups=load(m['group_catalog']['path'])['groups'];oldby={(g['source'],g['source_annotation_group_code']):g for g in oldgroups};history={g['stable_group_sha256']:g['global_subcluster_number'] for g in oldgroups}
    ar=Path(load(J/'outputs/joined_native_atlas_delivery_v1.json')['catalog']);ag=ar/'outputs/guarded_subset_identity_catalog_v1.json';assert sha(ag)==load(ar/'provenance/terminal_v1.json')['catalog_sha256']
    for g in load(ag)['groups']:history[g['stable_group_sha256']]=g['global_subcluster_number']
    expected=dict(history);currentkeys={g['stable_group_sha256'] for g in groups};nextnumber=max(expected.values())+1
    for k in sorted(currentkeys-set(expected)):expected[k]=nextnumber;nextnumber+=1
    assert sha(m['display_number_registry']['path'])==m['display_number_registry']['SHA256'] and load(m['display_number_registry']['path'])==expected
    sp=S/'outputs/TCR_changed_membership_native_index_v1.json';assert sha(sp)==load(S/'provenance/terminal_v1.json')['index_SHA256'];subsets={}
    for r in load(sp):
        p=S/'outputs'/r['metadata'];assert sha(p)==r['metadata_SHA256'];sm=load(p)
        for g in sm['groups']:subsets[(r['source'],g['stable_membership_SHA256'])]=(g,r['metadata_SHA256'])
    pi={r['source']:r for r in load(P/'outputs/TCR_boundary_proposal_index_v1.json')};totals=collections.Counter();expected_cross=[]
    for source,oldentry in oi.items():
        ne=ni[source]
        for entry,base in [(oldentry,oldbase),(ne,newbase)]:
            assert sha(base/entry['metadata'])==entry['metadata_SHA256'] and sha(base/entry['arrays'])==entry['arrays_SHA256']
        om=load(oldbase/oldentry['metadata']);nm=load(newbase/ne['metadata'])
        with np.load(oldbase/oldentry['arrays'],allow_pickle=False) as f:o={k:f[k] for k in f.files}
        with np.load(newbase/ne['arrays'],allow_pickle=False) as f:n={k:f[k] for k in f.files}
        assert set(o)==set(n)
        for k in o:
            if k not in ['consensus_major_code','consensus_celltype_code','annotation_review_code','annotation_group_code']:np.testing.assert_array_equal(o[k],n[k])
        pr=pi[source];pp=P/'outputs'/pr['metadata'];assert sha(pp)==pr['metadata_SHA256'];pm=load(pp);pa=P/'outputs'/pm['arrays'];assert sha(pa)==pm['arrays_SHA256']
        with np.load(pa,allow_pickle=False) as f:proposed=f['proposed_T'];np.testing.assert_array_equal(f['native_row'],o['native_row'])
        for field in ['consensus_major','consensus_celltype']:
            before=np.asarray(om['candidate_dictionaries'][field],object)[o[field+'_code']];after=np.asarray(nm['candidate_dictionaries'][field],object)[n[field+'_code']];want=before.copy();want[proposed]='T';np.testing.assert_array_equal(after,want);assert nm['major_counts' if field=='consensus_major' else 'celltype_counts']==dict(collections.Counter(want))
        review=o['annotation_review_code'].copy();review[proposed]=6;np.testing.assert_array_equal(n['annotation_review_code'],review)
        assert nm['annotation_review_counts']==dict(collections.Counter(map(str,review))) and nm['n_TCR_coarse_updates']==int(proposed.sum())
        assert nm['RNA_state_native_rows_and_source_gate_unchanged'] and nm['identity_state_and_row_arrays_unchanged']==(not proposed.any())
        expected_dict={k:list(v) for k,v in om['candidate_dictionaries'].items()}
        if proposed.any():
            for k in ['consensus_major','consensus_celltype']:
                if 'T' not in expected_dict[k]:expected_dict[k].append('T')
            expected_dict['annotation_review'].append('additional_TCR_core_reference_local_RNA_coarse_proposal_applied')
        assert nm['candidate_dictionaries']==expected_dict
        allowed={'previous_label_version','label_version','arrays','arrays_SHA256','candidate_dictionaries','n_source_annotation_groups','n_TCR_coarse_updates','major_counts','celltype_counts','annotation_review_counts','old_to_new_celltype_counts','group_catalog','identity_state_and_row_arrays_unchanged','RNA_state_native_rows_and_source_gate_unchanged'}
        for k,v in om.items():
            if k not in allowed:assert nm[k]==v,(source,k)
        sg=sorted([g for g in groups if g['source']==source],key=lambda g:g['source_annotation_group_code']);assert [g['source_annotation_group_code'] for g in sg]==list(range(len(sg)))
        for g in sg:
            code=g['source_annotation_group_code'];sel=n['annotation_group_code']==code;assert int(sel.sum())==g['n_scope_group_rows']
            assert g['native_rows_SHA256']==hashlib.sha256(np.asarray(n['native_row'][sel],'<i8').tobytes()).hexdigest()
            beforecodes=set(map(int,o['annotation_group_code'][sel]));assert len(beforecodes)==1;parent=beforecodes.pop();og=oldby[(source,parent)]
            assert g['global_subcluster_number']==expected[g['stable_group_sha256']] and g['global_subcluster_id']=='TME_annotation_group_'+str(g['global_subcluster_number'])
            if (source,g['stable_group_sha256']) in subsets:
                s,h=subsets[(source,g['stable_group_sha256'])];assert s['previous_group_SHA256']==og['stable_group_sha256'] and s['n_exact_subset']==int(sel.sum()) and s['native_rows_SHA256']==g['native_rows_SHA256']
                np.testing.assert_array_equal(sel,(o['annotation_group_code']==parent)&(proposed if s['RNA_coarse_proposal']=='T' else ~proposed))
                assert g['consensus_celltype']==s['RNA_coarse_proposal'] and g['major_distribution']=={g['consensus_celltype']:g['n_scope_group_rows']}
                choices=[v for v in s['own_full_native_gene_marker_candidates_not_DEG'] if re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9.-]*',v['gene']) and v['gene'].lower() not in ['nan','none','unknown']];best=choices[0] if choices else None
                assert g['marker']==(best['gene'] if best else None) and g['marker_origin']==(g['stable_group_sha256'] if best else 'unavailable') and g['marker_expression_evidence_metadata_SHA256']==h
                for field,k in [('marker_native_mean_CP10K','native_mean_CP10K'),('marker_detection_fraction','native_detection_fraction'),('marker_descriptive_log2_mean_ratio','descriptive_log2_mean_ratio')]:assert g[field]==(best[k] if best else None)
                relation='split'
            else:
                assert g==og and g['stable_group_sha256']==og['stable_group_sha256'];np.testing.assert_array_equal(sel,o['annotation_group_code']==parent);relation='same_membership'
            assert not g['parent_marker_borrowed'] and g['subcluster']==g['consensus_celltype']+'_'+(g['marker'] or 'marker_unresolved')+'_'+str(g['global_subcluster_number'])
            expected_cross.append({'source':source,'old_group_SHA256':og['stable_group_sha256'],'old_group_number':og['global_subcluster_number'],'old_display':og['subcluster'],'new_group_SHA256':g['stable_group_sha256'],'n_records':g['n_scope_group_rows'],'relation':relation,'new_celltype':g['consensus_celltype'],'new_group_number':g['global_subcluster_number'],'new_display':g['subcluster']})
        assert nm['n_source_annotation_groups']==len(sg);totals.update(n_sources=1,n_source_qualified_rows=len(n['native_row']),n_fit_candidates_not_final_QC=int(n['fit_pool_candidate_not_final_QC'].sum()),n_without_fitted_group=int((~n['fit_pool_candidate_not_final_QC']).sum()),n_TCR_coarse_updates=int(proposed.sum()),n_exact_scoped_groups=len(sg))
    xp=m['annotation_membership_crosswalk'];assert sha(xp['path'])==xp['SHA256'];actual=load(xp['path']);order=lambda r:(r['source'],r['new_group_SHA256']);assert sorted(actual,key=order)==sorted(expected_cross,key=order)
    assert sum(r['n_records'] for r in actual)==totals['n_fit_candidates_not_final_QC'];assert load(m['summary']['path'])['totals']==dict(totals)==load(U/'provenance/terminal_v1.json')['totals']
    save(root/'provenance/terminal_v1.json',{'status':'ALL_TME_TCR_WORKING_DELTAS_MEMBERSHIPS_NUMBERS_STATES_AND_OWN_MARKERS_ACCEPTED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'totals':dict(totals),'producer_manifest_SHA256':sha(mp),'producer_terminal_SHA256':sha(U/'provenance/terminal_v1.json'),'all_scope_deltas_and_crosswalks_checked':True,'final_module_input_ready':False,'wall_seconds':time.monotonic()-start})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
