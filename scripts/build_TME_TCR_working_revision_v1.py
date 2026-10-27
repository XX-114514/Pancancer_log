"""Apply audited TCR coarse proposals after exact-subset marker acceptance."""
import argparse,collections,datetime,hashlib,json,re,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def hist(v):return dict(collections.Counter(map(str,v)))

def main(root):
    start=time.monotonic();cfg=load(root/'config/spec_v1.json');S=Path(cfg['native_subsets']);A=Path(cfg['native_subset_audit']);sc=load(S/'config/spec_v1.json');P=Path(sc['producer']);pc=load(P/'config/spec_v1.json');J=Path(pc['joined'])
    for stage in [S,A,P,Path(sc['auditor'])]:assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(S/'provenance/terminal_v1.json')
    accepted=load(pc['working_acceptance']);oldmp=Path(accepted['producer_manifest']['path']);assert sha(oldmp)==accepted['producer_manifest']['SHA256'];old=load(oldmp)
    for name in ['index','group_catalog','summary']:assert sha(old[name]['path'])==old[name]['SHA256']
    sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(root/'scripts'))
    from atlas_identity_release_keys_v1 import assign_display_numbers,display_label
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    oldgroups=load(old['group_catalog']['path'])['groups'];history={g['stable_group_sha256']:g['global_subcluster_number'] for g in oldgroups}
    archived=Path(load(J/'outputs/joined_native_atlas_delivery_v1.json')['catalog']);ap=archived/'outputs/guarded_subset_identity_catalog_v1.json';assert sha(ap)==load(archived/'provenance/terminal_v1.json')['catalog_sha256']
    for g in load(ap)['groups']:
        if g['stable_group_sha256'] in history:assert history[g['stable_group_sha256']]==g['global_subcluster_number']
        history[g['stable_group_sha256']]=g['global_subcluster_number']
    subset_index=S/'outputs/TCR_changed_membership_native_index_v1.json';assert sha(subset_index)==load(S/'provenance/terminal_v1.json')['index_SHA256'];subsets={};subset_bindings={}
    for r in load(subset_index):
        p=S/'outputs'/r['metadata'];assert sha(p)==r['metadata_SHA256'];m=load(p);subsets[r['source']]=m;subset_bindings[r['source']]=r['metadata_SHA256']
    pi={r['source']:r for r in load(P/'outputs/TCR_boundary_proposal_index_v1.json')};base=Path(old['index']['path']).parent;catalog=[];metadata=[];cross=[];total=collections.Counter();sources=[]
    version='TME_source_qualified_working_annotation_20261027_v3_TCR_context'
    for item in load(old['index']['path']):
        storage_gate(root,pc);source=item['source'];mp=base/item['metadata'];arr=base/item['arrays'];assert sha(mp)==item['metadata_SHA256'] and sha(arr)==item['arrays_SHA256'];m=load(mp)
        with np.load(arr,allow_pickle=False) as f:q={k:f[k] for k in f.files}
        before=np.asarray(m['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code']];groups=sorted([dict(g) for g in oldgroups if g['source']==source],key=lambda g:g['source_annotation_group_code']);added=[];nupdates=0
        if source in subsets:
            r=pi[source];p=P/'outputs'/r['metadata'];assert sha(p)==r['metadata_SHA256'];pm=load(p);pp=P/'outputs'/pm['arrays'];assert sha(pp)==pm['arrays_SHA256']
            with np.load(pp,allow_pickle=False) as f:proposed=f['proposed_T'];np.testing.assert_array_equal(q['native_row'],f['native_row'])
            for field in ['consensus_major','consensus_celltype']:
                if 'T' not in m['candidate_dictionaries'][field]:m['candidate_dictionaries'][field].append('T')
                q[field+'_code'][proposed]=m['candidate_dictionaries'][field].index('T')
            assert len(m['candidate_dictionaries']['annotation_review'])==6;m['candidate_dictionaries']['annotation_review'].append('additional_TCR_core_reference_local_RNA_coarse_proposal_applied')
            q['annotation_review_code'][proposed]=6;nupdates=int(proposed.sum());nold=len(groups)
            pairs=collections.defaultdict(dict)
            for s in subsets[source]['groups']:pairs[s['previous_group_code']][s['RNA_coarse_proposal']]=s
            for parent,roles in sorted(pairs.items()):
                oldg=groups[parent];assert roles.keys()=={'T','unresolved'} and roles['T']['previous_group_SHA256']==oldg['stable_group_sha256'];oldsel=q['annotation_group_code']==parent
                for label in ['unresolved','T']:
                    s=roles[label];mask=oldsel&(proposed if label=='T' else ~proposed);assert int(mask.sum())==s['n_exact_subset'];g=dict(oldg)
                    for name in list(g):
                        if name=='marker' or name.startswith('marker_') or name in ['parent_marker_borrowed','borrow_reason','subcluster']:g['previous_working__'+name]=g.pop(name)
                    valid=[v for v in s['own_full_native_gene_marker_candidates_not_DEG'] if re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9.-]*',v['gene']) and v['gene'].lower() not in ['nan','none','unknown']];best=valid[0] if valid else None
                    code=parent if label=='unresolved' else nold+len(added)
                    g.update(source_annotation_group_code=code,stable_group_sha256=s['stable_membership_SHA256'],n_scope_group_rows=s['n_exact_subset'],consensus_celltype=label,major_distribution={label:s['n_exact_subset']},native_rows_SHA256=s['native_rows_SHA256'],
                        membership_identical_to_archive_parent=s['stable_membership_SHA256']==oldg['parent_archive_group_sha256'],previous_working_group_SHA256=oldg['stable_group_sha256'],previous_working_group_number=oldg['global_subcluster_number'],
                        marker=best['gene'] if best else None,marker_origin=s['stable_membership_SHA256'] if best else 'unavailable',marker_evidence_type=best['marker_evidence_type'] if best else 'no_accepted_current_exact_subset_marker',
                        marker_native_mean_CP10K=best['native_mean_CP10K'] if best else None,marker_detection_fraction=best['native_detection_fraction'] if best else None,marker_descriptive_log2_mean_ratio=best['descriptive_log2_mean_ratio'] if best else None,
                        marker_comparator_scope=s['comparator_scope'],marker_comparator_cell_count=s['n_comparator_cells'],marker_expression_evidence_metadata_SHA256=subset_bindings[source],marker_evidence_membership_matches_exact_group=True,
                        marker_selection_rule='accepted_changed_subset_all_native_gene_descriptive_ranking_not_DEG',marker_DEG_status='no_qualified_replicate_aware_DEG_supplied',marker_expression_unit_support_status='not_newly_estimated_no_patient_independence_claim',marker_small_group_status='small_group_descriptive_only' if s['n_exact_subset']<30 else 'descriptive_expression_not_inferential_DEG',parent_marker_borrowed=False,borrow_reason=None)
                    q['annotation_group_code'][mask]=code
                    if label=='T':added.append(g)
                    else:groups[parent]=g
                    cross.append({'source':source,'old_group_SHA256':oldg['stable_group_sha256'],'old_group_number':oldg['global_subcluster_number'],'old_display':oldg['subcluster'],'new_group_SHA256':g['stable_group_sha256'],'n_records':s['n_exact_subset'],'relation':'split','new_celltype':label})
            out=root/'outputs'/(source+'_TME_working_annotation_v3.npz')
            with out.open('xb') as f:np.savez_compressed(f,**q)
            with np.load(out,allow_pickle=False) as f:
                for k,v in q.items():np.testing.assert_array_equal(f[k],v)
            arrays_path=out
        else:arrays_path=arr.resolve()
        groups+=added
        changed_old={r['previous_working_group_SHA256'] for r in groups if 'previous_working_group_SHA256' in r}
        for g in groups:
            if g['stable_group_sha256'] not in changed_old and 'previous_working_group_SHA256' not in g:
                cross.append({'source':source,'old_group_SHA256':g['stable_group_sha256'],'old_group_number':g['global_subcluster_number'],'old_display':g['subcluster'],'new_group_SHA256':g['stable_group_sha256'],'n_records':g['n_scope_group_rows'],'relation':'same_membership','new_celltype':g['consensus_celltype']})
        catalog+=groups;major=np.asarray(m['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']];fine=np.asarray(m['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code']]
        m.update(previous_label_version=m['label_version'],label_version=version,arrays=str(arrays_path),arrays_SHA256=sha(arrays_path),n_source_annotation_groups=len(groups),n_TCR_coarse_updates=nupdates,major_counts=hist(major),celltype_counts=hist(fine),annotation_review_counts=hist(q['annotation_review_code']))
        m['identity_state_and_row_arrays_unchanged']=nupdates==0
        m['RNA_state_native_rows_and_source_gate_unchanged']=True
        m['old_to_new_celltype_counts']=[{'previous_working_celltype':a,'current_celltype':b,'n_records':n} for (a,b),n in sorted(collections.Counter(zip(before,fine)).items())]
        metadata.append(m);total.update(n_sources=1,n_source_qualified_rows=len(major),n_fit_candidates_not_final_QC=int(q['fit_pool_candidate_not_final_QC'].sum()),n_without_fitted_group=int((~q['fit_pool_candidate_not_final_QC']).sum()),n_TCR_coarse_updates=nupdates,n_exact_scoped_groups=len(groups));sources.append({'source':source,'major_counts':hist(major),'celltype_counts':hist(fine)})
    assigned=assign_display_numbers([g['stable_group_sha256'] for g in catalog],history)
    for g in catalog:g['global_subcluster_number']=assigned['current'][g['stable_group_sha256']];g['global_subcluster_id']='TME_annotation_group_'+str(g['global_subcluster_number']);g['subcluster']=display_label(g['consensus_celltype'],g['marker'],g['global_subcluster_number'])
    by={g['stable_group_sha256']:g for g in catalog}
    for r in cross:r.update(new_group_number=by[r['new_group_SHA256']]['global_subcluster_number'],new_display=by[r['new_group_SHA256']]['subcluster'])
    gp=root/'outputs/TME_working_group_catalog_v3.json';xp=root/'outputs/TME_annotation_membership_crosswalk_v3.json';rp=root/'outputs/TME_display_number_registry_v3.json';save(gp,{'groups':catalog});save(xp,cross);save(rp,assigned['registry']);index=[]
    for m in metadata:
        m['group_catalog']={'path':str(gp),'SHA256':sha(gp)};op=root/'outputs'/(m['source']+'_TME_working_annotation_v3.json');save(op,m);index.append({'source':m['source'],'metadata':op.name,'metadata_SHA256':sha(op),'arrays':m['arrays'],'arrays_SHA256':m['arrays_SHA256'],'n_source_qualified_rows':m['n_source_qualified_rows']})
    ip=root/'outputs/TME_working_annotation_index_v3.json';sp=root/'outputs/TME_working_annotation_summary_v3.json';save(ip,index);save(sp,{'totals':dict(total),'sources':sources,'final_module_input_ready':False})
    bindings=dict(old['bindings']);bindings.update({k:{'path':str(p),'SHA256':sha(p)} for k,p in {'previous_working_annotation':oldmp,'TCR_proposal_index':P/'outputs/TCR_boundary_proposal_index_v1.json','TCR_proposal_audit':Path(sc['auditor'])/'provenance/terminal_v1.json','changed_subset_index':subset_index,'changed_subset_audit':A/'provenance/terminal_v1.json'}.items()})
    manifest=dict(old,label_version=version,bindings=bindings,index={'path':str(ip),'SHA256':sha(ip)},group_catalog={'path':str(gp),'SHA256':sha(gp)},summary={'path':str(sp),'SHA256':sha(sp)},annotation_membership_crosswalk={'path':str(xp),'SHA256':sha(xp)},display_number_registry={'path':str(rp),'SHA256':sha(rp)},release_status='WORKING_TCR_COARSE_REVISION_INDEPENDENT_DELTA_AUDIT_REQUIRED_NOT_FINAL_TAXONOMY')
    manifest.pop('marker_display_crosswalk',None);mp=root/'outputs/TME_TCR_working_annotation_manifest_v3.json';save(mp,manifest)
    assert total['n_TCR_coarse_updates']==load(P/'provenance/terminal_v1.json')['totals']['n_proposed_T']
    save(root/'provenance/terminal_v1.json',{'status':'FULL_TME_TCR_COARSE_WORKING_REVISION_REQUIRES_INDEPENDENT_AUDIT','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(total),'manifest_SHA256':sha(mp),'new_CM_started':False,'final_module_input_ready':False,'wall_seconds':time.monotonic()-start})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
