"""Refresh display markers from accepted exact-membership native evidence.

Identity, cell membership, group numbers, states and underlying arrays are
unchanged. This is a working display release, not new biological or CM approval.
"""
import argparse,collections,datetime,hashlib,json,re,time
from pathlib import Path

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')

def main(root):
    start=time.monotonic();cfg=load(root/'config/spec_v1.json');E=Path(cfg['native_evidence']);A=Path(cfg['native_audit'])
    for stage in [E,A]:assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    at=load(A/'provenance/terminal_v1.json');assert at['producer_terminal_SHA256']==sha(E/'provenance/terminal_v1.json')
    accepted=load(cfg['previous_working_acceptance']);oldmp=Path(accepted['producer_manifest']['path']);assert sha(oldmp)==accepted['producer_manifest']['SHA256'];old=load(oldmp)
    for name in ['index','group_catalog']:assert sha(old[name]['path'])==old[name]['SHA256']
    oldgroups=load(old['group_catalog']['path'])['groups'];eindex=load(E/'outputs/TME_full_native_evidence_index_v1.json');groups=[];cross=[];counts=collections.Counter();evidence={}
    for item in eindex:
        p=E/'outputs'/item['metadata'];assert sha(p)==item['metadata_SHA256'];m=load(p)
        for r in m['groups']:evidence[(m['source'],r['source_annotation_group_code'])]=(r,item['metadata_SHA256'])
    for g in oldgroups:
        r,h=evidence[(g['source'],g['source_annotation_group_code'])];assert r['stable_group_sha256']==g['stable_group_sha256'] and r['n_exact_group_cells']==g['n_scope_group_rows']
        candidates=r['own_native_expression_marker_candidates']
        display_safe=lambda gene: bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9.-]*',gene)) and gene.lower() not in ['nan','none','unknown']
        best=next((v for v in candidates if display_safe(v['gene'])),None);new=dict(g)
        for name in ['marker','marker_origin','marker_evidence_type','parent_marker_borrowed','borrow_reason','subcluster','marker_selection_rule','marker_native_mean','marker_detection_fraction','marker_native_minimum']:
            if name in new:new['previous__'+name]=new.pop(name)
        new.update(marker=best['gene'] if best else None,marker_origin=g['stable_group_sha256'] if best else 'unavailable',
            marker_evidence_type=best['marker_evidence_type'] if best else 'no_accepted_current_exact_subset_marker',
            parent_marker_borrowed=False,borrow_reason=None,
            marker_selection_rule='accepted_current_exact_group_full_native_gene_descriptive_ranking_not_DEG',
            marker_native_mean_CP10K=best['native_mean_CP10K'] if best else None,
            marker_detection_fraction=best['native_detection_fraction'] if best else None,
            marker_descriptive_log2_mean_ratio=best['descriptive_log2_mean_ratio'] if best else None,
            marker_comparator_scope=r['comparator_scope'],marker_comparator_cell_count=r['n_comparator_cells'],
            marker_expression_evidence_metadata_SHA256=h,marker_evidence_membership_matches_exact_group=True,
            marker_display_excludes_unsafe_or_missing_symbol_tokens=True,
            marker_expression_unit_support_status='not_newly_estimated_for_current_marker_no_patient_independence_claim',
            marker_DEG_status='no_qualified_replicate_aware_DEG_input_supplied',
            marker_small_group_status='small_group_descriptive_only' if g['n_scope_group_rows']<30 else 'descriptive_expression_not_inferential_DEG')
        new['subcluster']=g['consensus_celltype']+'_'+(new['marker'] or 'marker_unresolved')+'_'+str(g['global_subcluster_number'])
        groups.append(new);counts.update(n_groups=1,n_previous_parent_borrowed=bool(g['parent_marker_borrowed']),n_with_current_own_marker=bool(best),n_marker_gene_changed=new['marker']!=g['marker'],n_display_changed=new['subcluster']!=g['subcluster'])
        cross.append({'source':g['source'],'stable_group_sha256':g['stable_group_sha256'],'global_subcluster_number':g['global_subcluster_number'],
            'old_display':g['subcluster'],'new_display':new['subcluster'],'old_marker':g['marker'],'new_marker':new['marker'],
            'n_cells_unchanged':g['n_scope_group_rows'],'identity_and_membership_and_number_unchanged':True,
            'relation':'same_membership_display_changed' if new['subcluster']!=g['subcluster'] else 'same_membership_display_unchanged'})
    gp=root/'outputs/TME_exact_own_marker_group_catalog_v2.json';xp=root/'outputs/TME_marker_display_crosswalk_v2.json';save(gp,{'groups':groups,'marker_only_working_revision':True,'biological_types_not_created_from_markers_or_numbers':True});save(xp,cross)
    version='TME_source_qualified_working_annotation_20261027_v2_exact_own_markers';index=[];base=Path(old['index']['path']).parent
    for item in load(old['index']['path']):
        p=base/item['metadata'];assert sha(p)==item['metadata_SHA256'];m=load(p);arrays=(base/item['arrays']).resolve();assert sha(arrays)==item['arrays_SHA256']
        m['previous_label_version']=m['label_version'];m['label_version']=version;m['arrays']=str(arrays);m['group_catalog']={'path':str(gp),'SHA256':sha(gp)}
        m['marker_display_refreshed_from_accepted_current_exact_native_evidence']=True;m['identity_state_and_row_arrays_unchanged']=True
        op=root/'outputs'/(item['source']+'_TME_scoped_annotation_v2.json');save(op,m)
        index.append(dict(item,metadata=op.name,metadata_SHA256=sha(op),arrays=str(arrays)))
    ip=root/'outputs/TME_scoped_annotation_index_v2.json';save(ip,index)
    bindings=dict(old['bindings']);bindings.update(previous_annotation_manifest={'path':str(oldmp),'SHA256':sha(oldmp)},accepted_native_evidence_terminal={'path':str(E/'provenance/terminal_v1.json'),'SHA256':sha(E/'provenance/terminal_v1.json')},accepted_native_evidence_audit_terminal={'path':str(A/'provenance/terminal_v1.json'),'SHA256':sha(A/'provenance/terminal_v1.json')})
    manifest=dict(old,bindings=bindings,index={'path':str(ip),'SHA256':sha(ip)},group_catalog={'path':str(gp),'SHA256':sha(gp)},label_version=version,
        marker_display_crosswalk={'path':str(xp),'SHA256':sha(xp)},release_status='WORKING_MARKER_REFRESH_INDEPENDENT_DISPLAY_AUDIT_REQUIRED_NOT_MODULE_FIT_READY')
    mp=root/'outputs/TME_exact_own_marker_working_manifest_v2.json';save(mp,manifest)
    save(root/'provenance/terminal_v1.json',{'status':'EXACT_CURRENT_SUBSET_OWN_MARKER_WORKING_DISPLAY_REFRESHED_INDEPENDENT_AUDIT_REQUIRED',
        'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'counts':dict(counts),'manifest_SHA256':sha(mp),'catalog_SHA256':sha(gp),'crosswalk_SHA256':sha(xp),
        'all_identity_and_state_arrays_reused_by_identical_SHA256':True,'no_cell_membership_or_numbers_or_identity_or_QC_or_CM_changed':True,
        'inferential_DEG_not_claimed':True,'wall_seconds':time.monotonic()-start})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
