"""Check all marker-only display changes against accepted native evidence."""
import argparse,collections,datetime,hashlib,json,re,time
from pathlib import Path

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')

def main(root):
    start=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    pc=load(P/'config/spec_v1.json');oldaccept=load(pc['previous_working_acceptance']);old=load(oldaccept['producer_manifest']['path']);E=Path(pc['native_evidence']);NA=Path(pc['native_audit'])
    assert load(NA/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    assert load(NA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(E/'provenance/terminal_v1.json')
    mp=P/'outputs/TME_exact_own_marker_working_manifest_v2.json';t=load(P/'provenance/terminal_v1.json');assert sha(mp)==t['manifest_SHA256'];new=load(mp)
    for b in new['bindings'].values():assert sha(b['path'])==b['SHA256']
    oldgroups=load(old['group_catalog']['path'])['groups'];newgroups=load(new['group_catalog']['path'])['groups'];assert len(oldgroups)==len(newgroups)
    oldby={(g['source'],g['source_annotation_group_code']):g for g in oldgroups};newby={(g['source'],g['source_annotation_group_code']):g for g in newgroups};assert set(oldby)==set(newby)
    cross=load(new['marker_display_crosswalk']['path']);assert sha(new['marker_display_crosswalk']['path'])==t['crosswalk_SHA256'];crossby={(g['source'],g['stable_group_sha256']):g for g in cross};assert len(crossby)==len(cross)==len(oldgroups)
    evid={}
    for entry in load(E/'outputs/TME_full_native_evidence_index_v1.json'):
        path=E/'outputs'/entry['metadata'];assert sha(path)==entry['metadata_SHA256'];m=load(path)
        for group in m['groups']:evid[(m['source'],group['source_annotation_group_code'])]=(group,entry['metadata_SHA256'])
    counts=collections.Counter();cells=0
    for key,g in oldby.items():
        n=newby[key];e,h=evid[key]
        for field in ['source','source_annotation_group_code','stable_group_sha256','parent_archive_group_sha256','parent_archive_display_number','parent_source_local_group_code','n_scope_group_rows','consensus_celltype','major_distribution','native_rows_SHA256','membership_identical_to_archive_parent','global_subcluster_number','global_subcluster_id']:assert n[field]==g[field],(key,field)
        candidates=e['own_native_expression_marker_candidates'];valid=[]
        for c in candidates:
            if re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9.-]*',c['gene']) and c['gene'].lower() not in ['nan','none','unknown']:valid.append(c)
        best=valid[0] if valid else None
        assert n['marker']==(best['gene'] if best else None) and n['parent_marker_borrowed'] is False and n['borrow_reason'] is None
        assert n['marker_origin']==(g['stable_group_sha256'] if best else 'unavailable')
        assert n['marker_evidence_type']==(best['marker_evidence_type'] if best else 'no_accepted_current_exact_subset_marker')
        assert n['marker_native_mean_CP10K']==(best['native_mean_CP10K'] if best else None)
        assert n['marker_detection_fraction']==(best['native_detection_fraction'] if best else None)
        assert n['marker_descriptive_log2_mean_ratio']==(best['descriptive_log2_mean_ratio'] if best else None)
        assert n['marker_expression_evidence_metadata_SHA256']==h and n['marker_evidence_membership_matches_exact_group']
        assert n['previous__marker']==g['marker'] and n['previous__subcluster']==g['subcluster']
        assert n['subcluster']==n['consensus_celltype']+'_'+(n['marker'] or 'marker_unresolved')+'_'+str(n['global_subcluster_number'])
        c=crossby[(g['source'],g['stable_group_sha256'])];assert c['old_display']==g['subcluster'] and c['new_display']==n['subcluster'] and c['n_cells_unchanged']==g['n_scope_group_rows'] and c['identity_and_membership_and_number_unchanged']
        counts.update(n_groups=1,n_previous_parent_borrowed=bool(g['parent_marker_borrowed']),n_with_current_own_marker=bool(best),n_marker_gene_changed=n['marker']!=g['marker'],n_display_changed=n['subcluster']!=g['subcluster']);cells+=g['n_scope_group_rows']
    oldindex={r['source']:r for r in load(old['index']['path'])};newindex={r['source']:r for r in load(new['index']['path'])};assert set(oldindex)==set(newindex)
    nrows=0
    for source,b in oldindex.items():
        a=newindex[source];assert a['arrays_SHA256']==b['arrays_SHA256'] and a['n_source_qualified_rows']==b['n_source_qualified_rows']
        oldarray=(Path(old['index']['path']).parent/b['arrays']).resolve();assert Path(a['arrays']).resolve()==oldarray and sha(oldarray)==a['arrays_SHA256']
        om=load(Path(old['index']['path']).parent/b['metadata']);nm=load(Path(new['index']['path']).parent/a['metadata'])
        for field in om:
            if field not in ['label_version','arrays','group_catalog']:assert nm[field]==om[field],(source,field)
        assert nm['previous_label_version']==om['label_version'] and nm['label_version']==new['label_version'] and nm['group_catalog']==new['group_catalog']
        assert nm['identity_state_and_row_arrays_unchanged'] is True;nrows+=a['n_source_qualified_rows']
    assert dict(counts)==t['counts']
    out={'status':'ALL_EXACT_MEMBERSHIPS_NUMBERS_IDENTITY_STATE_ARRAYS_AND_OWN_MARKER_DISPLAY_CROSSWALKS_ACCEPTED',
        'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'counts':dict(counts),'n_all_scope_records_preserved':nrows,'n_grouped_records_preserved':cells,
        'producer_manifest_SHA256':sha(mp),'producer_terminal_SHA256':sha(P/'provenance/terminal_v1.json'),
        'all_cell_arrays_verified_by_identical_SHA256':True,'no_new_biological_identity_or_QC_or_patient_or_CM_qualification':True,
        'inferential_DEG_not_claimed':True,'wall_seconds':time.monotonic()-start}
    op=root/'outputs/own_marker_display_refresh_audit_v1.json';save(op,out);save(root/'provenance/terminal_v1.json',dict(out,output_SHA256=sha(op)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
