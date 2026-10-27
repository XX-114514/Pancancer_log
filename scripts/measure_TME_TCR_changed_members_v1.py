"""All-native-gene evidence for every membership split by accepted TCR proposals.

Unchanged groups retain their already accepted full native evidence. This
incremental pass covers both proposed T and unresolved remainder in all affected
old groups; neither current consensus nor display numbers are changed here.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from TME_native_measurement_core_v1 import canonical_native,update_group_stats,expression_candidates

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')

def main(root):
    begin=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);A=Path(cfg['auditor']);pc=load(P/'config/spec_v1.json');E=Path(pc['native_evidence']);J=Path(pc['joined'])
    for s in [P,A]:assert load(s/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json')
    sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(root/'scripts'));from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
    from native_panel_aggregate_v1 import gather_csr_rows
    from atlas_identity_release_keys_v1 import membership_key
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    accepted=load(pc['working_acceptance']);w=load(accepted['producer_manifest']['path']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])}
    oldgroups={(g['source'],g['source_annotation_group_code']):g for g in load(w['group_catalog']['path'])['groups']}
    assert sha(w['index']['path'])==w['index']['SHA256'] and sha(w['group_catalog']['path'])==w['group_catalog']['SHA256']
    assert sha(P/'outputs/TCR_boundary_proposal_index_v1.json')==load(A/'provenance/terminal_v1.json')['producer_index_SHA256']
    assert sha(E/'outputs/TME_full_native_evidence_index_v1.json')==load(E/'provenance/terminal_v1.json')['index_SHA256']
    pi=load(P/'outputs/TCR_boundary_proposal_index_v1.json');ei={r['source']:r for r in load(E/'outputs/TME_full_native_evidence_index_v1.json')};native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest'])
    output=[];nrows=0;nT=0;ngroups=0
    for entry in pi:
        if entry['n_proposed_T']==0:continue
        storage_gate(root,pc);ds=entry['source'];mp=P/'outputs'/entry['metadata'];assert sha(mp)==entry['metadata_SHA256'];m=load(mp);pp=P/'outputs'/m['arrays'];assert sha(pp)==m['arrays_SHA256']
        with np.load(pp,allow_pickle=False) as f:proposed=f['proposed_T'];positions=f['source_native_position']
        wr=wi[ds];wa=base/wr['arrays'];wm=load(base/wr['metadata']);assert sha(wa)==wr['arrays_SHA256']
        with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','annotation_group_code','native_total_RNA']}
        np.testing.assert_array_equal(positions,q['source_native_position']);em=E/'outputs'/ei[ds]['metadata'];assert sha(em)==ei[ds]['metadata_SHA256'];e=load(em);genes=e['native_measured_gene_names'];cols=np.asarray(e['native_columns'],np.int64)
        parents=sorted(r['source_annotation_group_code'] for r in m['affected_previous_exact_groups']);take=np.flatnonzero(np.isin(q['annotation_group_code'],parents));old=q['annotation_group_code'][take];codes=np.searchsorted(parents,old)*2+proposed[take].astype(np.int64);ng=2*len(parents)
        stat={'n':np.zeros(ng,np.int64),**{k:np.zeros((ng,len(genes)),np.float64) for k in ['cp_sum','log_sum','detected']}}
        np.testing.assert_array_equal(native.rows(ds)[positions[take]],q['native_row'][take]);ids=native.native_ids(ds);assert text_axis_sha256(ids)==wm['native_ids_source_SHA256']
        key=native.legacy_config['counts_key'] if native.is_legacy(ds) else native.primary.partitions[ds]['counts_key'];csr=native.store(ds)[key]
        for low in range(0,len(take),2048):
            use=take[low:low+2048];raw=canonical_native(gather_csr_rows(csr,q['native_row'][use],max_span=8192));rna=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(rna,q['native_total_RNA'][use]);update_group_stats(stat,raw[:,cols].tocsr(),rna,codes[low:low+2048])
        assert (stat['n']>0).all();np.testing.assert_array_equal(stat['n'],np.bincount(codes,minlength=ng))
        oldstat=e['group_statistics'];assert sha(oldstat['path'])==oldstat['SHA256']
        with np.load(oldstat['path'],allow_pickle=False) as f:
            np.testing.assert_array_equal(f['native_columns'],cols)
            for j,parent in enumerate(parents):
                for k in ['n','detected']:np.testing.assert_array_equal(stat[k][j*2:j*2+2].sum(axis=0),f[k][parent])
                for k in ['cp_sum','log_sum']:np.testing.assert_allclose(stat[k][j*2:j*2+2].sum(axis=0),f[k][parent],rtol=2e-10,atol=1e-7)
        groups=[]
        for j in range(ng):
            selector=take[codes==j];parent=parents[j//2];oldgroup=oldgroups[(ds,parent)];label='T' if j%2 else 'unresolved';memberids=ids[positions[selector]].astype(str).tolist();stable=membership_key(ds,memberids)
            same_role=np.arange(ng)%2==j%2;same_role[j]=False
            if not same_role.any():same_role=np.arange(ng)!=j;scope='other_affected_membership_subsets_different_roles_fallback'
            else:scope='other_affected_membership_subsets_same_role_not_full_source_T_pool'
            nc=int(stat['n'][same_role].sum());othercp=stat['cp_sum'][same_role].sum(axis=0)/nc if nc else None;otherdet=stat['detected'][same_role].sum(axis=0)/nc if nc else None
            candidates=expression_candidates(genes,stat['cp_sum'][j]/stat['n'][j],stat['detected'][j]/stat['n'][j],othercp,otherdet,.1,20)
            groups.append({'measurement_subset_code':j,'RNA_coarse_proposal':label,'n_exact_subset':int(stat['n'][j]),'stable_membership_SHA256':stable,
                'native_rows_SHA256':hashlib.sha256(np.asarray(q['native_row'][selector],'<i8').tobytes()).hexdigest(),'previous_group_code':parent,'previous_group_number':oldgroup['global_subcluster_number'],'previous_group_SHA256':oldgroup['stable_group_sha256'],
                'own_full_native_gene_marker_candidates_not_DEG':candidates,'comparator_scope':scope,'n_comparator_cells':nc,'new_display_number_not_yet_assigned':True})
        op=root/'outputs'/(ds+'_changed_membership_full_gene_statistics_v1.npz')
        with op.open('xb') as f:np.savez_compressed(f,**stat,native_columns=cols,scope_position=take,native_row=q['native_row'][take],source_native_position=positions[take],measurement_subset_code=codes)
        with np.load(op,allow_pickle=False) as f:
            for k,v in stat.items():np.testing.assert_array_equal(f[k],v)
        meta={'source':ds,'n_affected_old_groups':len(parents),'n_new_exact_subsets':ng,'n_affected_records':len(take),'n_proposed_T':int(proposed[take].sum()),'native_measured_genes':genes,'groups':groups,
            'statistics':{'path':str(op),'SHA256':sha(op)},'proposal_metadata_SHA256':sha(mp),'working_arrays_SHA256':sha(wa),'old_full_native_evidence_metadata_SHA256':sha(em),
            'old_full_parent_statistics_recombined_for_every_gene':True,'independent_subgroup_native_replay_required':True,'current_annotation_not_applied':True}
        om=root/'outputs'/(ds+'_changed_membership_full_gene_evidence_v1.json');save(om,meta);output.append({'source':ds,'metadata':om.name,'metadata_SHA256':sha(om)});nrows+=len(take);nT+=meta['n_proposed_T'];ngroups+=ng;native._fields.clear()
    assert nT==load(P/'provenance/terminal_v1.json')['totals']['n_proposed_T']
    ip=root/'outputs/TCR_changed_membership_native_index_v1.json';save(ip,output)
    save(root/'provenance/terminal_v1.json',{'status':'ALL_TCR_CHANGED_MEMBERSHIP_FULL_GENE_STATISTICS_MEASURED_INDEPENDENT_REPLAY_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'n_affected_sources':len(output),'n_affected_records':nrows,'n_proposed_T':nT,'n_new_exact_subsets':ngroups,'index_SHA256':sha(ip),'wall_seconds':time.monotonic()-begin,'current_consensus_counts_states_and_CM_unchanged':True})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
