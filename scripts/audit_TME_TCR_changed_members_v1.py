"""Independent subgroup count replay plus accepted-parent subtraction proof."""
import argparse,datetime,hashlib,json,re,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def key(source,ids):
    assert len(ids)==len(set(ids));h=hashlib.sha256()
    for v in ['atlas_exact_membership_v1',source,*sorted(ids)]:
        b=v.encode();h.update(len(b).to_bytes(8,'big'));h.update(b)
    return h.hexdigest()

def main(root):
    begin=time.monotonic();S=Path(load(root/'config/spec_v1.json')['producer']);assert load(S/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    sc=load(S/'config/spec_v1.json');P=Path(sc['producer']);pc=load(P/'config/spec_v1.json');E=Path(pc['native_evidence']);J=Path(pc['joined']);NA=Path(pc['native_audit'])
    assert load(NA/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(NA/'scripts'))
    from native_archive_collection_v1 import NativeArchiveCollection
    from audit_full_TME_native_evidence_v1 import gather
    native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest'])
    w=load(load(pc['working_acceptance'])['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])}
    ei={r['source']:r for r in load(E/'outputs/TME_full_native_evidence_index_v1.json')};pi={r['source']:r for r in load(P/'outputs/TCR_boundary_proposal_index_v1.json')};ip=S/'outputs/TCR_changed_membership_native_index_v1.json';st=load(S/'provenance/terminal_v1.json');assert sha(ip)==st['index_SHA256']
    ns=nr=nt=ng=0
    for item in load(ip):
        ds=item['source'];mp=S/'outputs'/item['metadata'];assert sha(mp)==item['metadata_SHA256'];m=load(mp);sp=m['statistics'];assert sha(sp['path'])==sp['SHA256']
        with np.load(sp['path'],allow_pickle=False) as f:a={k:f[k] for k in f.files}
        em=E/'outputs'/ei[ds]['metadata'];assert sha(em)==ei[ds]['metadata_SHA256'];e=load(em);cols=np.asarray(e['native_columns']);genes=e['native_measured_gene_names'];assert genes==m['native_measured_genes'];np.testing.assert_array_equal(a['native_columns'],cols)
        wr=wi[ds];wa=base/wr['arrays'];assert sha(wa)==wr['arrays_SHA256'];
        with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','annotation_group_code','native_total_RNA']}
        pm=P/'outputs'/pi[ds]['metadata'];assert sha(pm)==pi[ds]['metadata_SHA256'];p=load(pm);pp=P/'outputs'/p['arrays'];assert sha(pp)==p['arrays_SHA256']
        with np.load(pp,allow_pickle=False) as f:proposed=f['proposed_T']
        parents=sorted(r['source_annotation_group_code'] for r in p['affected_previous_exact_groups']);take=np.flatnonzero(np.isin(q['annotation_group_code'],parents));codes=np.searchsorted(parents,q['annotation_group_code'][take])*2+proposed[take].astype(int)
        for name,v in [('scope_position',take),('measurement_subset_code',codes),('native_row',q['native_row'][take]),('source_native_position',q['source_native_position'][take])]:np.testing.assert_array_equal(a[name],v)
        n=len(m['groups']);assert n==2*len(parents);expected={'n':np.zeros(n,np.int64),**{k:np.zeros((n,len(cols))) for k in ['cp_sum','log_sum','detected']}}
        tpos=np.flatnonzero(proposed);csr=native.store(ds)[native.legacy_config['counts_key'] if native.is_legacy(ds) else native.primary.partitions[ds]['counts_key']]
        for lo in range(0,len(tpos),1024):
            use=tpos[lo:lo+1024];raw=gather(csr,q['native_row'][use]);rna=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(rna,q['native_total_RNA'][use]);which=np.searchsorted(parents,q['annotation_group_code'][use])*2+1;x=raw[:,cols].tocsr()
            for code in np.unique(which):
                sub=x[which==code];den=rna[which==code];cp=sub.multiply((10000/den)[:,None]).tocsr();logged=cp.copy();logged.data=np.log1p(logged.data)
                expected['n'][code]+=len(den)
                for name,values in [('cp_sum',cp),('log_sum',logged),('detected',sub>0)]:expected[name][code]+=np.asarray(values.sum(axis=0)).ravel()
        old=e['group_statistics'];assert sha(old['path'])==old['SHA256']
        with np.load(old['path'],allow_pickle=False) as f:
            for j,parent in enumerate(parents):
                for k in expected:expected[k][j*2]=f[k][parent]-expected[k][j*2+1]
        for k in ['n','detected']:np.testing.assert_array_equal(a[k],expected[k])
        for k in ['cp_sum','log_sum']:np.testing.assert_allclose(a[k],expected[k],rtol=2e-10,atol=1e-7)
        ids=native.native_ids(ds);assert m['proposal_metadata_SHA256']==sha(pm) and m['working_arrays_SHA256']==sha(wa) and m['old_full_native_evidence_metadata_SHA256']==sha(em)
        for j,g in enumerate(m['groups']):
            sel=take[codes==j];assert g['measurement_subset_code']==j and g['n_exact_subset']==len(sel)==a['n'][j]
            assert g['stable_membership_SHA256']==key(ds,ids[q['source_native_position'][sel]].astype(str).tolist())
            assert g['native_rows_SHA256']==hashlib.sha256(np.asarray(q['native_row'][sel],'<i8').tobytes()).hexdigest()
            assert g['RNA_coarse_proposal']==('T' if j%2 else 'unresolved')
            compare=np.arange(n)%2==j%2;compare[j]=False
            scope='other_affected_membership_subsets_same_role_not_full_source_T_pool'
            if not compare.any():compare=np.arange(n)!=j;scope='other_affected_membership_subsets_different_roles_fallback'
            count=int(a['n'][compare].sum());own=a['cp_sum'][j]/a['n'][j];det=a['detected'][j]/a['n'][j];other=a['cp_sum'][compare].sum(axis=0)/count if count else None;otherdet=a['detected'][compare].sum(axis=0)/count if count else None
            eligible=np.asarray([not re.match(r'^(MT-|RPL|RPS|ENSG|IG[HKL][VDJC])',v) for v in genes])&(det>=.1)&(own>0);effect=np.log2((own+.1)/(other+.1)) if count else None
            if count:
                rank=effect*np.maximum(det-otherdet,0);order=sorted(np.flatnonzero(eligible&(effect>0)&(det>otherdet)),key=lambda k:(-rank[k],genes[k]))[:20];kind='own_subset_expression_specificity_not_inferential_DEG'
            else:order=[];kind='own_subset_expression_without_comparator_not_inferential_DEG'
            if not order:
                order=sorted(np.flatnonzero(eligible),key=lambda k:(-own[k]*det[k],genes[k]))[:20]
                if count:kind='own_subset_expression_without_positive_specificity_not_inferential_DEG'
            actual=g['own_full_native_gene_marker_candidates_not_DEG'];assert [v['gene'] for v in actual]==[genes[k] for k in order]
            for k,c in zip(order,actual):assert c['native_mean_CP10K']==float(own[k]) and c['native_detection_fraction']==float(det[k]) and c['marker_evidence_type']==kind and c['descriptive_log2_mean_ratio']==(float(effect[k]) if count else None) and c['comparator_detection_fraction']==(float(otherdet[k]) if count else None)
            assert g['comparator_scope']==scope and g['n_comparator_cells']==count
        ns+=1;nr+=len(take);nt+=len(tpos);ng+=n;save(root/'outputs'/(ds+'_independent_changed_membership_native_audit_v1.json'),{'source':ds,'n_scope_records':len(take),'n_T_native_records_directly_reread':len(tpos),'n_exact_subsets':n,'all_native_gene_statistics_memberships_and_markers_checked':True});native._fields.clear()
    assert ns==st['n_affected_sources'] and nr==st['n_affected_records'] and nt==st['n_proposed_T'] and ng==st['n_new_exact_subsets']
    report={'status':'ALL_CHANGED_TCR_MEMBERSHIP_NATIVE_GENE_STATISTICS_AND_KEYS_ACCEPTED_NOT_APPLIED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'n_affected_sources':ns,'n_affected_records':nr,'n_proposed_T_directly_reread':nt,'n_new_exact_subsets':ng,'producer_terminal_SHA256':sha(S/'provenance/terminal_v1.json'),
        'independent_direct_T_counts_and_accepted_parent_subtraction_proof':True,'current_consensus_and_states_and_CM_unchanged':True,'wall_seconds':time.monotonic()-begin}
    save(root/'provenance/terminal_v1.json',report)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
