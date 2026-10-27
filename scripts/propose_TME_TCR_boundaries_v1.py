"""Full current-scope TCR-core evidence sidecars over accepted native amounts."""
import argparse,collections,datetime,gc,hashlib,json,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
import psutil
from TME_TCR_boundary_policy_v1 import CORE_GENES,FAMILIES,REASONS,evaluate

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')

def main(root):
    begin=time.monotonic();cfg=load(root/'config/spec_v1.json');E=Path(cfg['native_evidence']);A=Path(cfg['native_audit']);J=Path(cfg['joined']);F=Path(cfg['candidate_producer'])
    for b in cfg['bindings'].values():assert sha(b['path'])==b['SHA256']
    for stage in [E,A,Path(cfg['frozen_policy_replay'])]:assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(E/'provenance/terminal_v1.json')
    tests=subprocess.run([sys.executable,'-B',str(root/'tests/test_TME_TCR_boundary_policy_v1.py')],capture_output=True,text=True,timeout=60)
    (root/'provenance/test_stderr_v1.txt').write_text(tests.stderr)
    assert tests.returncode==0 and 'Ran 14 tests' in tests.stderr and 'skipped' not in tests.stderr.lower()
    sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(F/'scripts'))
    from scoped_low_reference_identity_v1 import scoped_low_identity
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    accepted=load(cfg['working_acceptance']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];w=load(accepted['producer_manifest']['path'])
    assert sha(w['index']['path'])==w['index']['SHA256'];wi={r['source']:r for r in load(w['index']['path'])};base=Path(w['index']['path']).parent
    ref=Path(load(F/'config/full_cell_candidate_spec_v1.json')['reference']);rip=ref/'outputs/full_archive_reference_index_v1.json';assert sha(rip)==load(ref/'provenance/terminal_v1.json')['index_sha256'];ri={r['source']:r for r in load(rip)}
    lmroot=Path(load(J/'outputs/joined_native_atlas_delivery_v1.json')['local_review']);li={r['source']:r for r in load(lmroot/'outputs/full_local_identity_review_index_v1.json')}
    models=load(ref/'outputs/reference_models_v1.json');lowT=np.asarray([scoped_low_identity(v)[0]=='T' for v in models['immune_low']['labels']],bool)
    records=[];totals=collections.Counter();reason_totals=collections.Counter();control_totals=collections.Counter()
    for entry in load(E/'outputs/TME_full_native_evidence_index_v1.json'):
        storage_gate(root,cfg);ds=entry['source'];ep=E/'outputs'/entry['metadata'];assert sha(ep)==entry['metadata_SHA256'];e=load(ep);wr=wi[ds]
        wp=base/wr['metadata'];wa=base/wr['arrays'];assert sha(wp)==wr['metadata_SHA256'] and sha(wa)==wr['arrays_SHA256'];wm=load(wp)
        cp=e['native_context'];assert sha(cp['path'])==cp['SHA256']
        if 4*e['n_current_scope_records']*len(e['context_genes'])*8+cfg['RAM_reserve_bytes']>psutil.virtual_memory().available:raise MemoryError('Fresh TCR context RAM gate failed')
        with np.load(cp['path'],allow_pickle=False) as f:c={k:f[k] for k in ['native_row','source_native_position','native_gene_values','measured','native_total_RNA']}
        with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','consensus_major_code','native_family_bitmask','reference_family_bitmask','fit_pool_candidate_not_final_QC','parent_source_local_group_code','annotation_group_code']}
        for k in ['native_row','source_native_position']:np.testing.assert_array_equal(c[k],q[k])
        genes=e['context_genes'];detected=c['native_gene_values']>0;measured=c['measured'].copy();restored=[]
        rr=ri[ds];rp=ref/'outputs'/rr['metadata'];ra=ref/'outputs'/rr['arrays'];assert sha(rp)==rr['metadata_sha256'] and sha(ra)==rr['arrays_sha256'];rm=load(rp);positions=q['source_native_position']
        with np.load(ra,allow_pickle=False) as f:
            np.testing.assert_array_equal(q['native_row'],f['native_row'][positions]);low=f['immune_low__unique_code'][positions];referenceT=np.zeros(len(low),bool);valid=low>=0;referenceT[valid]=lowT[low[valid]]
            # Restore only a uniquely measured alias absent from the new
            # canonical-symbol amount panel. Its evidence remains binary.
            for gene in sorted({g for p in CORE_GENES.values() for g in p}):
                if gene in genes and not measured[genes.index(gene)] and gene in rm['review_genes'] and rm['gene_measured'][rm['review_genes'].index(gene)]:
                    k=rm['review_genes'].index(gene);detected[:,genes.index(gene)]=((f['native_gene_detection_packbits'][positions,k//8]>>(k%8))&1).astype(bool);measured[genes.index(gene)]=True;restored.append(gene)
        lr=li[ds];lp=lmroot/'outputs'/lr['output'];assert sha(lp)==lr['output_sha256'];lm=load(lp);local={r['source_local_group_code']:r['local_major_proposal'] in ['T','T_NK'] for r in lm['groups']}
        localT=np.asarray([local.get(int(k),False) for k in q['parent_source_local_group_code']],bool)
        major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']];unresolved=major=='unresolved'
        result=evaluate(detected,genes,measured,q['native_family_bitmask'],q['reference_family_bitmask'],referenceT,localT,q['fit_pool_candidate_not_final_QC'],c['native_total_RNA']>0,unresolved)
        proposed=result['proposed_T'];reason=collections.Counter(REASONS[int(k)] for k in result['review_code'][unresolved]);assert sum(reason.values())==int(unresolved.sum())
        control={label:{'n_records':int((major==label).sum()),'n_with_TCR_core':int(((major==label)&result['TCR_core']).sum()),'n_with_TCR_and_competing_core':int(((major==label)&result['TCR_core']&(result['competing_core_mask']!=0)).sum())} for label in ['T','NK','Myeloid','Epithelial','unresolved']}
        # These are the exact old memberships affected by a proposal. Full
        # native marker statistics for changed memberships remain a next gate.
        affected=[]
        for code in np.unique(q['annotation_group_code'][proposed]):
            oldgroup=q['annotation_group_code']==code;take=oldgroup&proposed;values=c['native_gene_values'][take];rna=c['native_total_RNA'][take]
            affected.append({'source_annotation_group_code':int(code),'n_previous_exact_group':int(oldgroup.sum()),'n_proposed_T':int(take.sum()),
                'n_remainder':int(oldgroup.sum()-take.sum()),'own_proposed_subset_fixed_marker_expression_not_DEG':{
                    g:{'mean_native_count':float(values[:,genes.index(g)].mean()),'native_mean_CP10K':float((values[:,genes.index(g)]*10000/rna).mean()),'detection_fraction':float((values[:,genes.index(g)]>0).mean())} for g in ['CD3D','CD3E']},
                'changed_membership_full_gene_marker_gate_pending':int(take.sum())!=int(oldgroup.sum())})
        side={k:v for k,v in result.items() if isinstance(v,np.ndarray)};side.update(native_row=q['native_row'],source_native_position=positions,reference_low_T=referenceT,source_local_T_context=localT,baseline_unresolved=unresolved)
        npz=root/'outputs'/(ds+'_TCR_boundary_proposals_v1.npz')
        with npz.open('xb') as f:np.savez_compressed(f,**side)
        with np.load(npz,allow_pickle=False) as f:
            for k,v in side.items():np.testing.assert_array_equal(f[k],v)
        meta={'source':ds,'n_scope_records':len(major),'n_baseline_unresolved':int(unresolved.sum()),'n_proposed_T':int(proposed.sum()),'n_TCR_core_records':int(result['TCR_core'].sum()),
            'reason_counts_in_baseline_unresolved':dict(reason),'baseline_label_control_counts_not_accuracy':control,'competing_core_measurement_complete':result['competing_core_measurement_complete'],
            'binary_only_alias_restoration':restored,'arrays':npz.name,'arrays_SHA256':sha(npz),'affected_previous_exact_groups':affected,
            'inputs':{'working_metadata_SHA256':sha(wp),'working_arrays_SHA256':sha(wa),'native_context_metadata_SHA256':sha(ep),'reference_metadata_SHA256':sha(rp),'reference_arrays_SHA256':sha(ra),'local_review_SHA256':sha(lp)},
            'label_status':'additional_working_RNA_coarse_proposal_not_applied_to_current_consensus','not_orthogonal_truth_or_new_QC_malignancy_patient_or_CM_eligibility':True}
        mp=root/'outputs'/(ds+'_TCR_boundary_proposals_v1.json');save(mp,meta);records.append({'source':ds,'metadata':mp.name,'metadata_SHA256':sha(mp),'n_proposed_T':int(proposed.sum())})
        totals.update(n_sources=1,n_scope_records=len(major),n_baseline_unresolved=int(unresolved.sum()),n_proposed_T=int(proposed.sum()),n_affected_previous_groups=len(affected),n_split_previous_groups=sum(r['n_remainder']>0 for r in affected));reason_totals.update(reason)
        for label,v in control.items():
            for k,n in v.items():control_totals[(label,k)]+=n
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(totals)})
        del c,q,detected,side,result;gc.collect()
    assert totals['n_scope_records']==load(E/'provenance/terminal_v1.json')['totals']['n_current_scope_records'] and set(wi)=={r['source'] for r in records}
    ip=root/'outputs/TCR_boundary_proposal_index_v1.json';save(ip,records)
    summary={'totals':dict(totals),'reason_counts_in_baseline_unresolved':dict(reason_totals),'baseline_control_counts_not_accuracy':{label:{k:control_totals[(label,k)] for k in ['n_records','n_with_TCR_core','n_with_TCR_and_competing_core']} for label in ['T','NK','Myeloid','Epithelial','unresolved']},
        'proposal_not_applied_to_consensus':True,'new_CM_started':False,'final_module_input_ready':False}
    sp=root/'outputs/TCR_boundary_proposal_summary_v1.json';save(sp,summary)
    save(root/'provenance/terminal_v1.json',{'status':'FULL_CURRENT_SCOPE_ADDITIONAL_TCR_COARSE_PROPOSALS_REQUIRE_INDEPENDENT_AUDIT','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'totals':dict(totals),'index_SHA256':sha(ip),'summary_SHA256':sha(sp),'n_tests_actual':14,'n_tests_skipped':0,'wall_seconds':time.monotonic()-begin,
        'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'no_current_consensus_or_counts_or_states_or_CM_changed':True})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
