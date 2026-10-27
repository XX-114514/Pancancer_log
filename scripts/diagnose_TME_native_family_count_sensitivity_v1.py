"""Audit frozen native family support at prespecified count thresholds.

This is a diagnostic of count sensitivity, not a classifier, contamination
estimate, doublet call, or threshold optimization. No identities are changed.
"""
import argparse,collections,csv,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')

def main(root):
    begin=time.monotonic();cfg=load(root/'config/spec_v1.json');E=Path(cfg['native_evidence']);A=Path(cfg['native_audit']);R=Path(cfg['context_replay'])
    for stage in [E,A,R]:assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(E/'provenance/terminal_v1.json')
    replay=load(R/'provenance/terminal_v1.json');assert replay['totals']['n_frozen_family_mask_mismatches']==0
    accepted=load(cfg['working_acceptance']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];m=load(accepted['producer_manifest']['path'])
    assert sha(m['index']['path'])==m['index']['SHA256'];wi={r['source']:r for r in load(m['index']['path'])};base=Path(m['index']['path']).parent
    sys.path.insert(0,cfg['native_review_scripts']);from native_annotation_review_v1 import marker_calls
    sys.path.insert(0,cfg['routing_scripts']);from native_lineage_routing_v1 import FAMILIES,native_family_support
    thresholds=(1,2,3);counts=collections.Counter();transitions=collections.Counter();sources=[];rows=0;unknowns=0
    # Two families lack complete newly measured anchor/alias quantities.
    # Keep their audited support fixed, rather than inventing amounts or
    # dropping the families. This is explicitly a ten-family diagnostic.
    held_fixed=('Erythroid','Mural')
    pop=np.asarray([bin(i).count('1') for i in range(1<<len(FAMILIES))],np.uint8)
    for item in load(E/'outputs/TME_full_native_evidence_index_v1.json'):
        source=item['source'];p=E/'outputs'/item['metadata'];assert sha(p)==item['metadata_SHA256'];e=load(p);w=wi[source]
        ap=base/w['arrays'];mp=base/w['metadata'];assert sha(ap)==w['arrays_SHA256'] and sha(mp)==w['metadata_SHA256'];meta=load(mp)
        cp=e['native_context'];assert sha(cp['path'])==cp['SHA256']
        with np.load(cp['path'],allow_pickle=False) as f:x=f['native_gene_values'];native=f['native_row'];measured=f['measured'];rna=f['native_total_RNA']
        with np.load(ap,allow_pickle=False) as f:
            np.testing.assert_array_equal(native,f['native_row']);old=f['native_family_bitmask'];major=np.asarray(meta['candidate_dictionaries']['consensus_major'],object)[f['consensus_major_code']]
        assert np.isnan(x[:,~measured]).all() and np.isfinite(x[:,measured]).all();unresolved=major=='unresolved';masks={};local=collections.Counter()
        for threshold in thresholds:
            detected=x>=threshold;calls=marker_calls(detected,e['context_genes'],measured)
            by={g:k for k,g in enumerate(e['context_genes'])};anchors=[by[g] for g in ['CSF3R','FCGR3B'] if g in by and measured[by[g]]]
            calls['Neutrophil_related']&=detected[:,anchors].any(axis=1) if anchors else np.zeros(len(x),bool)
            support=native_family_support(calls)
            for family in held_fixed:support[family]=(old&(1<<FAMILIES.index(family)))!=0
            mask=np.zeros(len(x),np.uint16)
            for j,family in enumerate(FAMILIES):mask|=np.asarray(support[family],np.uint16)<<j
            masks[threshold]=mask
            if threshold==1:np.testing.assert_array_equal(mask,old)
            else:assert not np.any(mask & ~masks[1]),'Increasing native threshold created a family'
            for scope,selection in [('all_source_qualified_records',np.ones(len(x),bool)),('major_unresolved',unresolved)]:
                hist=np.bincount(pop[mask[selection]],minlength=len(FAMILIES)+1)
                assert hist.sum()==selection.sum()
                for n,v in enumerate(hist):
                    if v:counts[(source,scope,threshold,n)]+=int(v)
            local.update({'threshold_'+str(threshold)+'_unresolved_no_family':int((unresolved&(pop[mask]==0)).sum()),
                'threshold_'+str(threshold)+'_unresolved_single_family':int((unresolved&(pop[mask]==1)).sum()),
                'threshold_'+str(threshold)+'_unresolved_multiple_families':int((unresolved&(pop[mask]>1)).sum())})
        for threshold in [2,3]:
            # Retain the original pattern so an apparent simplification is
            # traceable. A one-family result is not a replacement identity.
            selection=unresolved&(pop[old]>1);pairs=np.stack([old[selection],masks[threshold][selection]],axis=1)
            patterns,ns=np.unique(pairs,axis=0,return_counts=True)
            for (oldmask,newmask),n in zip(patterns,ns):transitions[(source,threshold,int(oldmask),int(newmask))]+=int(n)
        sources.append({'source':source,'n_scope_records':len(x),'n_unresolved_records':int(unresolved.sum()),'n_context_genes_measured':int(measured.sum()),
            'n_zero_native_RNA':int((rna==0).sum()),'counts':dict(local),'no_identity_changes':True})
        rows+=len(x);unknowns+=int(unresolved.sum())
        save(root/'outputs'/(source+'_native_count_sensitivity_v1.json'),sources[-1])
    expected=load(E/'provenance/terminal_v1.json')['totals'];assert rows==expected['n_current_scope_records'] and len(sources)==expected['n_sources']
    def table(path,fieldnames,records):
        with path.open('x',newline='') as f:w=csv.DictWriter(f,fieldnames=fieldnames,delimiter='\t');w.writeheader();w.writerows(records)
    table(root/'outputs/native_family_count_histograms_v1.tsv',['source','scope','minimum_native_counts','n_supported_families','n_records'],
        [dict(zip(['source','scope','minimum_native_counts','n_supported_families','n_records'],(*key,value))) for key,value in sorted(counts.items())])
    names=lambda mask:'+'.join(f for j,f in enumerate(FAMILIES) if mask&(1<<j)) or 'no_detected_family'
    table(root/'outputs/unresolved_multiple_family_transitions_v1.tsv',['source','minimum_native_counts','original_native_families','threshold_native_families','n_records'],
        [{'source':s,'minimum_native_counts':t,'original_native_families':names(a),'threshold_native_families':names(b),'n_records':v} for (s,t,a,b),v in sorted(transitions.items())])
    totals={}
    for t in thresholds:
        totals[str(t)]={kind:sum(r['counts']['threshold_'+str(t)+'_unresolved_'+kind] for r in sources) for kind in ['no_family','single_family','multiple_families']}
        assert sum(totals[str(t)].values())==unknowns
    report={'status':'COMPLETE_FIXED_THRESHOLD_NATIVE_SUPPORT_DIAGNOSTIC_NO_LABEL_CHANGE','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'n_sources':len(sources),'n_scope_records':rows,'n_major_unresolved':unknowns,'prespecified_native_count_thresholds':list(thresholds),
        'unresolved_support_histograms':totals,'family_axis':list(FAMILIES),'threshold_one_reproduces_all_original_masks':True,
        'families_varied_with_new_native_counts':[f for f in FAMILIES if f not in held_fixed],
        'families_kept_fixed_at_audited_binary_support':list(held_fixed),
        'fixed_family_reason':'additional erythroid anchors and one mural gene alias lack new quantities; frozen support not fabricated amounts',
        'higher_threshold_is_not_new_identity_or_ambient_or_doublet_proof':True,'missing_features_not_called_negative':True,
        'depth_dependence_requires_source_and_RNA_review':True,'no_threshold_selected_by_yield_or_clinical_association':True,
        'identity_and_QC_and_CM_unchanged':True,'wall_seconds':time.monotonic()-begin}
    rp=root/'outputs/native_family_count_sensitivity_summary_v1.json';save(rp,report)
    save(root/'provenance/terminal_v1.json',dict(report,summary_SHA256=sha(rp)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
