"""Reconcile new native amounts with the complete frozen family policy.

The earlier family policy includes neutrophil and erythroid anchor filters.
Additional/alias-only genes absent from the new amount panel retain previously
audited binary detections here; no new native amounts are fabricated.
"""
import argparse,collections,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')

def main(root):
    begin=time.monotonic();cfg=load(root/'config/spec_v1.json');E=Path(cfg['native_evidence']);A=Path(cfg['native_audit']);F=Path(cfg['candidate_producer'])
    for stage in [E,A]:assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(E/'provenance/terminal_v1.json')
    accepted=load(cfg['working_acceptance']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];w=load(accepted['producer_manifest']['path'])
    assert sha(w['index']['path'])==w['index']['SHA256'];wi={r['source']:r for r in load(w['index']['path'])};base=Path(w['index']['path']).parent
    ref=Path(load(F/'config/full_cell_candidate_spec_v1.json')['reference']);rt=load(ref/'provenance/terminal_v1.json');ri=ref/'outputs/full_archive_reference_index_v1.json';assert sha(ri)==rt['index_sha256'];index={r['source']:r for r in load(ri)}
    sys.path.insert(0,cfg['native_review_scripts']);from native_annotation_review_v1 import marker_calls
    sys.path.insert(0,str(F/'scripts'));from native_lineage_routing_v1 import FAMILIES,native_family_support
    totals=collections.Counter();results=[];restored=collections.Counter()
    for item in load(E/'outputs/TME_full_native_evidence_index_v1.json'):
        source=item['source'];ep=E/'outputs'/item['metadata'];assert sha(ep)==item['metadata_SHA256'];e=load(ep);r=index[source];wr=wi[source]
        rp=ref/'outputs'/r['metadata'];ra=ref/'outputs'/r['arrays'];ap=base/wr['arrays'];assert sha(rp)==r['metadata_sha256'] and sha(ra)==r['arrays_sha256'] and sha(ap)==wr['arrays_SHA256'];rm=load(rp)
        c=e['native_context'];assert sha(c['path'])==c['SHA256']
        with np.load(c['path'],allow_pickle=False) as f:x=f['native_gene_values'];rows=f['native_row'];positions=f['source_native_position'];measured=f['measured']
        with np.load(ap,allow_pickle=False) as f:old=f['native_family_bitmask'];np.testing.assert_array_equal(rows,f['native_row'])
        with np.load(ra,allow_pickle=False) as f:
            np.testing.assert_array_equal(rows,f['native_row'][positions]);original=np.unpackbits(f['native_gene_detection_packbits'][positions],axis=1,bitorder='little')[:,:len(rm['review_genes'])].astype(bool)
            np.testing.assert_array_equal(old,f['native_family_bitmask'][positions])
        genes=rm['review_genes'];mask=np.asarray(rm['gene_measured'],bool);new=original.copy();context={g:j for j,g in enumerate(e['context_genes'])};fallback=[];checked=0
        for k,g in enumerate(genes):
            j=context.get(g)
            if j is not None and measured[j]:
                assert mask[k];np.testing.assert_array_equal(x[:,j]>0,original[:,k]);new[:,k]=x[:,j]>0;checked+=len(x)
            elif mask[k]:fallback.append(g);restored[g]+=1
        calls=marker_calls(new,genes,mask);by={g:k for k,g in enumerate(genes)}
        def any_gene(names):
            cols=[by[g] for g in names if g in by and mask[by[g]]]
            return new[:,cols].any(axis=1) if cols else np.zeros(len(x),bool)
        calls['Erythroid_related']&=any_gene(['AHSP','GYPA','GYPB','ALAS2','SLC4A1'])
        calls['Neutrophil_related']&=any_gene(['CSF3R','FCGR3B'])
        support=native_family_support(calls);actual=np.zeros(len(x),np.uint16)
        for j,family in enumerate(FAMILIES):actual|=np.asarray(support[family],np.uint16)<<j
        mismatches=int((actual!=old).sum());np.testing.assert_array_equal(actual,old)
        results.append({'source':source,'n_scope_records':len(x),'n_common_gene_detections_checked':checked,'additional_or_alias_only_genes_reusing_frozen_detection':fallback,'n_frozen_family_mask_mismatches':mismatches})
        totals.update(n_sources=1,n_scope_records=len(x),n_common_gene_detections_checked=checked,n_frozen_family_mask_mismatches=mismatches)
        save(root/'outputs'/(source+'_frozen_policy_replay_v1.json'),results[-1])
    expected=load(E/'provenance/terminal_v1.json')['totals'];assert totals['n_scope_records']==expected['n_current_scope_records'] and totals['n_sources']==expected['n_sources']
    report={'status':'ALL_COMMON_NATIVE_DETECTIONS_AND_COMPLETE_FROZEN_FAMILY_POLICY_RECONCILED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'totals':dict(totals),'sources':results,'additional_or_alias_only_gene_source_counts':dict(restored),
        'frozen_binary_fallback_does_not_supply_new_native_amounts':True,'higher_count_threshold_diagnostic_not_yet_qualified':True,
        'all_annotation_identity_and_state_and_QC_and_CM_unchanged':True,'wall_seconds':time.monotonic()-begin}
    op=root/'outputs/full_frozen_family_policy_replay_v1.json';save(op,report)
    save(root/'provenance/terminal_v1.json',dict(report,sources=None,output_SHA256=sha(op)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
