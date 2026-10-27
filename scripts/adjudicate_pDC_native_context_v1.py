"""Versioned RNA-context adjudication of a previously frozen disagreement set.

Only the recorded pDC disagreement frame can receive a proposal. Shared genes
alone do not veto it; specific competing programs and unmeasured context remain
review cases. This stage changes no archived labels or analysis eligibility.
"""
import argparse,collections,datetime,hashlib,json,resource,subprocess,time
from pathlib import Path
import numpy as np
from pDC_native_context_proposal_v1 import context_proposal

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def hist(values):
    a,b=np.unique(np.asarray(values).astype(str),return_counts=True)
    return dict(zip(a.tolist(),map(int,b)))

def main(root):
    start=time.monotonic();cfg=load(root/'config/spec_v1.json');F=Path(cfg['context']);A=Path(cfg['context_audit'])
    for b in cfg['bindings'].values():assert sha(b['path'])==b['SHA256']
    assert load(F/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    assert load(A/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    test=subprocess.run([cfg['python'],'-B',str(root/'tests/test_pDC_native_context_proposal_v1.py')],capture_output=True,text=True,timeout=60)
    (root/'provenance/test_stderr_v1.txt').write_text(test.stderr)
    assert test.returncode==0 and 'Ran 8 tests' in test.stderr
    totals=collections.Counter();index=[];summary=[]
    for item in load(F/'outputs/disputed_pDC_context_index_v1.json'):
        assert sha(item['path'])==item['SHA256'];m=load(item['path']);src=m['source'];npz=Path(item['path']).parent/m['arrays'];assert sha(npz)==m['arrays_SHA256']
        with np.load(npz,allow_pickle=False) as f:inp={k:f[k] for k in f.files}
        proposed=context_proposal(inp['native_gene_values'],m['genes'],inp['measured'])
        cases=inp['stratum_code']==0
        assert int(cases.sum())==m['n_all_disputed_cases_measured']
        arrays={k:inp[k][cases] for k in ['source_native_position','native_row','source_qualification_code','fit_pool_candidate','parent_source_local_group_code','guarded_source_subset_code','native_family_bitmask','reference_family_bitmask']}
        for k in ['disposition_code','additional_pDC_specific_state','competing_core_state']:arrays[k]=proposed[k][cases]
        p=root/'outputs'/(src+'_pDC_context_proposal_v1.npz')
        with p.open('xb') as f:np.savez_compressed(f,**arrays)
        with np.load(p,allow_pickle=False) as f:
            for k,v in arrays.items():np.testing.assert_array_equal(f[k],v)
        code=arrays['disposition_code'];eligible=arrays['source_qualification_code']==1;fit=arrays['fit_pool_candidate'].astype(bool)
        scopes={'all_disputed_cases':np.ones(len(code),bool),'source_qualified_cases':eligible,'source_qualified_fit_cases_not_final_QC':eligible&fit}
        counts={name:hist(code[mask]) for name,mask in scopes.items()}
        controls={m['frames'][j]['name']:hist(proposed['disposition_code'][inp['stratum_code']==j]) for j in [1,2]}
        out={'source':src,'arrays':p.name,'arrays_SHA256':sha(p),'n_case_rows':len(code),
             'native_context_metadata_SHA256':item['SHA256'],'native_context_arrays_SHA256':m['arrays_SHA256'],
             'disposition_order':proposed['disposition_order'],'competing_core_order':proposed['competing_core_order'],
             'scope_disposition_counts':counts,'bounded_controls_descriptive_not_identity_truth':controls,
             'proposed_major_for_supported_case':'Myeloid','proposed_celltype_for_supported_case':'pDC_related',
             'Myeloid_is_an_analysis_compartment_including_DCs_not_an_ontogeny_claim':True,
             'pDC_related_does_not_exclude_AXL_SIGLEC6_related_DC_states':True,
             'no_source_QC_or_malignancy_eligibility_change':True,'original_labels_preserved':True,
             'confidence':'RNA_context_supported_proposal_not_calibrated_probability_or_orthogonal_validation'}
        mp=p.with_suffix('.json');save(mp,out);index.append({'source':src,'path':str(mp),'SHA256':sha(mp)})
        summary.append({k:out[k] for k in ['source','n_case_rows','scope_disposition_counts','bounded_controls_descriptive_not_identity_truth']})
        totals.update(n_sources=1,n_case_rows=len(code),n_context_supported_cases=int((code==0).sum()),
                      n_source_qualified_context_supported_cases=int(((code==0)&eligible).sum()),
                      n_source_qualified_fit_context_supported_cases_not_final_QC=int(((code==0)&eligible&fit).sum()))
        save(root/'checkpoints'/(src+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_completed_sources':len(index)})
    assert totals['n_case_rows']==load(F/'outputs/disputed_pDC_context_summary_v1.json')['totals']['n_disputed_case_rows']
    ip=root/'outputs/pDC_context_proposal_index_v1.json';sp=root/'outputs/pDC_context_proposal_summary_v1.json'
    save(ip,index);save(sp,{'totals':dict(totals),'sources':summary,'no_archive_labels_overwritten':True,'rule_not_selected_by_CM_or_clinical_association':True})
    save(root/'provenance/terminal_v1.json',{'status':'ALL_FROZEN_DISPUTED_pDC_CONTEXT_PROPOSALS_MATERIALIZED_INDEPENDENT_AUDIT_REQUIRED',
         'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(totals),
         'index_SHA256':sha(ip),'summary_SHA256':sha(sp),'native_expression_reread':False,'old_labels_and_CM_unchanged':True,
         'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
