"""Independently replay every measured disagreement and bounded control row."""
import argparse,collections,datetime,hashlib,json,resource,time
from pathlib import Path
import numpy as np

PANELS={
 'B_specific':(['MS4A1','CD19','CD79B','PAX5','BANK1','CD22'],2),
 'T_receptor':(['CD3D','CD3E','CD247','TRAC','TRBC1','TRBC2'],2),
 'NK_specific':(['GNLY','KLRD1','NCR1','KLRF1'],2),
 'Endothelial':(['PECAM1','CDH5','CLDN5','VWF','EMCN'],2),
 'Fibroblast':(['COL1A1','COL1A2','DCN','LUM','PDGFRA','COL3A1'],3),
 'Glial':(['GFAP','AQP4','OLIG2','SOX10','PLP1'],2)}
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,ensure_ascii=False,indent=2);f.write('\n')
def truth(record,names,k):
    lower=sum(record.get(g) is not None and record[g]>=2 for g in names)
    upper=lower+sum(record.get(g) is None for g in names)
    return 1 if lower>=k else 0 if upper<k else 2
def conjunction(states):
    if 0 in states:return 0
    if 2 in states:return 2
    return 1
def hist(x):return {str(k):int(v) for k,v in sorted(collections.Counter(map(int,x)).items())}

def main(root):
    start=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);pc=load(P/'config/spec_v1.json');F=Path(pc['context'])
    assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    indexp=P/'outputs/pDC_context_proposal_index_v1.json';t=load(P/'provenance/terminal_v1.json');assert sha(indexp)==t['index_SHA256']
    source_index={r['source']:r for r in load(F/'outputs/disputed_pDC_context_index_v1.json')};totals=collections.Counter();n_controls=0
    for item in load(indexp):
        assert sha(item['path'])==item['SHA256'];m=load(item['path']);src=m['source'];original=source_index[src];assert sha(original['path'])==original['SHA256'];fm=load(original['path']);fp=Path(original['path']).parent/fm['arrays'];assert sha(fp)==fm['arrays_SHA256']
        with np.load(fp,allow_pickle=False) as f:x={k:f[k] for k in f.files}
        codes=[];specific=[];cores=[]
        for row in x['native_gene_values']:
            record={g:float(row[j]) if x['measured'][j] else None for j,g in enumerate(fm['genes'])}
            s=conjunction([truth(record,['LILRA4'],1),truth(record,['PLD4','TSPAN13'],1)])
            c=[truth(record,genes,k) for genes,k in PANELS.values()]
            c.append(conjunction([truth(record,['C1QA'],1),truth(record,['C1QB'],1),truth(record,['LYZ','LST1'],1)]))
            c.append(conjunction([truth(record,['SDC1','TNFRSF17'],1),truth(record,['MZB1','JCHAIN','DERL3','XBP1'],2)]))
            result=1 if 1 in c else 2 if s==2 or 2 in c else 0 if s==1 else 3
            codes.append(result);specific.append(s);cores.append(c)
        cases=x['stratum_code']==0;codes=np.asarray(codes,np.uint8);specific=np.asarray(specific,np.uint8);cores=np.asarray(cores,np.uint8)
        op=Path(item['path']).parent/m['arrays'];assert sha(op)==m['arrays_SHA256']
        with np.load(op,allow_pickle=False) as f:
            for k in ['source_native_position','native_row','source_qualification_code','fit_pool_candidate','parent_source_local_group_code','guarded_source_subset_code','native_family_bitmask','reference_family_bitmask']:np.testing.assert_array_equal(f[k],x[k][cases])
            np.testing.assert_array_equal(f['disposition_code'],codes[cases]);np.testing.assert_array_equal(f['additional_pDC_specific_state'],specific[cases]);np.testing.assert_array_equal(f['competing_core_state'],cores[cases])
        for j in [1,2]:assert m['bounded_controls_descriptive_not_identity_truth'][fm['frames'][j]['name']]==hist(codes[x['stratum_code']==j])
        source=x['source_qualification_code'][cases]==1;fit=x['fit_pool_candidate'][cases].astype(bool);c=codes[cases]
        assert m['scope_disposition_counts']=={'all_disputed_cases':hist(c),'source_qualified_cases':hist(c[source]),'source_qualified_fit_cases_not_final_QC':hist(c[source&fit])}
        totals.update(n_sources=1,n_case_rows=int(cases.sum()),n_context_supported_cases=int((c==0).sum()),n_source_qualified_context_supported_cases=int(((c==0)&source).sum()),n_source_qualified_fit_context_supported_cases_not_final_QC=int(((c==0)&source&fit).sum()))
        n_controls+=int((~cases).sum());save(root/'checkpoints'/(src+'.json'),{'n_sources_checked':totals['n_sources']})
    assert dict(totals)==t['totals']
    out={'status':'ALL_DISPUTED_pDC_CASES_AND_BOUNDED_CONTROLS_INDEPENDENTLY_REPLAYED_ACCEPTED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'totals':dict(totals),'n_controls_checked':n_controls,'producer_rule_helper_imported':False,
         'producer_terminal_SHA256':sha(P/'provenance/terminal_v1.json'),'all_case_arrays_and_controls_and_unknowns_checked':True,
         'not_orthogonal_biological_validation':True,'no_old_labels_or_QC_changed':True,
         'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024}
    p=root/'outputs/pDC_context_proposal_audit_v1.json';save(p,out);save(root/'provenance/terminal_v1.json',dict(out,output_SHA256=sha(p)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except BaseException:
        import traceback
        save(root/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
