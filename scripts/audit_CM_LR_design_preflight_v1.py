"""Independent QR reconstruction of every CM/LR design diagnostic, no inference."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,csv,gzip
import numpy as np,pandas as pd
from scipy.linalg import qr

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';L=P/'CM_edge_LR_extraction_v2';O=P/'CM_LR_model_design_preflight_v1';planpath=Path(os.environ['PANCANCER_CM_DESIGN_PLAN']);plan=json.loads(planpath.read_text());assert sha(__file__)==plan['code_sha256'][Path(__file__).name] and not (O/'independent_audit.json').exists();inv=json.loads((O/'inventory.json').read_text());assert inv['plan_sha256']==sha(planpath)
    for n,h in inv['outputs'].items():assert sha(O/n)==h
    assert sha(L/'inventory.json')==plan['inventory_sha256'] and sha(L/'independent_audit.json')==plan['audit_sha256'] and sha(L/'all314432_consensus_case_coverage.tsv.gz')==plan['case_table_sha256']
    li=json.loads((L/'inventory.json').read_text());groups={};cases={}
    with gzip.open(L/'all314432_consensus_case_coverage.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['data_coverage_gate_only']=='True':cases[r['case_id']]=r
    assert len(cases)==5377
    for gid in sorted({r['group_id'] for r in cases.values()}):
        p=L/'private_group_arrays'/(gid+'.npz');assert sha(p)==li['outputs'][str(p.relative_to(L))];groups[gid]=dict(np.load(p,allow_pickle=False))
    checked=0;seen=set();signatures={};names=['CM_activity','source_within_major','target_within_major','log_source_cells','log_target_cells'];counts={}
    with gzip.open(O/'all40491_model_design_diagnostics.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            case=cases[r['case_id']];b=groups[r['group_id']];cm=r['cm'];variant=r['variant'];key=(r['case_id'],cm,variant);assert key not in seen;seen.add(key);assert cm in case['CMs'].split(';')
            for field in ['group_id','pair_id','stratum','gse_id','cancerType','source','target','ligand','receptor']:assert r[field]==case[field]
            j=int(r['case_id'].rsplit('_LR',1)[1])-1;assert b['ligands'][j]==r['ligand'] and b['receptors'][j]==r['receptor'];mask=np.isfinite(b['scores'][:,j,0]);ids=np.flatnonzero(mask);blocks=b['patient_blocks'][mask];review=b['reviewed_patient'][mask];Y=np.log1p(b['scores'][mask,j,0]);ci=list(b['CMs']).index(cm);field='original_H' if variant=='original_H_library' else 'endpoint_excluded_activity';raw=np.array([[b[field][i,ci],b['source_within_major'][i],b['target_within_major'][i],np.log1p(b['source_n'][i]),np.log1p(b['target_n'][i])] for i in ids]);yy=Y;bb=blocks;rr=review
            if variant=='endpoint_excluded_block_mean_diagnostic':
                order=list(dict.fromkeys(blocks));raw=np.array([sum(raw[i] for i,k in enumerate(blocks) if k==block)/sum(blocks==block) for block in order]);yy=np.array([sum(Y[i] for i,k in enumerate(blocks) if k==block)/sum(blocks==block) for block in order]);rr=np.array([all(review[blocks==block]) for block in order]);bb=np.array(order)
            assert int(r['n_reviewed_rows'])==int(rr.sum()) and int(r['n_reviewed_source_patient_blocks'])==len(set(bb[rr]))
            sig=(r['group_id'],np.packbits(mask).tobytes(),cm,variant)
            if sig not in signatures:
                n=len(raw);G=len(set(bb));d=dict(n_rows=n,n_analysis_blocks=G)
                if np.isfinite(raw).all():
                    mean=np.array([sum(raw[:,i])/n for i in range(5)]);sd=np.array([np.sqrt(sum((raw[:,i]-mean[i])**2)/(n-1)) for i in range(5)]);const=sd<=1e-12;x=np.column_stack([np.ones(n)]+[(raw[:,i]-mean[i])/sd[i] if not const[i] else np.zeros(n) for i in range(5)]);Q,rx,_=qr(x,mode='economic',pivoting=True);tol=np.finfo(float).eps*max(x.shape)*abs(rx[0,0]);rank=int((abs(np.diag(rx))>tol).sum());Z=x[:,[0,2,3,4,5]];Qz,rz,_=qr(Z,mode='economic',pivoting=True);ztol=np.finfo(float).eps*max(Z.shape)*abs(rz[0,0]);zr=int((abs(np.diag(rz))>ztol).sum());v=x[:,1]-Qz[:,:zr]@(Qz[:,:zr].T@x[:,1]);resid=float(v@v);lev=float((Q[:,:rank]**2).sum(axis=1).max());cond=float(np.linalg.cond(x)) if rank==6 else None;unique=not const[0] and resid>1e-10*max(1.,float(x[:,1]@x[:,1]))
                    if n<20 or G<10:status='INSUFFICIENT20_ROWS_OR10_BLOCKS'
                    elif const[0]:status='CONSTANT_CM_ACTIVITY'
                    elif const[1:].any():status='CONSTANT_REQUIRED_NUISANCE_COLUMN'
                    elif rank<6:status='COLLINEAR_REQUIRED_DESIGN'
                    elif cond>=1e8:status='ILL_CONDITIONED_REQUIRED_DESIGN'
                    elif lev>.95:status='HIGH_INFLUENCE_DESIGN'
                    else:status='STRICT_DESIGN_SCREEN_PASSED_NOT_INFERENCE'
                    d.update(n_repeated_library_rows=n-G,design_columns=6,design_rank=rank,nuisance_span_rank=zr,constant_columns=';'.join(np.array(names)[const]),design_condition=cond,max_leverage=lev,CM_residual_sum_squares_given_all_nuisance_columns=resid,CM_coefficient_identifiable_in_full_nuisance_span=unique,status=status)
                else:d['status']='CM_ACTIVITY_UNAVAILABLE'
                signatures[sig]=d
            d=signatures[sig].copy();response=np.std(yy,ddof=1)>1e-12
            if not response:d['status']='NO_RESPONSE_VARIATION'
            assert r['response_variable']==str(bool(response)) and r['status']==d['status'],key
            for field in ['n_rows','n_analysis_blocks','n_repeated_library_rows','design_columns','design_rank','nuisance_span_rank']:
                if field in d:assert int(r[field])==d[field],(key,field)
            if 'constant_columns' in d:
                assert r['constant_columns']==d['constant_columns'] and r['CM_coefficient_identifiable_in_full_nuisance_span']==str(bool(d['CM_coefficient_identifiable_in_full_nuisance_span']))
                for field in ['design_condition','max_leverage','CM_residual_sum_squares_given_all_nuisance_columns']:
                    if d[field] is None:assert r[field]=='NA'
                    else:assert np.isclose(float(r[field]),d[field],rtol=1e-6,atol=1e-9),(key,field)
            counts.setdefault(variant,{});counts[variant][r['status']]=counts[variant].get(r['status'],0)+1;checked+=1
    assert checked==40491
    for v,summary in inv['summary'].items():assert counts[v]==summary['status_counts']
    result=dict(status='PASS_ALL40491_COMPLETE_DESIGNS_QR_NUISANCE_SPAN_LIBRARY_AND_BLOCK_MEAN_RECONSTRUCTION',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),diagnostic_rows_replayed=checked,distinct_QR_templates=len(signatures),new_coefficients=False,new_p_values=False,strict_primary_rank_policy_preserved=True,partial_identifiability_is_diagnostic_only=True,clinical_blocks_changed=False)
    (O/'independent_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
