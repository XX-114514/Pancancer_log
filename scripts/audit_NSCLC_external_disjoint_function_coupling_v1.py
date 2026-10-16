"""Independent sparse aggregation and full OLS for all functional couplings."""
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import rankdata

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from test_NSCLC_external_programme_coordination_v1 import read_inputs, QCOLS


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    src=root/'NSCLC_external_disjoint_function_coupling_20261016_v1'
    h=root/'HNSCC_receiver_followup_job_control_20261016_v2/describe_NSCLC_external_disjoint_function_coupling_v1'
    receipt=json.loads((h/'launch_receipt.json').read_text());deadline=time.monotonic()+3600
    while not (src/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists(),'Producer missing without terminal'
        f=p.read_text().rsplit(')',1)[1].split();assert f[0]!='Z' and f[19]==str(receipt['proc_start_ticks'])
        assert time.monotonic()<deadline,'Observation expired, do not infer termination'
        time.sleep(3)
    resource(root,1048576)
    terminal=json.loads((src/'private_terminal.json').read_text());assert terminal['source_sha256']==receipt['source_sha256']
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    assert terminal['main_coordination_results_or_clinical_values_loaded'] is False and terminal['new_P'] is False
    source,_,_,slots,cases=read_inputs(root)
    s=root/'NSCLC_external_programme_coordination_20261016_v1/private_all15600_primary_test_eligibility.tsv'
    scope=pd.read_csv(s,sep='\t');scope=scope[scope.test_index.ge(0)]
    assert sha(s)==terminal['scope_sha256']
    functional=root/'NSCLC_external_disjoint_functional_RNA_20261016_v1'
    assert sha(functional/'private_independent_disjoint_functional_RNA_audit.json')==terminal['functional_audit_sha256']
    with np.load(functional/'private_all6760_by50_disjoint_functional_RNA.npz',allow_pickle=False) as z:
        names=z['functions'].astype(str).tolist();values={k:z[k] for k in ['ULM_score','target_mean','nontarget_mean','expressed_target_genes']}
    table=pd.read_csv(src/'private_all_fixed_disjoint_function_couplings.tsv.gz',sep='\t')
    assert len(table)==len(scope)*50 and not table.duplicated(['record','reference']).any()
    table=table.set_index(['record','reference']);metadata=slots.set_index('slot')
    members=pd.read_csv(source/'private_same_library_pair_membership.tsv.gz',sep='\t')
    n=members.groupby('row_id').size().reindex(range(len(cases))).to_numpy(int)
    assert (n>0).all()
    aggregate={key:[] for key in values}
    for side in ['A','B']:
        matrix=csr_matrix((np.ones(len(members),np.int64),(members.row_id.to_numpy(int),members['unit_'+side].to_numpy(int))),shape=(len(cases),6760))
        assert matrix.max()==1
        for key,value in values.items():
            complete=(matrix@np.isfinite(value).astype(int))==n[:,None]
            expected=(matrix@np.nan_to_num(value))/n[:,None];expected[~complete]=np.nan
            aggregate[key].append(expected)
    maximum=0.;checked=0
    checkpoint=root/'HNSCC_receiver_followup_job_control_20261016_v2'/Path(__file__).stem/'private_audit_checkpoint.json'
    for item in scope.itertuples():
        meta=metadata.loc[item.slot];g=cases.iloc[int(meta.row_start):int(meta.row_stop)]
        if item.family=='reviewed_patients':g=g[g.reviewed_patient]
        frequencies=g.block_key.value_counts();g=g[g.block_key.map(frequencies).ge(2)]
        assert len(g)==item.active_units
        ids=g.row_id.to_numpy(int);blocks=g.block_key.to_numpy(str)
        quality=rankdata(g[QCOLS].to_numpy(float),axis=0);sd=quality.std(axis=0)
        quality=(quality-quality.mean(axis=0))/np.where(sd>0,sd,1)
        z=np.column_stack([np.ones(len(g)),*[(blocks==s).astype(float) for s in sorted(set(blocks))],quality])
        y=rankdata(g[['raw_A','raw_B']].to_numpy(float),axis=0)
        coefficient,_,rank,_=np.linalg.lstsq(z,y,rcond=1e-10);programme=y-z@coefficient
        assert rank==item.quality_rank and len(g)-rank==item.residual_df
        pn=np.linalg.norm(programme,axis=0)
        expected_rows=[{} for _ in range(50)]
        for key in ['ULM_score','target_mean','nontarget_mean']:
            rna=[];norm=[]
            for side in range(2):
                value=aggregate[key][side][ids];good=np.isfinite(value).all(axis=0)
                r=np.full_like(value,np.nan)
                ranks=rankdata(value[:,good],axis=0)
                beta=np.linalg.lstsq(z,ranks,rcond=1e-10)[0]
                r[:,good]=ranks-z@beta
                ns=np.linalg.norm(r,axis=0);ns[ns<1e-10]=np.nan
                rna.append(r);norm.append(ns)
            for j in range(50):
                combinations=[('programmeA_functionA',programme[:,0],rna[0][:,j],pn[0]*norm[0][j]),
                              ('programmeB_functionB',programme[:,1],rna[1][:,j],pn[1]*norm[1][j]),
                              ('programmeA_functionB',programme[:,0],rna[1][:,j],pn[0]*norm[1][j]),
                              ('programmeB_functionA',programme[:,1],rna[0][:,j],pn[1]*norm[0][j]),
                              ('functionA_functionB',rna[0][:,j],rna[1][:,j],norm[0][j]*norm[1][j])]
                for tag,a,b,den in combinations:expected_rows[j][key+'_'+tag+'_partial_rho']=float(np.clip(a@b/den,-1,1)) if np.isfinite(den) else np.nan
                expected_rows[j][key+'_same_cases_both_ends_available']=bool(np.isfinite(norm[0][j]) and np.isfinite(norm[1][j]))
        for j,name in enumerate(names):
            row=table.loc[(item.record,name)];assert row.active_units==len(g) and row.reviewed_units==g.reviewed_patient.sum()
            for side,label in enumerate(['A','B']):
                positive=aggregate['expressed_target_genes'][side][ids,j];target=aggregate['target_mean'][side][ids,j]
                good=np.isfinite(positive).all() and np.isfinite(target).all()
                expected_rows[j]['median_expressed_targets_'+label]=float(np.median(positive)) if good else np.nan
                expected_rows[j]['median_target_mean_'+label]=float(np.median(target)) if good else np.nan
                expected_rows[j]['low_foreground_support_'+label]=bool(good and (np.median(positive)<10 or np.median(target)<=0))
            for field,value in expected_rows[j].items():
                if isinstance(value,(bool,np.bool_)):assert row[field]==value,(item.record,name,field)
                else:
                    assert np.isclose(row[field],value,equal_nan=True,atol=1e-9,rtol=1e-9),(item.record,name,field)
                    if np.isfinite(value):maximum=max(maximum,abs(row[field]-value))
            checked+=1
        if item.test_index%10==0:checkpoint.write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),completed_record_index=int(item.test_index),total_records=len(scope)))+'\n')
    proof=src/'private_independent_functional_coupling_audit.json';assert not proof.exists()
    report=dict(status='PASS_ALL_ELIGIBLE_RECORDS50_FUNCTIONS_SPARSE_PATIENT_AGGREGATION_AND_FULL_OLS_COUPLINGS',
                finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(src/'private_terminal.json'),
                all_record_function_rows_checked=checked,maximum_error=maximum,main_coordination_results_or_clinical_values_loaded=False,new_P=False)
    proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
