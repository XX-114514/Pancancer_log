"""Independent OLS replay of all28 conditional descriptions and all LOO."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
import statsmodels.api as sm


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def independent_rho(x,y,D=None):
    a=rankdata(x).astype(float);b=rankdata(y).astype(float)
    if D is None:return float(np.corrcoef(a,b)[0,1])
    a=sm.OLS(a,D).fit().resid;b=sm.OLS(b,D).fit().resid
    return float(np.corrcoef(a,b)[0,1]) if np.linalg.norm(a)*np.linalg.norm(b)>1e-12 else np.nan


def equal(a,b):assert np.isclose(a,b,atol=1e-8,rtol=1e-8,equal_nan=True),(a,b)


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'DLG_LR_sender_receiver_decomposition_20261016_v1'
    receipt=R/'HNSCC_receiver_followup_job_control_20261016_v2/decompose_DLG_LR_sender_receiver_v1/launch_receipt.json'
    launch=json.loads(receipt.read_text());deadline=time.monotonic()+3600
    while not (O/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists()
        st=p.read_text().rsplit(')',1)[1].split();assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z'
        assert time.monotonic()<deadline;time.sleep(10)
    out=O/'private_independent_audit.json';assert not out.exists()
    terminal=json.loads((O/'private_terminal.json').read_text());path=O/'all28_private_sender_receiver_RNA_diagnostics.tsv'
    assert sha(path)==terminal['result_sha256']
    t=pd.read_csv(path,sep='\t');loo=pd.read_csv(O/'private_all_estimable_leave_one_block_out.tsv',sep='\t')
    assert len(t)==28 and t.record.is_unique and t.groupby(['context','source','target','MCP','ligand','receptor']).size().eq(2).all()
    assert not t.new_p_values.any() and t.post_selection_diagnostic.all()
    V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    original=pd.read_csv(V/'state_after_DIALOGUE_LR_response_20261016_v1/all_estimable_exploratory_LR_receiver_associations.tsv',sep='\t')
    keys=['context','source','target','MCP','ligand','receptor'];orig=original[original.BH_q.lt(.05)].set_index(keys)
    assert set(t.set_index(keys).index)==set(orig.index)
    count=0;maximum=0.
    for row in t.itertuples():
        z=np.load(O/'private_fixed_designs'/(row.record+'.npz'),allow_pickle=False)
        X=z['X'];Y=z['Y'];cov=z['nuisance'];D=np.column_stack([np.ones(len(X)),(cov-cov.mean(0))/cov.std(0)])
        assert len(set(z['blocks']))==len(X)==row.blocks
        assert row.design_columns==D.shape[1] and row.design_rank==np.linalg.matrix_rank(D)
        assert row.residual_df==len(X)-D.shape[1]-1;equal(row.design_condition,np.linalg.cond(D))
        baseD=D[:,:9];equal(independent_rho(X,Y,baseD),row.base_side_rho)
        v=independent_rho(X,Y,D);equal(v,row.conditional_rho);maximum=max(maximum,abs(v-row.conditional_rho))
        equal(independent_rho(z['interface'],Y,D),row.interface_conditional_on_opposite_rho)
        r=orig.loc[tuple(getattr(row,k) for k in keys)]
        equal(row.original_BH_q,r.BH_q);equal(independent_rho(z['interface'],Y,baseD),r.adjusted_masked_rho)
        tr=z['train'];te=z['test']
        if tr.sum()>=20 and te.sum()>=10 and np.linalg.matrix_rank(D[tr])==D.shape[1]:
            mx=sm.OLS(X[tr],D[tr]).fit();my=sm.OLS(Y[tr],D[tr]).fit()
            equal(independent_rho(X[te]-mx.predict(D[te]),Y[te]-my.predict(D[te])),row.internal_holdout_rho)
        else:assert pd.isna(row.internal_holdout_rho)
        values=[];f=loo[loo.record.eq(row.record)].set_index('removed_block_index');assert len(f)==len(X)
        for k in range(len(X)):
            use=np.arange(len(X))!=k
            v=independent_rho(X[use],Y[use],D[use]) if np.linalg.matrix_rank(D[use])==D.shape[1] else np.nan
            equal(v,f.loc[k,'rho']);count+=1
            if np.isfinite(v):values.append(v)
        assert len(values)==row.LOO_estimable and np.sum(np.sign(values)==np.sign(row.conditional_rho))==row.LOO_same_sign
        equal(min(values),row.LOO_minimum);equal(max(values),row.LOO_maximum)
    a=dict(status='PASS_ALL28_STATS_MODEL_OLS_POINTS_ALL_INTERNAL_HOLDOUT_ALL_LEAVE_ONE_BLOCK_OUT_DESCRIPTIONS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_result_sha256=sha(path),audit_code_sha256=sha(__file__),records=len(t),leave_one_out_records=count,maximum_point_error=maximum,
        no_new_p_values=True,post_selection_not_independent_validation=True)
    out.write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a))


if __name__=='__main__':main()
