"""Complete outcome-free source-count structure; no state/clinical validation."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib,json,os
import numpy as np
import pandas as pd
from scipy.stats import rankdata

REG=['4-1BB+ Treg','HLA-DR+ Treg','LAG3+ Treg','Prolif Treg']
OTHER=['Tn/Tmem','IFN-stim CD4T','GZMK+T','Tmem','Teff','Tstr','Tprolif']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def effect(y,c):
    d=c.copy()
    for j in [1,2]:d[:,j]=rankdata(d[:,j])
    d[:,1:3]-=d[:,1:3].mean(axis=0);d[:,1:3]/=np.maximum(d[:,1:3].std(axis=0),1e-15)
    rank=np.linalg.matrix_rank(d)
    if rank!=5 or len(y)-rank-1<10 or np.linalg.cond(d)>1e6:return np.full((4,7),np.nan)
    a=rankdata(y,axis=0);e=a-d@np.linalg.lstsq(d,a,rcond=1e-12)[0]
    n=np.linalg.norm(e,axis=0);r=e[:,:4].T@e[:,4:]/np.maximum(n[:4,None]*n[None,4:],1e-15)
    r[(n[:4,None]<1e-10)|(n[None,4:]<1e-10)]=np.nan
    # Separate QR projection checks every point/bootstrap/leave-one effect.
    q,_=np.linalg.qr(d,mode='reduced');f=a-q@(q.T@a);m=np.linalg.norm(f,axis=0)
    alt=f[:,:4].T@f[:,4:]/np.maximum(m[:4,None]*m[None,4:],1e-15)
    good=np.isfinite(r);np.testing.assert_allclose(r[good],alt[good],atol=1e-10,rtol=1e-10)
    if good.any(): assert np.max(np.abs(r[good]))<=1+1e-10
    return r


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);parent=Path(os.environ['PANCANCER_PARENT_OUTPUT']);out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent==Path('/dev/shm') and not out.exists()
    t=json.loads((parent/'private_terminal.json').read_text())
    assert t['status']=='COMPLETE_CM02_CONTEXT_AND_ALL_AUTHOR_STATE_COUNT_OPPORTUNITIES_NO_NEW_EFFECTS'
    previous=json.loads((parent/'plan_before_new_state_opportunities.json').read_text())
    p=root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    assert sha(p)==previous['author_join_SHA256']
    out.mkdir(mode=0o700);save=lambda n,v:(out/n).write_text(json.dumps(v,indent=2)+'\n')
    save('plan_before_any_new_count_association.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),parent_terminal_SHA256=sha(parent/'private_terminal.json'),source='GSE286827',module_context='2026-10-14 immune atlas v1 T_whole_TME_primary_tumor CM02 followup, not exact cell mapping',primary_phase='pre',secondary_phase='post_same_cohort_not_replication',all_positions=56,reg_author_states=REG,other_author_T_states=OTHER,main_denominators='reg within four author Treg states; other within seven remaining author T states; disjoint denominators',covariates=['log1p_author_Treg_cells','log1p_author_other_T_cells','treatment_DplusT','HPV_positive'],sensitivity_denominators=['all_author_T_cells','all_author_analyzed_cells'],bootstraps=2000,bootstrap_unit='author patient within treatment x HPV, no new independence review',minimum_valid_draws=1800,seed=2026102102,confidence_intervals='pointwise descriptive, no simultaneous coverage or multiplicity-controlled pair declarations',no_P_or_FDR_or_clinical_tests=True,outcome_columns_not_read=True,original_cells_and_formal_labels_unchanged=True))
    cols=['library_id','patient_registry','timepoint_registry','neo_type','HPV/p16','Cells analyzed']+REG+OTHER
    f=pd.read_csv(p,sep='\t',usecols=cols,dtype={'library_id':str,'patient_registry':str},keep_default_na=False)
    coverage={r['phase']:r for r in json.loads((parent/'private_all2_phase_scope.json').read_text())}
    rows=[];proof=[]
    for pi,phase in enumerate(['pre','post']):
        g=f[f.timepoint_registry.eq(phase)].copy();x=g[REG+OTHER].to_numpy(float);nr=x[:,:4].sum(axis=1);ne=x[:,4:].sum(axis=1)
        keep=(nr>=20)&(ne>=20)&g.neo_type.isin(['D','D+T']).to_numpy()&g['HPV/p16'].isin(['positive','negative']).to_numpy()&~g.patient_registry.duplicated(keep=False).to_numpy()
        g=g.loc[keep].sort_values('patient_registry');x=g[REG+OTHER].to_numpy(float);nr=x[:,:4].sum(axis=1);ne=x[:,4:].sum(axis=1)
        assert len(g)==coverage[phase]['eligible_for_count_structure']
        tx=g.neo_type.eq('D+T').to_numpy(int);hpv=g['HPV/p16'].eq('positive').to_numpy(int)
        c=np.column_stack([np.ones(len(g)),np.log1p(nr),np.log1p(ne),tx,hpv])
        y=np.column_stack([x[:,:4]/nr[:,None],x[:,4:]/ne[:,None]])
        main=effect(y,c);n=len(g)
        sensitivities={}
        for name,den in [('all_author_T_cells',nr+ne),('all_author_analyzed_cells',pd.to_numeric(g['Cells analyzed']).to_numpy(float))]:sensitivities[name]=effect(x/den[:,None],c)
        loo=np.array([effect(np.delete(y,i,0),np.delete(c,i,0)) for i in range(n)])
        rng=np.random.default_rng(2026102102+pi);strata=2*tx+hpv;blocks=[np.flatnonzero(strata==v) for v in np.unique(strata)]
        boot=[]
        for _ in range(2000):
            draw=np.concatenate([rng.choice(ix,len(ix),replace=True) for ix in blocks]);boot.append(effect(y[draw],c[draw]))
        boot=np.asarray(boot)
        for i,a in enumerate(REG):
            for j,b in enumerate(OTHER):
                draws=boot[:,i,j];draws=draws[np.isfinite(draws)];lv=loo[:,i,j];lv=lv[np.isfinite(lv)]
                ci=np.quantile(draws,[.025,.975]).tolist() if len(draws)>=1800 else [None,None]
                finite=np.isfinite(main[i,j]);r=float(main[i,j]) if finite else None
                rows.append(dict(source='GSE286827',phase=phase,reg_author_state=a,other_author_T_state=b,author_patients=n,partial_Spearman=r,pointwise_bootstrap_CI=ci,valid_bootstraps=len(draws),valid_leave_one=len(lv),leave_one_min=float(lv.min()) if len(lv) else None,leave_one_max=float(lv.max()) if len(lv) else None,leave_one_same_direction=int(np.sum(np.sign(lv)==np.sign(r))) if finite else None,sensitivity_partial_Spearman={k:float(v[i,j]) if np.isfinite(v[i,j]) else None for k,v in sensitivities.items()},no_P_or_multiplicity_controlled_pair_declaration=True,not_CD8_identity_CM02_mapping_or_clinical_validation=True))
        proof.append(dict(source='GSE286827',phase=phase,all_point_bootstrap_leave_one_and_denominator_sensitivity_effects_checked_by_QR=True,original_count_scope_exact=True,outcome_RNA_and_individual_rows_not_exported=True))
    assert len(rows)==56
    save('private_all56_source_count_combinations.json',rows);save('private_all_phase_numeric_proofs.json',proof)
    files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}
    assert sum(v['bytes'] for v in files.values())<1024**2
    terminal=dict(status='COMPLETE_ALL56_OUTCOME_FREE_AUTHOR_COUNT_COMBINATIONS_DESCRIPTIVE_ONLY',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),files=files)
    save('private_terminal.json',terminal);print(json.dumps(dict(status=terminal['status'],terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
