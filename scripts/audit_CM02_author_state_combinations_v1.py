"""Rebuild all count effects and resamples from the separate prior count table."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib,json,os
import numpy as np
import pandas as pd

REG=['4-1BB+ Treg','HLA-DR+ Treg','LAG3+ Treg','Prolif Treg']
OTHER=['Tn/Tmem','IFN-stim CD4T','GZMK+T','Tmem','Teff','Tstr','Tprolif']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def calculate(y,c):
    d=c.copy();d[:,1:3]=pd.DataFrame(c[:,1:3]).rank(method='average').to_numpy()
    d[:,1:3]-=d[:,1:3].mean(axis=0);d[:,1:3]/=np.maximum(d[:,1:3].std(axis=0),1e-15)
    u,s,_=np.linalg.svd(d,full_matrices=False);rank=sum(s>s.max()*max(d.shape)*np.finfo(float).eps)
    if rank!=5 or len(y)-rank-1<10 or s[0]/s[-1]>1e6:return np.full((4,7),np.nan)
    a=pd.DataFrame(y).rank(method='average').to_numpy();res=a-u@(u.T@a)
    sums=np.sum(res*res,axis=0);den=np.sqrt(sums[:4,None]*sums[None,4:])
    r=res[:,:4].T@res[:,4:]/np.maximum(den,1e-15)
    r[(sums[:4,None]<1e-20)|(sums[None,4:]<1e-20)]=np.nan
    return r


def compare(a,b):
    if b is None:assert not np.isfinite(a)
    else:assert np.isfinite(a) and abs(a-b)<1e-9


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);parent=Path(os.environ['PANCANCER_PARENT_OUTPUT']);out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent==Path('/dev/shm') and not out.exists()
    terminal=json.loads((parent/'private_terminal.json').read_text())
    for n,v in terminal['files'].items():assert sha(parent/n)==v['SHA256'] and (parent/n).stat().st_size==v['bytes']
    results=json.loads((parent/'private_all56_source_count_combinations.json').read_text())
    indexed={(r['phase'],r['reg_author_state'],r['other_author_T_state']):r for r in results}
    assert len(results)==len(indexed)==56
    meta=root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    fields=['library_id','patient_registry','timepoint_registry','neo_type','HPV/p16']
    f=pd.read_csv(meta,sep='\t',usecols=fields,dtype=str,keep_default_na=False)
    counts=pd.read_csv(root/'HNSCC_author_state_count_layer_preflight_20261016_v1/private_all57_original_author_counts.tsv',sep='\t',usecols=['library_id','Cells analyzed']+REG+OTHER,dtype={'library_id':str},keep_default_na=False)
    f=f.merge(counts,on='library_id',validate='one_to_one');checked=0
    for pi,phase in enumerate(['pre','post']):
        g=f[f.timepoint_registry.eq(phase)].copy();nr=g[REG].sum(axis=1);ne=g[OTHER].sum(axis=1)
        okay=nr.ge(20)&ne.ge(20)&g.neo_type.isin(['D','D+T'])&g['HPV/p16'].isin(['positive','negative'])&~g.patient_registry.duplicated(keep=False)
        g=g[okay].sort_values('patient_registry');nr=g[REG].sum(axis=1).to_numpy();ne=g[OTHER].sum(axis=1).to_numpy()
        x=g[REG+OTHER].to_numpy(float);y=np.column_stack([x[:,:4]/nr[:,None],x[:,4:]/ne[:,None]])
        tx=g.neo_type.eq('D+T').to_numpy(int);hpv=g['HPV/p16'].eq('positive').to_numpy(int);n=len(g)
        c=np.column_stack([np.ones(n),np.log1p(nr),np.log1p(ne),tx,hpv]);r=calculate(y,c)
        sensitivities={name:calculate(x/den[:,None],c) for name,den in [('all_author_T_cells',nr+ne),('all_author_analyzed_cells',g['Cells analyzed'].to_numpy(float))]}
        loo=np.asarray([calculate(np.delete(y,i,axis=0),np.delete(c,i,axis=0)) for i in range(n)])
        strata=2*tx+hpv;rng=np.random.default_rng(2026102102+pi);draws=[]
        for _ in range(2000):
            ix=np.concatenate([rng.choice(np.where(strata==s)[0],int(np.count_nonzero(strata==s)),replace=True) for s in sorted(set(strata))]);draws.append(calculate(y[ix],c[ix]))
        draws=np.asarray(draws)
        for i,a in enumerate(REG):
            for j,b in enumerate(OTHER):
                row=indexed[phase,a,b];assert row['author_patients']==n;compare(r[i,j],row['partial_Spearman'])
                z=draws[:,i,j];z=z[np.isfinite(z)];lv=loo[:,i,j];lv=lv[np.isfinite(lv)]
                assert row['valid_bootstraps']==len(z) and row['valid_leave_one']==len(lv)
                ci=np.quantile(z,[.025,.975]) if len(z)>=1800 else [np.nan,np.nan]
                for v,w in zip(ci,row['pointwise_bootstrap_CI']):compare(v,w)
                compare(lv.min() if len(lv) else np.nan,row['leave_one_min']);compare(lv.max() if len(lv) else np.nan,row['leave_one_max'])
                if np.isfinite(r[i,j]):assert row['leave_one_same_direction']==int(np.sum(np.sign(lv)==np.sign(r[i,j])))
                else:assert row['leave_one_same_direction'] is None
                for name,v in sensitivities.items():compare(v[i,j],row['sensitivity_partial_Spearman'][name])
                checked+=1
    assert checked==56
    result=dict(status='ACCEPTED_ALL56_COUNT_COMBINATIONS_ALL_BOOTSTRAPS_LOO_AND_DENOMINATOR_SENSITIVITIES',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),parent_terminal_SHA256=sha(parent/'private_terminal.json'),independent_prior_count_table_Pandas_ranks_and_SVD_projection=True,all56_positions_and_full_random_stream_replayed=True,pointwise_not_multiplicity_controlled_intervals=True,no_new_clinical_inference_or_identity_adoption=True,individual_rows_exported=False)
    out.mkdir(mode=0o700);(out/'private_independent_acceptance.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
