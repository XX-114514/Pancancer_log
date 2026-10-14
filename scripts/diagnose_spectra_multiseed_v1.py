"""Frozen full-grid diagnostics; aggregate outputs only, no scientific release."""
from pathlib import Path
from datetime import datetime, timezone
from itertools import combinations
import os, json, hashlib
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import rankdata

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608), b''): h.update(b)
    return h.hexdigest()

def correlation(x,y):
    a=rankdata(x); b=rankdata(y); a-=a.mean(); b-=b.mean()
    d=np.linalg.norm(a)*np.linalg.norm(b)
    return float(a@b/d) if d>0 else np.nan

def dummy(v):
    return pd.get_dummies(pd.Series(v,dtype=str),dtype=float).to_numpy()

def r2(y,d):
    beta=np.linalg.lstsq(d,y,rcond=None)[0]; e=y-d@beta
    den=np.sum((y-y.mean(axis=0))**2,axis=0)
    return np.divide(den-np.sum(e**2,axis=0),den,out=np.full(den.shape,np.nan),where=den>0)

def residual(y,d):
    return y-d@np.linalg.lstsq(d,y,rcond=None)[0]

def main():
    start=datetime.now(timezone.utc).isoformat()
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); P=Path(os.environ['PANCANCER_SPECTRA_DIAGNOSTIC_PLAN'])
    plan=json.loads(P.read_text()); assert sha(__file__)==plan['producer_sha256']
    S=R/'Spectra_state_discovery_20261014_v1'; inv=json.loads((S/'inventory.json').read_text())
    assert sha(S/'inventory.json')==plan['Spectra_inventory_sha256']
    for name,h in inv['outputs'].items(): assert sha(S/name)==h, name
    seeds=plan['seeds']; models=[]
    for seed in seeds:
        p=S/'private_seed_models'/('seed'+str(seed))
        with np.load(p/'factors_and_discovery_scores.npz',allow_pickle=False) as z:
            models.append({k:z[k].copy() for k in z.files})
    ref=models[0]; obs=pd.read_csv(S/'private_seed_models'/('seed'+str(seeds[0]))/'private_discovery_metadata.tsv',sep='\t',index_col=0,dtype={'node_id':str,'gse_id':str,'major':str})
    assert np.array_equal(obs.index.to_numpy(str),ref['cell_ids'])
    for m in models:
        for k in ['cell_ids','genes','factor_scope']: assert np.array_equal(m[k],ref[k])
        assert np.isfinite(m['factors']).all() and np.isfinite(m['cell_scores']).all()
    O=R/'Spectra_multiseed_diagnostics_20261014_v1'; assert not O.exists(); O.mkdir()
    top=[np.argsort(-m['factors'],axis=1,kind='stable')[:,:plan['top_genes']] for m in models]
    unit=[m['factors']/np.maximum(np.linalg.norm(m['factors'],axis=1)[:,None],1e-30) for m in models]
    maps={}; pairrows=[]; score_rows=[]; scope=ref['factor_scope']
    # Pairwise matching is restricted to the frozen biological mask, not score p values.
    for a,b in combinations(range(len(seeds)),2):
        mapping={}
        for s in ['global','T cell','Macrophage']:
            ix=np.where(scope==s)[0]; c=unit[a][ix]@unit[b][ix].T
            ri,ci=linear_sum_assignment(-c)
            for u,v in zip(ix[ri],ix[ci]): mapping[int(u)]=int(v)
        maps[a,b]=mapping
        for u,v in mapping.items():
            x=set(top[a][u]);y=set(top[b][v]);j=len(x&y)/len(x|y)
            pairrows.append(dict(seed_a=seeds[a],seed_b=seeds[b],factor_a=u,factor_b=v,scope=scope[u],cosine=float(unit[a][u]@unit[b][v]),top50_jaccard=j))
            majors=['T cell','Macrophage'] if scope[u]=='global' else [scope[u]]
            for major in majors:
                mask=obs.major.eq(major).to_numpy(); xa=models[a]['cell_scores'][mask,u];xb=models[b]['cell_scores'][mask,v]
                source=obs.gse_id.to_numpy()[mask];d=dummy(source)
                score_rows.append(dict(seed_a=seeds[a],seed_b=seeds[b],factor_a=u,factor_b=v,major=major,cells=int(mask.sum()),score_spearman=correlation(xa,xb),source_residual_spearman=correlation(residual(xa,d),residual(xb,d))))
    pd.DataFrame(pairrows).to_csv(O/'pairwise_factor_matching.tsv',sep='\t',index=False)
    pd.DataFrame(score_rows).to_csv(O/'pairwise_score_reproducibility.tsv',sep='\t',index=False)
    # Anchor is prespecified first seed, never the most favourable seed.
    diagnostics=[]; within=[]
    for major in ['T cell','Macrophage']:
        keep=obs.major.eq(major).to_numpy(); meta=obs.loc[keep]; ids=np.where((scope=='global')|(scope==major))[0]
        y=ref['cell_scores'][keep][:,ids]; ds=dummy(meta.gse_id.to_numpy())
        dn=dummy((meta.gse_id.astype(str)+'::'+meta.node_id.astype(str)).to_numpy())
        # Joint source/node categories contain the source-only model exactly.
        depth=np.column_stack([np.log1p(meta.full_RNA_total),np.log1p(meta.full_detected_genes)])
        depth=(depth-depth.mean(axis=0))/depth.std(axis=0)
        rs=r2(y,ds);rn=r2(y,dn);rd=r2(y,np.column_stack([dn,depth]));yr=residual(y,dn)
        for jj,j in enumerate(ids):
            diagnostics.append(dict(anchor_factor=int(j),scope=scope[j],major=major,cells=len(meta),source_r2=float(rs[jj]),source_node_r2=float(rn[jj]),node_increment_after_source=float(rn[jj]-rs[jj]),RNA_detected_increment_after_source_node=float(rd[jj]-rn[jj]),RNA_partial_rank_after_source_node=correlation(yr[:,jj],residual(depth[:,0],dn)),detected_partial_rank_after_source_node=correlation(yr[:,jj],residual(depth[:,1],dn))))
        for source in sorted(meta.gse_id.unique()):
            mm=meta.gse_id.eq(source).to_numpy();assert mm.sum()>=100
            for jj,j in enumerate(ids):
                values=y[mm,jj]
                within.append(dict(anchor_factor=int(j),scope=scope[j],major=major,gse_id=source,cells=int(mm.sum()),score_median=float(np.median(values)),score_q10=float(np.quantile(values,.1)),score_q90=float(np.quantile(values,.9)),RNA_spearman=correlation(values,depth[mm,0]),detected_spearman=correlation(values,depth[mm,1])))
    pd.DataFrame(diagnostics).to_csv(O/'source_identity_complexity_diagnostics.tsv',sep='\t',index=False)
    pd.DataFrame(within).to_csv(O/'within_source_major_distributions.tsv',sep='\t',index=False)
    pairs=pd.DataFrame(pairrows);scores=pd.DataFrame(score_rows);d=pd.DataFrame(diagnostics)
    dictionary=json.loads((S/'fixed_feature_dictionary.json').read_text())['reference']; candidates=[]
    for j in range(len(scope)):
        js=[j,maps[0,1][j],maps[0,2][j]]; cyc=maps[1,2][js[1]]==js[2]
        prows=[];srows=[]
        for a,b in combinations(range(3),2):
            prows.extend(pairs[(pairs.seed_a==seeds[a])&(pairs.seed_b==seeds[b])&(pairs.factor_a==js[a])].to_dict('records'))
            srows.extend(scores[(scores.seed_a==seeds[a])&(scores.seed_b==seeds[b])&(scores.factor_a==js[a])].to_dict('records'))
        cos=min(x['cosine'] for x in prows); jac=min(x['top50_jaccard'] for x in prows)
        sr=[x['source_residual_spearman'] for x in srows];rho=min(sr) if np.isfinite(sr).all() else np.nan
        stable=bool(cyc and cos>=plan['min_cosine'] and jac>=plan['min_top50_jaccard'] and np.isfinite(rho) and rho>=plan['min_source_residual_score_rho'])
        dg=d[d.anchor_factor.eq(j)]; complexity=bool((dg.RNA_partial_rank_after_source_node.abs()>=.5).any() or (dg.detected_partial_rank_after_source_node.abs()>=.5).any())
        # Similarity to the prior is a descriptor; never a confirmed state name.
        gs=set(ref['genes'][top[0][j]]); prior=[]
        for name,genes in dictionary[scope[j]].items():
            ss=set(genes);prior.append((len(gs&ss)/len(gs|ss),name,len(gs&ss)))
        best=max(prior) if prior else (0,'none',0)
        c=dict(candidate_id='SP_'+scope[j].replace(' ','_')+'_'+str(j+1).zfill(3),anchor_factor=j,scope=scope[j],seed2_factor=js[1],seed3_factor=js[2],cycle_consistent=bool(cyc),min_weight_cosine=cos,min_top50_jaccard=jac,min_source_residual_score_rho=rho,stable_by_frozen_rule=stable,complexity_flag=complexity,max_source_r2=float(dg.source_r2.max()),max_source_node_r2=float(dg.source_node_r2.max()),best_prior_descriptor=best[1],best_prior_top50_jaccard=best[0],best_prior_top50_overlap=best[2],top10_anchor_genes=';'.join(ref['genes'][top[0][j,:10]]))
        candidates.append(c)
    pd.DataFrame(candidates).to_csv(O/'all173_candidate_diagnostics.tsv',sep='\t',index=False)
    # Independent scalar formulas, incl. invariant, negative and undefined cases.
    checks=[]
    for row in [pairrows[0],min(pairrows,key=lambda x:x['cosine']),max(pairrows,key=lambda x:x['cosine'])]:
        a=seeds.index(row['seed_a']);b=seeds.index(row['seed_b']);u=row['factor_a'];v=row['factor_b'];x=models[a]['factors'][u].astype(float);y=models[b]['factors'][v].astype(float)
        direct=sum(float(xx)*float(yy) for xx,yy in zip(x,y))/(sum(float(xx)**2 for xx in x)**.5*sum(float(yy)**2 for yy in y)**.5)
        assert np.isclose(direct,row['cosine'],atol=1e-6)
        checks.append(dict(check='scalar_cosine',factor_a=u,factor_b=v,reported=row['cosine'],independent=direct,pass_check=True))
    assert np.isclose(correlation(np.arange(7),np.arange(7)),1) and np.isclose(correlation(np.arange(7),-np.arange(7)),-1) and np.isnan(correlation(np.ones(7),np.arange(7)))
    checks.append(dict(check='positive_negative_constant_rank',pass_check=True))
    (O/'independent_scalar_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    out=dict(status='COMPLETE_ALL173_DIAGNOSTICS_NOT_STATE_ACCEPTANCE',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),input_inventory_sha256=sha(S/'inventory.json'),cells=12000,seeds=seeds,candidates=173,pairwise_matches=len(pairrows),score_comparisons=len(score_rows),stable_count=sum(x['stable_by_frozen_rule'] for x in candidates),cycle_consistent_count=sum(x['cycle_consistent'] for x in candidates),complexity_flag_count=sum(x['complexity_flag'] for x in candidates),stable_by_scope={s:sum(x['stable_by_frozen_rule'] for x in candidates if x['scope']==s) for s in ['global','T cell','Macrophage']},formal_state_labels=False,UCell_complete=False,source_is_not_artifact_by_definition=True,cell_level_p_values_used=False,output_sha256={p.name:sha(p) for p in O.iterdir() if p.is_file()})
    (O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))

if __name__=='__main__':main()
