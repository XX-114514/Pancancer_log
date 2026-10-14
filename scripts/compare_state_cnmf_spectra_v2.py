"""All-K official cNMF diagnostics and exhaustive Spectra comparison."""
from pathlib import Path
from datetime import datetime,timezone
from itertools import combinations
import os,json,hashlib,gc
import numpy as np,pandas as pd,anndata as ad
from scipy.optimize import linear_sum_assignment
from scipy.stats import rankdata
from cnmf import cNMF

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def frame(p):
    with np.load(p,allow_pickle=True) as z:return pd.DataFrame(z['data'],index=z['index'],columns=z['columns'])
def unit(x):
    d=np.linalg.norm(x,axis=1,keepdims=True);return np.divide(x,d,out=np.full(x.shape,np.nan),where=d>0)
def ranks(x):
    y=rankdata(x,axis=0);y-=y.mean(axis=0);d=np.linalg.norm(y,axis=0);return np.divide(y,d,out=np.full(y.shape,np.nan),where=d>0)
def residual(x,g):
    y=x.astype(float).copy()
    for v in np.unique(g):
        keep=g==v;y[keep]-=y[keep].mean(axis=0)
    return y
def r2groups(x,g):
    den=np.sum((x-x.mean(axis=0))**2,axis=0);err=np.sum(residual(x,g)**2,axis=0)
    return np.divide(den-err,den,out=np.full(den.shape,np.nan),where=den>0)
def topgenes(weights,genes):return [set(genes[np.argsort(-x,kind='stable')[:50]]) for x in weights]
def jaccard(x,y):return len(x&y)/len(x|y)
def main():
    start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=Path(os.environ['PANCANCER_STATE_COMPARISON_PLAN']);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256']
    C=R/'cNMF_state_contrast_20261014_v1';S=R/'Spectra_state_discovery_20261014_v1';job=json.loads((R/'cNMF_state_training_20261014_v1.execution.json').read_text());assert job['exit_code']==0
    for role in ['worker_pid','child_pid']:assert not Path('/proc',str(job[role])).exists()
    inv=json.loads((C/'inventory.json').read_text());assert len(inv['summary'])==2 and all(x['all_requested_K_completed'] for x in inv['summary'])
    for n,h in inv['outputs'].items():assert sha(C/n)==h,n
    assert len(list(C.rglob('*spectra*k*iter*.npz')))==600
    assert sha(S/'inventory.json')==plan['Spectra_inventory_sha256']
    with np.load(S/'private_seed_models/seed20261014/factors_and_discovery_scores.npz',allow_pickle=False) as z:sp={k:z[k].copy() for k in z.files}
    obs=pd.read_csv(S/'private_seed_models/seed20261014/private_discovery_metadata.tsv',sep='\t',index_col=0,dtype={'node_id':str})
    assert np.array_equal(obs.index.to_numpy(str),sp['cell_ids'])
    O=R/'state_method_comparison_20261014_v2';assert not O.exists();O.mkdir()
    stats=[];runs=[];mods=[];cross=[];matches=[];within=[];all_models={}
    for major,name in [('T cell','T_cell'),('Macrophage','Macrophage')]:
        meta=obs.loc[obs.major.eq(major)];m=obs.major.eq(major).to_numpy();source=meta.gse_id.to_numpy();gn=(meta.gse_id.astype(str)+'::'+meta.node_id.astype(str)).to_numpy();model=cNMF(output_dir=str(C/name),name='cnmf');spids=np.where((sp['factor_scope']=='global')|(sp['factor_scope']==major))[0];sg=sp['genes'];sweights=sp['factors'][spids];stops=topgenes(sweights,sg);ss=sp['cell_scores'][m][:,spids];ssr=ranks(residual(ss,source));ssp=ranks(ss)
        for k in plan['K']:
            print(json.dumps(dict(major=major,K=k,phase='ALL_K_DIAGNOSTICS',UTC=datetime.now(timezone.utc).isoformat())),flush=True)
            base=C/name/'cnmf/cnmf_tmp';pre='cnmf.';spec=frame(base/(pre+f'spectra.k_{k}.dt_0_5.consensus.df.npz'));use=frame(base/(pre+f'usages.k_{k}.dt_0_5.consensus.df.npz'));coef=frame(base/(pre+f'gene_spectra_score.k_{k}.dt_0_5.df.npz'));tpm=frame(base/(pre+f'gene_spectra_tpm.k_{k}.dt_0_5.df.npz'));merged=frame(base/(pre+f'spectra.k_{k}.merged.df.npz'));density=frame(base/(pre+f'local_density_cache.k_{k}.merged.df.npz'))
            assert use.shape==(6000,k) and spec.shape==(k,2000) and coef.shape==(k,16625) and tpm.shape==(k,16625)
            assert use.index.equals(meta.index) and coef.index.equals(use.columns) and spec.index.equals(use.columns) and tpm.index.equals(coef.index)
            assert coef.columns.is_unique and np.isfinite(coef.values).all() and (use.values>=0).all() and np.isfinite(use.values).all()
            assert merged.shape==(100*k,2000) and merged.columns.equals(spec.columns) and density.index.equals(merged.index)
            # Official K-selection stats explicitly use unfiltered repetitions.
            ks=model.consensus(k=k,density_threshold=.5,local_neighborhood_size=.3,show_clustering=False,skip_density_and_return_after_stats=True)
            stats.append(dict(major=major,K=k,repeats=100,unfiltered_silhouette=float(ks.loc['silhouette','stats']),unfiltered_prediction_error=float(ks.loc['prediction_error','stats']),merged_factors=len(merged),density_retained_factors=int(density.iloc[:,0].lt(.5).sum()),density_retained_fraction=float(density.iloc[:,0].lt(.5).mean())))
            assert (use.values.sum(axis=1)>0).all()
            normalized=use.values/use.values.sum(axis=1,keepdims=True);cncoef=np.maximum(coef.values,0);genes=coef.columns.to_numpy(str);top=topgenes(cncoef,genes);toptpm=topgenes(tpm.values,genes);pos={g:i for i,g in enumerate(genes)};assert all(g in pos for g in sg);common=cncoef[:,[pos[g] for g in sg]]
            cos=unit(sweights)@unit(common).T;cs=ranks(normalized);cr=ranks(residual(normalized,source));ca=ranks(use.values)
            pooled=ssp.T@cs;src=ssr.T@cr;absolute=ssp.T@ca
            for i,sj in enumerate(spids):
                for jj,component in enumerate(use.columns):
                    matches.append(dict(major=major,K=k,cnmf_component=int(component),Spectra_anchor_factor=int(sj),Spectra_scope=sp['factor_scope'][sj],positive_coef_common6512_cosine=float(cos[i,jj]),top50_positive_coef_jaccard=jaccard(stops[i],top[jj]),top50_tpm_jaccard=jaccard(stops[i],toptpm[jj]),normalized_usage_score_rho=float(pooled[i,jj]),source_residual_normalized_usage_score_rho=float(src[i,jj]),absolute_usage_score_rho=float(absolute[i,jj]),gene_cosine_defined=bool(np.isfinite(cos[i,jj])),score_rank_defined=bool(np.isfinite(src[i,jj]))))
            # Hungarian one-to-one every original repetition, no missing-run acceptance.
            for repeat in range(100):
                keys=[f'iter{repeat}_topic{j+1}' for j in range(k)];assert all(v in merged.index for v in keys)
                x=unit(merged.loc[keys].values);c=x@unit(spec.values).T;ri,ci=linear_sum_assignment(-c)
                for i,j in zip(ri,ci):runs.append(dict(major=major,K=k,repeat=repeat,original_topic=i+1,consensus_component=int(spec.index[j]),weight_cosine=float(c[i,j]),density_retained=bool(density.loc[keys[i]].iloc[0]<.5)))
            rs=r2groups(normalized,source);rn=r2groups(normalized,gn);depth=np.column_stack([np.log1p(meta.full_RNA_total),np.log1p(meta.full_detected_genes)]);dr=ranks(residual(depth,gn));yr=ranks(residual(normalized,gn));dep=yr.T@dr
            ar=ranks(use.values).T@ranks(depth);nr=cs.T@ranks(depth)
            rr=pd.DataFrame(runs);rr=rr[(rr.major==major)&(rr.K==k)]
            for j,component in enumerate(use.columns):
                q=rr[rr.consensus_component.eq(component)];assert len(q)==100
                median=float(q.weight_cosine.median());p10=float(q.weight_cosine.quantile(.1));ret=float(q.density_retained.mean());stable=median>=plan['min_repeat_median_cosine'] and p10>=plan['min_repeat_p10_cosine'] and ret>=plan['min_repeat_density_retained_fraction']
                mods.append(dict(candidate_id=f'CN_{name}_K{k}_{int(component):02}',major=major,K=k,component=int(component),repeat_median_cosine=median,repeat_p10_cosine=p10,repeat_density_retained_fraction=ret,stable_repetition_rule=bool(stable),source_r2=float(rs[j]),source_node_r2=float(rn[j]),RNA_partial_rank_after_source_node=float(dep[j,0]),detected_partial_rank_after_source_node=float(dep[j,1]),absolute_usage_RNA_rho=float(ar[j,0]),normalized_usage_RNA_rho=float(nr[j,0]),positive_coef_genes=int((coef.iloc[j].values>0).sum()),top10_positive_coef_genes=';'.join(genes[np.argsort(-cncoef[j],kind='stable')[:10]]),top10_TPM_genes=';'.join(genes[np.argsort(-tpm.iloc[j].values,kind='stable')[:10]])))
                for g in np.unique(source):
                    keep=source==g;values=normalized[keep,j];within.append(dict(major=major,K=k,component=int(component),gse_id=g,cells=int(keep.sum()),normalized_usage_median=float(np.median(values)),normalized_usage_q10=float(np.quantile(values,.1)),normalized_usage_q90=float(np.quantile(values,.9))))
            all_models[major,k]=dict(coef=cncoef,genes=genes,normalized=normalized,components=use.columns.to_numpy(),top=top,source=source)
        for ka,kb in combinations(plan['K'],2):
            a=all_models[major,ka];b=all_models[major,kb];assert np.array_equal(a['genes'],b['genes']);c=unit(a['coef'])@unit(b['coef']).T;rho=ranks(residual(a['normalized'],source)).T@ranks(residual(b['normalized'],source))
            for i,u in enumerate(a['components']):
                for j,v in enumerate(b['components']):cross.append(dict(major=major,K_a=ka,component_a=int(u),K_b=kb,component_b=int(v),positive_coef_cosine=float(c[i,j]),top50_jaccard=jaccard(a['top'][i],b['top'][j]),source_residual_usage_rho=float(rho[i,j]),comparison_defined=bool(np.isfinite(c[i,j]) and np.isfinite(rho[i,j]))))
        del model;gc.collect()
    modules=pd.DataFrame(mods);cx=pd.DataFrame(cross)
    stable_map={(r.major,r.K,r.component):r.stable_repetition_rule for r in modules.itertuples()}
    for i,row in modules.iterrows():
        same=cx[cx.major.eq(row.major)];q=same[((same.K_a==row.K)&(same.component_a==row.component))|((same.K_b==row.K)&(same.component_b==row.component))];supported=q[(q.positive_coef_cosine>=plan['min_crossK_cosine'])&(q.top50_jaccard>=plan['min_crossK_top50_jaccard'])&(q.source_residual_usage_rho>=plan['min_crossK_source_residual_rho'])]
        count=0
        for p in supported.itertuples():
            other=(p.major,p.K_b,p.component_b) if p.K_a==row.K else (p.major,p.K_a,p.component_a)
            count+=bool(stable_map[other])
        modules.loc[i,'crossK_any_supported']=len(supported)>0;modules.loc[i,'crossK_supported']=count>0;modules.loc[i,'crossK_stable_support_pairs']=count
    for name,x in [('all90_cnmf_candidates.tsv',modules),('all9000_repetition_matches.tsv',pd.DataFrame(runs)),('all14355_Spectra_cnmf_comparisons.tsv',pd.DataFrame(matches)),('all1300_crossK_comparisons.tsv',cx),('all6_K_diagnostics.tsv',pd.DataFrame(stats)),('all360_source_distributions.tsv',pd.DataFrame(within))]:x.to_csv(O/name,sep='\t',index=False)
    assert len(mods)==90 and len(runs)==9000 and len(matches)==14355 and len(cross)==1300 and len(within)==360 and len(stats)==6
    # The official stats call must not mutate original fit/consensus outputs.
    for n,h in inv['outputs'].items():assert sha(C/n)==h,n
    out=dict(status='COMPLETE_ALL_K_METHOD_COMPARISON_PENDING_STATE_DICTIONARY',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),cnmf_inventory_sha256=sha(C/'inventory.json'),Spectra_inventory_sha256=sha(S/'inventory.json'),candidates=90,repeat_stable_candidates=int(modules.stable_repetition_rule.sum()),crossK_supported_candidates=int(modules.crossK_supported.sum()),repeat_and_crossK_candidates=int((modules.stable_repetition_rule&modules.crossK_supported).sum()),cNMF_stats_unfiltered=True,original_fit_outputs_unchanged=True,formal_state_labels=False,UCell_complete=False,clinical_or_LR_p_selected=False,output_sha256={p.name:sha(p) for p in O.iterdir() if p.is_file()});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
