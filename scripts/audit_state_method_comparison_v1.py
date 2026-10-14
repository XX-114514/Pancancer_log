"""Independent numeric readback; does not rerun training or import producer."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,os
import numpy as np,pandas as pd
from scipy.stats import spearmanr

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def f(p):
    with np.load(p,allow_pickle=True) as z:return pd.DataFrame(z['data'],index=z['index'],columns=z['columns'])
def cos(x,y):
    x=x.astype(float);y=y.astype(float);n=np.sqrt(np.dot(x,x)*np.dot(y,y));return np.dot(x,y)/n if n>0 else np.nan
def top(x,genes):return set(genes[np.argsort(-x,kind='stable')[:50]])
def rho_matrix(x,y):
    # pandas average ranks + pairwise Pearson differs from producer normalized-rank dot.
    v=pd.DataFrame(np.column_stack([x,y])).rank(axis=0,method='average');c=v.corr(method='pearson').to_numpy();return c[:x.shape[1],x.shape[1]:]
def demeangroup(x,g):
    return pd.DataFrame(x).sub(pd.DataFrame(x).groupby(pd.Series(g)).transform('mean')).to_numpy()
def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);D=R/'state_method_comparison_20261014_v2';C=R/'cNMF_state_contrast_20261014_v1';S=R/'Spectra_state_discovery_20261014_v1';inv=json.loads((D/'inventory.json').read_text())
    for n,h in inv['output_sha256'].items():assert sha(D/n)==h
    ci=json.loads((C/'inventory.json').read_text())
    for n,h in ci['outputs'].items():assert sha(C/n)==h
    assert sha(C/'inventory.json')==inv['cnmf_inventory_sha256']
    with np.load(S/'private_seed_models/seed20261014/factors_and_discovery_scores.npz',allow_pickle=False) as z:sp=dict(z)
    obs=pd.read_csv(S/'private_seed_models/seed20261014/private_discovery_metadata.tsv',sep='\t',index_col=0,dtype={'node_id':str});assert np.array_equal(obs.index.to_numpy(str),sp['cell_ids'])
    repeat=pd.read_csv(D/'all9000_repetition_matches.tsv',sep='\t');pairs=pd.read_csv(D/'all14355_Spectra_cnmf_comparisons.tsv',sep='\t');mods=pd.read_csv(D/'all90_cnmf_candidates.tsv',sep='\t');cross=pd.read_csv(D/'all1300_crossK_comparisons.tsv',sep='\t');stats=pd.read_csv(D/'all6_K_diagnostics.tsv',sep='\t');plan=json.loads(Path(os.environ['PANCANCER_STATE_COMPARISON_PLAN']).read_text());cache={};errors=[];score_errors=[];diag_errors=[];jchecks=0
    for major,name in [('T cell','T_cell'),('Macrophage','Macrophage')]:
        meta=obs[obs.major.eq(major)];g=meta.gse_id.to_numpy();scope=np.where((sp['factor_scope']=='global')|(sp['factor_scope']==major))[0];ss=sp['cell_scores'][obs.major.eq(major).to_numpy()][:,scope]
        for k in [10,15,20]:
            b=C/name/'cnmf/cnmf_tmp';use=f(b/f'cnmf.usages.k_{k}.dt_0_5.consensus.df.npz');coef=f(b/f'cnmf.gene_spectra_score.k_{k}.dt_0_5.df.npz');tp=f(b/f'cnmf.gene_spectra_tpm.k_{k}.dt_0_5.df.npz');spec=f(b/f'cnmf.spectra.k_{k}.dt_0_5.consensus.df.npz');merged=f(b/f'cnmf.spectra.k_{k}.merged.df.npz');den=f(b/f'cnmf.local_density_cache.k_{k}.merged.df.npz');assert use.index.equals(meta.index)
            q=repeat[(repeat.major==major)&(repeat.K==k)];assert len(q)==k*100
            for rep,qq in q.groupby('repeat'):
                assert len(qq)==k and qq.original_topic.is_unique and qq.consensus_component.is_unique
            for row in q.itertuples():
                key=f'iter{row.repeat}_topic{row.original_topic}';c=cos(merged.loc[key].values,spec.loc[row.consensus_component].values);assert np.isclose(c,row.weight_cosine,atol=1e-7);assert bool(den.loc[key].iloc[0]<.5)==row.density_retained;errors.append(abs(c-row.weight_cosine))
            norm=use.div(use.sum(axis=1),axis=0).to_numpy();raw=use.to_numpy();co=coef.clip(lower=0);sg=sp['genes'];common=co.loc[:,sg].to_numpy();gene=co.columns.to_numpy(str);cn_top={int(j):top(co.loc[j].values,gene) for j in co.index};tp_top={int(j):top(tp.loc[j].values,gene) for j in tp.index}
            r=rho_matrix(ss,norm);rr=rho_matrix(demeangroup(ss,g),demeangroup(norm,g));ar=rho_matrix(ss,raw);table=pairs[(pairs.major==major)&(pairs.K==k)];assert len(table)==len(scope)*k
            for row in table.itertuples():
                si=list(scope).index(row.Spectra_anchor_factor);j=list(co.index).index(row.cnmf_component);c=cos(sp['factors'][row.Spectra_anchor_factor],common[j]);st=top(sp['factors'][row.Spectra_anchor_factor],sg)
                for actual,reported in [(c,row.positive_coef_common6512_cosine),(r[si,j],row.normalized_usage_score_rho),(rr[si,j],row.source_residual_normalized_usage_score_rho),(ar[si,j],row.absolute_usage_score_rho)]:
                    assert np.isclose(actual,reported,atol=2e-6,equal_nan=True);score_errors.append(abs(actual-reported))
                for ct,reported in [(cn_top[row.cnmf_component],row.top50_positive_coef_jaccard),(tp_top[row.cnmf_component],row.top50_tpm_jaccard)]:assert abs(len(st&ct)/len(st|ct)-reported)<1e-12;jchecks+=1
                assert row.gene_cosine_defined==bool(np.isfinite(c)) and row.score_rank_defined==bool(np.isfinite(rr[si,j]))
            gn=(meta.gse_id.astype(str)+'::'+meta.node_id.astype(str)).to_numpy();depth=np.column_stack([np.log1p(meta.full_RNA_total),np.log1p(meta.full_detected_genes)]);rn=demeangroup(norm,gn);dr=demeangroup(depth,gn);dep=rho_matrix(rn,dr)
            for row in mods[(mods.major==major)&(mods.K==k)].itertuples():
                j=list(co.index).index(row.component);qq=q[q.consensus_component.eq(row.component)];med=float(qq.weight_cosine.median());p10=float(qq.weight_cosine.quantile(.1));ret=float(qq.density_retained.mean());stable=med>=plan['min_repeat_median_cosine'] and p10>=plan['min_repeat_p10_cosine'] and ret>=plan['min_repeat_density_retained_fraction'];assert stable==row.stable_repetition_rule
                assert np.allclose([med,p10,ret],[row.repeat_median_cosine,row.repeat_p10_cosine,row.repeat_density_retained_fraction],atol=1e-12)
                for groups,reported in [(g,row.source_r2),(gn,row.source_node_r2)]:
                    x=norm[:,j];a=1-np.sum(demeangroup(x[:,None],groups)**2)/np.sum((x-x.mean())**2);assert np.isclose(a,reported,atol=1e-7,equal_nan=True);diag_errors.append(abs(a-reported))
                for actual,reported in [(dep[j,0],row.RNA_partial_rank_after_source_node),(dep[j,1],row.detected_partial_rank_after_source_node)]:assert np.isclose(actual,reported,atol=2e-7,equal_nan=True);diag_errors.append(abs(actual-reported))
            cache[major,k]=dict(co=co,norm=norm,g=g)
    cross_ranks={}
    for row in cross.itertuples():
        a=cache[row.major,row.K_a];b=cache[row.major,row.K_b];assert a['co'].columns.equals(b['co'].columns);c=cos(a['co'].loc[row.component_a].values,b['co'].loc[row.component_b].values);ta=top(a['co'].loc[row.component_a].values,a['co'].columns.to_numpy(str));tb=top(b['co'].loc[row.component_b].values,b['co'].columns.to_numpy(str));j=len(ta&tb)/len(ta|tb);assert np.isclose(c,row.positive_coef_cosine,atol=1e-7,equal_nan=True) and np.isclose(j,row.top50_jaccard,atol=1e-12)
        key=(row.major,row.K_a,row.K_b)
        if key not in cross_ranks:cross_ranks[key]=rho_matrix(demeangroup(a['norm'],a['g']),demeangroup(b['norm'],b['g']))
        i=list(a['co'].index).index(row.component_a);jj=list(b['co'].index).index(row.component_b);rho=cross_ranks[key][i,jj];assert np.isclose(rho,row.source_residual_usage_rho,atol=2e-7,equal_nan=True);assert row.comparison_defined==bool(np.isfinite(c) and np.isfinite(rho))
    stable={(row.major,row.K,row.component):row.stable_repetition_rule for row in mods.itertuples()}
    for row in mods.itertuples():
        q=cross[cross.major.eq(row.major)];q=q[((q.K_a==row.K)&(q.component_a==row.component))|((q.K_b==row.K)&(q.component_b==row.component))];q=q[(q.positive_coef_cosine>=plan['min_crossK_cosine'])&(q.top50_jaccard>=plan['min_crossK_top50_jaccard'])&(q.source_residual_usage_rho>=plan['min_crossK_source_residual_rho'])];n=0
        for p in q.itertuples():n+=bool(stable[(p.major,p.K_b,p.component_b)] if p.K_a==row.K else stable[(p.major,p.K_a,p.component_a)])
        assert bool(n>0)==row.crossK_supported and bool(len(q)>0)==row.crossK_any_supported and n==row.crossK_stable_support_pairs
    assert len(stats)==6 and stats[['unfiltered_silhouette','unfiltered_prediction_error']].notna().all().all()
    O=R/'state_method_comparison_audit_20261014_v1';assert not O.exists();O.mkdir();out=dict(status='PASS_ALL_REPEAT_GENE_SCORE_RANK_IDENTITY_AND_CANDIDATE_RULE_READBACK',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(D/'inventory.json'),repeat_cosine_and_density_checks=len(errors),max_repeat_cosine_difference=max(errors),Spectra_cNMF_cosine_and_rank_statistics=len(score_errors),max_Spectra_cNMF_difference=float(np.nanmax(score_errors)),gene_overlap_checks=jchecks,source_node_complexity_statistics=len(diag_errors),max_diagnostic_difference=float(np.nanmax(diag_errors)),crossK_gene_checks=1300,candidate_rules_checked=90,original_fit_outputs_SHA_pass=True,official_K_statistics_independently_refit=False,biological_validation=False,formal_state_labels=False,UCell_complete=False);(O/'independent_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
