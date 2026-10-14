"""Independent all matching/rank/R2 readback. No producer imported."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
from scipy.stats import spearmanr

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def group_residual(x,g):
    y=x.copy().astype(float)
    for v in np.unique(g):
        m=g==v;y[m]-=y[m].mean(axis=0)
    return y
def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);D=R/'Spectra_multiseed_diagnostics_20261014_v1';S=R/'Spectra_state_discovery_20261014_v1';I=json.loads((D/'inventory.json').read_text())
    for name,h in I['output_sha256'].items(): assert sha(D/name)==h
    models=[dict(np.load(S/'private_seed_models'/('seed'+str(s))/'factors_and_discovery_scores.npz',allow_pickle=False)) for s in I['seeds']]
    obs=pd.read_csv(S/'private_seed_models/seed20261014/private_discovery_metadata.tsv',sep='\t',index_col=0,dtype={'node_id':str,'gse_id':str})
    assert np.array_equal(obs.index.to_numpy(str),models[0]['cell_ids'])
    pairs=pd.read_csv(D/'pairwise_factor_matching.tsv',sep='\t');sr=pd.read_csv(D/'pairwise_score_reproducibility.tsv',sep='\t');diag=pd.read_csv(D/'source_identity_complexity_diagnostics.tsv',sep='\t');can=pd.read_csv(D/'all173_candidate_diagnostics.tsv',sep='\t')
    dif=[];mapping={}
    for (sa,sb),p in pairs.groupby(['seed_a','seed_b']):
        assert p.factor_a.is_unique and p.factor_b.is_unique and len(p)==173
        a=I['seeds'].index(sa);b=I['seeds'].index(sb);mapping[a,b]=dict(zip(p.factor_a,p.factor_b))
        for row in p.itertuples():
            x=models[a]['factors'][row.factor_a].astype(float);y=models[b]['factors'][row.factor_b].astype(float)
            c=np.dot(x,y)/np.sqrt(np.dot(x,x)*np.dot(y,y))
            xx=set(sorted(range(len(x)),key=lambda v:(-x[v],v))[:50]);yy=set(sorted(range(len(y)),key=lambda v:(-y[v],v))[:50]);j=len(xx&yy)/len(xx|yy)
            assert abs(c-row.cosine)<1e-6 and abs(j-row.top50_jaccard)<1e-12
            assert models[a]['factor_scope'][row.factor_a]==models[b]['factor_scope'][row.factor_b]
            dif.append(abs(c-row.cosine))
    ranks=[]
    for row in sr.itertuples():
        a=I['seeds'].index(row.seed_a);b=I['seeds'].index(row.seed_b);m=obs.major.eq(row.major).to_numpy();assert m.sum()==row.cells
        x=models[a]['cell_scores'][m,row.factor_a].astype(float);y=models[b]['cell_scores'][m,row.factor_b].astype(float);g=obs.gse_id.to_numpy()[m]
        c=spearmanr(x,y).statistic;d=spearmanr(group_residual(x,g),group_residual(y,g)).statistic
        for actual,reported in [(c,row.score_spearman),(d,row.source_residual_spearman)]:
            assert np.isclose(actual,reported,atol=2e-6,equal_nan=True);ranks.append(abs(actual-reported))
    rdiff=[]
    for row in diag.itertuples():
        m=obs.major.eq(row.major).to_numpy();x=models[0]['cell_scores'][m,row.anchor_factor].astype(float);meta=obs.loc[m];g=meta.gse_id.to_numpy();gn=(meta.gse_id.astype(str)+'::'+meta.node_id.astype(str)).to_numpy();den=np.sum((x-x.mean())**2)
        for groups,reported in [(g,row.source_r2),(gn,row.source_node_r2)]:
            actual=1-np.sum(group_residual(x,groups)**2)/den;assert np.isclose(actual,reported,atol=2e-7,equal_nan=True);rdiff.append(abs(actual-reported))
    plan=json.loads(Path(os.environ['PANCANCER_SPECTRA_DIAGNOSTIC_PLAN']).read_text())
    for row in can.itertuples():
        j=row.anchor_factor; js=[j,mapping[0,1][j],mapping[0,2][j]];cyc=mapping[1,2][js[1]]==js[2];p=[];s=[]
        for a,b in [(0,1),(0,2),(1,2)]:
            p.extend(pairs[(pairs.seed_a==I['seeds'][a])&(pairs.seed_b==I['seeds'][b])&(pairs.factor_a==js[a])].to_dict('records'))
            s.extend(sr[(sr.seed_a==I['seeds'][a])&(sr.seed_b==I['seeds'][b])&(sr.factor_a==js[a])].to_dict('records'))
        rho=[x['source_residual_spearman'] for x in s]
        stable=bool(cyc and min(x['cosine'] for x in p)>=plan['min_cosine'] and min(x['top50_jaccard'] for x in p)>=plan['min_top50_jaccard'] and np.isfinite(rho).all() and min(rho)>=plan['min_source_residual_score_rho'])
        assert stable==row.stable_by_frozen_rule and cyc==row.cycle_consistent
    O=R/'Spectra_multiseed_diagnostics_audit_20261014_v1';assert not O.exists();O.mkdir()
    # A zero overlap returned a lexical prior name in v1; correct only new display table.
    corrected=can.copy();corrected.loc[corrected.best_prior_top50_overlap.eq(0),'best_prior_descriptor']='NO_PRIOR_OVERLAP_IN_TOP50'
    corrected['descriptor_is_confirmed_state']=False;corrected.to_csv(O/'all173_candidates_interpretation_safe.tsv',sep='\t',index=False)
    out=dict(status='PASS_ALL519_MATCHES_ALL957_SCORE_COMPARISONS_ALL_R2_AND173_STABILITY_RULES',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(D/'inventory.json'),cosine_checks=len(dif),max_cosine_difference=max(dif),score_statistics_checked=len(ranks),max_score_difference=float(np.nanmax(ranks)),source_node_r2_statistics_checked=len(rdiff),max_r2_difference=float(np.nanmax(rdiff)),stable_candidates=int(can.stable_by_frozen_rule.sum()),zero_prior_overlap_descriptor_correction=int(can.best_prior_top50_overlap.eq(0).sum()),biological_program_identity_validated=False,formal_state_labels=False,UCell_complete=False,corrected_table_sha256=sha(O/'all173_candidates_interpretation_safe.tsv'))
    (O/'independent_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
