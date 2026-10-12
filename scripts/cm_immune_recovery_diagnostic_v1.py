"""Frozen CM observability, fixed-W recovery and internal covariance checks.

No new identities, patient blocks, modules or causal claims. Private fixtures are
written outside public/; only cohort/module summaries may be published.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import rankdata

SEED, PERM, BOOT = 2026101209, 1999, 2000
SOURCES = ['GSE115978', 'GSE127465', 'GSE179994', 'GSE274229']
IMMUNE = {'B cell', 'T cell', 'T/NK lymphoid', 'Lymphoid cell',
          'Dendritic cell', 'Macrophage', 'Myeloid', 'Monocyte',
          'Neutrophil', 'Mast cell', 'Hematopoietic cell'}
COMPARTMENT = {'B cell': 'B', 'T cell': 'T', 'Macrophage': 'myeloid',
               'Dendritic cell': 'myeloid', 'Myeloid': 'myeloid',
               'Monocyte': 'myeloid', 'Neutrophil': 'myeloid',
               'Mast cell': 'myeloid', 'Endothelial cell': 'vascular',
               'Fibroblast': 'stromal', 'Mural cell': 'stromal',
               'Mesenchymal cell': 'stromal', 'Epithelial cell': 'epithelial',
               'Melanocytic cell': 'melanocytic', 'Glial cell': 'glial'}

def now(): return datetime.now(timezone.utc).isoformat()
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def seed(key): return (SEED + int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)) % 2**32
def save(p, obj): p.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
def bh(p):
    p = np.asarray(p); ix = np.argsort(p); q = np.empty(len(p))
    q[ix] = np.minimum(1, np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q

def kl(W, V, H):
    R = np.maximum(W @ H, 1e-15)
    return float(np.sum(R - V) + np.sum(V[V > 0] * np.log(V[V > 0]/R[V > 0])))

def project(W, V, key):
    """Multiplicative fixed-basis KL recovery, with explicit convergence data."""
    scale = W.sum(axis=0); A = W / scale
    rng = np.random.default_rng(seed(key)); best = None
    for attempt in range(3):
        H = np.maximum(V.sum(axis=0)[None, :]/W.shape[1], 1e-5) * rng.uniform(.7,1.3,(W.shape[1],V.shape[1]))
        old = np.inf
        for iteration in range(20000):
            H *= (A.T @ (V/np.maximum(A@H,1e-15))) / A.sum(axis=0)[:,None]
            H = np.maximum(H, 1e-15)
            if iteration % 20 == 0:
                loss = kl(A,V,H)
                grad = A.T @ (1-V/np.maximum(A@H,1e-15))
                kkt = float(np.max(np.abs(np.where(H>1e-7,grad,np.minimum(grad,0)))))
                if abs(old-loss) <= 1e-10*(1+abs(loss)) and kkt < 2e-5: break
                old = loss
        loss = kl(A,V,H)
        if best is None or loss < best[0]: best = (loss,H.copy(),iteration+1)
    loss,H,nit = best
    grad = A.T @ (1-V/np.maximum(A@H,1e-15))
    kkt = float(np.max(np.abs(np.where(H>1e-7,grad,np.minimum(grad,0)))))
    audits = []
    # SLSQP uses a separately written scalar objective/derivative.
    for j in range(min(3,V.shape[1])):
        v = V[:,j]
        def obj(h):
            z = np.maximum(np.sum(A*h[None,:],axis=1),1e-15)
            value = float(np.sum(z-v))
            nz = v > 0
            return value + float(np.sum(v[nz]*(np.log(v[nz])-np.log(z[nz]))))
        def jac(h):
            z = np.maximum(np.sum(A*h[None,:],axis=1),1e-15)
            return np.array([np.sum(A[:,k]*(1-v/z)) for k in range(A.shape[1])])
        fit = minimize(obj, np.full(A.shape[1],max(v.sum()/A.shape[1],.01)),
                       jac=jac, bounds=[(0,None)]*A.shape[1],method='SLSQP',
                       options={'maxiter':3000,'ftol':1e-11})
        gap = obj(H[:,j])-fit.fun
        audits.append(dict(fixture=j, success=bool(fit.success), objective_gap=float(gap),
                           production_objective=float(obj(H[:,j])), audit_objective=float(fit.fun)))
    passed = kkt < 1e-3 and all(abs(x['objective_gap']) < 1e-5*(1+abs(x['audit_objective'])) for x in audits)
    return H/scale[:,None], dict(loss=loss, iterations=nit, max_KKT=kkt,
                                numeric_fixture_pass=passed, independent_SLSQP=audits)

def mean_edge_rho(data, pairs):
    Z = rankdata(data, axis=0); Z -= Z.mean(axis=0)
    norm = np.sqrt(np.sum(Z*Z,axis=0))
    if np.any(norm <= 1e-12): return np.nan
    Z /= norm
    return float(np.mean([Z[:,i] @ Z[:,j] for i,j in pairs]))

def coherence(values, majors, pairs, key):
    Z = rankdata(values,axis=0); Z -= Z.mean(axis=0); Z /= np.sqrt(np.sum(Z*Z,axis=0))
    observed = float(np.mean([Z[:,i]@Z[:,j] for i,j in pairs]))
    rng = np.random.default_rng(seed(key)); groups = [np.flatnonzero(majors==m) for m in sorted(set(majors))]
    null = []
    for _ in range(PERM):
        Q = np.empty_like(Z)
        for ix in groups: Q[:,ix] = Z[rng.permutation(len(Z))][:,ix]
        null.append(float(np.mean([Q[:,i]@Q[:,j] for i,j in pairs])))
    p = (1+sum(v>=observed-1e-12 for v in null))/(PERM+1)
    boots = np.array([mean_edge_rho(values[rng.integers(len(values),size=len(values))],pairs) for _ in range(BOOT)])
    good = boots[np.isfinite(boots)]
    interval = np.quantile(good,[.025,.975]) if len(good)>=1800 else [np.nan,np.nan]
    # Independent scipy pairwise Spearman replay on observed and first32 nulls.
    from scipy.stats import spearmanr
    direct = np.mean([spearmanr(values[:,i],values[:,j]).statistic for i,j in pairs])
    assert abs(direct-observed)<1e-12
    return dict(mean_edge_rho=observed, permutation_p=p, bootstrap_valid=len(good),
                bootstrap_low=interval[0], bootstrap_high=interval[1]), np.array(null), boots

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--frozen-root',required=True)
    ap.add_argument('--review-root',required=True); ap.add_argument('--output',required=True)
    ap.add_argument('--plan',required=True); a = ap.parse_args()
    b,r,out = Path(a.frozen_root),Path(a.review_root),Path(a.output)
    assert not out.exists(); out.mkdir(); public=out/'public';public.mkdir();private=out/'private';private.mkdir()
    started=now(); plan=json.loads(Path(a.plan).read_text()); assert plan['sources']==SOURCES
    paths=dict(meta=b/'inputs/sample_metadata.tsv', counts=b/'inputs/sample_major_subcluster_counts.tsv',
               catalog=b/'downstream_review_v2/number_to_annotation_dictionary.tsv',
               W=b/'covarnet/tables/nmf_W_subcluster_by_cm.tsv', H=b/'covarnet/tables/nmf_H_cm_by_sample.tsv',
               F=b/'covarnet/tables/frequency_normalized.tsv',
               edges=b/'downstream_review_v2/corrected_CM_edges_source_and_marker_evidence.tsv',
               provenance=r/'immune_same_name_family_source_review_v1/patient_block_provenance_overlay_v1/private_all1290_patient_provenance_overlay.tsv')
    rawmeta=pd.read_csv(paths['meta'],sep='\t').set_index('sampleID');meta=rawmeta[rawmeta.eligible_min100_cells].copy()
    overlay=pd.read_csv(paths['provenance'],sep='\t').set_index('sampleID').reindex(meta.index)
    assert overlay.gse_id.notna().all() and len(meta)==1290
    for col in ['sample_class','treatment_phase','cancerType','effective_resampling_block','unit_kind']:
        if col in overlay: meta[col]=overlay[col]
    assert meta.effective_resampling_block.notna().all()
    cat=pd.read_csv(paths['catalog'],sep='\t').set_index('subCluster')
    counts=pd.read_csv(paths['counts'],sep='\t');counts=counts[counts.sampleID.isin(meta.index)]
    node=counts.pivot_table(index='sampleID',columns='subCluster',values='n_cells',aggfunc='sum',fill_value=0).reindex(index=meta.index,columns=cat.index,fill_value=0)
    major=counts.groupby(['sampleID','majorCluster']).n_cells.sum().unstack(fill_value=0).reindex(meta.index)
    assert np.array_equal(node.sum(axis=1).to_numpy(),meta.n_cells.to_numpy())
    denom=np.column_stack([major[cat.loc[n,'majorCluster']].to_numpy() for n in node.columns])
    freq=pd.DataFrame(np.divide(node.to_numpy(),denom,out=np.full(node.shape,np.nan),where=denom>0),index=node.index,columns=node.columns)
    W=pd.read_csv(paths['W'],sep='\t',index_col=0).reindex(cat.index)
    F=pd.read_csv(paths['F'],sep='\t',index_col=0).reindex(index=cat.index,columns=meta.index)
    H0=pd.read_csv(paths['H'],sep='\t',index_col=0).reindex(index=W.columns,columns=meta.index)
    assert W.shape==(172,9) and not F.isna().any().any()
    f0=freq.fillna(0).T.to_numpy(); low=f0.min(axis=1); span=f0.max(axis=1)-low
    recomputed=np.divide(f0-low[:,None],span[:,None],out=np.zeros_like(f0),where=span[:,None]>0)
    parity=float(np.max(abs(recomputed-F.to_numpy()))); assert parity < 1e-8
    edges=pd.read_csv(paths['edges'],sep='\t'); assert len(edges)==120
    cores={cm: sorted(set(edges.loc[edges.cm.eq(cm),'subCluster1'])|set(edges.loc[edges.cm.eq(cm),'subCluster2']), key=lambda n:-W.loc[n,cm]) for cm in W.columns}
    core_rows=[]
    for cm,nodes in cores.items():
        total=W.loc[nodes,cm].sum()
        for n in nodes:core_rows.append(dict(cm=cm,node=n,major=cat.loc[n,'majorCluster'],weight=float(W.loc[n,cm]),core_weight=float(W.loc[n,cm]/total),immune=cat.loc[n,'majorCluster'] in IMMUNE))
    pd.DataFrame(core_rows).to_csv(public/'fixed_CM_graph_core.tsv',sep='\t',index=False)
    cohort_rows=[]; results=[]; projections=[]; grids=[]
    for (gse,cancer,material,phase), md in meta[meta.gse_id.isin(SOURCES)].groupby(['gse_id','cancerType','sample_class','treatment_phase'],dropna=False,sort=True):
        key='|'.join(map(str,[gse,cancer,material,phase])); ix=md.index; nc=node.loc[ix].copy(); ma=major.loc[ix].copy()
        allowed=set(cat.majorCluster)
        if gse=='GSE179994':allowed={'T cell','T/NK lymphoid','Lymphoid cell'}
        allowed_nodes=cat.index[cat.majorCluster.isin(allowed)]
        nc.loc[:,~nc.columns.isin(allowed_nodes)]=0
        cohort_rows.append(dict(cohort=key,gse_id=gse,cancer=cancer,material=material,phase=phase,libraries=len(md),blocks=md.effective_resampling_block.nunique(),reviewed_patient_libraries=int(md.unit_kind.eq('reviewed_patient').sum()),cells=int(md.n_cells.sum()),discovery_overlap=True,whole_TME_eligible=gse!='GSE179994'))
        fx=freq.loc[ix].copy(); fx.loc[:,~fx.columns.isin(allowed_nodes)]=np.nan
        blockF=fx.groupby(md.effective_resampling_block).mean(); blockN=nc.groupby(md.effective_resampling_block).mean()
        # Use major sampling opportunity, not zeros for absent compartments.
        opportunity={m:bool((ma[m]>=20).mean()>=.5) and m in allowed for m in ma.columns}
        vary=np.ptp(F.loc[:,ix].to_numpy(),axis=1)>1e-12
        mask=vary & np.array([opportunity.get(cat.loc[n,'majorCluster'],False) for n in cat.index])
        AA=W.to_numpy()[mask]; rank=int(np.linalg.matrix_rank(AA)) if len(AA) else 0
        cond=float(np.linalg.cond(AA/np.maximum(AA.sum(axis=0),1e-300))) if len(AA)>=9 else np.inf
        recovered={}
        for norm in ['frozen_minmax','cohort_minmax']:
            pr=dict(cohort=key,normalization=norm,rows=int(mask.sum()),rank=rank,condition=cond if np.isfinite(cond) else np.nan,status='unidentifiable_masked_W')
            if rank==9 and cond<1e8:
                VV=F.loc[:,ix].to_numpy()[mask].copy()
                if norm=='cohort_minmax':
                    vv=f0[:,meta.index.get_indexer(ix)][mask];lo=vv.min(axis=1);sp=vv.max(axis=1)-lo;VV=(vv-lo[:,None])/sp[:,None]
                h,check=project(AA,VV,key+'|'+norm);recovered[norm]=h
                pr.update(status='fixed_W_numeric_checked' if check['numeric_fixture_pass'] else 'numeric_check_not_passed',loss=check['loss'],iterations=check['iterations'],max_KKT=check['max_KKT'],fixture_pass=check['numeric_fixture_pass'])
                save(private/('projection_audit_'+hashlib.sha256((key+norm).encode()).hexdigest()[:12]+'.json'),check)
                np.savez_compressed(private/('activities_'+hashlib.sha256((key+norm).encode()).hexdigest()[:12]+'.npz'),H=h,samples=np.asarray(ix,str),W=AA,V=VV)
            projections.append(pr)
        for cm,core in cores.items():
            cmaj=cat.loc[core,'majorCluster'];cw=W.loc[core,cm]/W.loc[core,cm].sum();cn=nc[core]
            present=cn>=20; weighted=present.to_numpy()@cw.to_numpy()
            ncore=present.sum(axis=1).to_numpy()
            nmajor=np.array([cmaj[present.iloc[i]].nunique() for i in range(len(md))])
            comp=np.array([len({COMPARTMENT[m] for m in cmaj[present.iloc[i]] if m in COMPARTMENT}) for i in range(len(md))])
            mincore=int(np.ceil(len(core)*.5));minmajor=min(3,cmaj.nunique());mincomp=min(2,len({COMPARTMENT[m] for m in cmaj if m in COMPARTMENT}))
            structural=(ncore>=mincore)&(nmajor>=minmajor)&(comp>=mincomp)
            sup=structural&(weighted>=.6); sup_blocks=md.loc[sup,'effective_resampling_block'].nunique()
            coverage=bool(sup.mean()>=.5 and sup_blocks>=10)
            rr=dict(cohort=key,cm=cm,n_core=len(core),core_immune_weight=float(cw[cmaj.isin(IMMUNE)].sum()),libraries=len(md),blocks=md.effective_resampling_block.nunique(),reviewed_patient_libraries=int(md.unit_kind.eq('reviewed_patient').sum()),median_core_ge20=float(np.median(ncore)),median_core_weight=float(np.median(weighted)),median_major=float(np.median(nmajor)),supported_libraries=int(sup.sum()),supported_blocks=int(sup_blocks),support_fraction=float(sup.mean()),coverage_gate=coverage,coherence_status='insufficient_capture',mean_edge_rho=np.nan,permutation_p=np.nan,bootstrap_low=np.nan,bootstrap_high=np.nan,independent_validation=False)
            for cut in [.4,.6,.8]:
                p=structural&(weighted>=cut)
                for fraction in [.5,.8]:grids.append(dict(cohort=key,cm=cm,weight_cut=cut,library_fraction_cut=fraction,supported_fraction=float(p.mean()),supported_blocks=int(md.loc[p,'effective_resampling_block'].nunique()),gate=bool(p.mean()>=fraction and md.loc[p,'effective_resampling_block'].nunique()>=10)))
            if 'frozen_minmax' in recovered:
                h=recovered['frozen_minmax'][W.columns.get_loc(cm)]
                rr['median_recovered_activity']=float(np.median(h));rr['score_positive_fraction']=float(np.mean(h>1e-6))
                rr['score_without_coverage']=bool(np.any(h>1e-6) and not coverage)
            em=edges[edges.cm.eq(cm)]
            # Fixed edges qualify by coverage/variation only; not direction or p.
            usable=[n for n in core if opportunity.get(cat.loc[n,'majorCluster'],False)]
            complete=blockF[usable].dropna() if usable else pd.DataFrame()
            if coverage and len(complete)>=10:
                bn=blockN.loc[complete.index,usable]; uv=[n for n in usable if (bn[n]>=20).mean()>=.5 and (bn[n]>=20).sum()>=10 and np.ptp(complete[n])>1e-12]
                pairs=[(x.subCluster1,x.subCluster2) for x in em.itertuples() if x.subCluster1 in uv and x.subCluster2 in uv and cat.loc[x.subCluster1,'majorCluster']!=cat.loc[x.subCluster2,'majorCluster']]
                pairs=sorted(set(tuple(sorted(p)) for p in pairs));used=sorted({n for pair in pairs for n in pair})
                rr.update(coherence_blocks=len(complete),eligible_edges=len(pairs))
                if len(pairs)>=3:
                    vv=complete[used].to_numpy(); pi=[(used.index(s),used.index(t)) for s,t in pairs]
                    stat,null,boot=coherence(vv,cat.loc[used,'majorCluster'].to_numpy(),pi,key+'|'+cm)
                    rr.update(stat,coherence_status='estimated_diagnostic')
                    np.savez_compressed(private/('coherence_'+hashlib.sha256((key+cm).encode()).hexdigest()[:12]+'.npz'),values=vv,pairs=np.asarray(pi),null=null,bootstrap=boot)
                else:rr['coherence_status']='too_few_variable_crossmajor_edges'
            elif coverage:rr['coherence_status']='too_few_complete_blocks'
            results.append(rr)
        print('COHORT_COMPLETE',gse,material,phase,len(md),'masked_rank',rank,flush=True)
    result=pd.DataFrame(results);q=bh(result.permutation_p.fillna(1).to_numpy());result['diagnostic_BH_q']=np.where(result.permutation_p.notna(),q,np.nan)
    result['internal_covariance_supported']=(result.coverage_gate & result.diagnostic_BH_q.lt(.05)&result.bootstrap_low.gt(0))
    for name,table in [('cohort_scope.tsv',cohort_rows),('all_CM_cohort_results.tsv',result),('all_projection_diagnostics.tsv',projections),('all_capture_threshold_sensitivities.tsv',grids)]:
        pd.DataFrame(table).fillna('NA').to_csv(public/name,sep='\t',index=False)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    pivot=result.pivot(index='cm',columns='cohort',values='median_core_weight');fig,ax=plt.subplots(figsize=(11,6),layout='constrained')
    im=ax.imshow(pivot.to_numpy(),vmin=0,vmax=1,cmap='YlGnBu',aspect='auto')
    ax.set_xticks(range(len(pivot.columns)),[x.replace('|','\n') for x in pivot.columns],fontsize=8)
    ax.set_yticks(range(len(pivot.index)),pivot.index);ax.set_title('Frozen CM core capture: median weight with >=20 labelled cells')
    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):ax.text(j,i,f'{pivot.iloc[i,j]:.2f}',ha='center',va='center',color='white' if pivot.iloc[i,j]>.6 else 'black')
    fig.colorbar(im,ax=ax,label='Observed core-weight fraction');fig.savefig(public/'CM_core_capture.png',dpi=180);plt.close(fig)
    summary=dict(status='COMPLETE_DIAGNOSTIC_INPUT_PARITY_AND_SELECTED_NUMERIC_CHECKS',started_UTC=started,finished_UTC=now(),sources=SOURCES,libraries=sum(x['libraries'] for x in cohort_rows),cohorts=len(cohort_rows),CM_tests=len(result),coverage_pass=int(result.coverage_gate.sum()),coherence_estimated=int(result.permutation_p.notna().sum()),internal_covariance_supported=int(result.internal_covariance_supported.sum()),coherence_status_counts=result.coherence_status.value_counts().to_dict(),frequency_reconstruction_max_difference=parity,projection_numeric_failures=sum(x['status']=='numeric_check_not_passed' for x in projections),input_sha256={k:sha(v) for k,v in paths.items()},plan_sha256=sha(a.plan),code_sha256=sha(__file__),independent_biological_validation=False,all_coherence_statistics_independent_replay=False,canonical_changed=False,formal_states=False,axis_validation_complete=False)
    save(public/'summary.json',summary);print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
