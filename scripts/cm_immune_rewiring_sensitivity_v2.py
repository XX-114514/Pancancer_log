"""All frozen-member covariance sensitivities; not new identity or validation."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,json,hashlib,itertools
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr
from statsmodels.stats.multitest import multipletests
from cm_immune_recovery_diagnostic_v1 import IMMUNE,COMPARTMENT,SOURCES,sha,save,bh,seed

def grouped_null(V,groups,pairs,key):
    rng=np.random.default_rng(seed(key));Z=rankdata(V,axis=0);Z-=Z.mean(axis=0);Z/=np.sqrt((Z*Z).sum(axis=0))
    obs=float(np.mean([Z[:,i]@Z[:,j] for i,j in pairs]));null=[];checks=[]
    for start in range(0,20000,500):
        Q=np.empty((500,*Z.shape));U=np.empty((500,*V.shape)) if start==0 else None
        for group in sorted(set(groups)):
            ix=np.flatnonzero(groups==group);per=np.argsort(rng.random((500,len(V))),axis=1)
            Q[:,:,ix]=Z[per][:,:,ix]
            if U is not None:U[:,:,ix]=V[per][:,:,ix]
        scores=np.mean(np.array([np.sum(Q[:,:,i]*Q[:,:,j],axis=1) for i,j in pairs]),axis=0);null.extend(scores)
        if start==0:
            for j in range(4):
                direct=float(np.mean([spearmanr(U[j,:,i],U[j,:,k]).statistic for i,k in pairs]));assert abs(direct-scores[j])<1e-12;checks.append(abs(direct-scores[j]))
    null=np.asarray(null);direct=float(np.mean([spearmanr(V[:,i],V[:,j]).statistic for i,j in pairs]));assert abs(direct-obs)<1e-12
    ii=rng.integers(len(V),size=(2000,len(V)));B=rankdata(V[ii],axis=1);B-=B.mean(axis=1,keepdims=True);sd=np.sqrt((B*B).sum(axis=1));valid=(sd>1e-12).all(axis=1)
    B[valid]/=sd[valid,None,:];bs=np.mean(np.array([np.sum(B[valid,:,i]*B[valid,:,j],axis=1) for i,j in pairs]),axis=0)
    lo,hi=np.quantile(bs,[.025,.975]) if valid.sum()>=1800 else [np.nan,np.nan]
    p=(1+np.sum(null>=obs-1e-12))/20001
    return dict(mean_rho=obs,permutation_p=float(p),bootstrap_low=lo,bootstrap_high=hi,bootstrap_valid=int(valid.sum()),direct_audit_max_diff=float(max(checks+[abs(direct-obs)]))),null,bs

def parent(m):
    if m in {'Macrophage','Monocyte','Myeloid','Dendritic cell','Neutrophil','Mast cell'}:return 'myeloid_analysis'
    if m in {'T cell','T/NK lymphoid','Lymphoid cell'}:return 'lymphoid_analysis'
    return m

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--frozen-root',required=True);ap.add_argument('--review-root',required=True);ap.add_argument('--output',required=True);ap.add_argument('--plan',required=True);a=ap.parse_args()
    b,r,out=Path(a.frozen_root),Path(a.review_root),Path(a.output);assert not out.exists();out.mkdir();p=out/'public';p.mkdir();private=out/'private';private.mkdir();started=datetime.now(timezone.utc).isoformat()
    cat=pd.read_csv(b/'downstream_review_v2/number_to_annotation_dictionary.tsv',sep='\t').set_index('subCluster')
    meta=pd.read_csv(r/'immune_same_name_family_source_review_v1/patient_block_provenance_overlay_v1/private_all1290_patient_provenance_overlay.tsv',sep='\t').set_index('sampleID')
    cnt=pd.read_csv(b/'inputs/sample_major_subcluster_counts.tsv',sep='\t');cnt=cnt[cnt.sampleID.isin(meta.index)]
    N=cnt.pivot_table(index='sampleID',columns='subCluster',values='n_cells',aggfunc='sum',fill_value=0).reindex(index=meta.index,columns=cat.index,fill_value=0)
    assert np.array_equal(N.sum(axis=1).to_numpy(),meta.n_cells.to_numpy())
    e=pd.read_csv(b/'downstream_review_v2/corrected_CM_edges_source_and_marker_evidence.tsv',sep='\t');cms=sorted(e.cm.unique())
    cores={cm:sorted(set(e.loc[e.cm.eq(cm),'subCluster1'])|set(e.loc[e.cm.eq(cm),'subCluster2'])) for cm in cms}
    rows=[];cache={};degen=[]
    for (gse,cancer,material,phase),md in meta[meta.gse_id.isin(SOURCES)].groupby(['gse_id','cancerType','sample_class','treatment_phase'],dropna=False,sort=True):
        key='|'.join(map(str,[gse,cancer,material,phase]));n=N.loc[md.index];allowed=set(cat.index)
        if gse=='GSE179994':allowed=set(cat.index[cat.majorCluster.isin({'T cell','T/NK lymphoid','Lymphoid cell'})])
        for mode in ['legacy_major','parent_analysis']:
            group=cat.majorCluster.map(parent) if mode=='parent_analysis' else cat.majorCluster
            den={m:n.loc[:,group.eq(m)].sum(axis=1) for m in group.unique()}
            fv=pd.DataFrame({node:n[node]/den[group[node]].replace(0,np.nan) for node in cat.index},index=n.index)
            opportunity={m:bool((total>=20).mean()>=.5) for m,total in den.items()}
            f=fv.groupby(md.effective_resampling_block).mean();nc=n.groupby(md.effective_resampling_block).mean()
            for cm in cms:
                for view in ['whole_core','immune_backbone']:
                    core=[node for node in cores[cm] if view=='whole_core' or cat.loc[node,'majorCluster'] in IMMUNE]
                    for threshold in [1,20]:
                        row=dict(cohort=key,cm=cm,view=view,denominator=mode,node_cell_cut=threshold,reference_members=len(core),status='insufficient_variable_members',mean_rho=np.nan,permutation_p=np.nan,bootstrap_low=np.nan,bootstrap_high=np.nan,bootstrap_valid=0,independent_validation=False)
                        measured=[node for node in core if node in allowed and opportunity[group[node]]]
                        complete=f[measured].dropna() if measured else pd.DataFrame()
                        variable=[node for node in measured if len(complete)>=10 and (nc.loc[complete.index,node]>=threshold).mean()>=.5 and (nc.loc[complete.index,node]>=threshold).sum()>=10 and np.ptp(complete[node])>1e-12]
                        row.update(complete_blocks=len(complete),variable_members=len(variable),unmeasured_members=len(core)-len(measured),constant_members=sum(np.ptp(complete[node])<=1e-12 for node in measured) if len(complete) else 0)
                        compartments={COMPARTMENT[cat.loc[node,'majorCluster']] for node in variable if cat.loc[node,'majorCluster'] in COMPARTMENT}
                        pair=[(s,t) for s,t in itertools.combinations(variable,2) if group[s]!=group[t]];row['candidate_pairs']=len(pair)
                        if len(compartments)<2:row['status']='fewer2_unambiguous_compartments'
                        elif len(pair)>=3:
                            used=sorted({node for pp in pair for node in pp});V=complete[used].to_numpy();gr=group.loc[used].to_numpy();pi=[(used.index(s),used.index(t)) for s,t in pair]
                            ident=hashlib.sha256((key+'|'+mode+'|'+';'.join(used)+'|'+str(pi)).encode()+V.tobytes()).hexdigest()
                            if ident not in cache:
                                stats,null,bs=grouped_null(V,gr,pi,key+'|'+ident);cache[ident]=stats
                                np.savez_compressed(private/(ident+'.npz'),values=V,groups=gr,pairs=np.asarray(pi),null=null,bootstrap=bs)
                            row.update(cache[ident],status='estimated_posthoc_diagnostic',statistical_key=ident)
                        rows.append(row)
                        if threshold==1:degen.append({k:row[k] for k in ['cohort','cm','view','denominator','reference_members','variable_members','constant_members','candidate_pairs','status']})
        print('REWIRING_COHORT_COMPLETE',gse,material,phase,flush=True)
    d=pd.DataFrame(rows);assert len(d)==432;rawp=d.permutation_p.fillna(1).to_numpy();q=bh(rawp);by=multipletests(rawp,method='fdr_by')[1];np.testing.assert_allclose(q,multipletests(rawp,method='fdr_bh')[1],atol=1e-14)
    d['diagnostic_BH']=np.where(d.permutation_p.notna(),q,np.nan);d['diagnostic_BY']=np.where(d.permutation_p.notna(),by,np.nan);d['positive_CI_BH_support']=d.diagnostic_BH.lt(.05)&d.bootstrap_low.gt(0)
    d.fillna('NA').to_csv(p/'all432_rewiring_sensitivities.tsv',sep='\t',index=False);pd.DataFrame(degen).to_csv(p/'all216_variable_member_diagnostics.tsv',sep='\t',index=False)
    summary=dict(status='COMPLETE_POSTHOC_REWIRING_DIAGNOSTIC_NOT_INDEPENDENT_VALIDATION',started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),configured_records=432,estimated_records=int(d.permutation_p.notna().sum()),unique_estimated_questions=len(cache),BH05=int(d.diagnostic_BH.lt(.05).sum()),BY05=int(d.diagnostic_BY.lt(.05).sum()),positive_CI_BH=int(d.positive_CI_BH_support.sum()),status_counts=d.status.value_counts().to_dict(),permutations=20000,bootstrap=2000,selected_numeric_audit_max_diff=max([x['direct_audit_max_diff'] for x in cache.values()],default=0),plan_sha256=sha(a.plan),code_sha256=sha(__file__),canonical_changed=False,formal_states=False,axis_validation_complete=False)
    save(p/'summary_v2.json',summary);print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
