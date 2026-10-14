"""Complete frozen T-programme screen, fixed receiver, patient-level models.
Deduplicate physical signatures; exploratory, not a new clinical validation.
"""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,zarr
from anndata._io.specs import read_elem
from scipy.stats import rankdata


def standard_rank(x):
    r=rankdata(x,axis=0);sd=np.std(r,axis=0)
    return np.divide(r-r.mean(axis=0),sd,out=np.zeros_like(r),where=sd>1e-12),sd>1e-12


def design(cov):
    r=rankdata(cov,axis=0) if cov.shape[1] else cov
    return np.column_stack((np.ones(len(cov)),r))


def bh(p):
    p=np.asarray(p,float);order=np.argsort(p);q=np.empty(len(p))
    q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q


def family(values,y,cov,seed):
    n=len(y);k=values.shape[1];N=design(cov);df=n-N.shape[1]-1
    defaults=[dict(status='rank_or_df_insufficient',beta=np.nan,t=np.nan,p=np.nan,lower=np.nan,upper=np.nan,valid_bootstrap=0,residual_df=df) for _ in range(k)]
    if n<20 or df<10 or np.linalg.matrix_rank(N)<N.shape[1]:return defaults
    X,variable=standard_rank(values);Y,validY=standard_rank(y)
    if not validY:return defaults
    Q,_=np.linalg.qr(N,mode='reduced');R=X-Q@(Q.T@X);ey=Y-Q@(Q.T@Y)
    norm=np.sum(R*R,axis=0);valid=variable&(norm>1e-10)
    beta=np.divide(R.T@ey,norm,out=np.zeros(k),where=norm>1e-10)
    sse=np.maximum(ey@ey-beta**2*norm,0);se=np.sqrt(np.divide(sse/df,norm,out=np.zeros(k),where=norm>1e-10))
    valid&=se>1e-12;t=np.divide(beta,se,out=np.zeros(k),where=se>1e-12)
    rng=np.random.default_rng(seed);extreme=np.zeros(k,int);reduced=Y-ey
    for lo in range(0,20000,250):
        size=min(250,20000-lo)
        perm=np.column_stack([reduced+ey[rng.permutation(n)] for _ in range(size)])
        py=perm-Q@(Q.T@perm);num=R.T@py
        pb=np.divide(num,norm[:,None],out=np.zeros((k,size)),where=norm[:,None]>1e-10)
        psse=np.maximum(np.sum(py*py,axis=0)[None,:]-pb**2*norm[:,None],0)
        pse=np.sqrt(np.divide(psse/df,norm[:,None],out=np.zeros_like(psse),where=norm[:,None]>1e-10))
        pt=np.divide(pb,pse,out=np.full_like(pb,np.inf),where=pse>1e-12)
        extreme+=np.sum(np.abs(pt)>=np.abs(t[:,None])-1e-12,axis=1)
    bootstrap=[[] for _ in range(k)]
    for _ in range(2000):
        ix=rng.integers(0,n,n);BN=design(cov[ix])
        if np.linalg.matrix_rank(BN)<BN.shape[1]:continue
        BX,vx=standard_rank(values[ix]);BY,vy=standard_rank(y[ix])
        if not vy:continue
        BQ,_=np.linalg.qr(BN,mode='reduced');BR=BX-BQ@(BQ.T@BX);BE=BY-BQ@(BQ.T@BY)
        bn=np.sum(BR*BR,axis=0);bb=np.divide(BR.T@BE,bn,out=np.zeros(k),where=bn>1e-10)
        for j in np.flatnonzero(vx&(bn>1e-10)):bootstrap[j].append(float(bb[j]))
    result=[]
    for j in range(k):
        row=defaults[j]
        if valid[j]:
            lo,hi=map(float,np.quantile(bootstrap[j],[.025,.975])) if len(bootstrap[j])>=1800 else (np.nan,np.nan)
            row.update(status='estimated_exploratory_patient_association',beta=float(beta[j]),t=float(t[j]),p=float((1+extreme[j])/20001),
              lower=lo,upper=hi,valid_bootstrap=len(bootstrap[j]),invalid_bootstrap=2000-len(bootstrap[j]),
              residual_predictor_variance=float(norm[j]/n))
        result.append(row)
    return result


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    inputs=root/'HNSCC_frozen_T_sender_states_20261014_v1'
    out=root/'HNSCC_frozen_T_sender_associations_20261014_v1';assert not out.exists();out.mkdir()
    terminal=json.loads((inputs/'private_terminal.json').read_text());assert terminal['candidate_ids']==51
    (out/'execution_contract_before_associations.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),
      producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      preparation_plan_sha256=hashlib.sha256((inputs/'plan_before_sender_scores_or_associations.json').read_bytes()).hexdigest(),
      models=['primary_RNA_capture','sensitivity_unadjusted','sensitivity_demographic','sensitivity_sender_gene_overlap_removed'],
      covariance_transform='all continuous nuisance and scores ranked within cohort; bootstrap recomputes ranks',
      exact_signature_dedup='original and overlap-removed physical signatures deduplicated separately; absent signatures retain frozen gene-set identity',
      inference='20000 shared-patient residual permutations per model;2000 shared-patient bootstraps; fixed seed',
      BH='primary admissible unique signatures; all-candidate diagnostic family also retained; sensitivities distinct families',
      receiver_chosen_after_clinical_results=True,not_independent_clinical_or_direct_LR_validation=True),indent=2)+'\n')
    md=pd.read_csv(inputs/'private_patient_receiver_covariates.tsv',sep='\t',index_col=0)
    scores=pd.read_csv(inputs/'private_patient_T_state_means.tsv',sep='\t',index_col=0).loc[md.index]
    orth=pd.read_csv(inputs/'private_patient_T_orthogonal_state_means.tsv',sep='\t',index_col=0).loc[md.index]
    candidates=pd.read_csv(inputs/'private_frozen_candidate_contract.tsv',sep='\t')
    dict_path=root/'state_programme_dictionary_20261014_v2/frozen_programme_dictionary.json'
    records={r['programme_id']:r for r in json.loads(dict_path.read_text())['programmes']}
    # Reconstruct overlap-removed physical groups only from the immutable source
    # axis and frozen dictionary; no expression/association based grouping.
    contract=json.loads((project/'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json').read_text())
    rec=next(r for r in contract['records'] if r['dataset']=='GSE286827')
    g=zarr.open_group(rec['path'],mode='r');var=read_elem(g[rec['var_key']]);order=list(read_elem(g['uns'])['source_order'])
    measured=np.asarray(g['varm/measured_in_source'])[:,order.index(rec['new_parent_dataset'])].astype(bool)
    lookup={}
    for i,(name,row) in enumerate(var.iterrows()):
        if not measured[i] or str(name).startswith('AMBIGUOUS::'):continue
        names={str(name),str(row.get('gene_symbol',''))}
        if row.get('feature_identity_status')=='stable_Ensembl_ID':names.update(str(row.get('source_gene_aliases','')).split('|'))
        for name in names:
            if name not in ['','nan','None']:lookup.setdefault(name,set()).add(i)
    receiver=pd.read_csv(root/'HNSCC_frozen_MCP_gene_axis_20261014_v2/Macnode18_private_frozen_gene_columns.tsv',sep='\t')
    excluded=set(receiver.loc[receiver.native_column.ge(0),'native_column'].astype(int))
    groups=[]
    for row in candidates.itertuples():
        ids=sorted({next(iter(lookup[gene])) for gene in records[row.programme_id]['primary_genes'] if len(lookup.get(gene,set()))==1}-excluded)
        basis=ids if row.orthogonal_measurable else sorted(records[row.programme_id]['primary_genes'])
        groups.append(hashlib.sha256(json.dumps(basis).encode()).hexdigest())
    candidates['orthogonal_group']=groups
    main_columns=['T_logRNA','T_loggenes','Mac_logRNA','Mac_loggenes','T_fraction','Mac_fraction','T_logcells','Mac_logcells']
    y=md.receiver_programme.to_numpy();main_cov=md[main_columns].to_numpy(float)
    models=[('primary_RNA_capture',scores,main_cov),('sensitivity_unadjusted',scores,np.empty((len(md),0))),
      ('sensitivity_demographic',scores,md[main_columns+['D_plus_T','HPV_positive','male','age']].to_numpy(float)),
      ('sensitivity_sender_gene_overlap_removed',orth,main_cov)]
    all_rows=[]
    for number,(name,frame,cov) in enumerate(models):
        group_col='orthogonal_group' if name.endswith('overlap_removed') else 'observed_signature_group'
        eligible=candidates.groupby(group_col).primary_admissible.any()
        representatives=candidates.sort_values('programme_id').drop_duplicates(group_col)
        existing=representatives[representatives.programme_id.isin(frame.columns)]
        values=frame[existing.programme_id].to_numpy()
        results=family(values,y,cov,20261014+number)
        by_group={group:r for group,r in zip(existing[group_col],results)}
        group_rows=[]
        for group,gp in candidates.groupby(group_col,sort=True):
            row=dict(model=name,signature_group=group,programme_ids='|'.join(sorted(gp.programme_id)),
              primary_admissible=bool(eligible.loc[group]),n_patients=len(md),nuisance_columns=cov.shape[1],
              foreign_RNA_flag=bool(gp.foreign_RNA_sensitive.any()),basic_expression_flag=bool(gp.basic_expression_review_flag.any()))
            row.update(by_group.get(group,dict(status='signature_not_measurable',beta=np.nan,t=np.nan,p=np.nan,lower=np.nan,upper=np.nan,valid_bootstrap=0)))
            group_rows.append(row)
        f=pd.DataFrame(group_rows);f['BH_all_candidate_unique_family']=np.where(f.p.notna(),bh(f.p.fillna(1)),np.nan)
        f['BH_admissible_unique_family']=np.nan;keep=f.primary_admissible.to_numpy()
        f.loc[keep,'BH_admissible_unique_family']=np.where(f.loc[keep,'p'].notna(),bh(f.loc[keep,'p'].fillna(1)),np.nan)
        all_rows.extend(f.to_dict('records'))
    results=pd.DataFrame(all_rows)
    results.to_csv(out/'all_private_sender_state_associations.tsv',sep='\t',index=False,na_rep='NA')
    candidates.to_csv(out/'private_signature_dedup_contract.tsv',sep='\t',index=False)
    report=dict(status='FULL_FIXED_T_SENDER_SCREEN_COMPLETE_EXPLORATORY_PRIVATE',finished_UTC=datetime.now(timezone.utc).isoformat(),
      patient_blocks=len(md),candidate_ids=51,association_rows=len(results),clinical_P_values_added=False,
      formal_states_or_direct_LR_verified=False,original_CM_labels_changed=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
