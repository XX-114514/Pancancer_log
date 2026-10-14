"""Descriptive nuisance/LOO diagnostics of fixed partial scores, no new p tests.
All scientific output stays server; primary source tests remain unchanged.
"""
import hashlib,json,os
from pathlib import Path
from datetime import datetime,timezone

import h5py,numpy as np,pandas as pd
from anndata._io.specs import read_elem
from scipy.stats import rankdata


def rng_seed(key):return 20261014+int(hashlib.sha256(key.encode()).hexdigest()[:8],16)


def partial(x,y,cov):
    X=np.column_stack((np.ones(len(x)),np.apply_along_axis(rankdata,0,cov))) if cov.shape[1] else np.ones((len(x),1))
    if len(x)-X.shape[1]<6 or np.linalg.matrix_rank(X)!=X.shape[1]:return np.nan
    x=rankdata(x);y=rankdata(y)
    ex=x-X@np.linalg.lstsq(X,x,rcond=None)[0];ey=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
    return float(np.corrcoef(ex,ey)[0,1]) if min(ex.std(),ey.std())>1e-12 else np.nan


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1'
    parent=crc/'private_existing_source_transfer_observed_part_v1';out=crc/'private_partial_source_covariates_v1'
    assert not out.exists();out.mkdir()
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      fixed_descriptive_variants=['unadjusted','both_fine_major_fractions','both_log_cell_counts','four_mean_cell_quality_covariates','joint_all8'],
      min_residual_df=6,leave_one_patient=True,no_new_p_values_or_BH=True,
      limitation='descriptive post-result nuisance diagnosis, does not replace frozen primary or earlier20-patient inference gate')
    (out/'plan_before_diagnostics.json').write_text(json.dumps(plan,indent=2)+'\n')
    score=pd.read_csv(parent/'private_patient_scores.tsv',sep='\t');rows=[];leave=[]
    labels=['CD8_T_cell_NKG7_9','Macrophage_IFI6_18'];types=['Tnode9','Macnode18']
    for source in ['GSE148071','GSE161529','GSE182109','GSE183904']:
        selected=pd.read_csv(parent/(source+'_private_selected_libraries.tsv'),sep='\t')
        ss=score[score.source.eq(source)&score.program.eq('MCP1')]
        if ss.empty:
            for variant in plan['fixed_descriptive_variants']:rows.append(dict(source=source,variant=variant,status='no_reviewed_scored_patient_pair',n_patients=0,rho=np.nan))
            continue
        paired=ss.pivot(index='analysis_block',columns='endpoint',values='score');assert paired.notna().all().all()
        with h5py.File(root/'state_singlecell_inputs_20261014_v1/private_native_sources'/(source+'.h5ad')) as f:obs=read_elem(f['obs'])
        major=obs.groupby(['sample_uid','major'],observed=True).size();quality={}
        for endpoint,label in zip(types,labels):
            candidates=obs[obs.sample_uid.isin(selected.sample_uid)&obs.subCluster.eq(label)]
            rng=np.random.default_rng(rng_seed(source+endpoint+'cells'))
            for block,group in candidates.groupby('analysis_block',observed=True,sort=True):
                ix=rng.choice(len(group),min(50,len(group)),replace=False);g=group.iloc[ix]
                quality[block,endpoint]=[float(np.log1p(g.full_RNA_total).mean()),float(np.log1p(g.full_detected_genes).mean())]
        lookup=selected.set_index('analysis_block');fractions=[];counts=[];qualities=[]
        for block in paired.index:
            s=lookup.loc[block];uid=s.sample_uid
            fractions.append([s[labels[0]]/major.loc[uid,'T cell'],s[labels[1]]/major.loc[uid,'Macrophage']])
            counts.append([np.log1p(s[labels[0]]),np.log1p(s[labels[1]])])
            qualities.append(quality[block,types[0]]+quality[block,types[1]])
        fractions=np.array(fractions);counts=np.array(counts);qualities=np.array(qualities)
        variants=[np.empty((len(paired),0)),fractions,counts,qualities,np.column_stack((fractions,counts,qualities))]
        x=paired[types[0]].to_numpy();y=paired[types[1]].to_numpy()
        for name,cov in zip(plan['fixed_descriptive_variants'],variants):
            effect=partial(x,y,cov);values=[]
            for i in range(len(paired)):
                mask=np.arange(len(paired))!=i;e=partial(x[mask],y[mask],cov[mask])
                leave.append(dict(source=source,variant=name,left_out_ordinal=i+1,rho=e))
                if np.isfinite(e):values.append(e)
            rows.append(dict(source=source,variant=name,status='descriptive_estimated' if np.isfinite(effect) else 'design_or_df_insufficient',
                             n_patients=len(paired),covariates=cov.shape[1],rho=effect,LOO_valid=len(values),
                             LOO_min=min(values) if values else np.nan,LOO_max=max(values) if values else np.nan))
    pd.DataFrame(rows).to_csv(out/'all20_private_descriptive_source_diagnostics.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(leave).to_csv(out/'private_leave_one_patient_diagnostics.tsv',sep='\t',index=False,na_rep='NA')
    (out/'private_terminal.json').write_text(json.dumps(dict(status='PARTIAL_SOURCE_DESCRIPTIVE_DIAGNOSTICS_COMPLETE_NO_NEW_TESTS',
      finished_UTC=datetime.now(timezone.utc).isoformat()),indent=2)+'\n')
    print('PARTIAL_SOURCE_DESCRIPTIVE_DIAGNOSTICS_COMPLETE_RESULTS_REMAIN_SERVER')


if __name__=='__main__':main()
