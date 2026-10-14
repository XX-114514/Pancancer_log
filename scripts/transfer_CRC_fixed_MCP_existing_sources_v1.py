"""Frozen CRC MCP transfer to all four existing non-CRC native sources.
Every scientific output remains server-side. No new discovery fitting.
"""
import hashlib
import json
import os
from pathlib import Path
from datetime import datetime, timezone

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def seed(key):
    return 20261014+int(hashlib.sha256(key.encode()).hexdigest()[:8],16)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1'
    out=crc/'private_existing_source_transfer_v1'
    assert not out.exists();out.mkdir()
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),
              sources=['GSE148071','GSE161529','GSE182109','GSE183904'],programs=['MCP1','MCP2'],
              labels=['CD8_T_cell_NKG7_9','Macrophage_IFI6_18'],endpoint_types=['Tnode9','Macnode18'],
              min_cells_each_endpoint=20,max_cells_each_patient_endpoint=50,
              sample_selection='lexicographic_first_eligible_library_per_reviewed_source_patient',
              min_patient_blocks=10,bootstrap=2000,permutations=20000,
              primary='all_frozen_weighted_genes_measured; native fullRNA logCP10K; CRC frozen normalization/quality',
              missing_genes='unmeasurable primary; no zero imputation',
              covariate_sensitivity='rank residualization on both fine_major fractions and log cell counts; min20 and residualdf10',
              BH_family='4sources*2programs=8 including unsupported slots internally p1',
              no_new_gene_selection_or_training=True,previously_used_CM_discovery_sources=True,
              not_independent_project_validation=True,not_direct_LR_or_clinical_validation=True,
              all_outputs_private=True)
    (out/'plan_before_results.json').write_text(json.dumps(plan,indent=2)+'\n')
    weights=crc/'private_transfer_fixed_weights_v1'
    rows=[];coverage=[];patient_scores=[];audit=[]
    for source in plan['sources']:
        native=root/'state_singlecell_inputs_20261014_v1/private_native_sources'/(source+'.h5ad')
        a=ad.read_h5ad(native);obs=a.obs;labels=plan['labels'];types=plan['endpoint_types']
        take=obs[obs.subCluster.isin(labels)]
        counts=take.groupby(['analysis_block','sample_uid','subCluster'],observed=True).size().unstack(fill_value=0)
        for label in labels:
            if label not in counts:counts[label]=0
        eligible=counts[(counts[labels]>=20).all(axis=1)].reset_index()
        selected=eligible.sort_values(['analysis_block','sample_uid']).drop_duplicates('analysis_block')
        review=obs.groupby('sample_uid',observed=True).patient_identity_reviewed.all()
        selected=selected[selected.sample_uid.map(review)]
        coverage.append(dict(source=source,eligible_libraries=len(eligible),reviewed_selected_patients=len(selected)))
        selected.to_csv(out/(source+'_private_selected_libraries.tsv'),sep='\t',index=False)
        major=obs.groupby(['sample_uid','major'],observed=True).size()
        score_by_program={p:{} for p in plan['programs']}
        for typ,label,parent in zip(types,labels,['T cell','Macrophage']):
            mask=obs.sample_uid.isin(selected.sample_uid)&obs.subCluster.eq(label)
            candidates=obs[mask];rng=np.random.default_rng(seed(source+typ+'cells'))
            positions=[]
            for _,group in candidates.groupby('analysis_block',observed=True,sort=True):
                ix=a.obs_names.get_indexer(group.index)
                positions.extend(rng.choice(ix,min(50,len(ix)),replace=False))
            positions=np.sort(positions).astype(int);md=obs.iloc[positions].copy()
            raw=a.X[positions].tocsr();total=np.asarray(raw.sum(1)).ravel()
            assert np.array_equal(total,md.full_RNA_total.to_numpy())
            qs=np.load(crc/'private_inputs'/(typ+'_frozen_projection.npz'))
            quality=(np.column_stack((np.log1p(total),np.log1p(md.full_detected_genes.to_numpy())))-qs['quality_mean'])/qs['quality_sd']
            for program in plan['programs']:
                path=weights/(typ+'_'+program+'_genes.tsv')
                if not path.exists():continue
                w=pd.read_csv(path,sep='\t');indices=a.var_names.get_indexer(w.gene)
                coverage.append(dict(source=source,endpoint=typ,program=program,
                                     frozen_genes=len(w),measured_genes=int((indices>=0).sum()),
                                     measured_absolute_weight_fraction=float(np.abs(w.weight[indices>=0]).sum()/np.abs(w.weight).sum())))
                if (indices<0).any():continue
                values=raw[:,indices].toarray().astype(float)
                values=np.log1p(values*10000/total[:,None])
                standardized=(values-w['mean'].to_numpy())/w.sd.to_numpy()
                qc=pd.read_csv(weights/(typ+'_'+program+'_quality.tsv'),sep='\t').set_index('term').value
                scores=standardized@w.weight.to_numpy()-qc['(Intercept)']-quality[:,0]*qc.cellQ-quality[:,1]*qc.geneQ
                assert np.isfinite(scores).all()
                frame=pd.DataFrame({'analysis_block':md.analysis_block.astype(str).to_numpy(),'score':scores})
                aggregated=frame.groupby('analysis_block',sort=True).score.mean()
                assert set(aggregated.index)==set(selected.analysis_block)
                score_by_program[program][typ]=aggregated
                for block,value in aggregated.items():
                    uid=selected.set_index('analysis_block').loc[block,'sample_uid']
                    patient_scores.append(dict(source=source,program=program,endpoint=typ,analysis_block=block,score=value))
                # Direct scalar projection check for deterministic first/last cells.
                if len(scores):
                    for index in sorted(set([0,len(scores)-1])):
                        direct=sum(float(w.weight.iloc[j])*float(standardized[index,j]) for j in range(len(w)))
                        direct-=qc['(Intercept)']+quality[index,0]*qc.cellQ+quality[index,1]*qc.geneQ
                        assert abs(direct-scores[index])<1e-7
                    audit.append(dict(source=source,endpoint=typ,program=program,scalar_parity=True))
        for program in plan['programs']:
            row=dict(source=source,program=program,status='unsupported_or_incomplete_gene_coverage',n_patients=len(selected),
                     rho=np.nan,lower=np.nan,upper=np.nan,p=np.nan,abundance_count_partial_rho=np.nan)
            both=score_by_program[program]
            if set(both)==set(types):
                paired=pd.concat([both[types[0]],both[types[1]]],axis=1);paired.columns=types
                assert paired.notna().all().all();n=len(paired)
                row['status']='patient_coverage_insufficient'
                if n>=10:
                    x=paired[types[0]].to_numpy();y=paired[types[1]].to_numpy()
                    if x.std()>1e-12 and y.std()>1e-12:
                        rho=float(spearmanr(x,y).statistic);rx=rankdata(x);ry=rankdata(y)
                        rx=(rx-rx.mean())/np.linalg.norm(rx-rx.mean());ry=(ry-ry.mean())/np.linalg.norm(ry-ry.mean())
                        rng=np.random.default_rng(seed(source+program+'statistics'));extreme=0
                        for _ in range(20000):extreme+=abs(rx@ry[rng.permutation(n)])>=abs(rho)-1e-12
                        boots=[]
                        for _ in range(2000):
                            ix=rng.integers(0,n,n);r=spearmanr(x[ix],y[ix]).statistic
                            if np.isfinite(r):boots.append(float(r))
                        row.update(status='estimated_transfer_consistency_only',rho=rho,p=(1+extreme)/20001,undefined_bootstrap=2000-len(boots))
                        if len(boots)>=1800:row['lower'],row['upper']=map(float,np.quantile(boots,[.025,.975]))
                        lookup=selected.set_index('analysis_block');cov=[]
                        for block in paired.index:
                            entry=lookup.loc[block];uid=entry.sample_uid
                            cov.append([entry[labels[0]]/major.loc[uid,'T cell'],entry[labels[1]]/major.loc[uid,'Macrophage'],
                                        np.log1p(entry[labels[0]]),np.log1p(entry[labels[1]])])
                        Z=np.column_stack((np.ones(n),np.apply_along_axis(rankdata,0,np.array(cov))))
                        if n>=20 and n-Z.shape[1]>=10 and np.linalg.matrix_rank(Z)==Z.shape[1]:
                            sx=rankdata(x);sy=rankdata(y);ex=sx-Z@np.linalg.lstsq(Z,sx,rcond=None)[0];ey=sy-Z@np.linalg.lstsq(Z,sy,rcond=None)[0]
                            if ex.std()>1e-12 and ey.std()>1e-12:row['abundance_count_partial_rho']=float(np.corrcoef(ex,ey)[0,1])
                    else:row['status']='constant_frozen_scores'
            rows.append(row)
    results=pd.DataFrame(rows);assert len(results)==8
    p=results.p.fillna(1).to_numpy();order=np.argsort(p);q=np.empty(8)
    q[order]=np.minimum(1,np.minimum.accumulate((p[order]*8/np.arange(1,9))[::-1])[::-1])
    results['BH_eight_slot_family']=np.where(results.p.notna(),q,np.nan)
    results.to_csv(out/'all8_private_source_consistency.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(coverage).to_csv(out/'private_capture_and_gene_coverage.tsv',sep='\t',index=False)
    pd.DataFrame(patient_scores).to_csv(out/'private_patient_scores.tsv',sep='\t',index=False)
    (out/'private_projection_scalar_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    (out/'private_terminal.json').write_text(json.dumps(dict(status='FROZEN_SOURCE_TRANSFER_COMPLETE_RESULTS_PRIVATE',
        finished_UTC=datetime.now(timezone.utc).isoformat(),no_new_expression_download=True,
        not_independent_project_validation=True,not_verified_interaction=True),indent=2)+'\n')
    print('FROZEN_SOURCE_TRANSFER_COMPLETE_RESULTS_REMAIN_SERVER')


if __name__=='__main__':
    main()
