"""Exact original signed programme/quality attribution; descriptive only.

All three existing macrophage subtype candidates, all five fixed algebraic
components and all50 nonoverlapping functions are kept. No new P/BH or refit.
"""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import numpy as np,pandas as pd
from diagnose_frozen_MCP_two_end_functions_v1 import hnscc_inputs
from diagnose_HNSCC_clinical_composition_v1 import beta
from diagnose_HNSCC_full50_orthogonal_functions_v1 import rho

COMPONENTS=['positive_weight_component','negative_weight_component','raw_gene_total','transferred_quality_correction','reconstructed_score']


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    source=root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1'
    crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1'
    out=root/'HNSCC_frozen_receiver_component_attribution_20261016_v1';assert not out.exists();out.mkdir()
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='All previously fixed C1QC/INFL/SPP1 subtype cell draws, patients, five algebraic components, all50 orthogonal Hallmarks.',
        components=COMPONENTS,attribution='score = positive weighted standardized genes + negative weighted standardized genes - original frozen intercept/quality coefficients. No baseline restandardization.',
        interpretation='Signed weights define centred contributions, not necessarily positive-valued expression states. Correlated component variance fractions are not causal explained fractions.',
        clinical='All five components x all existing subtype patients x arm/HPV, arm/HPV+nativeRNA/count, arm/HPV+age/sex; all LOO, no new P/BH.',
        function='Reuse original full50 scores with both original programme signatures excluded from reference AND rank background; unadjusted and native RNA/gene/count partial correlations, all LOO.',
        min_residual_df=10,all_unestimable_cases_retained=True,clinical_family_not_replaced=True,
        prior_candidates_not_promoted_by_this_diagnostic=True,no_new_expression_download=True)
    (out/'plan_before_component_diagnostics.json').write_text(json.dumps(plan,indent=2)+'\n')
    g,source_rows,cells,var,eligible,measured=hnscc_inputs(root,project)
    selected=pd.read_csv(source/'private_selected_cells.tsv',sep='\t')
    assert selected.cell.is_unique
    local=cells.index.get_indexer(selected.cell);assert (local>=0).all()
    md=cells.iloc[local]
    raw=ad.io.sparse_dataset(g['layers/counts'])[source_rows[local],:].tocsr().astype(float);raw.eliminate_zeros()
    total=np.asarray(raw.sum(1)).ravel();detected=np.diff(raw.indptr);assert (total>0).all()
    w=pd.read_csv(root/'HNSCC_frozen_MCP_gene_axis_20261014_v2/Macnode18_private_frozen_gene_columns.tsv',sep='\t')
    valid=w.native_column.ge(0);columns=w.loc[valid,'native_column'].to_numpy(int)
    assert measured[columns].all() and np.abs(w.loc[valid,'weight']).sum()/np.abs(w.weight).sum()>=.95
    values=np.log1p(raw[:,columns].toarray()*10000/total[:,None])
    z=(values-w.loc[valid,'mean'].to_numpy())/w.loc[valid,'sd'].to_numpy()
    weights=w.loc[valid,'weight'].to_numpy();positive=z[:,weights>0]@weights[weights>0];negative=z[:,weights<0]@weights[weights<0]
    frozen=np.load(crc/'private_inputs/Macnode18_frozen_projection.npz')
    coefficients=pd.read_csv(crc/'private_transfer_fixed_weights_v1/Macnode18_MCP1_quality.tsv',sep='\t').set_index('term').value
    Q=(np.column_stack((np.log1p(total),np.log1p(detected)))-frozen['quality_mean'])/frozen['quality_sd']
    correction=-coefficients['(Intercept)']-Q[:,0]*coefficients.cellQ-Q[:,1]*coefficients.geneQ
    frame=pd.DataFrame(dict(library_id=md.library_id.to_numpy(),subtype_candidate=selected.subtype_candidate.to_numpy(),
        positive_weight_component=positive,negative_weight_component=negative,raw_gene_total=positive+negative,
        transferred_quality_correction=correction,reconstructed_score=positive+negative+correction))
    means=frame.groupby(['library_id','subtype_candidate']).mean()
    original=pd.read_csv(source/'private_patient_subtype_scores.tsv',sep='\t').set_index(['library_id','subtype_candidate']).sort_index()
    means=means.loc[original.index]
    assert np.allclose(means.reconstructed_score,original.score,atol=1e-9,rtol=1e-9)
    assert np.allclose(means.raw_gene_total+means.transferred_quality_correction,means.reconstructed_score,atol=1e-10)
    author=pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',sep='\t')
    author=author[author.time.eq('Pre')].set_index('patient_registry')
    functions=pd.read_csv(root/'HNSCC_full50_orthogonal_function_diagnostics_20261014_v1/private_patient_subtype_function_means.tsv',sep='\t',index_col=[0,1]).loc[original.index]
    references=[c for c in functions if c.startswith('HALLMARK_')];assert len(references)==50
    clinical_rows=[];functional_rows=[];variance_rows=[]
    for subtype in sorted(original.index.get_level_values('subtype_candidate').unique()):
        old=original.xs(subtype,level='subtype_candidate').sort_index();m=means.xs(subtype,level='subtype_candidate').loc[old.index]
        f=functions.xs(subtype,level='subtype_candidate').loc[old.index]
        cm=author.loc[old.patient_id];y=cm['Tumor.regression'].to_numpy(float)
        arm=cm.neo_type.eq('D+T').to_numpy(float);HPV=cm['HPV/p16'].eq('positive').to_numpy(float)
        C=np.column_stack((np.ones(len(old)),arm,HPV))
        quality=np.column_stack((old.mean_logRNA,old.mean_log_genes,np.log1p(old.native_subtype_cells)))
        age=pd.to_numeric(cm.Age).to_numpy(float)
        male=cm.Gender.str.strip().str.lower().map({'m':1,'male':1,'f':0,'female':0}).to_numpy(float)
        parts=m[['positive_weight_component','negative_weight_component','transferred_quality_correction']].to_numpy()
        covariance=np.cov(parts,rowvar=False,ddof=0)
        assert np.isclose(covariance.sum(),np.var(m.reconstructed_score),atol=1e-9,rtol=1e-9)
        for j,name in enumerate(['positive_weight_component','negative_weight_component','transferred_quality_correction']):
            for k,other in enumerate(['positive_weight_component','negative_weight_component','transferred_quality_correction']):
                variance_rows.append(dict(subtype_candidate=subtype,component=name,other_component=other,covariance=float(covariance[j,k]),total_score_variance=float(np.var(m.reconstructed_score))))
        for component in COMPONENTS:
            x=m[component].to_numpy()
            for variant,N,cov in [('arm_HPV',C,np.empty((len(old),0))),('arm_HPV_native_quality_count',C,quality),('arm_HPV_age_sex',np.column_stack((C,male)),age[:,None])]:
                effect=beta(x,y,N,cov);leave=[]
                for i in range(len(x)):
                    mask=np.arange(len(x))!=i;value=beta(x[mask],y[mask],N[mask],cov[mask])
                    if np.isfinite(value):leave.append(value)
                clinical_rows.append(dict(subtype_candidate=subtype,component=component,variant=variant,n_patients=len(x),beta=effect,
                    LOO_valid=len(leave),LOO_min=min(leave) if leave else np.nan,LOO_max=max(leave) if leave else np.nan,
                    status='descriptive_estimated' if np.isfinite(effect) else 'rank_df_or_variation_insufficient'))
            for reference in references:
                yy=f[reference].to_numpy()
                for variant,cov in [('unadjusted',np.empty((len(x),0))),('native_quality_count',quality)]:
                    effect=rho(x,yy,cov);leave=[]
                    for i in range(len(x)):
                        mask=np.arange(len(x))!=i;value=rho(x[mask],yy[mask],cov[mask])
                        if np.isfinite(value):leave.append(value)
                    functional_rows.append(dict(subtype_candidate=subtype,component=component,reference=reference,variant=variant,
                        n_patients=len(x),rho=effect,LOO_valid=len(leave),LOO_min=min(leave) if leave else np.nan,LOO_max=max(leave) if leave else np.nan,
                        status='descriptive_estimated' if np.isfinite(effect) else 'rank_df_or_variation_insufficient'))
    assert len(clinical_rows)==45 and len(functional_rows)==1500 and len(variance_rows)==27
    means.to_csv(out/'private_patient_signed_and_quality_components.tsv',sep='\t')
    pd.DataFrame(clinical_rows).to_csv(out/'all45_private_component_clinical_descriptors.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(functional_rows).to_csv(out/'all1500_private_component_function_descriptors.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(variance_rows).to_csv(out/'private_component_covariance_decomposition.tsv',sep='\t',index=False)
    report=dict(status='EXACT_FROZEN_SIGNED_QUALITY_ATTRIBUTION_ALL3_SUBTYPES_COMPLETE_NO_NEW_P',finished_UTC=datetime.now(timezone.utc).isoformat(),
        original_patient_score_replay_max_error=float(np.max(np.abs(means.reconstructed_score-original.score))),
        clinical_descriptors=45,function_descriptors=1500,clinical_P_BH_changed=False,
        component_table_sha256=hashlib.sha256((out/'private_patient_signed_and_quality_components.tsv').read_bytes()).hexdigest())
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    c=pd.DataFrame(clinical_rows);print(c[c.subtype_candidate.eq('Macro_C1QC')].to_json(orient='records'))
    f=pd.DataFrame(functional_rows);print(f[f.subtype_candidate.eq('Macro_C1QC')&f.variant.eq('native_quality_count')&f.reference.str.contains('INTERFERON')].to_json(orient='records'))


if __name__=='__main__':main()
