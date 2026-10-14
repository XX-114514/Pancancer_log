"""Frozen programme variation within all macrophage subtype candidates.
Post-result descriptive follow-up, no new P/BH or programme learning.
"""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import h5py,numpy as np,pandas as pd,zarr
from anndata._io.specs import read_elem
from scipy.stats import spearmanr
from diagnose_HNSCC_clinical_composition_v1 import beta


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    inputs=root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1'
    crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1'
    out=root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1';assert not out.exists();out.mkdir()
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      scope='all observed Macro_ multigene subtype candidates, original baseline clinical population',
      no_outcome_based_subtype_selection=True,min_cells_each_patient_subtype=20,max_cells_each_patient_subtype=50,
      random_seed=20261014,min_residual_df=10,
      fixed_variants=['clinical_arm_HPV','clinical_arm_HPV_age','clinical_arm_HPV_native_quality_count'],
      no_new_P_or_BH=True,not_formal_state_annotation=True,no_gene_weights_refitted=True,
      fixed_measurement='previous >=95% measured absolute weight; raw source CP10K with CRC quality coefficients')
    (out/'plan_before_diagnostics.json').write_text(json.dumps(plan,indent=2)+'\n')
    baseline=pd.read_csv(inputs/'private_eligible_libraries_clinical.tsv',sep='\t')
    baseline=baseline[baseline.timepoint.eq('pre')].set_index('library_id');assert baseline.patient_id.is_unique
    contract=json.loads((atlas/'config/source_contract_v2_private.json').read_text())
    rec=next(r for r in contract['records'] if r['dataset']=='GSE286827')
    g=h5py.File(rec['path'],'r') if str(rec['path']).endswith('.h5ad') else zarr.open_group(rec['path'],mode='r')
    native=read_elem(g['obs']);source_rows=np.load(rec['row_positions_ref'])
    names=pd.Index('GSE286827::'+native.index[source_rows].astype(str))
    chunks=[]
    for chunk in pd.read_csv(atlas/'tables/all_core_annotation_private.tsv.gz',sep='\t',index_col=0,dtype=str,keep_default_na=False,chunksize=200000):
        take=chunk[chunk.dataset.eq('GSE286827')&chunk.library_id.isin(baseline.index)&chunk.major.eq('Myeloid')]
        if len(take):chunks.append(take[['library_id','subcluster']])
    cells=pd.concat(chunks)
    provisional=pd.read_csv(atlas/'tables/subcluster_catalog_private.tsv',sep='\t')
    refined=pd.read_csv(atlas/'tables/subcluster_catalog_refined_private.tsv',sep='\t')
    bridge=provisional[['subcluster','node_id']].merge(refined[['node_id','subtype_candidate','subtype_evidence_status']],on='node_id',validate='one_to_one').set_index('subcluster')
    cells['subtype']=cells.subcluster.map(bridge.subtype_candidate)
    cells['evidence']=cells.subcluster.map(bridge.subtype_evidence_status)
    cells=cells[cells.subtype.str.startswith('Macro_')&cells.evidence.eq('multigene_subtype_candidate')]
    assert cells.index.isin(names).all()
    w=pd.read_csv(root/'HNSCC_frozen_MCP_gene_axis_20261014_v2/Macnode18_private_frozen_gene_columns.tsv',sep='\t')
    valid=w.native_column.ge(0).to_numpy()
    assert np.abs(w.loc[valid,'weight']).sum()/np.abs(w.weight).sum()>=.95
    qc=pd.read_csv(crc/'private_transfer_fixed_weights_v1/Macnode18_MCP1_quality.tsv',sep='\t').set_index('term').value
    quality=np.load(crc/'private_inputs/Macnode18_frozen_projection.npz')
    matrix=ad.io.sparse_dataset(g[rec['counts_key']])
    author=pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',sep='\t')
    age=author[author.time.eq('Pre')].set_index('patient_registry').Age
    all_means=[];all_diagnostics=[];selected_cells=[];coverage=[]
    for subtype,group in cells.groupby('subtype',sort=True):
        rng=np.random.default_rng(20261014+int(hashlib.sha256(subtype.encode()).hexdigest()[:8],16))
        indices=[];counts={}
        for library,sub in group.groupby('library_id',sort=True):
            counts[library]=len(sub)
            if len(sub)<20:continue
            positions=names.get_indexer(sub.index);indices.extend(rng.choice(positions,min(50,len(sub)),replace=False))
        indices=np.sort(np.asarray(indices,int))
        coverage.append(dict(subtype_candidate=subtype,baseline_patients_available=len(baseline),patients_with_any_cells=len(counts),patients_at_least20_cells=len({cells.loc[names[i],'library_id'] for i in indices})))
        if not len(indices):continue
        chosen=cells.loc[names[indices]]
        raw=matrix[source_rows[indices],:].tocsr();total=np.asarray(raw.sum(1)).ravel();detected=np.asarray((raw>0).sum(1)).ravel()
        assert (total>0).all() and np.isfinite(raw.data).all() and (raw.data>=0).all()
        values=np.log1p(raw[:,w.loc[valid,'native_column'].to_numpy(int)].toarray()*10000/total[:,None])
        z=(values-w.loc[valid,'mean'].to_numpy())/w.loc[valid,'sd'].to_numpy()
        Q=(np.column_stack((np.log1p(total),np.log1p(detected)))-quality['quality_mean'])/quality['quality_sd']
        score=z@w.loc[valid,'weight'].to_numpy()-qc['(Intercept)']-Q[:,0]*qc.cellQ-Q[:,1]*qc.geneQ
        assert np.isfinite(score).all()
        # Independent scalar check for fixed first and last cells.
        for i in sorted(set([0,len(score)-1])):
            direct=sum(float(z[i,j])*float(w.loc[valid,'weight'].iloc[j]) for j in range(valid.sum()))-qc['(Intercept)']-Q[i,0]*qc.cellQ-Q[i,1]*qc.geneQ
            assert abs(direct-score[i])<1e-7
        frame=pd.DataFrame(dict(library_id=chosen.library_id.to_numpy(),score=score,mean_logRNA=np.log1p(total),mean_log_genes=np.log1p(detected)))
        means=frame.groupby('library_id').mean();md=baseline.loc[means.index]
        means['patient_id']=md.patient_id;means['subtype_candidate']=subtype
        means['native_subtype_cells']=pd.Series(counts).loc[means.index]
        all_means.extend(means.reset_index().to_dict('records'))
        selected_cells.extend(dict(cell=names[i],subtype_candidate=subtype) for i in indices)
        clinical=np.column_stack((np.ones(len(means)),md.treatment_x.eq('D+T').to_numpy(float),md.HPV.eq('positive').to_numpy(float)))
        x=means.score.to_numpy();y=md.regression.to_numpy(float)
        variants=[np.empty((len(means),0)),age.loc[md.patient_id].to_numpy(float)[:,None],
                  np.column_stack((means.mean_logRNA,means.mean_log_genes,np.log1p(means.native_subtype_cells)))]
        for variant,cov in zip(plan['fixed_variants'],variants):
            point=beta(x,y,clinical,cov);loo=[]
            for i in range(len(x)):
                mask=np.arange(len(x))!=i;value=beta(x[mask],y[mask],clinical[mask],cov[mask])
                if np.isfinite(value):loo.append(value)
            all_diagnostics.append(dict(subtype_candidate=subtype,variant=variant,n_patients=len(means),
              beta=point,status='descriptive_estimated' if np.isfinite(point) else 'rank_or_residual_df_insufficient',
              patient_score_SD=float(x.std()),unadjusted_rho=float(spearmanr(x,y).statistic) if min(x.std(),y.std())>0 else np.nan,
              LOO_valid=len(loo),LOO_min=min(loo) if loo else np.nan,LOO_max=max(loo) if loo else np.nan))
    pd.DataFrame(all_means).to_csv(out/'private_patient_subtype_scores.tsv',sep='\t',index=False)
    pd.DataFrame(selected_cells).to_csv(out/'private_selected_cells.tsv',sep='\t',index=False)
    pd.DataFrame(coverage).to_csv(out/'aggregate_subtype_capture.tsv',sep='\t',index=False)
    pd.DataFrame(all_diagnostics).to_csv(out/'all_private_subtype_descriptive_diagnostics.tsv',sep='\t',index=False,na_rep='NA')
    report=dict(status='ALL_MACRO_SUBTYPE_FIXED_PROGRAMME_DIAGNOSTICS_COMPLETE_NO_NEW_P',
      finished_UTC=datetime.now(timezone.utc).isoformat(),coverage=coverage,diagnostics=all_diagnostics,
      gene_weights_changed=False,original_labels_changed=False,formal_state_verified=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
