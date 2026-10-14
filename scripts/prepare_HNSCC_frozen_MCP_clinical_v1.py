"""Freeze clinical definitions and project unchanged CRC weights privately.

The author binary discrepancies stay unresolved. Continuous pathology is
consistent across visits; the >=50% label is explicitly derived, not a repair.
New-atlas multigene subtype candidates are not the original 172 identities.
"""
import gzip
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import anndata as ad
import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.io import mmread
from scipy.stats import spearmanr


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1'
    table_root=root/'HNSCC_author_clinical_table_user_20261014_v1'
    table=table_root/'private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    axis=root/'HNSCC_frozen_MCP_gene_axis_20261014_v2'
    capture=root/'HNSCC_identity_capture_review_20261014_v2/all57_private_cell_capture.tsv'
    out=root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1'
    assert not out.exists();out.mkdir()
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),
      input_sha256={p.name:sha(p) for p in [table,capture,axis/'private_terminal.json']},
      endpoints=['Tnode9','Macnode18','equal_mean_CRC_training_z'],program_slots=['MCP1','MCP2'],
      primary_outcome='continuous author Tumor.regression, confirmed identical across patient visits',
      binary_secondary='explicit paper-defined derived regression>=50%; do not overwrite author binary',
      author_binary_sensitivity='only patients with identical binary across all visits and agreement with pathology threshold',
      primary_assay='baseline only; both multigene-candidate CD8 and macrophage>=20; primary pre/post paired cohort separate',
      original_172_identity_mapping=False,formal_state_labels=False,
      identity_rule='major T, subtype starts CD8, or major Myeloid, subtype starts Macro_; multigene_subtype_candidate',
      raw_counts='preserved native counts; complete measured source-axis RNA denominator and detected genes; not ambient-corrected X',
      cell_selection='fixed seed, up to50 cells per library and endpoint, independent of expression/outcome',
      score_measurement='>=95% absolute fixed weight; omit unmeasured centered contributions, no gene/weight/quality refit',
      calibration='same source-measurable frozen genes against full CRC scores in discovery and holdout; descriptive only',
      primary_model='patient baseline rank outcome ~ standardized rank score + treatment D+T + HPV/p16 positive',
      model_gate='>=20 patients, >=10 residual df, full design rank; no automatic covariate deletion',
      primary_tests=6,primary_BH='3 score endpoints * 2 original slots; unsupported slots internally p1',
      primary_permutation='20000 Freedman-Lane residual permutations under fixed treatment/HPV nuisance model',
      bootstrap='2000 patient blocks, invalid rank/constant resamples recorded; >=1800 valid for CI',
      secondary='unadjusted continuous Spearman and rank-biserial binary difference; separate 6-slot families',
      binary_min_group_patients=5,
      original_author_label_conflicts_resolved=False,post_scores_not_response_prediction=True,
      response_association_not_direct_interaction_or_causal_therapy_effect=True,
      random_seed=20261014,all_outputs_private=True)
    (out/'plan_before_scores_or_effects.json').write_text(json.dumps(plan,indent=2)+'\n')
    author=pd.read_csv(table,sep='\t',keep_default_na=False)
    patient=[]
    for key,group in author.groupby('patient_registry',sort=True):
        for col in ['patient.tag','Tumor.regression','neo_type','HPV/p16']:
            assert group[col].nunique()==1,'Clinical covariate differs across visits: '+col
        regression=float(group['Tumor.regression'].iloc[0]);assert 0<=regression<=100
        raw_binary=group['Partial.pathologic.response'].astype(float).astype(int)
        consistent=raw_binary.nunique()==1
        concordant=consistent and raw_binary.iloc[0]==int(regression>=50)
        patient.append(dict(patient_id=key,regression=regression,derived_response=int(regression>=50),
          treatment=group.neo_type.iloc[0],HPV=group['HPV/p16'].iloc[0],
          author_binary_consistent=consistent,author_binary_concordant=concordant,
          author_binary=int(raw_binary.iloc[0]) if concordant else np.nan))
    patient=pd.DataFrame(patient);assert patient.patient_id.is_unique
    patient.to_csv(out/'private_patient_clinical_definition.tsv',sep='\t',index=False)
    contract=json.loads((atlas/'config/source_contract_v2_private.json').read_text())
    rec=next(r for r in contract['records'] if r['dataset']=='GSE286827')
    g=h5py.File(rec['path'],'r') if str(rec['path']).endswith('.h5ad') else zarr.open_group(rec['path'],mode='r')
    native_obs=read_elem(g['obs']);source_rows=np.load(rec['row_positions_ref'])
    source_names=pd.Index('GSE286827::'+native_obs.index[source_rows].astype(str))
    assert source_names.is_unique
    columns=['dataset','library_id','patient_id','timepoint','treatment','major','subcluster']
    annotation=[]
    for chunk in pd.read_csv(atlas/'tables/all_core_annotation_private.tsv.gz',sep='\t',index_col=0,keep_default_na=False,chunksize=200000):
        take=chunk[chunk.dataset.eq('GSE286827')]
        if len(take):annotation.append(take[columns])
    cells=pd.concat(annotation);assert len(cells)==len(source_names) and set(cells.index)==set(source_names)
    cells=cells.loc[source_names]
    provisional=pd.read_csv(atlas/'tables/subcluster_catalog_private.tsv',sep='\t',keep_default_na=False)
    refined=pd.read_csv(atlas/'tables/subcluster_catalog_refined_private.tsv',sep='\t',keep_default_na=False)
    assert provisional.node_id.is_unique and refined.node_id.is_unique
    bridge=provisional[['subcluster','node_id']].merge(refined[['node_id','subtype_candidate','subtype_evidence_status']],on='node_id',validate='one_to_one').set_index('subcluster')
    cells['subtype_candidate']=cells.subcluster.map(bridge.subtype_candidate)
    cells['subtype_evidence_status']=cells.subcluster.map(bridge.subtype_evidence_status)
    assert cells.subtype_evidence_status.notna().all()
    unit=pd.read_csv(capture,sep='\t',keep_default_na=False)
    eligible=unit[(unit.CD8_candidate_cells>=20)&(unit.Macrophage_candidate_cells>=20)&unit.timepoint.isin(['pre','post'])]
    assert not eligible.duplicated(['patient_id','timepoint']).any()
    eligible=eligible.merge(patient,on='patient_id',validate='many_to_one')
    assert len(eligible)==len(unit[(unit.CD8_candidate_cells>=20)&(unit.Macrophage_candidate_cells>=20)&unit.timepoint.isin(['pre','post'])])
    eligible.to_csv(out/'private_eligible_libraries_clinical.tsv',sep='\t',index=False)
    matrix=ad.io.sparse_dataset(g[rec['counts_key']])
    order=list(read_elem(g['uns'])['source_order'])
    measured=np.asarray(g['varm/measured_in_source'])[:,order.index(rec['new_parent_dataset'])].astype(bool)
    scores=[];calibration=[];scalar=[]
    for endpoint,major,start in [('Tnode9','T','CD8'),('Macnode18','Myeloid','Macro_')]:
        candidates=cells[cells.library_id.isin(eligible.library_id)&cells.major.eq(major)&cells.subtype_candidate.str.startswith(start)&cells.subtype_evidence_status.eq('multigene_subtype_candidate')]
        rng=np.random.default_rng(20261014+(0 if endpoint=='Tnode9' else 1))
        local=[]
        for _,group in candidates.groupby('library_id',sort=True):
            pos=source_names.get_indexer(group.index);assert (pos>=0).all()
            local.extend(rng.choice(pos,min(50,len(pos)),replace=False))
        local=np.sort(np.asarray(local,dtype=int));md=cells.iloc[local]
        rows=source_rows[local];raw=matrix[rows,:].tocsr()
        assert np.isfinite(raw.data).all() and (raw.data>=0).all() and np.allclose(raw.data,np.rint(raw.data))
        assert raw[:,~measured].nnz==0
        total=np.asarray(raw.sum(1)).ravel();detected=np.asarray((raw>0).sum(1)).ravel()
        assert (total>0).all()
        w=pd.read_csv(axis/(endpoint+'_private_frozen_gene_columns.tsv'),sep='\t')
        valid=w.native_column.ge(0).to_numpy();fraction=float(np.abs(w.loc[valid,'weight']).sum()/np.abs(w.weight).sum())
        assert fraction>=.95
        qc=pd.read_csv(crc/'private_transfer_fixed_weights_v1'/(endpoint+'_MCP1_quality.tsv'),sep='\t').set_index('term').value
        quality=np.load(crc/'private_inputs'/(endpoint+'_frozen_projection.npz'))
        q=(np.column_stack((np.log1p(total),np.log1p(detected)))-quality['quality_mean'])/quality['quality_sd']
        values=np.log1p(raw[:,w.loc[valid,'native_column'].to_numpy(dtype=int)].toarray()*10000/total[:,None])
        z=(values-w.loc[valid,'mean'].to_numpy())/w.loc[valid,'sd'].to_numpy()
        score=z@w.loc[valid,'weight'].to_numpy()-qc['(Intercept)']-q[:,0]*qc.cellQ-q[:,1]*qc.geneQ
        assert np.isfinite(score).all()
        for ix in sorted(set([0,len(score)-1])):
            direct=sum(float(z[ix,j])*float(w.loc[valid,'weight'].iloc[j]) for j in range(valid.sum()))
            direct-=qc['(Intercept)']+q[ix,0]*qc.cellQ+q[ix,1]*qc.geneQ
            assert abs(direct-score[ix])<1e-7
        frame=pd.DataFrame(dict(library_id=md.library_id.to_numpy(),patient_id=md.patient_id.to_numpy(),timepoint=md.timepoint.to_numpy(),score=score))
        aggregated=frame.groupby(['library_id','patient_id','timepoint'],sort=True).score.mean().reset_index()
        assert set(aggregated.library_id)==set(eligible.library_id)
        # Original full CRC discovery patient scale, never source refitted.
        frozen_scale=None
        for partition in ['discovery','holdout']:
            prefix=crc/'private_inputs'/(endpoint+'_'+partition)
            genes=pd.read_csv(str(prefix)+'_genes.tsv',header=None)[0]
            ix=pd.Index(genes).get_indexer(w.gene);assert (ix>=0).all()
            with gzip.open(str(prefix)+'_tpm.mtx.gz','rb') as handle:x=mmread(handle).tocsr()[ix].toarray().T
            train_md=pd.read_csv(str(prefix)+'_metadata.tsv',sep='\t')
            standardized=(x-w['mean'].to_numpy())/w.sd.to_numpy()
            correction=qc['(Intercept)']+train_md.cellQ.to_numpy()*qc.cellQ+train_md.geneQ.to_numpy()*qc.geneQ
            full=standardized@w.weight.to_numpy()-correction
            partial=standardized[:,valid]@w.loc[valid,'weight'].to_numpy()-correction
            means=pd.DataFrame(dict(patient=train_md['sample'],full=full,partial=partial)).groupby('patient').mean()
            calibration.append(dict(endpoint=endpoint,partition=partition,patient_rho=float(spearmanr(means.full,means.partial).statistic),cell_rho=float(spearmanr(full,partial).statistic)))
            if partition=='discovery':frozen_scale=(float(means.full.mean()),float(means.full.std(ddof=0)))
        assert frozen_scale[1]>0
        aggregated['CRC_training_z']=(aggregated.score-frozen_scale[0])/frozen_scale[1]
        aggregated['endpoint']=endpoint;aggregated['program']='MCP1';scores.extend(aggregated.to_dict('records'))
        scalar.append(dict(endpoint=endpoint,n_cells=len(md),n_libraries=len(aggregated),scalar_parity=True,absolute_weight_retained=fraction))
    pd.DataFrame(scores).to_csv(out/'private_patient_library_frozen_scores.tsv',sep='\t',index=False)
    pd.DataFrame(calibration).to_csv(out/'aggregate_source_axis_CRC_calibration.tsv',sep='\t',index=False)
    baseline=eligible[eligible.timepoint.eq('pre')]
    terminal=dict(status='CLINICAL_DEFINITIONS_AND_FIXED_PARTIAL_SCORES_PREPARED_NO_EFFECT_TESTS',
      finished_UTC=datetime.now(timezone.utc).isoformat(),baseline_patients=baseline.patient_id.nunique(),
      baseline_derived_response_counts=baseline.derived_response.value_counts().to_dict(),
      baseline_HPV_categories=baseline.HPV.value_counts().to_dict(),
      baseline_author_binary_concordant_patients=int(baseline.author_binary_concordant.sum()),
      projection_audit=scalar,calibration=calibration,clinical_P_values_computed=False,
      original_labels_changed=False,author_binary_overwritten=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(terminal))


if __name__=='__main__':main()
