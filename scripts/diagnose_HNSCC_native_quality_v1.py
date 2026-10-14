"""Native RNA and age diagnostics, descriptive only; unchanged sample draw."""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from diagnose_HNSCC_clinical_composition_v1 import beta


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    inputs=root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1'
    out=root/'HNSCC_native_quality_descriptive_20261014_v1';assert not out.exists();out.mkdir()
    (out/'plan_before_diagnostics.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),
      scope='all3 clinical score endpoints; same deterministic cell sample; baseline only',
      variants=['clinical_plus_four_native_RNA_quality_means','clinical_plus_baseline_age','clinical_plus_joint_RNA_quality_age'],
      no_new_p_values_or_BH=True,min_residual_df=10,no_new_cell_draw=True),indent=2)+'\n')
    md=pd.read_csv(inputs/'private_eligible_libraries_clinical.tsv',sep='\t')
    baseline=md[md.timepoint.eq('pre')].set_index('patient_id')
    scores=pd.read_csv(inputs/'private_patient_library_frozen_scores.tsv',sep='\t')
    s=scores[scores.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='score')
    z=scores[scores.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='CRC_training_z')
    s['equal_mean_CRC_training_z']=z.mean(axis=1);baseline=baseline.loc[s.index]
    contract=json.loads((atlas/'config/source_contract_v2_private.json').read_text())
    rec=next(r for r in contract['records'] if r['dataset']=='GSE286827')
    g=h5py.File(rec['path'],'r') if str(rec['path']).endswith('.h5ad') else zarr.open_group(rec['path'],mode='r')
    native=read_elem(g['obs']);source_rows=np.load(rec['row_positions_ref'])
    names=pd.Index('GSE286827::'+native.index[source_rows].astype(str))
    chunks=[]
    for chunk in pd.read_csv(atlas/'tables/all_core_annotation_private.tsv.gz',sep='\t',index_col=0,dtype=str,keep_default_na=False,chunksize=200000):
        take=chunk[chunk.dataset.eq('GSE286827')]
        if len(take):chunks.append(take[['library_id','major','subcluster']])
    cells=pd.concat(chunks).loc[names];assert len(cells)==len(names)
    provisional=pd.read_csv(atlas/'tables/subcluster_catalog_private.tsv',sep='\t')
    refined=pd.read_csv(atlas/'tables/subcluster_catalog_refined_private.tsv',sep='\t')
    bridge=provisional[['subcluster','node_id']].merge(refined[['node_id','subtype_candidate','subtype_evidence_status']],on='node_id',validate='one_to_one').set_index('subcluster')
    cells['subtype']=cells.subcluster.map(bridge.subtype_candidate)
    cells['evidence']=cells.subcluster.map(bridge.subtype_evidence_status)
    matrix=ad.io.sparse_dataset(g[rec['counts_key']]);quality=[]
    for endpoint,major,start,offset in [('Tnode9','T','CD8',0),('Macnode18','Myeloid','Macro_',1)]:
        candidate=cells[cells.library_id.isin(md.library_id)&cells.major.eq(major)&cells.subtype.str.startswith(start)&cells.evidence.eq('multigene_subtype_candidate')]
        rng=np.random.default_rng(20261014+offset);local=[]
        for library,group in candidate.groupby('library_id',sort=True):
            ix=names.get_indexer(group.index);chosen=rng.choice(ix,min(50,len(ix)),replace=False)
            if library in set(baseline.library_id):local.extend(chosen)
        local=np.sort(np.asarray(local,int));cm=cells.iloc[local]
        raw=matrix[source_rows[local],:].tocsr();total=np.asarray(raw.sum(1)).ravel();detected=np.asarray((raw>0).sum(1)).ravel()
        assert (total>0).all() and (detected>0).all() and (raw.data>=0).all()
        frame=pd.DataFrame(dict(library_id=cm.library_id.to_numpy(),RNA=np.log1p(total),genes=np.log1p(detected)))
        q=frame.groupby('library_id').mean();q.columns=[endpoint+'_'+col for col in q.columns]
        quality.append(q)
    q=pd.concat(quality,axis=1).loc[baseline.library_id].to_numpy();assert np.isfinite(q).all()
    author=pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',sep='\t')
    pre=author[author.time.eq('Pre')].set_index('patient_registry');assert pre.index.is_unique
    age=pd.to_numeric(pre.loc[baseline.index,'Age'],errors='raise').to_numpy(float)[:,None]
    assert np.isfinite(age).all()
    clinical=np.column_stack((np.ones(len(baseline)),baseline.treatment_x.eq('D+T').to_numpy(float),baseline.HPV.eq('positive').to_numpy(float)))
    y=baseline.regression.to_numpy(float);rows=[]
    for endpoint in s.columns:
        x=s[endpoint].to_numpy()
        for name,cov in [('clinical_plus_four_native_RNA_quality_means',q),('clinical_plus_baseline_age',age),('clinical_plus_joint_RNA_quality_age',np.column_stack((q,age)))]:
            point=beta(x,y,clinical,cov);leave=[]
            for i in range(len(x)):
                mask=np.arange(len(x))!=i;v=beta(x[mask],y[mask],clinical[mask],cov[mask])
                if np.isfinite(v):leave.append(v)
            rows.append(dict(endpoint=endpoint,variant=name,n_patients=len(x),beta=point,
              status='descriptive_estimated' if np.isfinite(point) else 'rank_or_df_insufficient',
              LOO_valid=len(leave),LOO_min=min(leave) if leave else np.nan,LOO_max=max(leave) if leave else np.nan))
    result=pd.DataFrame(rows);assert len(result)==9
    result.to_csv(out/'all9_private_native_quality_diagnostics.tsv',sep='\t',index=False,na_rep='NA')
    (out/'private_terminal.json').write_text(json.dumps(dict(status='NATIVE_QUALITY_AND_AGE_DIAGNOSTICS_COMPLETE_NO_NEW_P',finished_UTC=datetime.now(timezone.utc).isoformat()),indent=2)+'\n')
    print(result.to_json(orient='records'))


if __name__=='__main__':main()
