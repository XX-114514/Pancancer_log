"""Count existing new-atlas HNSCC identities, exact GEO titles and sampling.
No score/response imputation; subtype candidates are not frozen172 identities.
"""
import json,os,re
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd


def main():
    project=Path(os.environ['PANCANCER_PROJECT_ROOT']);root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    out=root/'HNSCC_identity_capture_review_20261014_v2';assert not out.exists();out.mkdir()
    path=atlas/'tables/all_core_annotation_private.tsv.gz'
    cols=['dataset','library_id','sampleID','patient_id','timepoint','treatment','major','subcluster']
    header=pd.read_csv(path,sep='\t',nrows=0).columns
    assert set(cols).issubset(header), 'Expected completed annotation contract'
    chunks=[]
    for chunk in pd.read_csv(path,sep='\t',usecols=cols,chunksize=200000,keep_default_na=False):
        take=chunk[chunk.dataset.eq('GSE286827')]
        if len(take):chunks.append(take)
    cells=pd.concat(chunks,ignore_index=True);assert len(cells)==375019
    catalog=pd.read_csv(atlas/'tables/subcluster_catalog_refined_private.tsv',sep='\t',keep_default_na=False)
    assert catalog.subcluster.is_unique
    provisional=pd.read_csv(atlas/'tables/subcluster_catalog_private.tsv',sep='\t',keep_default_na=False)
    assert provisional.subcluster.is_unique and provisional.node_id.is_unique and catalog.node_id.is_unique
    assert set(provisional.node_id)==set(catalog.node_id)
    cells=cells.merge(provisional[['subcluster','node_id']],on='subcluster',how='left',validate='many_to_one')
    assert cells.node_id.notna().all()
    cells=cells.merge(catalog[['node_id','subtype_candidate','subtype_evidence_status']],on='node_id',how='left',validate='many_to_one')
    assert cells.subtype_evidence_status.notna().all()
    units=[]
    for library,g in cells.groupby('library_id',sort=True):
        for c in ['patient_id','timepoint','treatment']:assert g[c].nunique()==1
        md=g.iloc[0]
        cd8=g.major.eq('T')&g.subtype_candidate.str.startswith('CD8')&g.subtype_evidence_status.eq('multigene_subtype_candidate')
        macro=g.major.eq('Myeloid')&g.subtype_candidate.str.startswith('Macro_')&g.subtype_evidence_status.eq('multigene_subtype_candidate')
        units.append(dict(library_id=library,patient_id=md.patient_id,timepoint=md.timepoint,treatment=md.treatment,
          T_cells=int(g.major.eq('T').sum()),Myeloid_cells=int(g.major.eq('Myeloid').sum()),
          CD8_candidate_cells=int(cd8.sum()),Macrophage_candidate_cells=int(macro.sum())))
    units=pd.DataFrame(units);assert len(units)==57
    collection=project/'Project_v5/runs/20261010_clinical_scrna_expansion_v1'
    geo=json.loads((collection/'metadata/GSE286827/official_geo.json').read_text())
    by_gsm={g['gsm_id']:g for g in geo['samples']};assert len(by_gsm)==len(geo['samples'])
    checks=[]
    for row in units.itertuples():
        assert row.library_id in by_gsm
        title=by_gsm[row.library_id]['sample_title']
        match=re.fullmatch(r'(P\d+)_(pre|post)(\d*)_GEX',title);assert match
        checks.append(dict(library_id=row.library_id,patient_title_token=match[1],title_timepoint=match[2]+match[3],
                           patient_registry=row.patient_id,timepoint_registry=row.timepoint))
        assert row.timepoint==match[2]+match[3]
    exact=pd.DataFrame(checks)
    assert exact.groupby('patient_title_token').patient_registry.nunique().max()==1
    assert exact.groupby('patient_registry').patient_title_token.nunique().max()==1
    units.to_csv(out/'all57_private_cell_capture.tsv',sep='\t',index=False)
    exact.to_csv(out/'private_GEO_patient_title_crosscheck.tsv',sep='\t',index=False)
    rows=[]
    for endpoint,s,t in [('broad_T_Myeloid','T_cells','Myeloid_cells'),('multigene_candidate_CD8_Macrophage','CD8_candidate_cells','Macrophage_candidate_cells')]:
        valid=units[(units[s]>=20)&(units[t]>=20)]
        for phase in ['pre','post','post2']:
            q=valid[valid.timepoint.eq(phase)]
            rows.append(dict(endpoint_scope=endpoint,phase=phase,min_cells=20,libraries=len(q),source_patient_tokens=q.patient_id.nunique(),
                             treatment_counts=q.treatment.value_counts().to_dict(),old172_identity_mapping_accepted=False))
        pairs=set(valid.loc[valid.timepoint.eq('pre'),'patient_id'])&set(valid.loc[valid.timepoint.eq('post'),'patient_id'])
        rows.append(dict(endpoint_scope=endpoint,phase='paired_pre_primary_post',min_cells=20,source_patient_tokens=len(pairs),
                         old172_identity_mapping_accepted=False))
    pd.DataFrame(rows).to_csv(out/'aggregate_HNSCC_identity_capture_opportunities.tsv',sep='\t',index=False)
    (out/'private_terminal.json').write_text(json.dumps(dict(status='HNSCC_CAPTURE_AND_ALL57_GEO_PATIENT_TIME_KEYS_AUDITED',
      observed_UTC=datetime.now(timezone.utc).isoformat(),direct_old172_identity_mapping=False,response_labels_missing=True,
      clinical_or_treatment_MCP_tests_executed=False,clinical_value_verified=False),indent=2)+'\n')
    print('HNSCC_CAPTURE_AND_GEO_PATIENT_TIME_KEYS_AUDITED_NO_NEW_TESTS')


if __name__=='__main__':main()
