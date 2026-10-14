"""Existing spatial panels and clinical input audit; no expression scoring.
Private rows/gene identities stay on server. Output only execution status.
"""
import hashlib,json,os,re
from pathlib import Path
from datetime import datetime,timezone
import h5py,numpy as np,pandas as pd,zarr
from anndata._io.specs import read_elem


def frame(path,key):
    if str(path).endswith('.h5ad'):
        with h5py.File(path) as f:return read_elem(f[key])
    return read_elem(zarr.open_group(str(path),mode='r')[key])


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1';out=root/'frozen_MCP_validation_input_review_20261014_v1'
    assert not out.exists();out.mkdir()
    (out/'plan_before_results.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),
      spatial='all existing repaired GSE287472/GSE299786 spatial assay panels, exact unambiguous aliases',
      clinical='all9 immune sources captured registry; GSE286827 official fields review',
      no_expression_scoring=True,no_new_clinical_P_values=True,results_private=True),indent=2)+'\n')
    manifest_path=root/'immune_same_name_family_source_review_v1/spatial_native_input_repair_v2/analysis_input_manifest.tsv'
    manifest=pd.read_csv(manifest_path,sep='\t')
    spatial=[]
    for row in manifest[manifest.dataset_id.isin(['GSE287472','GSE299786'])].itertuples():
        var=frame(row.analysis_path,'var');lookup={}
        for i,(name,r) in enumerate(var.iterrows()):
            names={str(name)}
            if 'gene_symbol' in var:names.add(str(r.gene_symbol))
            for name in names:
                if name not in {'','nan','None'}:lookup.setdefault(name,set()).add(i)
        for endpoint in ['Tnode9','Macnode18']:
            w=pd.read_csv(crc/'private_transfer_fixed_weights_v1'/(endpoint+'_MCP1_genes.tsv'),sep='\t')
            measured=np.array([len(lookup.get(g,set()))==1 for g in w.gene]);fraction=float(np.abs(w.weight[measured]).sum()/np.abs(w.weight).sum())
            spatial.append(dict(dataset=row.dataset_id,assay_private=row.sample_id,endpoint=endpoint,
              panel_genes=len(var),frozen_genes=len(w),measured_genes=int(measured.sum()),absolute_weight_measurable=fraction,
              ambiguous_genes=sum(len(lookup.get(g,set()))>1 for g in w.gene),
              assay_RNA_denominator='targeted_panel_not_fullRNA',
              original_CP10K_quality_coefficients_compatible=False,
              patient_map_reviewed=str(row.independent_patient_mapping_reviewed)))
    pd.DataFrame(spatial).to_csv(out/'private_spatial_assay_weight_coverage.tsv',sep='\t',index=False)
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    contract=json.loads((atlas/'config/source_contract_v2_private.json').read_text())
    columns=['dataset','library_id','patient_id','patient_identity_verified','timepoint','treatment','tissue','capture_scope',
             'clinical_validation_eligible','native_response_endpoint','clinical_endpoint_definition','source_author_label']
    header=pd.read_csv(contract['all_rows_ref'],sep='\t',nrows=0).columns
    registry=pd.read_csv(contract['all_rows_ref'],sep='\t',usecols=[c for c in columns if c in header],keep_default_na=False)
    cohorts=[]
    production=set(contract['production_roles'])
    for source in ['GSE286827','GSE256326','GSE294109','GSE292046','GSE305078','GSE169246','GSE176021','GSE316195','CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea']:
        md=registry[registry.dataset.eq(source)];unique=md.drop_duplicates('library_id')
        per=[]
        for library,g in md.groupby('library_id',observed=True):
            for c in ['patient_id','timepoint','treatment','native_response_endpoint']:
                if c in g:assert g[c].nunique()==1,'Inconsistent per-library metadata: '+c
            per.append(g.iloc[0])
        units=pd.DataFrame(per)
        cohorts.append(dict(source=source,libraries=len(units),cells=len(md),
          nonmissing_patient_tokens=int(units.patient_id.ne('').sum()),
          unique_patient_tokens=int(units.loc[units.patient_id.ne(''),'patient_id'].nunique()),
          timepoint_counts=units.timepoint.value_counts().to_dict(),
          treatment_counts=units.treatment.value_counts().to_dict(),
          response_counts=units.native_response_endpoint.value_counts().to_dict(),
          clinical_eligible_libraries=int(units.clinical_validation_eligible.astype(str).str.lower().isin(['true','1']).sum())))
    (out/'private_all9_clinical_registry_scope.json').write_text(json.dumps(cohorts,indent=2)+'\n')
    collection=project/'Project_v5/runs/20261010_clinical_scrna_expansion_v1'
    geo=json.loads((collection/'metadata/GSE286827/official_geo.json').read_text());characteristics={}
    for sample in geo['samples']:
        raw=json.loads(sample['raw_characteristics_json'])
        if isinstance(raw,dict):
            for key in raw:characteristics[key]=characteristics.get(key,0)+1
    (out/'HNSCC_official_field_names.json').write_text(json.dumps(characteristics,indent=2)+'\n')
    rec=next(r for r in contract['records'] if r['dataset']=='GSE286827')
    obs=frame(rec['path'],'obs');selected=obs.iloc[np.load(rec['row_positions_ref'])]
    candidate=[c for c in selected if re.search('response|pathol|regress|treat|time|author|cell.*type|hpv',c,re.I)]
    fields=[]
    for column in candidate:
        value=selected[column].astype(str)
        fields.append(dict(field=column,nonmissing=int((~value.isin(['','nan','None','NA','__SOURCE_VALUE_MISSING__'])).sum()),distinct_values=value.nunique()))
    pd.DataFrame(fields).to_csv(out/'HNSCC_native_candidate_field_coverage.tsv',sep='\t',index=False)
    (out/'private_terminal.json').write_text(json.dumps(dict(status='EXISTING_SPATIAL_AND_CLINICAL_INPUT_REVIEW_COMPLETE_NO_NEW_MODELS',
      finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      no_new_expression_download=True,clinical_MCP_validation_complete=False),indent=2)+'\n')
    print('EXISTING_SPATIAL_AND_CLINICAL_INPUT_REVIEW_COMPLETE_RESULTS_REMAIN_SERVER')


if __name__=='__main__':main()
