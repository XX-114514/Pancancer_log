"""Outcome metadata qualification only; exact original-cell bridge, no scores."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

TARGET='CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea'
DOI='10.1016/j.ccell.2023.08.015'


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    collection=root.parents[1]/'runs/20261010_clinical_scrna_expansion_v1';cm=collection/'metadata/collections/d9a9760c-e43d-404d-86aa-447c9a1e9120.json';public=json.loads(cm.read_text());assert public['doi']==DOI
    matches=[x for x in public['datasets'] if x['dataset_id']==TARGET.removeprefix('CXG_').replace('_','-')];assert len(matches)==1 and matches[0]['cell_count']==417167
    cvpath=collection/'manifests'/f'{TARGET}_native_conversion.json';conversion=json.loads(cvpath.read_text());rawpath=Path(conversion['source_path']);assert rawpath.exists()
    prepath=collection/'manifests'/f'{TARGET}_public_preflight.json';pre=json.loads(prepath.read_text());assert pre['source_doi']==DOI and pre['source_sha256']==conversion['source_sha256']
    acceptedpath=collection/'manifests'/f'accepted_{TARGET}.json';accepted=json.loads(acceptedpath.read_text());assert accepted['source_sha256']==conversion['source_sha256']
    out=root/'CAR_product_RNA_MCP_clinical_metadata_qualification_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,scope='All original public infusion product metadata and exact atlas original-cell IDs. No programme scores, expression matrix or outcome association accessed. Preserve author response3m R/NR verbatim; no inferred label for unknown. Separate author binary response from ORR.',source_DOI=DOI,official_collection_sha256=sha(cm),original_count_conversion_manifest_sha256=sha(cvpath),public_source_preflight_sha256=sha(prepath),accepted_original_source_manifest_sha256=sha(acceptedpath),original_H5AD_previously_accepted_sha256=conversion['source_sha256'],original_H5AD_live_size_bytes=rawpath.stat().st_size,reference_endpoint='Original article PMID37738975 states3monthPET/CT followup. Do not infer binary definition from abstract; verify exact published response3m column and joint PET coding.',no_new_expression_or_download_or_effects=True)
    (out/'plan_before_exact_original_metadata_bridge.json').write_text(json.dumps(plan,indent=2)+'\n')
    fields=['sample_id','donor_id','response3m','3mo PET/CT','cell_type']
    with h5py.File(rawpath,'r') as h:
        obs=h['obs'];index=np.asarray(read_elem(obs[str(obs.attrs.get('_index','_index'))])).astype(str)
        source=pd.DataFrame({f:np.asarray(read_elem(obs[f])).astype(str) for f in fields},index=pd.Index(index,name='original_cell_id'))
    assert len(source)==417167 and source.index.is_unique
    z=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r');get=lambda field:np.asarray(read_elem(z['obs'][field])).astype(str)
    mask=get('dataset')==TARGET;ids=get('original_cell_id')[mask];positions=source.index.get_indexer(ids);assert (positions>=0).all() and pd.Index(ids).is_unique
    current=pd.DataFrame({field:get(field)[mask] for field in ['library_id','patient_block','source_response','source_author_label']})
    current['original_cell_id']=ids
    for field in fields:current['original_'+field]=source[field].to_numpy(str)[positions]
    assert current.original_donor_id.nunique()==59 and current.library_id.nunique()==59 and current.patient_block.nunique()==59
    for key in ['original_donor_id','library_id','patient_block']:
        for field in ['original_donor_id','library_id','patient_block','original_sample_id','original_response3m','original_3mo PET/CT','source_response']:
            assert current.groupby(key)[field].nunique(dropna=False).eq(1).all()
    source_known=current.original_response3m.str.strip().str.upper().isin(['R','NR']);current_known=current.source_response.str.strip().str.upper().isin(['R','NR']);assert source_known.equals(current_known)
    assert current.loc[source_known,'source_response'].str.strip().str.upper().equals(current.loc[source_known,'original_response3m'].str.strip().str.upper())
    grouped=current.groupby('library_id',sort=True).first();rows=[]
    for library,r in grouped.iterrows():
        response=r.original_response3m.strip().upper();pet=r['original_3mo PET/CT'].strip().upper();known=response in ['R','NR']
        rows.append(dict(library_id=library,patient_block=r.patient_block,original_donor_id=r.original_donor_id,original_sample_id=r.original_sample_id,author_response3m=r.original_response3m,author_PET_CT_3month=r['original_3mo PET/CT'],author_binary_response=response if known else 'UNREVIEWED_OR_UNKNOWN',clinical_binary_eligible=known,binary_value=int(response=='R') if known else np.nan,clinical_endpoint='author_response3m_binary_3month_PET_CT',retained_atlas_cells=int(current.library_id.eq(library).sum()),exact_original_cell_bridge=True,not_ORR_recode=True))
    table=pd.DataFrame(rows);assert len(table)==59 and table.clinical_binary_eligible.sum()==35
    # Author binary status is CR versus non-CR among these labelled products.
    # In particular, PR is not silently recoded as a positive ORR endpoint.
    valid=table[table.clinical_binary_eligible]
    joint=valid.groupby(['author_binary_response','author_PET_CT_3month']).size().reset_index(name='products')
    assert set(valid.loc[valid.author_binary_response.eq('R'),'author_PET_CT_3month'].str.upper())=={'CR'}
    assert set(valid.loc[valid.author_binary_response.eq('NR'),'author_PET_CT_3month'].str.upper()).issubset({'PR','SD','PD'})
    profile=current.groupby(['library_id','source_author_label']).size().reset_index(name='cells')
    table.to_csv(out/'private_all59_exact_author_endpoint_and_patient_bridge.tsv',sep='\t',index=False);joint.to_csv(out/'private_source_author_binary_and_PET_joint_counts.tsv',sep='\t',index=False);profile.to_csv(out/'private_all59_original_author_label_profiles.tsv',sep='\t',index=False)
    # Bind the live metadata payload without re-reading any expression block.
    payload=hashlib.sha256()
    for field,values in [('index',source.index.to_numpy(str))]+[(f,source[f].to_numpy(str)) for f in fields]:
        payload.update(field.encode()+b'\0')
        for value in values:payload.update(value.encode()+b'\0')
    summary=dict(dataset=TARGET,study_DOI=DOI,related_GEO='GSE241783',study_alias_same_publication=True,exact_GEO_sample_accession_bridge_not_yet_independently_verified=True,source_products=59,retained_atlas_products=59,original_source_cells=417167,exact_atlas_cell_bridges=len(ids),binary_endpoint_products=35,author_R_products=int(valid.author_binary_response.eq('R').sum()),author_NR_products=int(valid.author_binary_response.eq('NR').sum()),unreviewed_products=int((~table.clinical_binary_eligible).sum()),author_binary_PET_joint=joint.to_dict(orient='records'),endpoint='Author response3m binary at3monthPET/CT. Observed R=CR; NR contains partial/progressive disease. This is not an ORR CR+PR recode.',product_context='Pre-infusion manufactured CAR-T product, not tumour ICI microenvironment. Source already present in project; not untouched project-independent validation.',live_original_obs_payload_sha256=payload.hexdigest(),scores_or_new_association_read=False,original_Zarr_or_identity_modified=False)
    (out/'private_source_endpoint_qualification_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<4*1024**2
    terminal=dict(status='COMPLETE_EXACT_AUTHOR3MONTH_BINARY_METADATA_BRIDGE_NO_SCORES_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,products=59,exact_cell_bridges=len(ids),clinical_products=35,individual_rows_or_scores_exported=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'),exact_cell_bridges=len(ids))))


if __name__=='__main__':main()
