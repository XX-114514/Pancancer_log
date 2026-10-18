"""Independent original metadata joins and all product endpoint records."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2)
    producer=root/'CAR_product_RNA_MCP_clinical_metadata_qualification_20261019_v2';assert sha(producer/'private_terminal.json')=='00a13975d441239e4844b17359124e96b73d127e2f04d747fbdfc6bf1929cdd9'
    term=json.loads((producer/'private_terminal.json').read_text())
    for name,v in term['files'].items():assert sha(producer/name)==v['sha256']
    target='CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea';collection=root.parents[1]/'runs/20261010_clinical_scrna_expansion_v1'
    conversion=json.loads((collection/'manifests'/f'{target}_native_conversion.json').read_text());accepted=json.loads((collection/'manifests'/f'accepted_{target}.json').read_text());cm=collection/'metadata/collections/d9a9760c-e43d-404d-86aa-447c9a1e9120.json';assert sha(cm)==accepted['official_source_metadata_sha256'];assert json.loads(cm.read_text())['doi']=='10.1016/j.ccell.2023.08.015'
    fields=['sample_id','donor_id','response3m','3mo PET/CT','cell_type']
    with h5py.File(conversion['source_path'],'r') as h:
        obs=h['obs'];key=str(obs.attrs['_index']);source=pd.DataFrame({f:np.asarray(read_elem(obs[f])).astype(str) for f in fields},index=pd.Index(np.asarray(read_elem(obs[key])).astype(str),name='original_cell_id'))
    qc=zarr.open_group(accepted['output'],mode='r');countsource=str(read_elem(qc['uns/raw_counts_source']))
    with h5py.File(countsource,'r') as h:
        obs=h['obs'];assert source.index.equals(pd.Index(np.asarray(read_elem(obs[str(obs.attrs['_index'])])).astype(str),name='original_cell_id'))
        for f in fields:assert np.array_equal(source[f],np.asarray(read_elem(obs[f])).astype(str))
    z=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r');get=lambda f:np.asarray(read_elem(z['obs'][f])).astype(str);take=get('dataset')==target
    atlas=pd.DataFrame({f:get(f)[take] for f in ['original_cell_id','library_id','patient_block','source_response','source_author_label']}).set_index('original_cell_id');assert atlas.index.is_unique and source.index.is_unique
    merged=atlas.join(source,how='left',validate='one_to_one');assert not merged[fields].isna().any().any() and len(merged)==396593
    # Original author labels are a distinct metadata provenance check, never
    # a formal cell-state call from the learned programme.
    missing={'','nan','None','<NA>','NA','__SOURCE_VALUE_MISSING__','unknown','Unknown'}
    normalize=lambda a:a.astype(str).where(~a.astype(str).isin(missing),'__DECLARED_MISSING__')
    assert normalize(merged.source_author_label).equals(normalize(merged.cell_type))
    raw_author_differences=int((merged.source_author_label!=merged.cell_type).sum())
    table=pd.read_csv(producer/'private_all59_exact_author_endpoint_and_patient_bridge.tsv',sep='\t',keep_default_na=False).set_index('library_id');assert len(table)==59
    for lib,g in merged.groupby('library_id'):
        r=table.loc[lib]
        for f in ['patient_block','donor_id','sample_id','response3m','3mo PET/CT','source_response']:assert g[f].nunique(dropna=False)==1
        original=g.iloc[0];binary=original.response3m.strip().upper();known=binary in ['R','NR']
        assert r.patient_block==original.patient_block and r.original_donor_id==original.donor_id and r.original_sample_id==original.sample_id and r.author_response3m==original.response3m and r.author_PET_CT_3month==original['3mo PET/CT']
        assert r.author_binary_response==(binary if known else 'UNREVIEWED_OR_UNKNOWN') and r.clinical_binary_eligible==known and r.retained_atlas_cells==len(g) and r.exact_original_cell_bridge and r.not_ORR_recode
        if known:assert float(r.binary_value)==float(binary=='R') and original.source_response.strip().upper()==binary
        else:assert r.binary_value=='' and original.source_response.strip().upper() not in ['R','NR']
        assert r.clinical_endpoint=='author_response3m_binary_3month_PET_CT'
    assert table.original_donor_id.is_unique and table.patient_block.is_unique
    valid=table[table.clinical_binary_eligible];assert len(valid)==35 and valid.author_binary_response.eq('R').sum()==15 and valid.author_binary_response.eq('NR').sum()==20
    joint=valid.groupby(['author_binary_response','author_PET_CT_3month']).size().sort_index();saved=pd.read_csv(producer/'private_source_author_binary_and_PET_joint_counts.tsv',sep='\t').set_index(['author_binary_response','author_PET_CT_3month']).products.sort_index();pd.testing.assert_series_equal(joint,saved,check_names=False)
    assert set(valid[valid.author_binary_response.eq('R')].author_PET_CT_3month.str.upper())=={'CR'}
    assert set(valid[valid.author_binary_response.eq('NR')].author_PET_CT_3month.str.upper())=={'PR','PD'}
    profile=merged.groupby(['library_id','source_author_label']).size().sort_index();saved=pd.read_csv(producer/'private_all59_original_author_label_profiles.tsv',sep='\t',keep_default_na=False).set_index(['library_id','source_author_label']).cells.sort_index();pd.testing.assert_series_equal(profile,saved,check_names=False)
    summary=json.loads((producer/'private_source_endpoint_qualification_summary.json').read_text());assert summary['exact_atlas_cell_bridges']==len(merged) and summary['binary_endpoint_products']==35 and summary['unreviewed_products']==24 and not summary['scores_or_new_association_read']
    out=root/'CAR_product_RNA_MCP_clinical_metadata_audit_20261019_v3';assert not out.exists();out.mkdir(mode=0o700)
    receipt=dict(status='ACCEPTED_ALL_ORIGINAL_METADATA_IDENTITIES_EXACT_CELL_BRIDGES_AND_AUTHOR3MONTH_BINARY_ENDPOINT_NO_SCORES',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(producer/'private_terminal.json'),resource_guard=guard,original_metadata_identities=417167,original_metadata_fields=5,exact_atlas_cell_bridges=len(merged),products=59,author_binary_endpoint_products=35,missing_endpoint_products=24,clinical_endpoint='author_response3m_binary_3month_PET_CT_observed_CR_vs_nonCR_not_ORR',author_labels_exact_original=raw_author_differences==0,raw_author_label_different_cells=raw_author_differences,author_labels_match_after_declared_missing_normalization=True,GEO_accession_material_bridge_independently_accepted=False,programme_scores_or_effects_read=False,authority_Zarr_modified=False)
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
