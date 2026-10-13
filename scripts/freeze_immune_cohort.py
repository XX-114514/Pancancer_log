#!/usr/bin/env python3
"""Freeze a native-count immune cohort after independent additional-QC audit.

Compiled, selectively retained exports stay supplemental. Row-level clinical
keys and matrix locations are private; public outputs contain aggregates only.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import anndata as ad
from anndata._io.specs import read_elem,write_elem
import numpy as np
import pandas as pd
import zarr
from merge_public_preprocessed_zarr import CSRAppender,csr_read,adopted_counts
from immune_atlas_io import save_json_new


def main(args):
    config=json.loads(args.contract.read_text())
    proof=json.loads(args.audit.read_text())
    if not proof['all_checks_pass'] or proof['n_datasets']!=1:
        raise ValueError('Additional source requires a complete independent QC audit')
    entry=proof['datasets'][0]
    if entry['dataset']!='GSE305078':raise ValueError('Unexpected additional study')
    directory=args.prepared/'GSE305078'
    summary=json.loads((directory/'dataset_summary.json').read_text())
    if summary['status']!='success' or summary.get('pilot_only'):raise ValueError('Additional QC not complete')
    destination=args.scratch/'GSE305078_prepared_native_union.zarr'
    if destination.exists():raise ValueError('Existing preparation attempt preserved')
    frames=[];order=[];var=None;checked={r['library_token']:r for r in entry['libraries']}
    for row in sorted(summary['libraries'],key=lambda r:r['library_token']):
        if row['status']!='success' or not checked[row['library_token']]['all_checks_pass']:
            raise ValueError('A library failed independent audit')
        p=directory/row['library_token']/'prepared.zarr'
        group=zarr.open_group(str(p),mode='r');local_var=read_elem(group['var'])
        if var is not None and not var.index.equals(local_var.index):raise ValueError('Native gene order drift')
        var=local_var
        obs=read_elem(group['obs']);obs.index=pd.Index(obs.original_cell_id.astype(str),name='original_cell_id_key')
        frames.append(obs);order.append((p,row))
    native_obs=pd.concat(frames)
    if not native_obs.index.is_unique:raise ValueError('Additional original cell IDs repeat')
    destination.mkdir(mode=0o700)
    g=zarr.open_group(str(destination),mode='w');g.attrs.update({'encoding-type':'anndata','encoding-version':'0.1.0'})
    layers=g.create_group('layers');layers.attrs.update({'encoding-type':'dict','encoding-version':'0.1.0'})
    raw_out=CSRAppender(layers,'counts','i4',len(var));removal_out=CSRAppender(layers,'ambient_removed_counts','f8',len(var))
    xout=CSRAppender(g,'X','f4',len(var));adopted=0
    for path,row in order:
        z=zarr.open_group(str(path),mode='r');obs=read_elem(z['obs'])
        for lo in range(0,len(obs),2048):
            hi=min(lo+2048,len(obs))
            raw=csr_read(z['layers/counts'],lo,hi);estimate=csr_read(z['layers/counts_decontx'],lo,hi)
            removed,x,total,selected=adopted_counts(raw,estimate,obs.iloc[lo:hi],row,'GSE305078')
            raw_out.append(raw);removal_out.append(removed);xout.append(x);adopted+=int(selected.sum())
    write_elem(g,'obs',native_obs);write_elem(g,'var',var)
    for key in ['obsm','obsp','varm','varp']:write_elem(g,key,{})
    write_elem(g,'uns',dict(native_counts=True,independent_QC_audit_sha256=hashlib.sha256(args.audit.read_bytes()).hexdigest()))
    zarr.consolidate_metadata(str(destination))
    obs=pd.read_csv(config['all_rows_ref'],sep='\t',dtype=str,keep_default_na=False,index_col=0)
    old305=obs.loc[obs.dataset.eq('GSE305078')].copy()
    old305=old305.loc['GSE305078::'+native_obs.index.astype(str)].copy()
    old305['source_row_position']=np.arange(len(native_obs))
    old305['prior_QC_status']='native_sample_QC_DD_ambient_assessment_independently_audited'
    old305['analysis_role']='core'
    old305['capture_scope']='whole_tumor_patient_mapping_pending'
    for key in ['residual_doublet_flag','ambient_high_loss_review','ambient_status','original_native_row_position']:
        old305[key]=native_obs[key].to_numpy()
    obs=pd.concat([obs.loc[~obs.dataset.eq('GSE305078')],old305])
    for record in config['records']:
        if record['role']=='core_source_review':
            record['role']='supplemental_retained_export'
            record['CM_composition_eligible']=False
            record['composition_constraint']='author_donor_celltype_timepoint_min20_filter_not_population_census'
            obs.loc[obs.dataset.eq(record['dataset']),'analysis_role']=record['role']
        if record['dataset']=='GSE305078':
            selector=args.run/'manifests/GSE305078_postQC_source_rows.npy';np.save(selector,np.arange(len(native_obs)))
            record.update(path=str(destination),counts_key='layers/counts',var_key='var',
                          row_positions_ref=str(selector),n_cells=len(native_obs),role='core',
                          ambient_removed_key='layers/ambient_removed_counts',
                          prior_QC_status='native_sample_QC_DD_ambient_assessment_independently_audited')
    sidecar=pd.read_csv(config['clinical_sidecar_ref'],sep='\t',dtype=str,keep_default_na=False)
    if sidecar.duplicated(['dataset','library_id']).any():raise ValueError('Clinical library mapping ambiguous')
    obs['patient_identity_verified']=False
    obs['clinical_validation_eligible']=False
    obs['CM_composition_eligible']=obs.analysis_role.isin(['core','CAR_product'])
    obs['ex_vivo_stimulation']=False
    obs['clinical_endpoint_definition']='missing_or_not_frozen'
    obs['native_response_endpoint']=''
    for name in ['GSE169246','GSE176021','GSE316195']:
        mask=obs.dataset.eq(name);frame=sidecar.loc[sidecar.dataset.eq(name)].set_index('library_id')
        libraries=obs.loc[mask,'library_id']
        if not libraries.isin(frame.index).all():raise ValueError('Clinical library coverage incomplete')
        rows=frame.loc[libraries].copy();rows.index=obs.index[mask]
        verified=rows.patient_key_status.isin(['title_patient_and_published_clinical_table_joined','raw_GEO_patient_id_verified'])
        if name=='GSE176021':
            parsed=rows.sample_title.str.extract(r'^([^_]+)_([^_]+)_([0-9]+)$')
            verified=parsed[0].notna()&parsed[0].eq(rows.patient_id)&rows.patient_id.ne('')
            stimulation=parsed[1].fillna('').str.startswith(('Flu','MANA'))
            obs.loc[mask,'ex_vivo_stimulation']=stimulation.to_numpy()
            obs.loc[rows.index[stimulation],'CM_composition_eligible']=False
            obs.loc[rows.index[stimulation],'capture_scope']='ex_vivo_antigen_stimulated_lymphocytes'
        for col in ['patient_id','tissue','timepoint','treatment','response_source']:
            obs.loc[mask,col]=rows[col].to_numpy()
        if name=='GSE176021':
            longitudinal=parsed[1].fillna('').str.fullmatch(r'[WM][0-9]+')&rows.tissue.str.contains('blood',case=False)
            obs.loc[rows.index[longitudinal],'timepoint']=parsed.loc[longitudinal,1]
        obs.loc[mask,'treatment_arm']=rows.treatment_arm_normalized.to_numpy()
        obs.loc[mask,'patient_identity_verified']=verified.to_numpy()
        obs.loc[mask,'patient_id_source_status']=np.where(verified,'exact_original_library_official_metadata_verified','pending')
        obs.loc[mask,'native_response_endpoint']=rows.response.to_numpy()
        obs.loc[mask,'clinical_endpoint_definition']=('author_reported_clinical_response_not_cross_study_RECIST' if name!='GSE176021'
                                                     else 'GEO_reported_MPR_not_RECIST')
        if name=='GSE169246':
            eligible=verified&rows.response_evaluable.eq('True')&rows.treatment_arm_normalized.isin(['ICI_plus_chemo','chemo'])
            obs.loc[mask,'clinical_validation_eligible']=eligible.to_numpy()
    if 'patient_identity_status' in obs:
        linked=obs.patient_identity_status.eq('source_linked')&obs.patient_id.ne('')
        obs.loc[linked,'patient_identity_verified']=True
    car=obs.analysis_role.eq('CAR_product')
    if 'response_definition' in obs:
        obs.loc[car,'clinical_endpoint_definition']=obs.loc[car,'response_definition']
    obs.loc[car,'native_response_endpoint']=obs.loc[car,'source_response']
    obs.loc[car,'clinical_validation_eligible']=(obs.loc[car,'patient_identity_verified'].astype(bool)&obs.loc[car,'native_response_endpoint'].isin(['R','NR']))
    obs['patient_block']=np.where(obs.patient_identity_verified.astype(bool),obs.dataset+'::'+obs.patient_id,'')
    obs['independent_validation_eligible']=False
    obs['independence_gate']='P3_original_study_overlap_not_seen_but_library_alias_and_prior_analysis_use_review_pending'
    obs['sample_unit']=obs.sampleID
    # Merge documented same-patient technical libraries only when every key is verified.
    known=obs.patient_identity_verified.astype(bool)&obs.tissue.ne('')
    obs.loc[known,'sample_unit']=obs.loc[known,'patient_block']+'::'+obs.loc[known,'tissue']+'::'+obs.loc[known,'timepoint']+'::'+obs.loc[known,'capture_scope']
    for col in obs:
        if obs[col].dtype==object:obs[col]=obs[col].fillna('').astype(str)
    path=args.run/'manifests/cohort_cell_registry_v2_private.tsv.gz'
    obs.to_csv(path,sep='\t')
    config['all_rows_ref']=str(path);config['schema']='immune_atlas_source_contract_v2'
    config['additional_QC_audit_sha256']=hashlib.sha256(args.audit.read_bytes()).hexdigest()
    config['production_roles']=['core','CAR_product']
    save_json_new(args.run/'config/source_contract_v2_private.json',config)
    model=obs.loc[obs.analysis_role.isin(config['production_roles'])]
    report=dict(observed_utc=datetime.now(timezone.utc).isoformat(),n_model_cells=len(model),
                n_sources=model.dataset.nunique(),n_libraries=model.sampleID.nunique(),
                source_cell_counts=model.dataset.value_counts().to_dict(),
                additional_QC_cells=len(native_obs),additional_adopted_DecontX_cells=adopted,
                n_patient_verified_cells=int(model.patient_identity_verified.astype(bool).sum()),
                n_native_clinical_endpoint_eligible_cells=int(model.clinical_validation_eligible.astype(bool).sum()),
                n_ex_vivo_stimulated_cells=int(model.ex_vivo_stimulation.astype(bool).sum()),
                supplementary_ICB_census_excluded_cells=322008,
                independent_validation_frozen=False,authoritative_release=False,
                clinical_endpoint_scope='study_specific_no_response_harmonization',
                source_review='compiled_min20_celltype_filtered_exports_and_baseline_overlap_not_CM_discovery_inputs')
    save_json_new(args.run/'manifests/cohort_freeze_v2_public.json',report)
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['contract','audit','prepared','scratch','run']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
