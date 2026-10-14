#!/usr/bin/env python3
"""Prepare composition scopes with measured major denominators, no NA-to-zero."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
from anndata._io.specs import read_elem
import numpy as np
import pandas as pd
import zarr
from immune_atlas_io import save_json_new


def truth(values):return values.astype(str).str.lower().isin(['true','1'])


def material_classes(obs):
    """Conservative source-material grouping; generic LN is not a metastasis."""
    tissue=obs.tissue.astype('string').fillna('').astype(str).str.lower().str.strip()
    capture=obs.capture_scope.astype('string').fillna('').astype(str)
    material=pd.Series('unknown_material',index=obs.index)
    whole=capture.str.contains('whole_tumor|serial_tumor',regex=True)
    material.loc[whole|tissue.isin(['tumor','glioblastoma tumor','primary tumor','primary tumor (multiple regions)'])]='primary_tumor'
    material.loc[tissue.str.contains('metastas',regex=False)]='metastatic_tumor'
    material.loc[tissue.isin(['ln','lymph_node','lymph node'])]='lymph_node_involvement_unknown'
    material.loc[tissue.str.contains('normal|adjacent',regex=True)]='normal_adjacent'
    material.loc[tissue.str.contains('blood|pbmc',regex=True)]='blood'
    material.loc[obs.analysis_role.astype(str).eq('CAR_product')]='infusion_product'
    return material


def prepare(args):
    obs=read_elem(zarr.open_group(str(args.scratch/'immune_refined.zarr'),mode='r')['obs'])
    eligible=truth(obs.CM_composition_eligible)&~truth(obs.annotation_exclude_from_primary)&~truth(obs.ex_vivo_stimulation)&~obs.major.astype(str).isin(['Unresolved','T_NK_boundary'])
    car=obs.analysis_role.astype(str).eq('CAR_product')
    capture=obs.capture_scope.astype(str)
    scopes=dict(whole_TME=eligible&~car&capture.str.contains('whole_tumor|serial_tumor|whole_tumor_patient'),
                immune_enriched=eligible&~car&capture.str.contains('CD45_enriched|lymphocyte_enriched|T_enriched'),
                CAR_product=eligible&car,
                T_lineage_biopsy=eligible&~car&obs.major.astype(str).eq('T'))
    mode=getattr(args,'scope_mode','broad')
    material=material_classes(obs)
    if mode=='material':
        scopes={}
        capture_families={
            'whole_TME':~car&capture.str.contains('whole_tumor|serial_tumor'),
            'immune_enriched':~car&capture.str.contains('CD45_enriched|lymphocyte_enriched|T_enriched'),
            'mixed_enriched_nuclei':~car&capture.eq('mixed_tumor_nuclei_NFAT_enrichment')}
        for family,family_mask in capture_families.items():
            for label in ['primary_tumor','metastatic_tumor','blood','normal_adjacent','lymph_node_involvement_unknown','unknown_material']:
                selected=eligible&family_mask&material.eq(label)
                if not selected.any():continue
                scopes[family+'_'+label]=selected
                if family!='mixed_enriched_nuclei':
                    scopes['T_'+family+'_'+label]=selected&obs.major.astype(str).eq('T')
    reports=[]
    for name,mask in scopes.items():
        frame=obs.loc[mask].copy();folder=args.run/'CM'/name/'inputs';folder.mkdir(parents=True,exist_ok=False,mode=0o700)
        counts=pd.crosstab(frame.subcluster.astype(str),frame.sample_unit.astype(str))
        major_map=frame[['subcluster','major']].astype(str).drop_duplicates().set_index('subcluster').major
        if major_map.index.duplicated().any():raise ValueError('Subcluster crosses major identities')
        major_counts=counts.groupby(major_map.reindex(counts.index)).sum()
        valid_units=counts.columns[counts.sum(axis=0)>=100]
        major_coverage=(major_counts[valid_units]>=20).mean(axis=1)
        major_keep=major_coverage.index[major_coverage>=.8]
        if len(major_keep):valid_units=valid_units[(major_counts.loc[major_keep,valid_units]>=20).all(axis=0)]
        else:valid_units=[]
        counts=counts.loc[major_map.reindex(counts.index).isin(major_keep),valid_units]
        retained_nodes=(counts.sum(axis=1)>=100)&((counts>0).sum(axis=1)>=4)
        counts=counts.loc[retained_nodes]
        major_map=major_map.reindex(counts.index)
        denominators=major_counts.reindex(index=major_map.to_numpy(),columns=counts.columns).to_numpy()
        if denominators.size and (denominators<20).any():raise ValueError('Unmeasured major denominator entered CM')
        frequency=counts/denominators if denominators.size else counts.astype(float)
        finite_variation=np.isfinite(frequency).all(axis=1)&frequency.var(axis=1).gt(1e-8)
        frequency=frequency.loc[finite_variation];major_map=major_map.reindex(frequency.index)
        metadata=[]
        for unit,group in frame.groupby('sample_unit',observed=True):
            if unit not in frequency.columns:continue
            if group.dataset.nunique()!=1:raise ValueError('Composition unit crosses original studies')
            patients=group.patient_block.astype(str).unique()
            verified=bool(len(patients)==1 and patients[0])
            metadata.append(dict(sampleID=str(unit),dataset=str(group.dataset.iloc[0]),
                patient_block=patients[0] if verified else '',
                independent_unit=patients[0] if verified else 'unverified_library::'+str(unit),
                patient_verified=verified,n_cells=len(group),
                cancer_type='|'.join(sorted(group.cancer_type.astype(str).unique())),
                capture_scope='|'.join(sorted(group.capture_scope.astype(str).unique())),
                material_class='|'.join(sorted(material.loc[group.index].unique())),
                original_tissue='|'.join(sorted(group.tissue.astype('string').fillna('').astype(str).unique())),
                n_original_libraries=group.sampleID.nunique(),
                timepoint='|'.join(sorted(group.timepoint.astype(str).unique()))))
        meta=pd.DataFrame(metadata,columns=['sampleID','dataset','patient_block','independent_unit','patient_verified','n_cells','cancer_type','capture_scope','material_class','original_tissue','n_original_libraries','timepoint'])
        if mode=='material' and meta.material_class.str.contains('|',regex=False).any():
            raise ValueError('Composition library crosses material classes')
        export=frequency.copy();export.insert(0,'majorCluster',major_map)
        export.to_csv(folder/'frequency_private.tsv',sep='\t',index_label='subCluster')
        meta.to_csv(folder/'sample_metadata_private.tsv',sep='\t',index=False)
        major_coverage.to_csv(folder/'major_denominator_coverage_private.tsv',sep='\t')
        report=dict(scope=name,n_input_cells=len(frame),n_complete_denominator_units=frequency.shape[1],
                n_nodes=frequency.shape[0],n_major_denominators=len(major_keep),
                major_denominators=list(major_keep),n_sources=meta.dataset.nunique(),
                n_verified_patient_blocks=meta.loc[meta.patient_verified,'patient_block'].nunique(),
                n_unverified_library_units=int((~meta.patient_verified).sum()),
                n_independent_blocks=meta.independent_unit.nunique(),
                measured_major_min_cells=20,total_unit_min_cells=100,major_unit_coverage_min=.8,
                ready=bool(frequency.shape[0]>=4 and frequency.shape[1]>=12 and meta.independent_unit.nunique()>=8),
                interpretation='discovery_descriptive;unknown_patient_units_exploratory;no_clinical_response_tests',
                scope_mode=mode,material_classes=sorted(material.loc[frame.index].unique()),
                material_mixed=len(material.loc[frame.index].unique())>1,
                composition_capture_bias_not_removed=True,
                longitudinal_patient_mean_is_not_response_or_timepoint_test=True,
                ex_vivo_excluded=True,selectively_retained_compiled_exports_excluded=True)
        reports.append(report);save_json_new(folder/'input_scope_public.json',report)
    name='CM_material_preparation_public.json' if mode=='material' else 'CM_preparation_public.json'
    save_json_new(args.run/'manifests'/name,dict(observed_utc=datetime.now(timezone.utc).isoformat(),scopes=reports,
                   fixed_P3_9CM_projection_completed=False,new_CM_discovery_not_validation_of_retrained_W=True))
    print(json.dumps(reports),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--scratch',type=Path,required=True)
    p.add_argument('--scope-mode',choices=['broad','material'],default='broad')
    prepare(p.parse_args())
