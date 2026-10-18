"""Schema only for every assay with complete unchanged CM02 arms and pan-T.

No expression effects, cell-state assignments, patient rows or new downloads.
"""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,4*1024**2)
    coverage=root/'fixed_clinical_panels_existing_spatial_axes_20261018_v1'
    source=coverage/'private_all28_functional_axis_assay_coverage.tsv'
    terminal=json.loads((coverage/'private_terminal.json').read_text());assert sha(source)==terminal['files'][source.name]['sha256']
    d=pd.read_csv(source,sep='\t');use=d[d.axis.eq('pCM02_regulatory_cytotoxic')&d.both_arm_sets_complete&d.pan_T_complete]
    state=Path(json.loads((root/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    axespath=state/'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json'
    axes=json.loads(axespath.read_text())
    manifest=root/'immune_same_name_family_source_review_v1/spatial_native_input_repair_v2/analysis_input_manifest.tsv'
    meta=pd.read_csv(manifest,sep='\t',keep_default_na=False)
    out=root/'CM02_existing_targeted_spatial_schema_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,coverage_sha256=sha(source),axes_sha256=sha(axespath),manifest_sha256=sha(manifest),
        eligibility='All original assays with both unchanged CM02 fixed arm gene sets and pan-T fully measurable, selected before any new expression or spatial effect. Incomplete IFN/other panels explicitly retained as limitations; no model controls silently removed.',
        scope='Input SHA, HDF5 axis/storage/metadata column names, coordinate dimensions/finiteness/uniqueness and low-cardinality grouping FIELD NAMES only. No obs values or patient identifiers read.',
        caveat='TMA/section/FOV/core/patient units must be distinguished before spatial inference; gene coverage is not phenotype/capture validation. Old additional9 counts/QC not assumed passed.')
    (out/'plan_before_targeted_input_schema.json').write_text(json.dumps(plan,indent=2)+'\n')
    rows=[]
    for item in use.itertuples():
        m=meta[meta.sample_id.eq(item.sample_id)];assert len(m)==1
        p=Path(m.iloc[0].analysis_path);assert sha(p)==item.input_sha256
        axis=next(r for r in axes if r['sample_id']==item.sample_id);assert axis['input_sha256']==item.input_sha256
        with h5py.File(p,'r') as h:
            x=h['X'];shape=list(x.shape) if isinstance(x,h5py.Dataset) else list(map(int,x.attrs['shape']))
            obs=h['obs'];var=h['var'];obsfields=list(obs.keys());layers=list(h['layers'].keys()) if 'layers' in h else []
            coordinates={}
            if 'obsm' in h:
                for k in h['obsm']:
                    a=h['obsm'][k]
                    if k in ['spatial','X_spatial'] and isinstance(a,h5py.Dataset):
                        v=np.asarray(a);assert v.ndim==2 and len(v)==shape[0] and np.isfinite(v).all()
                        coordinates[k]=dict(shape=list(v.shape),unique_xy=int(len(np.unique(v[:,:2],axis=0))),minimum=v.min(axis=0).tolist(),maximum=v.max(axis=0).tolist(),physical_unit='NOT_ASSUMED_FROM_MAGNITUDE')
            rows.append(dict(dataset_id=item.dataset_id,sample_id=item.sample_id,input_sha256=item.input_sha256,analysis_file=str(p),shape=shape,X_storage='dense' if isinstance(x,h5py.Dataset) else str(x.attrs.get('encoding-type')),X_dtype=str(x.dtype) if isinstance(x,h5py.Dataset) else str(x['data'].dtype),layers=layers,raw_present='raw/X' in h,obs_columns=obsfields,var_columns=list(var.keys()),uns_keys=list(h['uns'].keys()) if 'uns' in h else [],obsm_keys=list(h['obsm'].keys()) if 'obsm' in h else [],coordinates=coordinates,possible_grouping_fields=[c for c in obsfields if any(s in c.lower() for s in ['core','tma','fov','region','roi','patient','sample','donor'])],manifest_field_names=list(meta.columns)))
    (out/'private_all_eligible_targeted_input_schema.json').write_text(json.dumps(rows,indent=2)+'\n')
    files=[p for p in out.iterdir() if p.is_file()]
    t=dict(status='COMPLETE_ALL_COVERAGE_ELIGIBLE_TARGETED_SCHEMA_ONLY',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),assays=len(rows),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},expression_QC_or_spatial_or_patient_validation=False)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],completed_UTC=t['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))

if __name__=='__main__':main()
