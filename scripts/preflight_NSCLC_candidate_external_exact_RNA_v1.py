"""Metadata-only exact-RNA routes for every original internal candidate/source.

No external expression, association effects, clinical outcomes or downloads.
"""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

from NSCLC_external_lineage_inputs_v2 import load
from external_functional_coverage_inputs_v1 import read_inputs
from describe_NSCLC_native_interface_mechanistic_followup_v1 import read_scope, SCOPE_SHA
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); budget=resource(root,8*1024**2)
    proof=root/'NSCLC_native_interface_mechanistic_followup_20261018_v1/private_independent_mechanistic_description_audit.json'
    assert sha(proof)=='98a2fd9abf5c5e125dbc2f33616032b87a1cce033899d225f102ba25b44052ca'
    scope,candidates=read_scope(dict(root=root))
    out=root/'NSCLC_candidate_external_exact_RNA_preflight_20261018_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),scope_sha256=SCOPE_SHA,
        mechanistic_audit_sha256=sha(proof),budget=budget,
        coverage='All original candidates, all existing61 source-arms and6760 frozen role/library units. Do not remove weak biological candidates or sources.',
        metadata_only=True,new_expression_downloads_effects_or_clinical_outcomes=False,
        desired_measurement='True mean per-cell log1pCP10K ligand; same-cell complete curated receptor at count>=1 and>=2 sensitivity. Original full physical RNA denominator.',
        limits='Mapping and readable raw containers do not prove actual expression, receptor co-detection, comparable cell identity, valid model or replication.')
    (out/'plan_before_any_external_exact_RNA_read.json').write_text(json.dumps(plan,indent=2)+'\n')
    data=load();lib=data['lib'];axes=read_inputs(root)['source_axes']
    assert len(axes)==61 and len(lib)==1690
    desired=sorted(set(scope['ligand_specs'])|{g for s in scope['ligand_specs'].values() for c in s['complexes'] for g in c})
    marker_root=root/'NSCLC_native_MFA_external_role_marker_profiles_20261016_v2'
    manifest=json.loads((marker_root/'private_terminal.json').read_text())
    marker=marker_root/'private_fixed_marker_pool_sums.npz'
    assert sha(marker)==manifest['files'][marker.name]['sha256']
    with np.load(marker,allow_pickle=False) as z:
        existing_markers=set(z['markers'].astype(str));counts=z['n_cells']
    assert len(counts)==len(lib)*4
    oldbase=root.parents[2]/'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1/by_gse'
    sources=[];genes=[];routes=[]
    containers={}
    for arm,path in [('frozen_original_annotation',data['old_zarr']),('archived_immune_atlas',data['atlas'])]:
        z=zarr.open_group(str(path),mode='r');raw=z['X'] if arm=='frozen_original_annotation' else z['layers']['counts']
        assert raw.attrs.get('encoding-type')=='csr_matrix'
        nr,nf=map(int,raw.attrs['shape'])
        assert raw['indptr'].shape==(nr+1,) and raw['indices'].shape==raw['data'].shape
        containers[arm]=dict(rows=nr,features=nf,stored_nonzero_coordinates=int(raw['data'].shape[0]),
            route='Original canonical global rows plus frozen roles' if arm=='frozen_original_annotation' else 'Original archive rows plus frozen roles',
            raw_count_values_read=False)
    for (arm,dataset),axis in sorted(axes.items()):
        mapping={g:(next(iter(axis['lookup'][g])) if g in axis['lookup'] and len(axis['lookup'][g])==1 else -1) for g in desired}
        seen={p for p in mapping.values() if p>=0}
        for g,p in list(mapping.items()):
            if p>=0 and sum(v==p for v in mapping.values())>1:mapping[g]=-1
        libraries=lib.index[lib.arm.eq(arm)&lib.dataset.eq(dataset)].to_numpy(int)
        sources.append(dict(arm=arm,dataset=dataset,libraries=len(libraries),role_units=len(libraries)*4,
            unique_measured_requested_genes=sum(v>=0 for v in mapping.values()),requested_genes=len(desired),
            pool_identity='Frozen subtype extension, not exact native full-author views',
            B_pool='Frozen original B-major only; plasma not silently added' if arm=='frozen_original_annotation' else 'Frozen archive B and Plasma'))
        if arm=='frozen_original_annotation':
            base=oldbase/dataset;m=json.loads((base/'manifest.json').read_text())
            assert m['status']=='COMPLETE_EXACT_MEMBER_EXPRESSION_EVIDENCE'
            for name in ['units.tsv','genes.tsv','sum_log1p_CP10K.npz','positive_counts.npz']:
                assert (base/name).is_file() and name in m['output_sha256']
            oldgenes=pd.read_csv(base/'genes.tsv',sep='\t',keep_default_na=False)
            assert oldgenes.gene.is_unique and oldgenes.core_column.is_unique
        for gene,position in mapping.items():
            genes.append(dict(arm=arm,dataset=dataset,gene=gene,source_column=position,
                status='UNIQUE_MEASURED_MAPPING' if position>=0 else 'MISSING_OR_AMBIGUOUS_MAPPING',
                existing44_marker_cache_has_gene=gene in existing_markers,
                exact_mean_route='Existing full-gene native-unit cache' if arm=='frozen_original_annotation' else 'Existing raw count stream',
                actual_same_cell_receptor_route='Existing raw count stream required'))
        for c in candidates.itertuples():
            spec=scope['ligand_specs'][c.ligand]
            complete=[v for v in spec['complexes'] if all(mapping[g]>=0 for g in v)]
            routes.append(dict(record=c.record,arm=arm,dataset=dataset,ligand_unique_measured=mapping[c.ligand]>=0,
                complete_functional_complexes_measured=len(complete),functional_complexes_requested=len(spec['complexes']),
                measurement_opportunity=bool(mapping[c.ligand]>=0 and complete),
                ligand_in44_marker_cache=c.ligand in existing_markers,
                exact_joint_counts_available_in_pseudobulk=False,
                in_native_model_dataset=dataset=='GSE243013'))
    artifacts={'private_all61_exact_RNA_source_routes.tsv':pd.DataFrame(sources),
        'private_all_requested_gene_mappings.tsv':pd.DataFrame(genes),
        'private_all_candidate_source_routes.tsv':pd.DataFrame(routes)}
    for name,table in artifacts.items():table.to_csv(out/name,sep='\t',index=False)
    terminal=dict(status='COMPLETE_METADATA_ONLY_EXACT_RNA_ROUTES_NOT_EXPRESSION_OR_MODEL_QUALIFICATION',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),scope_sha256=SCOPE_SHA,
        libraries=len(lib),source_arms=len(axes),frozen_role_units=len(counts),input_hashes=data['hashes'],containers=containers,
        files={name:dict(sha256=sha(out/name),bytes=(out/name).stat().st_size) for name in artifacts},
        storage_projection='Persistent selected-panel per-library sums and joint counts are small; precise raw-read worker/RAM and independent audit contract required before expression extraction.',
        no_new_expression_effects_or_downloads=True)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],completed_UTC=terminal['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
