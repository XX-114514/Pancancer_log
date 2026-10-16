"""I/O and frozen pool identities for the reviewed HNSCC longitudinal cohort."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd,zarr
from anndata._io.specs import read_elem
from external_functional_coverage_inputs_v1 import read_inputs

POOLS=['candidate::CD8','candidate::Macrophage','candidate::Treg','extension::B_plus_Plasma']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def inputs(root):
    r=Path(root);capture=r/'complete_external_lineage_opportunities_20261016_v2';gene=r/'complete_external_functional_gene_coverage_20261016_v1';readiness=r/'complete_external_programme_readiness_20261016_v1';paths={}
    for name,p in [('capture_audit',capture/'private_independent_opportunity_audit.json'),('gene_audit',gene/'private_independent_gene_coverage_audit.json')]:
        assert json.loads(p.read_text())['status'].startswith('PASS_ALL');paths[name]=p
    paths.update(readiness_terminal=readiness/'private_terminal.json',pool_definitions=capture/'private_all_pool_node_definitions.tsv',pool_counts=capture/'private_all_library_pool_counts.tsv.gz',metadata=capture/'private_all_library_latest_metadata.tsv',coverage=gene/'private_all48800_functional_native_coverage.tsv',clinical_definition=r/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_patient_clinical_definition.tsv',official=r/'fullgene_receiver_target_code_v1/collectri_reference_v1/official_decoupler192_ULM_source.txt')
    rt=json.loads(paths['readiness_terminal'].read_text());assert rt['status'].startswith('COMPLETE_ALL_ORIGINAL_CONTEXT_PROGRAMME_SOURCE_ROLES')
    meta=pd.read_csv(paths['metadata'],sep='\t',keep_default_na=False,dtype={'library_id':str,'patient_block':str});meta=meta[meta.arm.eq('archived_immune_atlas')&meta.dataset.eq('GSE286827')].sort_values('library_id').reset_index(drop=True);assert len(meta)==57 and meta.library_id.is_unique and meta.reviewed_patient.all()
    definitions=pd.read_csv(paths['pool_definitions'],sep='\t',dtype=str);definitions=definitions[definitions.arm.eq('archived_immune_atlas')&definitions.pool.isin(POOLS)];node_to_pool={}
    for row in definitions.itertuples():assert row.node_id not in node_to_pool;node_to_pool[row.node_id]=row.pool
    assert set(node_to_pool.values())==set(POOLS)
    counts=pd.read_csv(paths['pool_counts'],sep='\t',dtype={'library_id':str});counts=counts[counts.arm.eq('archived_immune_atlas')&counts.dataset.eq('GSE286827')&counts.pool.isin(POOLS)].set_index(['library_id','pool']).n_cells
    axis=[]
    for lib in meta.itertuples():
        for pool in POOLS:
            n=int(counts.loc[(lib.library_id,pool)]);axis.append(dict(pool_index=len(axis),library_id=lib.library_id,pool=pool,patient_block=lib.patient_block,timepoint=lib.timepoint,n_cells=n,fraction_of_all_source_library_cells=n/lib.annotated_scope_cells,treatment_arm=lib.treatment_arm,tissue=lib.tissue,reviewed_patient=True,captured_ge20=n>=20,new_identity_or_clinical_acceptance=False))
    axis=pd.DataFrame(axis);assert len(axis)==228
    atlas=r.parents[2]/'Project_v5/runs/20261013_immunotherapy_atlas_v1';z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r');fields={k:np.asarray(read_elem(z['obs'][k])).astype(str) for k in ['dataset','library_id','node_id']};assert len(fields['dataset'])==2670057
    lookup={(x.library_id,x.pool):int(x.pool_index) for x in axis.itertuples()};mapping=np.full(len(fields['dataset']),-1,np.int32)
    for i in np.flatnonzero(fields['dataset']=='GSE286827'):
        pool=node_to_pool.get(fields['node_id'][i])
        if pool is not None:mapping[i]=lookup[(fields['library_id'][i],pool)]
    assert np.array_equal(np.bincount(mapping[mapping>=0],minlength=len(axis)),axis.n_cells)
    source_totals=pd.Series(fields['library_id'][fields['dataset']=='GSE286827']).value_counts();assert np.array_equal(meta.library_id.map(source_totals),meta.annotated_scope_cells)
    data=read_inputs(r);physical=data['source_axes'][('archived_immune_atlas','GSE286827')];coverage=pd.read_csv(paths['coverage'],sep='\t');coverage=coverage[coverage.arm.eq('archived_immune_atlas')&coverage.dataset.eq('GSE286827')];assert len(coverage)==800 and coverage.technical_measurement_gate.all()
    measured=np.asarray(read_elem(z['varm']['measured_in_source']),bool);order=np.asarray(read_elem(z['uns']['source_order'])).astype(str).tolist();native=measured[:,order.index('GSE286827')];symbols=np.asarray(read_elem(z['var']['gene_symbol'])).astype(str)
    hashes={k:sha(p) for k,p in paths.items()};hashes.update(raw_pool_mapping=hashlib.sha256(mapping.tobytes()).hexdigest(),atlas_release=sha(atlas/'manifests/immune_authority_release_public.json'))
    return dict(axis=axis,mapping=mapping,zarr=z,physical=physical,parameters=data['parameters'],masks=data['masks'],references=data['references'],coverage=coverage,native_mask=native,symbols=symbols,input_sha256=hashes,official_path=paths['official'],coverage_input_sha256=data['input_sha256'])
