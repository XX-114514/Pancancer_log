"""Metadata-only inputs for full frozen annotation and archived immune atlas arms."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd,zarr
from anndata._io.specs import read_elem

STRATA=['cancer_label','tissue','sample_class','timepoint','treatment_arm','capture_scope','ex_vivo','official_stage']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def atom(values):
    v={str(x) for x in values}
    return next(iter(v)) if len(v)==1 else 'CONFLICTING_WITHIN_LIBRARY_METADATA'


def load(root):
    """Read inherited metadata only; source axes and matrices are never changed."""
    r=Path(root);f=r/'all_CM_state_native_input_preflight_20261014_v1';iv=json.loads((f/'inventory.json').read_text());proof=r/'all_CM_state_row_preflight_audit_20261014_v1/independent_audit.json'
    assert json.loads(proof.read_text())['status'].startswith('PASS_ALL4532933') and json.loads(proof.read_text())['producer_inventory_sha256']==sha(f/'inventory.json')
    unit_path=f/'private_all_frozen_sample_node_units.tsv.gz';scope_path=r/'immune_same_name_family_source_review_v1/patient_block_provenance_overlay_v1/private_all1290_patient_provenance_overlay.tsv'
    assert sha(unit_path)==iv['private_unit_table_sha256'] and sha(scope_path)==iv['input_sha256']['scope']
    units=pd.read_csv(unit_path,sep='\t',keep_default_na=False,dtype={'gse_id':str,'sample_uid':str,'node_id':str});scope=pd.read_csv(scope_path,sep='\t',keep_default_na=False,dtype={'sampleID':str,'gse_id':str,'effective_patient_block':str,'effective_resampling_block':str});assert len(units)==35567 and scope.sampleID.is_unique and len(scope)==1290
    assert set(units.sample_uid)==set(scope.sampleID) and units.n_cells.sum()==iv['selected_cells']
    scope=scope.set_index('sampleID');assert units.gse_id.eq(units.sample_uid.map(scope.gse_id)).all()
    assert units.unit_kind.eq(units.sample_uid.map(scope.unit_kind)).all() and units.effective_resampling_block.eq(units.sample_uid.map(scope.effective_resampling_block)).all()
    old=units.rename(columns={'gse_id':'dataset','sample_uid':'library_id'})[['dataset','library_id','node_id','major','n_cells']].copy();old.insert(0,'arm','frozen_original_annotation')
    libraries=[]
    for key,s in scope.iterrows():
        reviewed=s.unit_kind=='reviewed_patient';patient=str(s.effective_patient_block)
        assert not reviewed or patient not in ['', 'unknown','None','nan']
        libraries.append(dict(arm='frozen_original_annotation',dataset=s.gse_id,library_id=key,patient_block=patient if reviewed else 'unreviewed_library::'+key,reviewed_patient=reviewed,resampling_block=str(s.effective_resampling_block),
            patient_evidence=s.effective_patient_evidence_status,cancer_label=s.cancerType,tissue=s.tissue,sample_class=s.sample_class,timepoint=s.treatment_phase,treatment_arm='not_reviewed',capture_scope=s.sorting_enrichment_status,ex_vivo='not_reviewed',official_stage='not_reviewed',
            recorded_clinical_endpoint='reviewed_MMR_only' if str(s.reviewed_mmrstatus) not in ['', 'unknown','not_reviewed','None','nan'] else 'not_reviewed',clinical_endpoint_review='inherited_scope_not_new_endpoint_acceptance',clinical_eligible_cell_count=0,latest_sidecar='frozen_patient_provenance_overlay',new_clinical_model_accepted=False,exact_NSCLC_endpoint_bridge=False))
    project=r.parents[2];ar=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1';release=ar/'manifests/immune_authority_release_public.json';rel=json.loads(release.read_text());assert rel['all_checks_pass'] and rel['technical_authority_released']
    z=zarr.open_group(str(ar/'outputs/immune_validation_authority_v1.zarr'),mode='r');assert tuple(z['layers/counts'].attrs['shape'])==(2670057,49563)
    fields=['dataset','library_id','node_id','major','patient_id','patient_block','patient_identity_verified','timepoint','tissue','treatment_arm','capture_scope','ex_vivo_stimulation','cancer_type','native_response_endpoint','clinical_endpoint_definition_status','clinical_validation_eligible']
    md=pd.DataFrame({k:pd.Categorical(np.asarray(read_elem(z['obs'][k])).astype(str)) for k in fields});assert len(md)==2670057
    text=read_elem(z['uns']['subcluster_catalog_json']);text=text.item() if isinstance(text,np.ndarray) else text;catalog=pd.DataFrame(json.loads(text));assert len(catalog)==270 and catalog.node_id.is_unique
    assert md.node_id.isin(catalog.node_id).all();cm=catalog.set_index('node_id').major;assert md.major.astype(str).eq(md.node_id.astype(str).map(cm)).all()
    new=md.groupby(['dataset','library_id','node_id','major'],observed=True,sort=True).size().reset_index(name='n_cells');new=new.astype({'dataset':str,'library_id':str,'node_id':str,'major':str});new.insert(0,'arm','archived_immune_atlas')
    H=r/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv';h=pd.read_csv(H,sep='\t',keep_default_na=False,dtype={'library_id':str,'patient_registry':str,'timepoint_registry':str}).set_index('library_id');assert len(h)==57 and h.index.is_unique
    clinical=r/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_patient_clinical_definition.tsv';hc=pd.read_csv(clinical,sep='\t',keep_default_na=False,dtype={'patient_id':str}).set_index('patient_id');assert hc.index.is_unique and set(h.patient_registry)==set(hc.index)
    S=r/'GSE305078_official_sampling_sidecar_20261016_v1/private_exact_official_sampling_sidecar.tsv';sampling=pd.read_csv(S,sep='\t',keep_default_na=False,dtype={'library_id':str}).set_index('library_id');assert len(sampling)==32 and sampling.index.is_unique
    assert sampling.exact_official_key_matched.astype(str).eq('True').all()
    Hkeys=set();Skeys=set()
    for (source,key),g in md.groupby(['dataset','library_id'],observed=True,sort=True):
        raw={name:atom(g[name]) for name in fields if name not in ['dataset','library_id','node_id','major','clinical_validation_eligible']};reviewed=raw['patient_identity_verified']=='True';patient=raw['patient_block'] if reviewed else 'unreviewed_library::'+key
        assert not reviewed or patient not in ['', 'unknown','None','nan','CONFLICTING_WITHIN_LIBRARY_METADATA']
        row=dict(arm='archived_immune_atlas',dataset=str(source),library_id=str(key),patient_block=patient,reviewed_patient=reviewed,resampling_block=patient,patient_evidence='inherited_atlas_patient_review' if reviewed else 'unreviewed_sample_fallback',
            cancer_label=raw['cancer_type'],tissue=raw['tissue'],sample_class='author_atlas_tissue_only',timepoint=raw['timepoint'],treatment_arm=raw['treatment_arm'],capture_scope=raw['capture_scope'],ex_vivo=raw['ex_vivo_stimulation'],official_stage='not_reviewed',recorded_clinical_endpoint=raw['native_response_endpoint'],clinical_endpoint_review=raw['clinical_endpoint_definition_status'],
            clinical_eligible_cell_count=int(g.clinical_validation_eligible.astype(str).eq('True').sum()),latest_sidecar='inherited_atlas_fields',new_clinical_model_accepted=False,exact_NSCLC_endpoint_bridge=False)
        if source=='GSE286827':
            assert key in h.index;v=h.loc[key];Hkeys.add(key);assert raw['patient_id']==str(v.patient_registry) and raw['timepoint']==str(v.timepoint_registry)
            assert str(v.patient_registry) in hc.index
            row.update(patient_block='GSE286827::author_patient::'+str(v.patient_registry),resampling_block='GSE286827::author_patient::'+str(v.patient_registry),reviewed_patient=True,patient_evidence='reviewed_author_table_exact_library_join',timepoint=str(v.timepoint_registry),treatment_arm=str(v.neo_type),recorded_clinical_endpoint='continuous_author_Tumor.regression',clinical_endpoint_review='continuous_pathology_reviewed_author_binary_conflicts_preserved',latest_sidecar='author_Table_S1_and_frozen_continuous_definition')
        if source=='GSE305078':
            assert key in sampling.index;v=sampling.loc[key];Skeys.add(key);assert not reviewed and str(v.patient_identity_newly_verified)=='False'
            row.update(timepoint=str(v.author_sampling_phase),tissue=str(v.official_tissue),official_stage=str(v.official_stage),latest_sidecar='official_GEO_double_identifier_sampling_only')
        libraries.append(row)
    assert Hkeys==set(h.index) and Skeys==set(sampling.index)
    lib=pd.DataFrame(libraries);counts=pd.concat([old,new],ignore_index=True);assert not lib.duplicated(['arm','dataset','library_id']).any()
    sizes=counts.groupby(['arm','dataset','library_id'],sort=True).n_cells.sum();assert len(sizes)==len(lib)==1690
    lib['annotated_scope_cells']=sizes.reindex(pd.MultiIndex.from_frame(lib[['arm','dataset','library_id']])).to_numpy(int)
    raw_major={arm:sorted(g.major.unique()) for arm,g in counts.groupby('arm')};assert all(len(x)==18 for x in raw_major.values())
    pool_nodes={}
    for arm in raw_major:
        cats=counts[counts.arm.eq(arm)][['node_id','major']].drop_duplicates();assert cats.node_id.is_unique
        pool_nodes[arm]={'major::'+major:set(cats.loc[cats.major.eq(major),'node_id']) for major in raw_major[arm]}
    pools=pool_nodes['archived_immune_atlas'];candidate=catalog.subtype_evidence_status.eq('multigene_subtype_candidate')
    pools.update({'candidate::CD8':set(catalog.loc[catalog.major.eq('T')&catalog.subtype_candidate.eq('CD8_Effector')&candidate,'node_id']),
        'candidate::Treg':set(catalog.loc[catalog.major.eq('T')&catalog.subtype_candidate.eq('Treg')&candidate,'node_id']),
        'candidate::Macrophage':set(catalog.loc[catalog.major.eq('Myeloid')&catalog.subtype_candidate.str.startswith('Macro_')&candidate,'node_id']),
        'extension::B_plus_Plasma':set(catalog.loc[catalog.major.isin(['B','Plasma']),'node_id'])})
    assert all(len(v)>0 for v in pools.values())
    shared=set(counts[counts.arm.eq('frozen_original_annotation')].dataset)&set(counts[counts.arm.eq('archived_immune_atlas')].dataset)
    lib['accession_present_in_both_arms']=lib.dataset.isin(shared);lib['cross_source_patient_overlap_reviewed']=False;lib['untouched_independent_validation_accepted']=False
    hashes=dict(baseline_inventory=sha(f/'inventory.json'),baseline_audit=sha(proof),baseline_units=sha(unit_path),patient_overlay=sha(scope_path),atlas_release=sha(release),HNSCC_table=sha(H),HNSCC_continuous_definition=sha(clinical),TNBC_official_sidecar=sha(S),catalog_text=hashlib.sha256(str(text).encode()).hexdigest())
    return dict(libraries=lib,counts=counts,raw_majors=raw_major,pool_nodes=pool_nodes,catalog=catalog,input_sha256=hashes)
