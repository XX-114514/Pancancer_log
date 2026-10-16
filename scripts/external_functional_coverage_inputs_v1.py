"""Frozen reference and physical native gene axes, no expression or endpoint effects."""
import hashlib,json
from collections import defaultdict
from pathlib import Path
import numpy as np,pandas as pd,zarr
from anndata._io.specs import read_elem


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_inputs(root):
    r=Path(root);pre=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1';par=r/'NSCLC_context_frozen_parameters_20261016_v1';parameter_path=par/'private_all53_frozen_gene_parameters.tsv'
    pt=json.loads((par/'private_terminal.json').read_text());assert sha(parameter_path)==pt['files'][parameter_path.name]
    parameters=pd.read_csv(parameter_path,sep='\t');scope_path=r/'NSCLC_context_native_programme_projection_20261016_v1/private_all56_projection_scope.tsv';scope=pd.read_csv(scope_path,sep='\t',keep_default_na=False);assert len(scope)==56 and parameters.groupby(['context','type','program']).ngroups==53
    cov_path=pre/'private_all1200_functional_scope.tsv';old=pd.read_csv(cov_path,sep='\t');masks=old[['mask_id','context','program']].drop_duplicates().sort_values('mask_id');assert len(masks)==16
    refroot=r/'CRC_fine_pair_DIALOGUE_code_20261014_v1/hallmark50_public_reference_v3';prov=refroot/'reference_provenance.json';reference={};ref_hash={}
    for record in json.loads(prov.read_text())['records']:
        p=refroot/(record['name']+'.json');assert sha(p)==record['sha256'];reference[record['name']]=set(json.loads(p.read_text())[record['name']]['geneSymbols']);ref_hash[record['name']]=sha(p)
    assert len(reference)==50
    registry_path=r/'complete_external_lineage_opportunities_20261016_v2/private_all_library_latest_metadata.tsv';registry=pd.read_csv(registry_path,sep='\t',keep_default_na=False,dtype={'dataset':str});entries=registry[['arm','dataset']].drop_duplicates().sort_values(['arm','dataset']);assert len(entries)==61
    legacy_path=r/'all_CM_state_native_input_preflight_20261014_v1/all_source_native_scope.tsv';legacy=pd.read_csv(legacy_path,sep='\t').set_index('gse_id');genedir=r.parents[2]/'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1/by_gse';sources={};source_hash={}
    for row in entries[entries.arm.eq('frozen_original_annotation')].itertuples():
        p=genedir/row.dataset/'genes.tsv';assert sha(p)==legacy.loc[row.dataset,'native_gene_axis_sha256'];g=pd.read_csv(p,sep='\t');assert g.gene.is_unique and g.core_column.is_unique and len(g)==legacy.loc[row.dataset,'native_genes']
        names=g.gene.astype(str).tolist();lookup={name:{i} for i,name in enumerate(names) if name not in ['', 'nan','None'] and not name.startswith('AMBIGUOUS::')};sources[(row.arm,row.dataset)]=dict(lookup=lookup,native=set(c for ss in lookup.values() for c in ss),physical_features=len(names),mapping='frozen_native_symbol_axis_exact_no_new_alias_inference');source_hash[row.dataset]=sha(p)
    ar=r.parents[2]/'Project_v5/runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr';z=zarr.open_group(str(ar),mode='r');order=np.asarray(read_elem(z['uns']['source_order'])).astype(str).tolist();measured=np.asarray(read_elem(z['varm']['measured_in_source']),bool)
    columns={field:np.asarray(read_elem(z['var'][field])).astype(str) for field in ['gene_symbol','source_gene_aliases','canonical_gene_identity','feature_identity_status']};N=len(columns['gene_symbol']);assert measured.shape==(N,len(order))
    for row in entries[entries.arm.eq('archived_immune_atlas')].itertuples():
        si=order.index(row.dataset);native=set();lookup=defaultdict(set)
        for i in np.flatnonzero(measured[:,si]):
            if columns['canonical_gene_identity'][i].startswith('AMBIGUOUS::') or columns['feature_identity_status'][i]=='ambiguous_source_feature':continue
            native.add(int(i));names={columns['gene_symbol'][i],columns['canonical_gene_identity'][i]}
            if columns['feature_identity_status'][i]=='stable_Ensembl_ID':names.update(columns['source_gene_aliases'][i].split('|'))
            for name in names-{'','nan','None','__SOURCE_VALUE_MISSING__'}:lookup[name].add(int(i))
        sources[(row.arm,row.dataset)]=dict(lookup=dict(lookup),native=native,physical_features=int(measured[:,si].sum()),mapping='native_unambiguous_physical_column_exact_symbol_canonical_or_existing_stable_Ensembl_alias')
    ns=r/'NSCLC_full_native_gene_units_20261016_v1/private_native_gene_axis.tsv';ng=pd.read_csv(ns,sep='\t').gene.astype(str).tolist();assert len(ng)==len(set(ng))
    calibration=dict(lookup={name:{i} for i,name in enumerate(ng)},native=set(range(len(ng))),physical_features=len(ng),mapping='original_NSCLC_exact_physical_symbol_calibration')
    hashes=dict(parameters=sha(parameter_path),programme_scope=sha(scope_path),original_function_scope=sha(cov_path),reference_provenance=sha(prov),source_registry=sha(registry_path),legacy_source_scope=sha(legacy_path),NSCLC_gene_axis=sha(ns),atlas_source_order=hashlib.sha256(json.dumps(order).encode()).hexdigest(),atlas_measured_mask=hashlib.sha256(measured.tobytes()).hexdigest())
    for k,v in columns.items():hashes['atlas_'+k]=hashlib.sha256(json.dumps(v.tolist(),ensure_ascii=False,separators=(',',':')).encode()).hexdigest()
    return dict(parameters=parameters,programmes=scope,masks=masks,references=reference,source_axes=sources,NSCLC_calibration=calibration,original_scope=old,input_sha256=hashes,reference_sha256=ref_hash,legacy_gene_axis_sha256=source_hash)
