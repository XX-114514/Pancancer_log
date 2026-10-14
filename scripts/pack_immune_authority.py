#!/usr/bin/env python3
"""Full-native-gene immune atlas, independent whole-row audit, final-volume copy.

Temporary processing never becomes the authority. Clinical missingness and
independence limitations remain explicit in the released metadata.
"""
import argparse
from datetime import datetime,timezone
import gc,hashlib,json,shutil,subprocess
from pathlib import Path
import anndata as ad
from anndata._io.specs import read_elem,write_elem
import numpy as np
import pandas as pd
from scipy import sparse
import zarr
from immune_atlas_io import load_contract,normalized,projected,save_json_new
from merge_public_preprocessed_zarr import CSRAppender,csr_read

RESERVE=500_000_000_000


def capacity(path,group):
    physical=shutil.disk_usage(path).free
    mount=subprocess.check_output(['df','--output=source',str(path)],text=True).splitlines()[-1].strip()
    result=subprocess.run(['quota','-g','-v',group],capture_output=True,text=True,check=True)
    lines=result.stdout.splitlines()
    for i,line in enumerate(lines):
        if mount not in line:continue
        suffix=line.split(mount,1)[1].strip()
        values=(suffix if suffix else lines[i+1].strip()).split()
        used=int(values[0].rstrip('*'));limits=[int(v) for v in values[1:3] if int(v)>0]
        available=(min(limits)-used)*1024 if limits else physical
        return min(physical,max(0,available))
    raise RuntimeError('Cannot establish effective final-volume group quota')


def pack(args):
    config,sources,var,identities,diagnostics=load_contract(args.contract,['core','CAR_product'])
    model=ad.read_zarr(args.scratch/'immune_refined.zarr')
    cnv=pd.read_csv(args.CNV_results or args.run/'tables/all_core_CNV_private.tsv.gz',sep='\t',index_col=0,keep_default_na=False)
    if not cnv.index.is_unique or not cnv.index.equals(model.obs.index):raise ValueError('CNV cell mapping differs')
    if not np.array_equal(cnv.sampleID.astype(str),model.obs.sampleID.astype(str)):raise ValueError('CNV sample mapping differs')
    obs=model.obs.copy()
    for key in cnv:
        if key=='sampleID':continue
        if key=='cnv_score':obs[key]=pd.to_numeric(cnv[key],errors='coerce').astype(np.float32)
        else:obs[key]=cnv[key]
    target=args.scratch/'immune_fullgene_authority_candidate_v1.zarr'
    if target.exists():raise ValueError('Existing fullgene candidate preserved')
    target.mkdir(mode=0o700)
    g=zarr.open_group(str(target),mode='w');g.attrs.update({'encoding-type':'anndata','encoding-version':'0.1.0'})
    layers=g.create_group('layers');layers.attrs.update({'encoding-type':'dict','encoding-version':'0.1.0'})
    arrays=dict(counts=CSRAppender(layers,'counts','i4',len(var)),
                ambient_removed_counts=CSRAppender(layers,'ambient_removed_counts','f8',len(var)),
                X=CSRAppender(g,'X','f4',len(var)))
    rows=0;denominators=np.zeros(len(obs),dtype=np.float64)
    for source in sources:
        local=obs.loc[obs.dataset.astype(str).eq(source.record['dataset'])]
        expected=source.record['dataset']+'::'+source.obs.index[source.rows].astype(str)
        if not np.array_equal(local.index.astype(str),expected):raise ValueError('Fullgene source ID/order mismatch')
        for lo in range(0,len(local),2048):
            ids=np.arange(lo,min(lo+2048,len(local)))
            raw,chosen,total=source.analysis(ids)
            raw=projected(raw.astype(np.int64),source.global_columns,len(var))
            chosen=projected(chosen,source.global_columns,len(var))
            removal=(raw.astype(np.float64)-chosen).tocsr();removal.eliminate_zeros()
            if removal.nnz and removal.data.min()<-1e-7:raise ValueError('Fullgene correction added counts')
            removal.data=np.maximum(removal.data,0);removal.eliminate_zeros()
            for key,block in [('counts',raw),('ambient_removed_counts',removal),('X',normalized(chosen,total))]:arrays[key].append(block)
            denominators[rows:rows+len(ids)]=total;rows+=len(ids)
        print(json.dumps(dict(phase='fullgene_pack_source',dataset=source.record['dataset'],n_output_cells=rows)),flush=True)
    if rows!=len(obs):raise ValueError('Fullgene rows incomplete')
    obs['native_analysis_denominator']=denominators
    obs['data_authority']='immune_native_count_archive_v1_technical_release'
    obs['identity_scope']='native_multigene_major_descriptive_subclusters_fine_review_pending'
    obs['malignancy_scope']='per_original_library_CNV_single_method_not_dual_method_gold'
    var['in_modeling_feature_space']=var.index.isin(model.var_names)
    var['highly_variable']=False
    var.loc[model.var_names,'highly_variable']=model.var.highly_variable.to_numpy()
    measured=np.column_stack([var.index.isin(s.canonical_genes) for s in sources])
    write_elem(g,'obs',obs);write_elem(g,'var',var)
    write_elem(g,'varm',dict(measured_in_source=measured))
    write_elem(g,'varp',{})
    write_elem(g,'obsm',dict(model.obsm));write_elem(g,'obsp',dict(model.obsp))
    uns=dict(model.uns)
    catalog=pd.read_csv(args.run/'tables/subcluster_catalog_refined_private.tsv',sep='\t',keep_default_na=False)
    if set(catalog.subcluster)!=set(obs.subcluster.astype(str).unique()):raise ValueError('Final annotation catalog differs')
    uns['subcluster_catalog_json']=catalog.to_json(orient='records')
    CM_reports=[];CM_memberships={}
    CM_attempts={}
    if args.CM_discovery_index:
        index=json.loads(args.CM_discovery_index.read_text())
        if not index['all_ready_scopes_completed']:raise ValueError('CM discovery index incomplete')
        CM_attempts={p['scope']:p['discovery_attempt'] for p in index['results']}
    for manifest in ['CM_preparation_public.json','CM_material_preparation_public.json']:
        path=args.run/'manifests'/manifest
        if not path.exists():continue
        for scope in json.loads(path.read_text())['scopes']:
            if not scope['ready']:continue
            attempt=CM_attempts.get(scope['scope'],args.CM_discovery_attempt)
            if not __import__('re').fullmatch(r'discovery_v[0-9]+',attempt):raise ValueError('Invalid CM artifact attempt')
            folder=args.run/'CM'/scope['scope']/attempt
            evidence=json.loads((folder/'CM_public.json').read_text())
            if evidence['scope']!=scope['scope']:raise ValueError('Final CM scope differs')
            CM_reports.append(dict(input_scope=scope,discovery=evidence))
            try:
                table=pd.read_csv(folder/'each_node.tsv',sep='\t')
                members=table.astype(object).where(pd.notna(table),None).to_dict(orient='records')
            except pd.errors.EmptyDataError:
                if evidence['n_network_nodes']!=0:raise
                members=[]
            CM_memberships[scope['scope']]=members
    if not CM_reports:raise ValueError('No completed CM evidence to archive')
    uns['CM_scope_reports_json']=json.dumps(CM_reports)
    uns['CM_node_membership_json']=json.dumps(CM_memberships,allow_nan=False)
    uns.update(source_order=np.asarray([s.record['dataset'] for s in sources],dtype=str),
               feature_identity_diagnostics=diagnostics,raw_counts_layer='counts',
               chosen_count_formula='layers/counts - layers/ambient_removed_counts',
               X_semantics='log1p(CP10K);native_complete-library_chosen-count_denominator',
               authority_scope='technical_count_metadata_archive;clinical_eligibility_per_cell;independent_validation_not_frozen',
               independent_validation_frozen=False,all_cells_retained=True,
               CM_discovery_scope='candidate_new_modules;not_fixed_P3_reference_validation',
               source_contract_sha256=hashlib.sha256(args.contract.read_bytes()).hexdigest())
    if 'PCs' in model.varm:
        pcs=np.zeros((len(var),model.varm['PCs'].shape[1]),dtype=np.float32)
        pcs[var.index.get_indexer(model.var_names)]=model.varm['PCs']
        write_elem(g['varm'],'PCs',pcs)
    write_elem(g,'uns',uns);zarr.consolidate_metadata(str(target))
    del model;gc.collect()
    # Independent pass re-reads every output row and every original source count.
    checks=dict(all_counts_exact=True,all_correction_exact=True,all_normalization_exact=True,
                cell_source_order_exact=True,measured_feature_mask_exact=True,
                all_counts_nonnegative_integer=True,all_X_finite=True,all_coordinates_finite=True)
    out=zarr.open_group(str(target),mode='r');outobs=read_elem(out['obs']);outvar=read_elem(out['var'])
    checks['cell_source_order_exact']=outobs.index.equals(obs.index) and outvar.index.equals(var.index)
    checks['measured_feature_mask_exact']=np.array_equal(np.asarray(out['varm/measured_in_source']),measured)
    for name in ['X_pca','X_pca_harmony','X_umap']:
        matrix=out['obsm'][name]
        checks['all_coordinates_finite']&=bool(matrix.shape[0]==len(obs) and all(np.isfinite(matrix[lo:lo+10000]).all() for lo in range(0,len(obs),10000)))
    offset=0;max_error=0.
    for source in sources:
        for lo in range(0,len(source.rows),2048):
            ids=np.arange(lo,min(lo+2048,len(source.rows)));hi=offset+len(ids)
            raw,chosen,total=source.analysis(ids)
            expected_raw=projected(raw.astype(np.int64),source.global_columns,len(var))
            expected_chosen=projected(chosen,source.global_columns,len(var))
            actual=csr_read(out['layers/counts'],offset,hi)
            diff=actual.astype(np.int64)-expected_raw;diff.eliminate_zeros();checks['all_counts_exact']&=diff.nnz==0
            checks['all_counts_nonnegative_integer']&=bool(np.isfinite(actual.data).all() and (actual.data>=0).all() and np.array_equal(actual.data,np.rint(actual.data)))
            removal=csr_read(out['layers/ambient_removed_counts'],offset,hi)
            diff=(actual.astype(np.float64)-removal-expected_chosen).tocsr();diff.eliminate_zeros()
            checks['all_correction_exact']&=bool(not diff.nnz or abs(diff.data).max()<=1e-7)
            x=csr_read(out['X'],offset,hi);expected_x=normalized(expected_chosen,total)
            diff=x-expected_x;diff.eliminate_zeros();error=float(abs(diff.data).max()) if diff.nnz else 0.
            max_error=max(max_error,error);checks['all_normalization_exact']&=error<=1e-6
            checks['all_X_finite']&=bool(np.isfinite(x.data).all())
            offset=hi
        print(json.dumps(dict(phase='independent_fullgene_source_audit',dataset=source.record['dataset'],n_rows_audited=offset)),flush=True)
    if offset!=len(obs) or not all(checks.values()):raise RuntimeError('Independent fullgene audit failed: '+json.dumps(checks))
    report=dict(observed_utc=datetime.now(timezone.utc).isoformat(),n_cells=len(obs),n_features=len(var),
                n_sources=len(sources),all_checks_pass=True,checks={k:bool(v) for k,v in checks.items()},
                max_X_absolute_error=max_error,full_native_counts_and_all_cells=True,
                annotation_major_counts=obs.major.astype(str).value_counts().to_dict(),
                malignancy_status_counts=obs.malignant_status.astype(str).value_counts().to_dict(),
                independent_validation_frozen=False,clinical_scope='per_cell_native_endpoint_eligibility')
    save_json_new(args.run/'manifests/independent_immune_fullgene_audit_public.json',report)
    allocated=sum(p.stat().st_blocks*512 for p in target.rglob('*') if p.is_file())
    free=capacity(args.run,args.quota_group)
    budget=dict(observed_utc=datetime.now(timezone.utc).isoformat(),actual_candidate_allocated_bytes=allocated,
                final_effective_free_bytes=free,reserve_bytes=RESERVE,additional_copy_allowance_bytes=1000000000,
                fits=free>=allocated+RESERVE+1000000000)
    save_json_new(args.run/'manifests/immune_final_storage_budget_public.json',budget)
    if not budget['fits']:raise RuntimeError('Final result exceeds effective quota with agreed reserve; audited temporary candidate preserved')
    destination=args.run/'outputs/immune_validation_authority_v1.zarr'
    destination.parent.mkdir(exist_ok=True,mode=0o700)
    if destination.exists():raise ValueError('Existing final authority preserved')
    shutil.copytree(target,destination)
    files=[]
    for src in sorted(target.rglob('*')):
        if not src.is_file():continue
        relative=src.relative_to(target);dst=destination/relative
        a=hashlib.sha256(src.read_bytes()).hexdigest();b=hashlib.sha256(dst.read_bytes()).hexdigest()
        if a!=b:raise ValueError('Final copy hash mismatch')
        files.append(dict(relative_path=str(relative),sha256=a,size=src.stat().st_size))
    save_json_new(args.run/'manifests/final_Zarr_file_checksums_private.json',dict(files=files))
    final=zarr.open_group(str(destination),mode='r')
    if final['X'].attrs['shape']!=[len(obs),len(var)] or not read_elem(final['obs']).index.equals(obs.index):raise ValueError('Final Zarr readback differs')
    save_json_new(args.run/'manifests/immune_authority_release_public.json',dict(report,observed_utc=datetime.now(timezone.utc).isoformat(),
       technical_authority_released=True,final_copy_all_file_checksums_verified=True,
       final_storage='original_project_volume',temporary_processing_volume_is_not_authority=True,
       final_allocated_bytes=sum(p.stat().st_blocks*512 for p in destination.rglob('*') if p.is_file()),
       final_effective_free_bytes=capacity(args.run,args.quota_group),source_final_AC=False,
       clinical_eligibility_is_source_specific=True,malignancy_gold_standard=False,
       new_CM_is_discovery_not_external_fixed_reference_validation=True))
    print(json.dumps(dict(phase='immune_technical_authority_release',n_cells=len(obs),n_features=len(var),all_checks_pass=True)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['contract','run','scratch']:p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--quota-group',required=True)
    p.add_argument('--CM-discovery-attempt',default='discovery_v1')
    p.add_argument('--CM-discovery-index',type=Path)
    p.add_argument('--CNV-results',type=Path)
    pack(p.parse_args())
