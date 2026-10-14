#!/usr/bin/env python3
"""Per-original-library native-count CNV with conservative coverage semantics."""
import argparse
from datetime import datetime,timezone
import cProfile,gc,hashlib,importlib.util,logging,json,sys,time,traceback
from pathlib import Path
import anndata as ad
from anndata._io.specs import read_elem
import numpy as np
import pandas as pd
from scipy import sparse
import zarr
from immune_atlas_io import load_contract,save_json_new


def native_cnv(args):
    sys.path.insert(0,str(args.core_dir))
    spec=importlib.util.spec_from_file_location('immune_corrected_malignant_core',args.corrected_core)
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    import scanpy as sc
    sc.settings.n_jobs=4
    logger=logging.getLogger('immune_native_CNV');logger.setLevel(logging.INFO)
    logger.addHandler(logging.StreamHandler())
    config,sources,var,identities,diagnostics=load_contract(args.contract,['core','CAR_product'])
    model=zarr.open_group(str(args.scratch/'immune_refined.zarr'),mode='r')
    obs=read_elem(model['obs'])
    coords=pd.read_csv(args.gene_annotation,sep=r'\s+',header=None,names=['gene_symbol','chromosome','start','end'])
    coords=coords.loc[coords.chromosome.str.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2])')].drop_duplicates('gene_symbol')
    coord_path=args.run/'config'/('CNV_autosomal_coordinate_contract_'+args.attempt_version+'.tsv')
    if not coord_path.exists():raise ValueError('Coordinator must freeze coordinate contract before workers')
    frozen=pd.read_csv(coord_path,sep='\t',header=None,names=list(coords.columns))
    if not frozen.reset_index(drop=True).equals(coords.reset_index(drop=True)):raise ValueError('Coordinate contract differs')
    output=[];summary_rows=[];library_number=0;profiled=False
    for source in sources:
        source_obs=obs.loc[obs.dataset.astype(str).eq(source.record['dataset'])]
        expected=source.record['dataset']+'::'+source.obs.index[source.rows].astype(str)
        if not np.array_equal(source_obs.index.astype(str),expected):raise ValueError('CNV raw-row identity mismatch')
        symbols=np.asarray(source.symbols,dtype=str)
        indexed=pd.Index(symbols).isin(coords.gene_symbol)
        coordinate_set=set(coords.gene_symbol.tolist())
        direct=np.fromiter((g in coordinate_set for g in symbols),dtype=bool,count=len(symbols))
        if not np.array_equal(indexed,direct):raise ValueError('Coordinate membership semantics differ')
        keep=(~pd.Index(symbols).duplicated(keep=False))&indexed
        for library,positions in source_obs.groupby('library_id',observed=True,sort=False).indices.items():
            assigned=library_number%args.shard_count==args.shard_index;library_number+=1
            if not assigned:continue
            positions=np.asarray(positions);frame=source_obs.iloc[positions].copy()
            token=source.record['dataset']+'_'+hashlib.sha256(str(library).encode()).hexdigest()[:16]
            out=args.scratch/('CNV_'+args.attempt_version)/token;out.mkdir(parents=True,exist_ok=False,mode=0o700)
            if args.reuse_attempt_version:
                old=args.scratch/('CNV_'+args.reuse_attempt_version)/token
                method_path=old/'method_summary_private.json';cells_path=old/'cell_CNV_results_private.tsv.gz'
                if method_path.exists() and cells_path.exists():
                    old_method=json.loads(method_path.read_text())
                    accepted=old_method['status']=='success' and old_method.get('integer_counts_preserved') and old_method.get('all_expression_entries_independently_checked') and old_method.get('normalization_max_absolute_error',1.)<=1e-6 and old_method.get('CNV_expression')=='sample_native_raw_counts_log1p_CP10K_full_native_denominator'
                    if accepted:
                        result=pd.read_csv(cells_path,sep='\t',index_col=0,keep_default_na=False)
                        if not result.index.equals(frame.index) or not np.array_equal(result.sampleID.astype(str),frame.sampleID.astype(str)):raise ValueError('Completed correct-expression reuse mapping differs')
                        def truth(v):return v.astype(str).str.lower().isin(['true','1'])
                        row=dict(dataset=source.record['dataset'],library_token=token,status='success',n_cells=len(result),
                            n_candidates=int(truth(result.tumor_lineage_candidate).sum()),n_references=int(truth(result.tumor_reference_cell).sum()),
                            n_evaluable=int(truth(result.cnv_evaluable).sum()),n_CNV_supported_malignant=int(result.malignant_status.eq('CNV_supported_malignant').sum()),detail='completed_correct_log_expression_reused')
                        save_json_new(out/'reused_correct_CNV_private.json',dict(original_result=str(cells_path),
                            original_result_sha256=hashlib.sha256(cells_path.read_bytes()).hexdigest(),original_method_sha256=hashlib.sha256(method_path.read_bytes()).hexdigest(),
                            original_core_patch_sha256=hashlib.sha256((args.run/'manifests/CNV_logspace_core_patch_public.json').read_bytes()).hexdigest(),
                            native_gene_membership_independently_equal=True,cell_and_original_library_IDs_exact=True))
                        summary_rows.append(row);output.append(result)
                        print(json.dumps({k:v for k,v in row.items() if k not in ['library_token','detail']}),flush=True)
                        continue
            result=pd.DataFrame(index=frame.index)
            result['sampleID']=frame.sampleID.astype(str)
            result['malignant_status']='not_evaluated'
            result['cnv_score']=np.nan
            result['cnv_evaluable']=False
            result['cnv_cluster']=''
            result['infercnv_positive']=False
            result['tumor_lineage_candidate']=False
            result['tumor_reference_cell']=False
            result['CNV_method_status']='pending'
            result['CNV_native_total_counts']=np.nan
            method=dict(status='skipped',detail='capture_scope_not_tumor_CNVevaluable',candidate_cells=0,reference_cells=0)
            # Products and T-enriched samples provide no valid tumor-cell census.
            if source.record['role']!='CAR_product' and not frame.capture_scope.astype(str).str.contains('T_enriched|lymphocyte_enriched|CD45_enriched|ex_vivo').all():
                major=frame.major.astype(str)
                author=frame.source_author_label.astype(str).str.lower()
                author_tumor=author.str.contains(r'tumou?r|malignan|cancer|epithelial',regex=True)
                residual=frame.annotation_exclude_from_primary.astype(str).str.lower().isin(['true','1'])
                brain_neural=source.record['dataset']=='GSE292046'
                neural_candidates=major.isin(['Neuron','Oligodendrocyte']) if brain_neural else np.zeros(len(frame),dtype=bool)
                candidate=(major.isin(['Epithelial','Glial','Melanocytic','Unresolved'])|neural_candidates|author_tumor)&~residual
                reference=major.isin(['T','NK','B','Plasma','Myeloid','Mast','Neutrophil','Endothelial','Fibroblast','Pericyte_SMC'])&~author_tumor&~candidate&~residual
                result['tumor_lineage_candidate']=candidate.to_numpy();result['tumor_reference_cell']=reference.to_numpy()
                result['CNV_candidate_basis']=np.where(major.eq('Unresolved'),'unresolved_in_tumor_capture_not_identity_assignment',
                    np.where(neural_candidates,'neural_glial_state_in_GBM_not_automatically_normal','marker_or_author_tumor_lineage_candidate'))
                profile=None
                try:
                    started=time.perf_counter()
                    counts=source.raw(positions)
                    total=np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
                    if not np.isfinite(total).all() or (total<=0).any():raise ValueError('Invalid native library RNA denominators')
                    result['CNV_native_total_counts']=total
                    if keep.sum()<1000:raise ValueError('Insufficient unambiguous autosomal gene coordinate coverage')
                    local=ad.AnnData(sparse.csr_matrix(counts[:,keep],dtype=np.int32),
                                     obs=frame[['major','subcluster']].copy(),var=pd.DataFrame(index=pd.Index(symbols[keep],name='gene_symbol')))
                    local.layers['counts']=local.X.copy()
                    expression=(sparse.diags(10000./total)@local.layers['counts']).tocsr()
                    expression.data=np.log1p(expression.data)
                    local.layers['cnv_log1p_CP10K']=expression.astype(np.float32)
                    expected=local.layers['counts'].astype(np.float64).multiply((10000./total)[:,None]).tocsr()
                    expected.data=np.log1p(expected.data)
                    difference=local.layers['cnv_log1p_CP10K'].astype(np.float64)-expected
                    difference.eliminate_zeros()
                    max_error=float(np.abs(difference.data).max()) if difference.nnz else 0.
                    if max_error>1e-6:raise ValueError('Sample log normalization independent arithmetic differs')
                    del expected,difference
                    local.obs['CNV_native_total_counts']=total
                    if not np.isfinite(expression.data).all() or (expression.data<0).any():raise ValueError('Invalid sample log expression')
                    local.obs['CNV_role']=np.where(candidate,'candidate',np.where(reference,'reference','excluded'))
                    local.obs['batch']=str(library)
                    cfg=dict(candidate_category_values=['candidate'],reference_category_values=['reference'],
                             candidate_label_terms=[],reference_label_terms=[],exclude_label_terms=[],
                             min_candidate_cells=50,min_reference_cells=50,window_size=250,
                             reference_quantile=.99,random_seed=20261013,
                             copykat={'enabled':False},cluster_min_infercnv_fraction=.60,
                             cluster_min_cells=30,cluster_mean_burden_multiplier=1.05,
                             infercnv_expression_layer='cnv_log1p_CP10K',infercnv_n_jobs=2)
                    raw_digest=hashlib.sha256(local.layers['counts'].data.tobytes()).hexdigest()
                    profile=cProfile.Profile() if not profiled else None
                    if profile:profile.enable()
                    method=core.run_malignant_inference(local,label_col='subcluster',category_col='CNV_role',
                           table_dir=out/'tables',figure_dir=out/'figures',logger=logger,
                           gene_annotation_path=coord_path,cfg=cfg,batch_col='batch',run_label=token)
                    if profile:
                        profile.disable();profile.dump_stats(out/'core_profile_private.pstats');profiled=True
                    if hashlib.sha256(local.layers['counts'].data.tobytes()).hexdigest()!=raw_digest:raise ValueError('CNV core altered raw counts')
                    method.update(CNV_expression='sample_native_raw_counts_log1p_CP10K_full_native_denominator',
                        integer_counts_preserved=True,normalization_max_absolute_error=max_error,
                        all_expression_entries_independently_checked=True,infercnv_process_workers=2,elapsed_seconds=time.perf_counter()-started)
                    if method['status']=='success':
                        score=pd.to_numeric(local.obs.cnv_score,errors='coerce')
                        if not np.isfinite(score.to_numpy()).all() or score.max()<=0:
                            raise ValueError('CNV returned nonfinite or completely uninformative signal')
                        evaluated=candidate&np.isfinite(score.to_numpy())
                        result.loc[evaluated,'cnv_evaluable']=True
                        result.loc[evaluated,'cnv_score']=score.loc[evaluated]
                        labels=local.obs.loc[evaluated,'malignant_label'].astype(str)
                        result.loc[evaluated,'cnv_cluster']=local.obs.loc[evaluated,'cnv_cluster'].astype(str)
                        result.loc[evaluated,'infercnv_positive']=local.obs.loc[evaluated,'infercnv_positive'].astype(bool)
                        # A single-method negative is not proof of normal identity.
                        result.loc[evaluated,'malignant_status']=labels.replace({'malignant':'CNV_supported_malignant',
                            'uncertain':'CNV_uncertain','non_malignant':'CNV_negative_candidate'})
                        result.loc[reference,'malignant_status']='lineage_reference_not_tumor_candidate'
                        result.loc[reference,'cnv_score']=score.loc[reference]
                    del local,counts;gc.collect()
                except Exception as exc:
                    if profile:profile.disable();profile.dump_stats(out/'failed_core_profile_private.pstats')
                    method=dict(status='failed',detail=str(exc),candidate_cells=int(candidate.sum()),reference_cells=int(reference.sum()))
                    (out/'traceback_private.txt').write_text(traceback.format_exc())
            result['CNV_method_status']=method['status']
            result['CNV_input']='sample_native_raw_counts_to_log1p_CP10K_autosomal_genes_same_original_library'
            result['CNV_evidence_tier']='single_infercnvpy_cluster_support_not_dual_method_gold'
            # Guard the historical core default: never retain a skipped normal/0 call.
            if method['status']!='success':
                result['malignant_status']='not_evaluated_'+method['status']
                result['cnv_score']=np.nan;result['cnv_evaluable']=False
            result.to_csv(out/'cell_CNV_results_private.tsv.gz',sep='\t')
            save_json_new(out/'method_summary_private.json',method)
            row=dict(dataset=source.record['dataset'],library_token=token,status=method['status'],
                     n_cells=len(result),n_candidates=int(result.tumor_lineage_candidate.sum()),
                     n_references=int(result.tumor_reference_cell.sum()),
                     n_evaluable=int(result.cnv_evaluable.sum()),n_CNV_supported_malignant=int(result.malignant_status.eq('CNV_supported_malignant').sum()),
                     detail=method.get('detail',''))
            summary_rows.append(row);output.append(result)
            print(json.dumps({k:v for k,v in row.items() if k not in ['library_token','detail']}),flush=True)
    merged=pd.concat(output)
    expected=obs.index[obs.index.isin(merged.index)];merged=merged.reindex(expected)
    if not merged.index.is_unique or merged.sampleID.isna().any():raise ValueError('CNV shard cell mapping incomplete')
    suffix=args.attempt_version+'_shard_'+str(args.shard_index)
    destination=args.run/'tables'/('all_core_CNV_'+suffix+'_private.tsv.gz')
    if destination.exists():raise ValueError('Existing shard preserved')
    merged.to_csv(destination,sep='\t')
    pd.DataFrame(summary_rows).to_csv(args.run/'tables'/('CNV_library_method_index_'+suffix+'_private.tsv'),sep='\t',index=False)
    public=dict(observed_utc=datetime.now(timezone.utc).isoformat(),n_cells=len(merged),n_libraries=len(summary_rows),
        shard_index=args.shard_index,shard_count=args.shard_count,attempt_version=args.attempt_version,
        library_method_status_counts=pd.Series([r['status'] for r in summary_rows]).value_counts().to_dict(),
        cell_malignant_status_counts=merged.malignant_status.value_counts().to_dict(),
        original_gene_annotation_sha256=hashlib.sha256(args.gene_annotation.read_bytes()).hexdigest(),
        core_sha256=hashlib.sha256(args.corrected_core.read_bytes()).hexdigest(),
        raw_count_input=True,same_original_library_reference=True,skips_never_normal=True,
        inference_input='sample_native_raw_counts_to_log1p_CP10K',raw_counts_preserved=True,
        CNV_evidence_tier='single_infercnvpy_cluster_support',n_method_failed=sum(r['status']=='failed' for r in summary_rows),
        authoritative_release=False)
    save_json_new(args.run/'manifests'/('CNV_'+suffix+'_public.json'),public)
    if public['n_method_failed']:raise RuntimeError('CNV method failures require review before final release')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['contract','run','scratch','core-dir','gene-annotation']:p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--corrected-core',type=Path,required=True)
    p.add_argument('--attempt-version',default='logspace_v1')
    p.add_argument('--shard-index',type=int,default=0);p.add_argument('--shard-count',type=int,default=1)
    p.add_argument('--reuse-attempt-version')
    args=p.parse_args()
    if args.shard_count<1 or not 0<=args.shard_index<args.shard_count:raise ValueError('Invalid original-library shard')
    native_cnv(args)
