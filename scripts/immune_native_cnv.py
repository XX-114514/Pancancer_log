#!/usr/bin/env python3
"""Per-original-library native-count CNV with conservative coverage semantics."""
import argparse
from datetime import datetime,timezone
import gc,hashlib,importlib,logging,json,sys,traceback
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
    core=importlib.import_module('malignant_inference')
    import scanpy as sc
    sc.settings.n_jobs=4
    logger=logging.getLogger('immune_native_CNV');logger.setLevel(logging.INFO)
    logger.addHandler(logging.StreamHandler())
    config,sources,var,identities,diagnostics=load_contract(args.contract,['core','CAR_product'])
    model=zarr.open_group(str(args.scratch/'immune_refined.zarr'),mode='r')
    obs=read_elem(model['obs'])
    coords=pd.read_csv(args.gene_annotation,sep=r'\s+',header=None,names=['gene_symbol','chromosome','start','end'])
    coords=coords.loc[coords.chromosome.str.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2])')].drop_duplicates('gene_symbol')
    coord_path=args.run/'config/CNV_autosomal_coordinate_contract_v1.tsv'
    if coord_path.exists():raise ValueError('Existing coordinate attempt preserved')
    coords.to_csv(coord_path,sep='\t',index=False,header=False)
    output=[];summary_rows=[]
    for source in sources:
        source_obs=obs.loc[obs.dataset.astype(str).eq(source.record['dataset'])]
        expected=source.record['dataset']+'::'+source.obs.index[source.rows].astype(str)
        if not np.array_equal(source_obs.index.astype(str),expected):raise ValueError('CNV raw-row identity mismatch')
        for library,positions in source_obs.groupby('library_id',observed=True,sort=False).indices.items():
            positions=np.asarray(positions);frame=source_obs.iloc[positions].copy()
            token=source.record['dataset']+'_'+hashlib.sha256(str(library).encode()).hexdigest()[:16]
            out=args.scratch/'CNV'/token;out.mkdir(parents=True,exist_ok=False,mode=0o700)
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
            method=dict(status='skipped',detail='capture_scope_not_tumor_CNVevaluable',candidate_cells=0,reference_cells=0)
            # Products and T-enriched samples provide no valid tumor-cell census.
            if source.record['role']!='CAR_product' and not frame.capture_scope.astype(str).str.contains('T_enriched|lymphocyte_enriched|CD45_enriched').all():
                major=frame.major.astype(str)
                author=frame.source_author_label.astype(str).str.lower()
                author_tumor=author.str.contains(r'tumou?r|malignan|cancer|epithelial',regex=True)
                residual=frame.annotation_exclude_from_primary.astype(str).str.lower().isin(['true','1'])
                candidate=(major.isin(['Epithelial','Glial','Melanocytic'])|author_tumor)&~residual
                reference=major.isin(['T','NK','B','Plasma','Myeloid','Mast','Neutrophil','Endothelial','Fibroblast','Pericyte_SMC'])&~author_tumor&~candidate&~residual
                result['tumor_lineage_candidate']=candidate.to_numpy();result['tumor_reference_cell']=reference.to_numpy()
                try:
                    counts=source.raw(positions)
                    symbols=np.asarray(source.symbols,dtype=str)
                    keep=(~pd.Index(symbols).duplicated(keep=False))&np.isin(symbols,coords.gene_symbol)
                    if keep.sum()<1000:raise ValueError('Insufficient unambiguous autosomal gene coordinate coverage')
                    local=ad.AnnData(sparse.csr_matrix(counts[:,keep],dtype=np.int32),
                                     obs=frame[['major','subcluster']].copy(),var=pd.DataFrame(index=pd.Index(symbols[keep],name='gene_symbol')))
                    local.layers['counts']=local.X.copy()
                    local.obs['CNV_role']=np.where(candidate,'candidate',np.where(reference,'reference','excluded'))
                    local.obs['batch']=str(library)
                    cfg=dict(candidate_category_values=['candidate'],reference_category_values=['reference'],
                             candidate_label_terms=[],reference_label_terms=[],exclude_label_terms=[],
                             min_candidate_cells=50,min_reference_cells=50,window_size=250,
                             reference_quantile=.99,random_seed=20261013,
                             copykat={'enabled':False},cluster_min_infercnv_fraction=.60,
                             cluster_min_cells=30,cluster_mean_burden_multiplier=1.05)
                    method=core.run_malignant_inference(local,label_col='subcluster',category_col='CNV_role',
                           table_dir=out/'tables',figure_dir=out/'figures',logger=logger,
                           gene_annotation_path=coord_path,cfg=cfg,batch_col='batch',run_label=token)
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
                    method=dict(status='failed',detail=str(exc),candidate_cells=int(candidate.sum()),reference_cells=int(reference.sum()))
                    (out/'traceback_private.txt').write_text(traceback.format_exc())
            result['CNV_method_status']=method['status']
            result['CNV_input']='native_integer_counts_autosomal_genes_same_original_library'
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
    merged=pd.concat(output).reindex(obs.index)
    if not merged.index.is_unique or merged.sampleID.isna().any() or len(merged)!=len(obs):raise ValueError('CNV full-cell mapping incomplete')
    merged.to_csv(args.run/'tables/all_core_CNV_private.tsv.gz',sep='\t')
    pd.DataFrame(summary_rows).to_csv(args.run/'tables/CNV_library_method_index_private.tsv',sep='\t',index=False)
    public=dict(observed_utc=datetime.now(timezone.utc).isoformat(),n_cells=len(obs),n_libraries=len(summary_rows),
        library_method_status_counts=pd.Series([r['status'] for r in summary_rows]).value_counts().to_dict(),
        cell_malignant_status_counts=merged.malignant_status.value_counts().to_dict(),
        original_gene_annotation_sha256=hashlib.sha256(args.gene_annotation.read_bytes()).hexdigest(),
        core_sha256=hashlib.sha256((args.core_dir/'malignant_inference.py').read_bytes()).hexdigest(),
        raw_count_input=True,same_original_library_reference=True,skips_never_normal=True,
        CNV_evidence_tier='single_infercnvpy_cluster_support',n_method_failed=sum(r['status']=='failed' for r in summary_rows),
        authoritative_release=False)
    save_json_new(args.run/'manifests/CNV_public.json',public)
    if public['n_method_failed']:raise RuntimeError('CNV method failures require review before final release')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['contract','run','scratch','core-dir','gene-annotation']:p.add_argument('--'+key,type=Path,required=True)
    native_cnv(p.parse_args())
