#!/usr/bin/env python3
"""Full-native-gene descriptive expression contrasts for every new subcluster.

No p values or DEG claims are manufactured. Sparse expression statistics use
native denominators and source feature masks; unmeasured genes are missing.
"""
import argparse
from datetime import datetime,timezone
import gc,json,re
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from immune_atlas_io import load_contract,normalized,projected,save_json_new


def allowed(gene):
    return bool(gene and gene not in {'MALAT1','NEAT1','ACTB','ACTG1','B2M','FOS','FOSB','JUN','JUNB','JUND','EGR1','HLA-A','HLA-B','HLA-C'}
                and not re.match(r'^(MT-|RPL|RPS|HSP|ENSG|ERCC-|IG[HKL][VJCD])',gene))


def refine(args):
    model=ad.read_zarr(args.scratch/'immune_annotated.zarr')
    config,sources,var,identities,diagnostics=load_contract(args.contract,['core','CAR_product'])
    catalog=pd.read_csv(args.run/'tables/subcluster_catalog_private.tsv',sep='\t',dtype={'node_id':str})
    nodes=catalog.subcluster.astype(str).tolist();node_lookup={v:i for i,v in enumerate(nodes)}
    labels=np.array([node_lookup[x] for x in model.obs.subcluster.astype(str)],dtype=np.int32)
    n=len(nodes);width=len(var)
    sums=np.zeros((n,width),dtype=np.float64);positive=np.zeros((n,width),dtype=np.int64)
    measured=np.zeros((n,width),dtype=np.int64);study_stats=[];offset=0
    for source in sources:
        local_obs=model.obs.loc[model.obs.dataset.astype(str).eq(source.record['dataset'])]
        expected=source.record['dataset']+'::'+source.obs.index[source.rows].astype(str)
        if not np.array_equal(local_obs.index.astype(str),expected):raise ValueError('Native statistics order mismatch')
        groups=labels[offset:offset+len(local_obs)]
        local_sum=np.zeros((n,len(source.var)),dtype=np.float64)
        local_positive=np.zeros((n,len(source.var)),dtype=np.int64)
        for lo in range(0,len(local_obs),2048):
            hi=min(lo+2048,len(local_obs));raw,chosen,total=source.analysis(np.arange(lo,hi))
            cp=(sparse.diags(10000./total)@chosen).tocsr()
            membership=sparse.csr_matrix((np.ones(hi-lo),(groups[lo:hi],np.arange(hi-lo))),shape=(n,hi-lo))
            local_sum+=(membership@cp).toarray()
            mask=chosen.copy();mask.data=(mask.data>0).astype(np.int32);mask.eliminate_zeros()
            local_positive+=(membership@mask).toarray().astype(np.int64)
            canonical_mask=projected(chosen,source.global_columns,width)
            canonical_mask.data=(canonical_mask.data>0).astype(np.int32);canonical_mask.eliminate_zeros()
            positive+=(membership@canonical_mask).toarray().astype(np.int64)
        # Sum aliases first for means; positive fractions are bounded later and
        # ambiguous/repeated native identities are excluded from marker naming.
        sums+=projected(sparse.csr_matrix(local_sum),source.global_columns,width).toarray()
        counts=np.bincount(groups,minlength=n)
        cols=np.unique(source.global_columns)
        measured[:,cols]+=counts[:,None]
        study_stats.append((source.record['dataset'],source.global_columns.copy(),counts,local_sum,local_positive))
        offset+=len(local_obs)
        print(json.dumps(dict(phase='fullgene_source_statistics',dataset=source.record['dataset'],n_cells=len(local_obs),n_genes=len(source.var))),flush=True)
    if offset!=model.n_obs:raise ValueError('All-cell statistics coverage mismatch')
    # Alias-positive summation cannot be used as detection when aliases duplicate.
    repeated=set()
    for source in sources:
        counts=pd.Series(source.canonical_genes).value_counts()
        repeated.update(counts.index[counts>1])
    symbols=var.gene_symbol.astype(str)
    unique_symbol=~symbols.duplicated(keep=False).to_numpy()
    usable=np.array([allowed(x) for x in symbols])&unique_symbol&~var.index.isin(repeated)&var.feature_identity_status.ne('ambiguous_source_feature').to_numpy()
    rows=[];mapping={};full_evidence=[]
    major_values=catalog.major.astype(str).to_numpy()
    major_totals={m:(sums[major_values==m].sum(axis=0),positive[major_values==m].sum(axis=0),measured[major_values==m].sum(axis=0)) for m in np.unique(major_values)}
    for i,node in enumerate(nodes):
        major=catalog.iloc[i].major
        total_sum,total_positive,total_measured=major_totals[major]
        outside_n=total_measured-measured[i]
        pct=np.divide(positive[i],measured[i],out=np.full(width,np.nan),where=measured[i]>0)
        mean=np.divide(sums[i],measured[i],out=np.full(width,np.nan),where=measured[i]>0)
        outside_pct=np.divide(total_positive-positive[i],outside_n,out=np.full(width,np.nan),where=outside_n>0)
        outside_mean=np.divide(total_sum-sums[i],outside_n,out=np.full(width,np.nan),where=outside_n>0)
        log_ratio=np.log2((mean+0.01)/(outside_mean+0.01));delta=pct-outside_pct
        valid=usable&(measured[i]>=30)&(pct>=0.2)
        specific=valid&(outside_n>=30)&(delta>=0.1)&(log_ratio>=0.5)
        candidates=np.flatnonzero(specific)
        role='own_group_fullnative_expression_contrast'
        if not len(candidates):candidates=np.flatnonzero(valid);role='own_group_fullnative_representative_non_specific'
        ranking=sorted(candidates,key=lambda j:(-float(delta[j]) if np.isfinite(delta[j]) else 0.,-float(log_ratio[j]) if np.isfinite(log_ratio[j]) else 0.,-float(mean[j]),str(symbols.iloc[j])))
        shortlist=ranking[:20]
        support={j:0 for j in shortlist};coverage={j:0 for j in shortlist}
        for name,columns,counts,local_sum,local_positive in study_stats:
            if counts[i]<30:continue
            lookup={int(v):k for k,v in enumerate(columns)}
            others=(major_values==major)&(np.arange(n)!=i)
            out_count=counts[others].sum()
            for j in shortlist:
                if j not in lookup:continue
                k=lookup[j];coverage[j]+=1
                local_pct=local_positive[i,k]/counts[i]
                local_out=local_positive[others,k].sum()/out_count if out_count>=30 else np.nan
                local_mean=local_sum[i,k]/counts[i]
                local_outmean=local_sum[others,k].sum()/out_count if out_count>=30 else np.nan
                if local_pct>=.2 and out_count>=30 and local_pct-local_out>=.1 and np.log2((local_mean+.01)/(local_outmean+.01))>=.5:support[j]+=1
        if shortlist:
            best=sorted(shortlist,key=lambda j:(-support[j],ranking.index(j)))[0]
            marker=str(symbols.iloc[best]);identity=str(var.index[best])
        else:best=None;marker='NA';identity=''
        display=f'{major}_{marker}_{int(catalog.iloc[i].number)}';mapping[node]=display
        row=catalog.iloc[i].to_dict();row.update(provisional_subcluster=node,subcluster=display,marker=marker,
             marker_gene_identity=identity,marker_role=role if best is not None else 'no_valid_native_marker',
             marker_positive_fraction=float(pct[best]) if best is not None else None,
             marker_mean_CP10K=float(mean[best]) if best is not None else None,
             marker_source_coverage=coverage.get(best,0),marker_source_contrast_support=support.get(best,0),
             marker_origin='no_local_DEG_exists_for_new_clusters;full_native_group_expression',
             inference_scope='descriptive_no_DEG_pvalue',identity_promoted=False)
        rows.append(row)
        for j in shortlist:
            full_evidence.append(dict(node_id=row['node_id'],gene=str(symbols.iloc[j]),gene_identity=str(var.index[j]),
                positive_fraction=float(pct[j]),mean_CP10K=float(mean[j]),
                within_major_delta_pct=float(delta[j]) if np.isfinite(delta[j]) else None,
                within_major_log2ratio=float(log_ratio[j]) if np.isfinite(log_ratio[j]) else None,
                n_measured_cells=int(measured[i,j]),n_measured_sources=coverage[j],n_source_contrasts_supported=support[j]))
    model.obs['subcluster_provisional']=model.obs.subcluster.astype(str)
    model.obs['subcluster']=pd.Categorical(model.obs.subcluster_provisional.map(mapping))
    node_map={r['subcluster']:r['node_id'] for r in rows}
    model.obs['node_id']=pd.Categorical(model.obs.subcluster.astype(str).map(node_map))
    model.obs['majorCluster']=model.obs.major
    model.obs['subCluster']=model.obs.subcluster
    model.uns['marker_naming_method']='own_fullnative_expression_contrast_or_own_representative;one_symbol_number;SHA_separate;no_DEG_significance_claim'
    model.uns['native_fullgene_statistics_scope']=dict(n_cells=model.n_obs,n_features=width,all_cells_processed=True)
    out=args.scratch/'immune_refined.zarr'
    if out.exists():raise ValueError('Existing marker refinement preserved')
    model.write_zarr(out)
    pd.DataFrame(rows).to_csv(args.run/'tables/subcluster_catalog_refined_private.tsv',sep='\t',index=False)
    pd.DataFrame(full_evidence).to_csv(args.run/'tables/fullnative_marker_contrasts_private.tsv',sep='\t',index=False)
    stats=args.scratch/'fullgene_statistics_v1';stats.mkdir(mode=0o700)
    np.save(stats/'sum_CP10K.npy',sums);np.save(stats/'positive_counts.npy',positive);np.save(stats/'measured_counts.npy',measured)
    var.to_csv(stats/'genes_private.tsv',sep='\t')
    pd.DataFrame(rows).to_csv(stats/'nodes_private.tsv',sep='\t',index=False)
    save_json_new(args.run/'manifests/marker_refinement_public.json',dict(observed_utc=datetime.now(timezone.utc).isoformat(),
        n_cells=model.n_obs,n_nodes=n,n_native_union_genes=width,all_cell_native_expression_statistics=True,
        marker_origin_counts=pd.Series([r['marker_role'] for r in rows]).value_counts().to_dict(),
        one_marker_per_display=True,SHA_in_separate_node_id=True,formal_DEG_tests=False,authoritative_release=False))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['contract','run','scratch']:p.add_argument('--'+key,type=Path,required=True)
    refine(p.parse_args())
