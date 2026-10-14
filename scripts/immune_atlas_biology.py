#!/usr/bin/env python3
"""Native RNA, full-cell integration and conservative immune-atlas annotations.

PCA fitting uses balanced landmarks; all eligible cells retain native marker
measurements and receive PCA, Harmony, neighborhood and cluster coordinates.
The modeling object is temporary and never substitutes for the final fullgene
count archive. Patient outcomes are not used for clustering or annotation.
"""
import argparse
from datetime import datetime, timezone
import gc
import hashlib
import json
import os
from pathlib import Path
import time
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
import harmonypy

from immune_atlas_io import load_contract, read_registry, normalized, projected, save_json_new

SEED=20261013


def event(phase, **values):
    print(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),phase=phase,**values)),flush=True)


def panels(path):
    major=json.loads((path/'major_markers.json').read_text())
    major['NK']=[g for g in major['NK'] if g!='TRDC']
    major.update(Neutrophil=['CSF3R','FCGR3B','FPR1','CXCR2','S100A8','S100A9'],
                 Glial=['GFAP','AQP4','SLC1A3','OLIG1','OLIG2','SOX2'],
                 Oligodendrocyte=['MBP','PLP1','MOG','MAG','MOBP'],
                 Neuron=['SNAP25','SYT1','RBFOX3','STMN2','SLC17A7'],
                 Melanocytic=['PMEL','MLANA','TYR','DCT','MITF'],
                 Erythroid=['GYPA','GYPB','SLC4A1','AHSP','HBB','HBA1'])
    subtype=json.loads((path/'lineage_markers.json').read_text())
    # Shared activation and cytotoxic programs are state evidence, not CD4/CD8 truth.
    subtype['T'].pop('CD4_Naive',None)
    subtype['T'].pop('CD8_Naive',None)
    subtype['T'].pop('CD4_Memory',None)
    subtype['T'].pop('Cycling_T',None)
    subtype['T'].pop('CD8_Cytotoxic',None)
    subtype['T'].pop('T_Exhausted',None)
    subtype['T']['Naive_memory_program']=['CCR7','TCF7','LEF1','IL7R','MAL','LTB']
    subtype['NK'].pop('NKT',None)
    subtype['NK']['NK_FCGR3A']=[g for g in subtype['NK']['NK_FCGR3A'] if g!='TRAC']
    for values in subtype.values():
        for label in values:values[label]=list(dict.fromkeys(g for g in values[label] if g not in ['MALAT1','FOS']))
    states=dict(naive_memory=['TCF7','LEF1','CCR7','IL7R','MAL'],
                cytotoxic=['PRF1','GZMB','GZMH','NKG7','CCL5','CST7'],
                dysfunction=['PDCD1','HAVCR2','LAG3','TOX','ENTPD1','CXCL13'],
                interferon=['ISG15','IFIT1','IFIT3','MX1','OAS1','STAT1'],
                cycling=['MKI67','TOP2A','TYMS','UBE2C','STMN1'])
    return major,subtype,states


def capture_marker_indices(var, genes):
    lookup={}
    for gene in genes:
        positions=np.flatnonzero(var.gene_symbol.astype(str).eq(gene).to_numpy())
        if len(positions)==1 and var.iloc[positions[0]].feature_identity_status!='ambiguous_source_feature':
            lookup[gene]=int(positions[0])
    return lookup


def build_model(args):
    config,sources,var,identities,diagnostics=load_contract(args.contract,['core','CAR_product'])
    obs=read_registry(config,sources)
    major,subtype,states=panels(args.panel_dir)
    marker_names=sorted(set(g for genes in major.values() for g in genes)|
                        set(g for lineage in subtype.values() for genes in lineage.values() for g in genes)|
                        set(g for genes in states.values() for g in genes))
    marker_lookup=capture_marker_indices(var,marker_names)
    common=set.intersection(*(set(s.canonical_genes) for s in sources))
    common=sorted(g for g in common if not g.startswith('AMBIGUOUS::'))
    if len(common)<2000:raise ValueError('Shared native feature space is too small for a joint model')
    common_lookup={gene:i for i,gene in enumerate(common)}
    width=len(common)
    marker_rows=sum(min(600,len(s.rows)) for s in sources) if args.phase=='pilot' else len(obs)
    marker_matrix=np.zeros((marker_rows,len(marker_names)),dtype=np.float32)
    marker_available=np.zeros((len(sources),len(marker_names)),dtype=bool)
    expression=[];counts=[];start=0;metrics=[];registry_parts=[]
    for number,source in enumerate(sources):
        local=obs.loc[obs.dataset.eq(source.record['dataset'])].copy()
        selected=np.arange(len(local))
        if args.phase=='pilot':
            selected=np.sort(np.random.default_rng(SEED+number).choice(len(local),min(600,len(local)),replace=False))
            local=local.iloc[selected].copy()
        source_rows=[];source_raw=[]
        common_pos=np.array([common_lookup.get(g,-1) for g in source.canonical_genes])
        valid=common_pos>=0
        global_to_marker={int(marker_lookup[name]):j for j,name in enumerate(marker_names) if name in marker_lookup}
        marker_destination=np.array([global_to_marker.get(int(column),-1) for column in source.global_columns])
        marker_present=marker_destination>=0
        marker_available[number]=np.isin(np.arange(len(marker_names)),marker_destination[marker_present])
        all_totals=[]
        for lo in range(0,len(selected),2048):
            ids=selected[lo:lo+2048]
            raw,chosen,total=source.analysis(ids)
            if raw.nnz and raw.data.max()>np.iinfo(np.int32).max:raise ValueError('Count exceeds int32')
            projected_raw=projected(raw[:,valid].astype(np.int64),common_pos[valid],width)
            if projected_raw.nnz and projected_raw.data.max()>np.iinfo(np.int32).max:raise ValueError('Alias-summed count exceeds int32')
            chosen_common=projected(chosen[:,valid],common_pos[valid],width)
            source_raw.append(projected_raw.astype(np.int32))
            source_rows.append(normalized(chosen_common,total))
            # Marker expression always uses the complete native library denominator.
            marker_counts=projected(chosen[:,marker_present],marker_destination[marker_present],len(marker_names)).toarray()
            marker_matrix[start+lo:start+lo+len(ids)]=np.log1p(marker_counts*(10000./total[:,None])).astype(np.float32)
            all_totals.extend(total.tolist())
            if lo%20480==0:event('native_RNA_block',dataset=source.record['dataset'],processed=lo+len(ids),total=len(selected))
        local['native_analysis_denominator']=np.asarray(all_totals)
        local['source_record_number']=number
        registry_parts.append(local)
        expression.append(sparse.vstack(source_rows,format='csr'))
        counts.append(sparse.vstack(source_raw,format='csr'))
        metrics.append(dict(dataset=source.record['dataset'],n_cells=len(local),n_native_features=len(source.var),
                            n_shared_features=width,n_native_marker_genes=int(marker_available[number].sum())))
        start+=len(selected);gc.collect()
    obs=pd.concat(registry_parts)
    marker_matrix=marker_matrix[:len(obs)]
    model=ad.AnnData(sparse.vstack(expression,format='csr'),obs=obs,
                    var=var.loc[common].copy())
    model.layers['counts']=sparse.vstack(counts,format='csr')
    del expression,counts;gc.collect()
    model.obsm['X_native_marker']=marker_matrix
    model.uns['native_marker_names']=np.asarray(marker_names,dtype=str)
    model.uns['native_marker_available']=marker_available
    model.uns['native_source_order']=np.asarray([s.record['dataset'] for s in sources],dtype=str)
    model.uns['model_genes_scope']='shared-native-feature modeling view; upstream fullgene counts preserved'
    model.uns['source_roles_frozen_not_outcomes']=True
    model.uns['gene_identity_diagnostics']=diagnostics
    model.uns['normalization']='native chosen counts / native full-library denominator *10000; log1p'
    model.uns['assay_mask_semantics']='structural or retained-export absence not a biological zero'
    model.obs['major_marker_proposal'],model.obs['major_marker_margin']=broad_proposals(model,major)
    model.obs['annotation_exclude_from_primary']=model.obs.get('residual_doublet_flag',pd.Series('',index=model.obs.index)).astype(str).str.lower().isin(['true','1'])
    event('full_native_model_matrix',n_cells=model.n_obs,n_shared_features=width)
    rng=np.random.default_rng(SEED)
    groups=model.obs.loc[~model.obs.annotation_exclude_from_primary].groupby('dataset',observed=True).indices
    per=min(50000,max(1,250000//len(groups)))
    eligible_positions=np.flatnonzero(~model.obs.annotation_exclude_from_primary.to_numpy())
    landmarks=np.sort(np.concatenate([eligible_positions[rng.choice(v,min(per,len(v)),replace=False)] for v in groups.values()]))
    fit=model[landmarks].copy()
    sc.pp.highly_variable_genes(fit,n_top_genes=3000,flavor='seurat_v3',batch_key='dataset',layer='counts')
    hvg=fit.var.highly_variable.to_numpy();model.var['highly_variable']=hvg
    if hvg.sum()<500:raise ValueError('Too few usable HVGs')
    sc.tl.pca(fit,n_comps=50,mask_var='highly_variable',svd_solver='arpack',random_state=SEED)
    basis=fit.varm['PCs'][hvg].astype(np.float32)
    center=np.asarray(fit.X[:,hvg].mean(axis=0)).ravel()@basis
    pca=np.empty((model.n_obs,50),dtype=np.float32)
    for lo in range(0,model.n_obs,4096):pca[lo:lo+4096]=model.X[lo:lo+4096,hvg]@basis-center
    model.obsm['X_pca']=pca
    model.varm['PCs']=fit.varm['PCs'].astype(np.float32)
    model.uns['pca']=fit.uns['pca']
    model.uns['PCA_fitting_scope']=dict(n_landmarks=len(landmarks),max_per_dataset=per,all_cells_projected=True)
    del fit;gc.collect()
    sample=np.sort(rng.choice(model.n_obs,min(10000,model.n_obs),replace=False))
    for phase in ['before','after']:
        if phase=='after':
            event('Harmony_full_cells_start',n_cells=model.n_obs,n_pcs=50,batch_key='dataset')
            result=harmonypy.run_harmony(model.obsm['X_pca'],model.obs,'dataset',theta=1.,max_iter_harmony=20,random_state=SEED)
            corrected=np.asarray(result.Z_corr)
            if corrected.shape==model.obsm['X_pca'].T.shape:corrected=corrected.T
            if corrected.shape!=model.obsm['X_pca'].shape or not np.isfinite(corrected).all():raise ValueError('Harmony result invalid')
            model.obsm['X_pca_harmony']=corrected.astype(np.float32)
            del result,corrected;gc.collect()
        metric=integration_metrics(model,sample,'X_pca' if phase=='before' else 'X_pca_harmony')
        metrics.append(dict(phase=phase,**metric))
    sc.settings.n_jobs=min(16,os.cpu_count() or 1)
    event('full_neighbor_graph_start',n_cells=model.n_obs)
    sc.pp.neighbors(model,n_neighbors=30,use_rep='X_pca_harmony',random_state=SEED)
    if model.n_obs<=100000:
        event('full_UMAP_start',n_cells=model.n_obs,n_epochs=100)
        sc.tl.umap(model,maxiter=100,min_dist=0.3,random_state=SEED)
        model.uns['UMAP_fitting_scope']='all_cells'
    else:
        # Large-cohort skill: fit a bounded view and retain coordinates for all cells.
        import umap
        take=landmarks if len(landmarks)<=100000 else np.sort(rng.choice(landmarks,100000,replace=False))
        event('balanced_UMAP_start',n_fit_cells=len(take),n_total_cells=model.n_obs)
        embedding=umap.UMAP(n_neighbors=30,min_dist=0.3,n_components=2,n_epochs=100,
                            random_state=SEED,transform_seed=SEED,n_jobs=1)
        fitted=embedding.fit_transform(model.obsm['X_pca_harmony'][take])
        coords=np.empty((model.n_obs,2),dtype=np.float32)
        coords[take]=fitted
        remaining=np.setdiff1d(np.arange(model.n_obs),take)
        for lo in range(0,len(remaining),10000):
            rows=remaining[lo:lo+10000]
            coords[rows]=embedding.transform(model.obsm['X_pca_harmony'][rows])
            if lo%100000==0:event('UMAP_full_projection',n_projected=lo+len(rows),n_remaining=len(remaining))
        if not np.isfinite(coords).all():raise ValueError('Nonfinite projected UMAP')
        model.obsm['X_umap']=coords
        model.uns['UMAP_fitting_scope']=dict(n_balanced_fit_cells=len(take),n_transformed_cells=len(remaining),all_cells_retained=True)
    sc.tl.leiden(model,resolution=0.1,key_added='coarse_cluster',flavor='igraph',n_iterations=2,directed=False,random_state=SEED)
    model.uns['integration_metrics_public']=json.dumps(metrics)
    target=args.scratch/('immune_model_'+args.phase+'.zarr')
    if target.exists():raise ValueError('Existing model attempt preserved')
    model.write_zarr(target)
    save_json_new(args.run/'manifests'/('integration_'+args.phase+'_public.json'),dict(
        observed_utc=datetime.now(timezone.utc).isoformat(),n_cells=model.n_obs,n_shared_features=model.n_vars,
        n_HVG=int(hvg.sum()),n_PCA_landmarks=len(landmarks),all_cells_PCA_Harmony_graph_and_UMAP_coordinates=True,
        UMAP_fitting_scope=model.uns['UMAP_fitting_scope'],
        n_coarse_clusters=model.obs.coarse_cluster.nunique(),pilot_only=args.phase=='pilot',
        metrics=metrics,authoritative_release=False,clinical_validation_frozen=False))
    event('integration_complete',n_cells=model.n_obs,pilot_only=args.phase=='pilot')


def panel_scores(model,panel):
    names=list(model.uns['native_marker_names']);matrix=model.obsm['X_native_marker']
    sources=pd.to_numeric(model.obs.source_record_number).to_numpy(dtype=int)
    available=np.asarray(model.uns['native_marker_available'])
    values=np.zeros((model.n_obs,len(panel)),dtype=np.float32)
    anchors=np.zeros_like(values,dtype=np.int16)
    for i,(label,genes) in enumerate(panel.items()):
        columns=[names.index(g) for g in genes if g in names]
        if not columns:continue
        measured=available[:,columns][sources]
        denominator=measured.sum(axis=1)
        expr=matrix[:,columns]
        values[:,i]=np.divide((expr*measured).sum(axis=1),denominator,out=np.zeros(model.n_obs,dtype=np.float32),where=denominator>=2)
        anchors[:,i]=((expr>0.5)&measured).sum(axis=1)
    return values,anchors


def broad_proposals(model,major):
    score,anchors=panel_scores(model,major)
    order=np.argsort(score,axis=1)
    best=order[:,-1];second=order[:,-2]
    rows=np.arange(len(score));margin=score[rows,best]-score[rows,second]
    labels=np.array(list(major),dtype=object)[best]
    accepted=(score[rows,best]>=0.25)&(anchors[rows,best]>=2)&(margin>=0.08)
    labels[~accepted]='Unresolved'
    # Keep related-lineage boundary evidence rather than forcing T/NK identity.
    names=list(major)
    if 'T' in names and 'NK' in names:
        marker_names=list(model.uns['native_marker_names'])
        columns=[marker_names.index(g) for g in ['CD3D','CD3E','TRAC'] if g in marker_names]
        coherent_T=(np.asarray(model.obsm['X_native_marker'][:,columns])>0.5).sum(axis=1)>=2
        labels[coherent_T&np.isin(labels,['T','NK','Unresolved','T_NK_boundary'])]='T'
        boundary=(~accepted)&(anchors[:,names.index('T')]>=2)&(anchors[:,names.index('NK')]>=2)
        labels[boundary&~coherent_T]='T_NK_boundary'
    return labels,margin.astype(np.float32)


def integration_metrics(model,sample,basis):
    from sklearn.neighbors import NearestNeighbors
    coords=model.obsm[basis][sample]
    indices=NearestNeighbors(n_neighbors=min(31,len(sample))).fit(coords).kneighbors(coords,return_distance=False)[:,1:]
    batch=model.obs.dataset.astype(str).to_numpy()[sample]
    major=model.obs.major_marker_proposal.astype(str).to_numpy()[sample]
    entropy=[];purity=[]
    for i,row in enumerate(indices):
        counts=np.unique(batch[row],return_counts=True)[1];p=counts/counts.sum()
        entropy.append(float(-(p*np.log(p)).sum()))
        if major[i] not in ['Unresolved','T_NK_boundary']:purity.append(float(np.mean(major[row]==major[i])))
    return dict(n_evaluation_cells=len(sample),mean_batch_neighbor_entropy=float(np.mean(entropy)),
                known_marker_major_neighbor_purity=float(np.mean(purity)) if purity else None,
                biology_labels='native_multigene_proposals_not_ground_truth',
                clinical_outcomes_used=False)


def lineage_hvg(fit, name):
    """Reuse historical robust sparse-lineage fallback with an explicit audit."""
    sizes=fit.obs.dataset.astype(str).value_counts()
    minimum=int(sizes.min())
    fallback=fit.n_obs<5000 or minimum<200
    reason='small_lineage_or_source_batch_below_200'
    if not fallback:
        try:
            sc.pp.highly_variable_genes(fit,n_top_genes=3000,flavor='seurat_v3',batch_key='dataset',layer='counts')
            return dict(major=name,used='batch_aware_seurat_v3',input='native_integer_counts',
                        n_fit_cells=fit.n_obs,min_source_cells=minimum,n_HVG=int(fit.var.highly_variable.sum()))
        except ValueError as exc:
            if not any(word in str(exc).lower() for word in ['singular','loess','extrapolation','svddc']):raise
            reason='observed_seurat_v3_numerical_failure: '+str(exc)
    matrix=fit.X
    means=np.asarray(matrix.mean(axis=0)).ravel().astype(np.float64)
    mean_square=np.asarray(matrix.power(2).mean(axis=0)).ravel().astype(np.float64)
    variance=np.maximum(mean_square-means*means,0.)
    dispersion=variance/np.maximum(means,1e-12)
    log_dispersion=np.log1p(dispersion)
    eligible=np.flatnonzero((means>0)&(variance>0)&np.isfinite(log_dispersion))
    if len(eligible)<3:raise ValueError('Too few variable native-lineage genes')
    bins=pd.qcut(means[eligible],q=min(20,len(np.unique(means[eligible]))),labels=False,duplicates='drop')
    bins=np.asarray(bins,dtype=float)
    if np.isnan(bins).all():bins=np.zeros(len(eligible))
    scores=np.full(fit.n_vars,-np.inf)
    for bin_id in np.unique(bins[np.isfinite(bins)]):
        positions=eligible[bins==bin_id];values=log_dispersion[positions]
        scale=max(float(values.std()),1e-8)
        scores[positions]=(values-values.mean())/scale
    ordered=eligible[np.argsort(-scores[eligible],kind='stable')]
    chosen=ordered[:min(3000,len(eligible))]
    mask=np.zeros(fit.n_vars,dtype=bool);mask[chosen]=True
    fit.var['highly_variable']=mask
    report=dict(major=name,used='historical_robust_dispersion_qcut',input='native_log1p_CP10K_X',
                n_fit_cells=fit.n_obs,min_source_cells=minimum,n_HVG=len(chosen),reason=reason,
                native_counts_unchanged=True,batch_correction_after_HVG='Harmony_on_PCA_by_source')
    event('lineage_HVG_fallback',**report)
    return report


def annotate(args):
    target=args.scratch/'immune_model_integrate.zarr'
    model=ad.read_zarr(target)
    major,subtype,states=panels(args.panel_dir)
    model.obs['major']=model.obs.major_marker_proposal.astype(str)
    # Official protocol retains tumor/stroma with a broad NFAT gate and adds
    # a stringent enriched fraction. It is not an immune-only sample.
    mixed=model.obs.dataset.astype(str).eq('GSE292046')
    model.obs['capture_scope']=model.obs.capture_scope.astype(str)
    model.obs.loc[mixed,'capture_scope']='mixed_tumor_nuclei_NFAT_enrichment'
    model.obs['capture_bias_status']='source_capture_scope_preserved'
    model.obs.loc[mixed,'capture_bias_status']='broad_gate_plus_T_enriched_nuclei_composition_biased'
    score,_=panel_scores(model,states)
    for i,name in enumerate(states):model.obs['state_score_'+name]=score[:,i]
    cell_groups=np.full(model.n_obs,'',dtype=object);catalog=[];number=0;hvg_reports=[]
    checkpoints=args.scratch/'lineage_annotation_checkpoints_v5';checkpoints.mkdir(exist_ok=True,mode=0o700)
    sc.settings.n_jobs=min(16,int(os.environ.get('NUMBA_NUM_THREADS','16')))
    for name in sorted(model.obs.major.astype(str).unique()):
        mask=model.obs.major.astype(str).eq(name).to_numpy();positions=np.flatnonzero(mask)
        if len(positions)>=100:
            local=model[positions].copy()
            local.obsp.clear()
            rng=np.random.default_rng(SEED)
            eligible=~local.obs.annotation_exclude_from_primary.astype(str).str.lower().isin(['true','1']).to_numpy()
            eligible_positions=np.flatnonzero(eligible)
            groups_local=local.obs.iloc[eligible_positions].groupby('dataset',observed=True).indices
            per_source=max(1,100000//len(groups_local))
            fit_indices=np.sort(np.concatenate([eligible_positions[rng.choice(v,min(len(v),per_source),replace=False)] for v in groups_local.values()]))
            fit=local[fit_indices].copy()
            event('lineage_HVG_start',major=name,n_fit_cells=len(fit))
            report=lineage_hvg(fit,name);hvg_reports.append(report)
            save_json_new(checkpoints/(name+'_HVG_public.json'),report)
            sc.tl.pca(fit,n_comps=min(30,len(fit)-1,int(fit.var.highly_variable.sum())-1),mask_var='highly_variable',svd_solver='arpack',random_state=SEED)
            hvg=fit.var.highly_variable.to_numpy();basis=fit.varm['PCs'][hvg]
            center=np.asarray(fit.X[:,hvg].mean(axis=0)).ravel()@basis
            projected_pca=np.vstack([local.X[lo:lo+4096,hvg]@basis-center for lo in range(0,len(local),4096)])
            local.obsm['X_lineage_pca']=projected_pca.astype(np.float32)
            del fit;gc.collect()
            if local.obs.dataset.nunique()>1:
                result=harmonypy.run_harmony(projected_pca,local.obs,'dataset',theta=1.,max_iter_harmony=10,random_state=SEED)
                z=np.asarray(result.Z_corr)
                if z.shape==projected_pca.T.shape:z=z.T
                if z.shape!=projected_pca.shape or not np.isfinite(z).all():raise ValueError('Invalid lineage Harmony')
                local.obsm['X_lineage_harmony']=z.astype(np.float32)
                del result,z
            else:local.obsm['X_lineage_harmony']=local.obsm['X_lineage_pca']
            event('lineage_specific_neighbors',major=name,n_cells=len(local),n_balanced_PCA_fit_cells=len(fit_indices))
            sc.pp.neighbors(local,n_neighbors=min(30,len(local)-1),use_rep='X_lineage_harmony',random_state=SEED)
            sc.tl.leiden(local,resolution=0.6,key_added='within_major',flavor='igraph',n_iterations=2,directed=False,random_state=SEED)
            groups=local.obs.within_major.astype(str).to_numpy()
            np.save(checkpoints/(name+'_clusters_private.npy'),groups)
            save_json_new(checkpoints/(name+'_complete_public.json'),dict(major=name,n_cells=len(local),n_clusters=len(np.unique(groups)),all_cells_projected=True))
        else:groups=np.repeat('0',len(positions));local=model[positions].copy()
        fine_panels=subtype.get(name,{})
        fine_score,_=panel_scores(local,fine_panels) if fine_panels else (np.zeros((len(local),0)),None)
        marker_names=list(model.uns['native_marker_names'])
        availability=np.asarray(model.uns['native_marker_available'],dtype=bool)
        source_number=pd.to_numeric(model.obs.source_record_number).to_numpy(dtype=int)
        major_source_counts=np.bincount(source_number[positions],minlength=len(availability))
        major_measured=major_source_counts@availability
        major_markers=model.obsm['X_native_marker'][positions]
        major_positive=(major_markers>0).sum(axis=0)
        major_cp=np.expm1(major_markers).sum(axis=0,dtype=np.float64)
        for group in sorted(np.unique(groups),key=lambda x:(len(x),x)):
            inside=groups==group;global_positions=positions[inside];number+=1
            subtype_label=name;fine_status='major_only';best_fine=''
            if fine_panels:
                mean=np.mean(fine_score[inside],axis=0);order=np.argsort(mean)
                best_fine=list(fine_panels)[order[-1]]
                if len(order)>1 and mean[order[-1]]>=0.5 and mean[order[-1]]-mean[order[-2]]>=0.15:
                    # Candidate subtype is a reviewable evidence field; naming stays descriptive.
                    fine_status='multigene_subtype_candidate'
            native=model.obsm['X_native_marker'][global_positions]
            group_measured=np.bincount(source_number[global_positions],minlength=len(availability))@availability
            group_positive=(native>0).sum(axis=0)
            group_cp=np.expm1(native).sum(axis=0,dtype=np.float64)
            pct=np.divide(group_positive,group_measured,out=np.full(len(marker_names),np.nan),where=group_measured>0)
            mean=np.divide(native.sum(axis=0,dtype=np.float64),group_measured,out=np.full(len(marker_names),np.nan),where=group_measured>0)
            outside_n=major_measured-group_measured
            outside_pct=np.divide(major_positive-group_positive,outside_n,out=np.full(len(marker_names),np.nan),where=outside_n>0)
            cp_in=np.divide(group_cp,group_measured,out=np.full(len(marker_names),np.nan),where=group_measured>0)
            cp_out=np.divide(major_cp-group_cp,outside_n,out=np.full(len(marker_names),np.nan),where=outside_n>0)
            ratio=np.log2((cp_in+0.01)/(cp_out+0.01))
            genes=major.get(name,[])+[g for panel in fine_panels.values() for g in panel]
            candidates=[]
            for gene in sorted(set(genes)):
                if gene not in marker_names or gene in ['MALAT1','FOS','ACTB','B2M']:continue
                j=marker_names.index(gene)
                if pct[j]>=0.2:
                    specific=bool(np.isfinite(outside_pct[j]) and pct[j]-outside_pct[j]>=0.1 and ratio[j]>=0.5)
                    candidates.append((specific,float(mean[j]),float(pct[j]),gene))
            choice=sorted(candidates,key=lambda x:(not x[0],-x[1],-x[2],x[3]))[0] if candidates else None
            marker=choice[3] if choice else 'NA'
            group_id=hashlib.sha256((name+'::'+group+'::'+str(SEED)).encode()).hexdigest()[:16]
            display=f'{subtype_label}_{marker}_{number}'
            cell_groups[global_positions]=display
            catalog.append(dict(node_id=group_id,number=number,major=name,within_major_cluster=group,
                                subcluster=display,marker=marker,marker_role='own_group_expression_contrast' if choice and choice[0] else 'own_group_representative_non_specific',
                                subtype_candidate=best_fine,subtype_evidence_status=fine_status,n_cells=len(global_positions),
                                n_libraries=model.obs.iloc[global_positions].sampleID.nunique(),
                                n_sources=model.obs.iloc[global_positions].dataset.nunique(),identity_promoted=False))
        del local,major_markers;gc.collect()
    model.obs['subcluster']=pd.Categorical(cell_groups)
    model.obs['annotation_status']='multigene_major_proposal_descriptive_subcluster_fine_review_pending'
    model.uns['lineage_HVG_audit']=json.dumps(hvg_reports)
    model.obs['malignant_status']='not_yet_evaluated_per_native_library'
    destination=args.scratch/'immune_annotated.zarr'
    if destination.exists():raise ValueError('Existing annotation attempt preserved')
    model.write_zarr(destination)
    pd.DataFrame(catalog).to_csv(args.run/'tables/subcluster_catalog_private.tsv',sep='\t',index=False)
    model.obs.to_csv(args.run/'tables/all_core_annotation_private.tsv.gz',sep='\t')
    save_json_new(args.run/'manifests/annotation_public.json',dict(observed_utc=datetime.now(timezone.utc).isoformat(),
        n_cells=model.n_obs,n_nodes=len(catalog),major_counts=model.obs.major.astype(str).value_counts().to_dict(),
        state_axes=list(states),all_raw_counts_preserved_in_upstream=True,display_marker_count=1,
        lineage_HVG_audit=hvg_reports,fine_identity_promoted=False,CNV_completed=False,authoritative_release=False))
    event('annotation_complete',n_cells=model.n_obs,n_nodes=len(catalog))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--contract',type=Path,required=True);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--scratch',type=Path,required=True);p.add_argument('--panel-dir',type=Path,required=True)
    p.add_argument('--phase',choices=['pilot','integrate','annotate'],required=True)
    args=p.parse_args();args.scratch.mkdir(parents=True,exist_ok=True,mode=0o700)
    if args.phase=='annotate':annotate(args)
    else:build_model(args)
