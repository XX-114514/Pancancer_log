#!/usr/bin/env python3
"""Append-only library preprocessing; real data stay on the execution server.

Reuses audited QC/Scrublet, fits an independent doublet method on pre-Scrublet
cells, and keeps original counts separate from estimated decontaminated counts.
No source AC or cross-study integration is granted by this script.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

CONTEXT = {}


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8*1024**2), b''):
            h.update(b)
    return h.hexdigest()


def save(path, value):
    with open(path, 'x') as f:
        json.dump(value, f, indent=2, allow_nan=False)


def event(dataset, phase, **kwargs):
    print(json.dumps(dict(utc=now(), dataset=dataset, phase=phase, **kwargs)), flush=True)


def external(command, log, dataset, library):
    with log.open('x') as f:
        p = subprocess.Popen(command, stdout=f, stderr=subprocess.STDOUT)
        while p.poll() is None:
            event(dataset, 'external_heartbeat', library=library, pid=p.pid,
                  log_bytes=log.stat().st_size)
            time.sleep(20)
    if p.returncode:
        raise RuntimeError('External method exit %d; see private log' % p.returncode)


def marker_review(raw, corrected, groups, genes):
    import numpy as np
    import pandas as pd
    sets = {'T_NK': ['CD3D','CD3E','TRAC','NKG7','GNLY'],
            'B_plasma': ['CD79A','CD79B','MS4A1','MZB1','JCHAIN'],
            'myeloid': ['LYZ','LST1','TYROBP','FCER1G'],
            'epithelial': ['EPCAM','KRT8','KRT18','KRT19'],
            'endothelial': ['PECAM1','VWF','CLDN5'],
            'fibroblast': ['COL1A1','COL1A2','DCN','LUM'],
            'mast': ['TPSAB1','TPSB2','CPA3'],
            'pericyte': ['RGS5','PDGFRB','CSPG4']}
    gene_to_pos = {g:i for i,g in enumerate(genes)}
    records = []
    # Cluster profiles are diagnostics, never final biological annotation.
    for cluster in sorted(set(groups)):
        pos = np.flatnonzero(groups == cluster)
        for lineage, markers in sets.items():
            for g in markers:
                if g not in gene_to_pos:
                    continue
                j = gene_to_pos[g]; a = raw[pos,j]; b = corrected[pos,j]
                before = float(a.sum()); after = float(b.sum())
                records.append(dict(cluster=str(cluster), diagnostic_lineage=lineage,
                                    marker=g, n_cells=len(pos), raw_UMI=before,
                                    corrected_UMI=after,
                                    raw_detected_fraction=float(a.getnnz()/len(pos)),
                                    corrected_detected_fraction=float((b>0.1).sum()/len(pos)),
                                    UMI_retained_fraction=after/before if before else 1.0))
    return pd.DataFrame(records)


def broad_groups(b, dataset):
    """Merge coarse clusters by marker programs; keep unknowns explicitly."""
    import numpy as np
    from scipy import sparse
    sets = {'T_NK': ['CD3D','CD3E','TRAC','NKG7','GNLY'],
            'B_plasma': ['CD79A','CD79B','MS4A1','MZB1','JCHAIN'],
            'myeloid': ['LYZ','LST1','TYROBP','FCER1G'],
            'epithelial': ['EPCAM','KRT8','KRT18','KRT19'],
            'endothelial': ['PECAM1','VWF','CLDN5'],
            'fibroblast': ['COL1A1','COL1A2','DCN','LUM'],
            'mast': ['TPSAB1','TPSB2','CPA3'],
            'pericyte': ['RGS5','PDGFRB','CSPG4']}
    labels = np.full(b.n_obs, 'unknown', dtype=object)
    if dataset.startswith('CXG_5b9b89d2'):
        t = b.obs['cell_type'].astype(str)
        labels[t.str.contains('T cell|thymocyte', case=False).to_numpy()] = 'T_NK'
        labels[t.str.contains('monocyte', case=False).to_numpy()] = 'myeloid'
        return labels
    genes = {g:i for i,g in enumerate(b.var_names)}
    for cluster in sorted(set(b.obs['preliminary_cluster'].astype(str))):
        positions = np.flatnonzero(b.obs['preliminary_cluster'].astype(str).to_numpy()==cluster)
        mean = np.asarray(b.X[positions].mean(axis=0)).ravel()
        scores = {}
        for k,v in sets.items():
            ids = [genes[g] for g in v if g in genes]
            if len(ids) < 2:
                continue
            values = mean[ids]
            scores[k] = float(np.mean(values)) if np.count_nonzero(values>0.15)>=2 else 0.0
        ranking = sorted(scores.items(), key=lambda z:z[1], reverse=True)
        if ranking and ranking[0][1]>0.25 and ranking[0][1]>1.5*ranking[1][1]:
            labels[positions] = ranking[0][0]
    return labels


def library_job(a, x_pre, pre_ledger, positions, library, dataset, destination,
                rscript, rworker, seed):
    import anndata as ad
    import doubletdetection
    import numpy as np
    import pandas as pd
    from scipy import sparse, io
    import scanpy as sc
    started = now(); destination.mkdir(parents=True, exist_ok=False)
    result = dict(dataset=dataset, library_token=destination.name, status='running',
                  start_time=started, n_input_QC=len(positions), n_pre_Scrublet=x_pre.shape[0])
    params = dict(n_iters=10, clustering_algorithm='louvain', standard_scaling=True,
                  pseudocount=0.1, n_jobs=1, random_state=seed,
                  n_components=min(30, x_pre.shape[0]-2), n_top_var_genes=10000,
                  p_thresh=1e-16, voter_thresh=0.5, max_remove_fraction=0.08)
    save(destination/'params.json', dict(doubletdetection=params,
         decontX=dict(seed=seed,maxIter=500,delta=[10,10],estimateDelta=True,
                      background='none_available',z='marker_merged_broad_clusters'),
         normalization=dict(target_sum=10000,log1p=True,HVG=2000),
         counts_semantics='counts=immutable_original; counts_decontx=fractional_estimate'))
    try:
        event(dataset, 'doublet_fit', library=destination.name, n_cells=x_pre.shape[0])
        fit_params = {k:v for k,v in params.items() if k not in ['p_thresh','voter_thresh','max_remove_fraction']}
        nonzero = x_pre.getnnz(axis=0)>0
        doublet_evaluable=x_pre.shape[0]>=500
        if doublet_evaluable:
            clf = doubletdetection.BoostClassifier(**fit_params)
            calls = clf.fit(x_pre[:,nonzero]).predict(p_thresh=params['p_thresh'], voter_thresh=params['voter_thresh'])
            scores = clf.doublet_score()
        else:
            calls=np.zeros(x_pre.shape[0],dtype=bool);scores=np.full(x_pre.shape[0],np.nan)
        if (doublet_evaluable and not np.isfinite(scores).all()) or not np.isfinite(calls).all():
            raise ValueError('Doublet scores or calls nonfinite')
        calls = calls.astype(bool)
        ledger = pre_ledger.copy()
        ledger['doubletdetection_score'] = scores
        ledger['doubletdetection_predicted'] = calls
        old_calls = ledger.doublet_status.eq('predicted_doublet').to_numpy()
        result.update(n_DD_predicted=int(calls.sum()), n_prior_Scrublet_predicted=int(old_calls.sum()),
                      n_method_overlap=int(np.count_nonzero(calls & old_calls)))
        result['independent_doublet_status']='computed' if doublet_evaluable else 'skipped_lt500_flags_only'
        result['prior_doublet_method']='published_author_flag_or_none' if 'prior_doublet_method' in ledger else 'audited_Scrublet'
        cap = int(np.floor(len(calls)*0.08))
        union = calls | old_calls
        overcap = int(union.sum())>cap
        # Preserve all old QC decisions; remove additional highest-scoring candidates
        # only within the residual cap. Remaining candidates keep flags and block release.
        added = np.flatnonzero(calls & ~old_calls)
        allowance = max(0, cap-int(old_calls.sum()))
        selected_added = added[np.argsort(scores[added])[::-1][:allowance]]
        removed = old_calls.copy(); removed[selected_added] = True
        ledger['additional_doublet_removed'] = False
        ledger.loc[ledger.index[selected_added], 'additional_doublet_removed'] = True
        ledger.to_csv(destination/'doublet_scores_private.tsv.gz',sep='\t',index=False)
        lookup = ledger.set_index('original_cell_id')
        b = a[positions].copy()
        ids = b.obs.original_cell_id.astype(str)
        if not ids.is_unique or not ids.isin(lookup.index).all():
            raise ValueError('Doublet-to-QC cell mapping not exact')
        mapped = lookup.loc[ids]
        b.obs['doubletdetection_score'] = mapped.doubletdetection_score.to_numpy()
        b.obs['doubletdetection_predicted'] = mapped.doubletdetection_predicted.to_numpy()
        b.obs['additional_doublet_removed'] = mapped.additional_doublet_removed.to_numpy()
        result['n_additional_removed'] = int(b.obs.additional_doublet_removed.sum())
        b = b[~b.obs.additional_doublet_removed.to_numpy()].copy()
        raw = sparse.csr_matrix(b.layers['counts']).copy()
        gene_strings=b.var_names.astype(str)
        for metric,mask in [('ribo',gene_strings.str.upper().str.match(r'^RP[SL]')),
                            ('hb',gene_strings.str.upper().isin(['HBA1','HBA2','HBB','HBD','HBG1','HBG2','HBE1','HBZ'])),
                            ('ercc',gene_strings.str.upper().str.startswith('ERCC-'))]:
            mask=np.asarray(mask); totals=np.asarray(raw.sum(axis=1)).ravel()
            amounts=np.asarray(raw[:,mask].sum(axis=1)).ravel() if mask.any() else np.zeros(b.n_obs)
            values=np.divide(amounts*100,totals,out=np.zeros(b.n_obs),where=totals>0)
            b.obs['pct_counts_'+metric]=values;b.obs[metric+'_metric_evaluable']=bool(mask.any() and np.any(amounts>0))
            median=np.median(values);mad=np.median(np.abs(values-median))
            b.obs[metric+'_diagnostic_outlier_flag']=values>median+3*mad
        b.X = raw.astype(np.float32).copy()
        sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b)
        if b.n_obs<100:
            b.layers['counts']=raw;b.layers['counts_decontx']=raw.astype(np.float32)
            b.obs['ambient_status']='not_identifiable_small_library'
            b.obs['decontX_contamination']=np.nan
            b.obs['preintegration_release_eligible']=False
            b.obs['preintegration_cell_eligible']=False
            b.uns['preprocess_v1']=dict(normalized_from='counts',source_final_AC=False,
                       ambient_background='missing_empty_droplets',embedding='skipped_small_library')
            b.write_zarr(destination/'prepared.zarr')
            result.update(status='success',end_time=now(),n_output=b.n_obs,n_genes=b.n_vars,
                ambient_status='not_identifiable_small_library',ambient_validation_pass=False,
                n_residual_doublet_flag=0,preintegration_release_eligible=False,
                n_credible_broad_lineages=0,broad_diagnostic_counts={},
                raw_UMI=int(raw.sum()),corrected_UMI=float(raw.sum()),
                median_removed_UMI_fraction=0.0,p95_removed_UMI_fraction=0.0)
            save(destination/'summary.json',result);return result
        sc.pp.highly_variable_genes(b,n_top_genes=min(2000,b.n_vars-1),flavor='seurat')
        hvg = b.var.highly_variable.to_numpy()
        if hvg.sum()<500:
            raise ValueError('Fewer than 500 HVGs; gene/count diagnosis required')
        temp = ad.AnnData(b.X[:,hvg].copy())
        pc = min(15 if b.n_obs<=5000 else 20 if b.n_obs<=10000 else 30 if b.n_obs<=40000 else 40,
                 temp.n_vars-1,temp.n_obs-2)
        sc.pp.scale(temp,max_value=10,zero_center=False)
        sc.tl.pca(temp,n_comps=pc,zero_center=True,svd_solver='arpack',random_state=seed)
        b.obsm['X_pca_raw_review'] = temp.obsm['X_pca'].copy()
        sc.pp.neighbors(b, n_neighbors=15,n_pcs=pc,use_rep='X_pca_raw_review',random_state=seed)
        sc.tl.leiden(b,resolution=0.3,random_state=seed,key_added='preliminary_cluster',
                     flavor='igraph',n_iterations=2,directed=False)
        groups = broad_groups(b,dataset)
        b.obs['ambient_diagnostic_group'] = groups.astype(str)
        known = pd.Series(groups)[pd.Series(groups).ne('unknown')].value_counts()
        credible = known[known.ge(max(50,int(b.n_obs*0.01)))]
        identifiable = len(credible)>=2
        result['n_credible_broad_lineages'] = len(credible)
        result['broad_diagnostic_counts'] = pd.Series(groups).value_counts().astype(int).to_dict()
        if identifiable:
            event(dataset,'ambient_fit',library=destination.name,n_cells=b.n_obs)
            io.mmwrite(str(destination/'input_counts_private.mtx'),raw.T)
            pd.DataFrame(dict(stable_cell_id=b.obs_names.astype(str),
                              broad_group=groups.astype(str))).to_csv(destination/'cells_private.tsv',sep='\t',index=False)
            pd.Series(b.var_names.astype(str)).to_csv(destination/'genes_private.tsv',index=False,header=False)
            external([rscript,str(rworker),str(destination),str(seed)],destination/'decontX.log',dataset,destination.name)
            dx = io.mmread(str(destination/'corrected_counts_private.mtx')).T.tocsr().astype(np.float32)
            cell_result = pd.read_csv(destination/'decontX_cells_private.tsv',sep='\t',keep_default_na=False)
            if not np.array_equal(cell_result.stable_cell_id.astype(str).to_numpy(), b.obs_names.to_numpy()):
                raise ValueError('R correction cell order mismatch')
            b.obs['decontX_contamination'] = cell_result.contamination.to_numpy(dtype=float)
            b.obs['decontX_cluster'] = cell_result.cluster.astype(str).to_numpy()
            result['ambient_fit'] = json.loads((destination/'decontX_summary.json').read_text())
            status = 'DecontX_filtered_only_estimated'
        else:
            dx = raw.astype(np.float32).copy()
            b.obs['decontX_contamination'] = np.nan
            b.obs['decontX_cluster'] = groups.astype(str)
            status = 'not_identifiable_without_empty_droplets'
        if dx.shape!=raw.shape or np.any(~np.isfinite(dx.data)) or np.any(dx.data<0):
            raise ValueError('Invalid corrected shape/values')
        diff = dx-raw
        if diff.nnz and np.max(diff.data)>1e-4:
            raise ValueError('Correction invents counts')
        b.layers['counts'] = raw
        b.layers['counts_decontx'] = dx
        before = np.asarray(raw.sum(axis=1)).ravel(); after = np.asarray(dx.sum(axis=1)).ravel()
        retained = np.divide(after,before,out=np.zeros(len(after)),where=before>0)
        b.obs['ambient_UMI_removed_fraction'] = 1-retained
        b.obs['ambient_high_loss_review'] = retained<0.5
        b.obs['ambient_status'] = status
        b.obs['residual_doublet_flag'] = b.obs.predicted_doublet | b.obs.doubletdetection_predicted
        review = marker_review(raw,dx,groups,b.var_names)
        review.to_csv(destination/'marker_preservation_diagnostics.tsv',sep='\t',index=False)
        # Dominant native-marker losses trigger review, not automatic removal of cells.
        native = review[(review.cluster==review.diagnostic_lineage) &
                        review.raw_detected_fraction.ge(0.25) & review.raw_UMI.ge(100)]
        marker_loss = bool(native.UMI_retained_fraction.lt(0.5).any())
        high_loss = bool(np.median(1-retained)>0.3 or np.mean(retained<0.5)>0.2)
        residual = int(b.obs.residual_doublet_flag.sum())
        convergence_pass=bool(result.get('ambient_fit',{}).get('converged',False)) if identifiable else False
        correction_pass = identifiable and convergence_pass and not marker_loss and not high_loss
        result.update(ambient_status=status, ambient_validation_pass=bool(correction_pass),
                      native_marker_loss_review=marker_loss,high_UMI_loss_review=high_loss,
                      ambient_convergence_pass=convergence_pass,
                      doublet_overcap_review=overcap,n_residual_doublet_flag=residual,
                      n_output=b.n_obs,n_genes=b.n_vars,
                      raw_UMI=int(raw.sum()),corrected_UMI=float(dx.sum(dtype=np.float64)),
                      median_removed_UMI_fraction=float(np.median(1-retained)),
                      p95_removed_UMI_fraction=float(np.quantile(1-retained,0.95)))
        # A review failure never promotes correction to the selected analysis counts.
        selected = dx if correction_pass else raw.astype(np.float32)
        b.X = selected.copy(); sc.pp.normalize_total(b,target_sum=1e4); sc.pp.log1p(b)
        sc.pp.highly_variable_genes(b,n_top_genes=min(2000,b.n_vars-1),flavor='seurat')
        b.var['feature_measured_in_source'] = True
        b.obs['preintegration_release_eligible'] = bool(correction_pass and doublet_evaluable and not overcap and residual==0)
        b.obs['preintegration_cell_eligible'] = b.obs.preintegration_release_eligible & ~b.obs.ambient_high_loss_review
        b.uns['preprocess_v1'] = dict(dataset=dataset,original_library_id=str(library),
            seed=seed,source_QC_ref=str(a.uns.get('raw_counts_source','')),
            normalized_from='counts_decontx' if correction_pass else 'counts',
            raw_counts_immutable=True,ambient_background='missing_empty_droplets',
            source_final_AC=False,cross_study_integration_performed=False,
            HVG_semantics='per_library_preparation_only; select_batch_aware_integration_HVG_later')
        b.write_zarr(destination/'prepared.zarr')
        # Diagnostic UMAP is library-local, never an integrated embedding.
        sc.tl.umap(b,random_state=seed)
        import zarr
        from anndata._io.specs import write_elem
        write_elem(zarr.open_group(str(destination/'prepared.zarr'),mode='a')['obsm'],
                   'X_umap_raw_review',b.obsm['X_umap'])
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1,3,figsize=(13,4))
        ax[0].scatter(before,b.obs.n_genes,s=1,alpha=.15);ax[0].set_xlabel('Original UMI');ax[0].set_ylabel('Detected genes')
        ax[1].hist(1-retained,bins=40);ax[1].set_xlabel('Estimated removed UMI fraction')
        coords=b.obsm['X_umap']; ax[2].scatter(coords[:,0],coords[:,1],s=1,c=1-retained,cmap='viridis')
        ax[2].set_title('Raw-count library-local UMAP'); fig.tight_layout();fig.savefig(destination/'QC_ambient_review.png',dpi=150);plt.close(fig)
        fig,axes=plt.subplots(1,2,figsize=(11,4))
        axes[0].violinplot([np.log1p(before),np.log1p(b.obs.n_genes.to_numpy())],showmedians=True)
        axes[0].set_xticks([1,2],['log1p original UMI','log1p original genes'])
        axes[1].violinplot([b.obs['pct_counts_'+k].to_numpy() for k in ['mt','ribo','hb']],showmedians=True)
        axes[1].set_xticks([1,2,3],['mt %','ribo %','Hb %']);fig.tight_layout()
        fig.savefig(destination/'QC_violin.png',dpi=150);plt.close(fig)
        result.update(status='success',end_time=now(),
                      preintegration_release_eligible=bool(correction_pass and doublet_evaluable and not overcap and residual==0),
                      n_preintegration_cell_eligible=int(b.obs.preintegration_cell_eligible.sum()),
                      output_sha_semantics='directory independently audited; input/ref hashes in manifest')
        save(destination/'summary.json',result)
        event(dataset,'library_complete',library=destination.name,n_output=b.n_obs,ambient=status,
              release_eligible=result['preintegration_release_eligible'])
    except Exception as exc:
        result.update(status='failed',end_time=now(),error=repr(exc))
        (destination/'traceback.txt').write_text(traceback.format_exc())
        save(destination/'summary.json',result)
        event(dataset,'library_failed',library=destination.name,error=repr(exc))
    return result


def context_job(item):
    import numpy as np
    i,lib,pos=item; c=CONTEXT; args=c['args']
    token='%03d_%s'%(i,hashlib.sha256(str(lib).encode()).hexdigest()[:12])
    prepos=np.flatnonzero(c['good'] & c['ledger'].library_id.astype(str).eq(str(lib)).to_numpy())
    return library_job(c['a'],c['x'][prepos].tocsr(),c['ledger'].iloc[prepos].copy(),np.asarray(pos),
                       lib,args.dataset,c['dest']/token,args.rscript,Path(args.rworker),20261011+i)


def native_contract(args):
    """QC an audited native shard; reuse explicitly published author PASS decisions."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    a=ad.read_zarr(args.native_object)
    x=sparse.csr_matrix(a.layers['counts'] if 'counts' in a.layers else a.X)
    if not a.var_names.is_unique or np.any(~np.isfinite(x.data)) or np.any(x.data<0) or np.any(x.data!=np.rint(x.data)):
        raise ValueError('Native counts/features invalid')
    total=np.asarray(x.sum(axis=1)).ravel(); ng=x.getnnz(axis=1)
    mt=np.asarray(a.var_names.str.upper().str.startswith('MT-'))
    mt_count=np.asarray(x[:,mt].sum(axis=1)).ravel() if mt.any() else np.zeros(a.n_obs)
    pct=np.divide(mt_count*100,total,out=np.zeros(a.n_obs),where=total>0)
    author_pass='pass' in a.obs and a.obs['pass'].astype(str).str.lower().eq('true').all()
    floor=(ng>=200)&(total>=500) if not author_pass else (total>0)
    keep=floor.copy();reasons=np.where(floor,'','low_genes_or_counts').astype(object)
    thresholds=[]
    groups=a.obs.groupby(a.obs.original_library_id.astype(str),sort=False).indices
    for lib,pos in groups.items():
        pos=np.asarray(pos);fp=pos[floor[pos]]
        med=float(np.median(pct[fp])) if len(fp) else 0
        mad=float(np.median(np.abs(pct[fp]-med))) if len(fp) else 0
        limit=min(25,max(15,med+3*mad));evaluable=bool(mt.any() and np.any(mt_count[pos]>0))
        # Author PASS is reused explicitly: no additional generic floor/mt deletions.
        if evaluable and not author_pass:
            bad=pos[pct[pos]>limit];keep[bad]=False;reasons[bad]='high_mitochondrial_fraction'
        thresholds.append(dict(library_token=hashlib.sha256(str(lib).encode()).hexdigest()[:12],
              n_raw=len(pos),n_floor=len(fp),mt_evaluable=evaluable,mt_threshold=limit if evaluable else None,
              author_PASS_reused=bool(author_pass)))
    prior=np.zeros(a.n_obs,dtype=bool)
    if 'doublet' in a.obs:
        prior=a.obs.doublet.astype(str).str.lower().isin(['true','1']).to_numpy()
    calls=prior&keep;reasons[calls]='sample_doublet';keep[calls]=False
    original=a.obs.original_cell_id.astype(str).to_numpy()
    ledger=pd.DataFrame(dict(original_cell_id=original,library_id=a.obs.original_library_id.astype(str).to_numpy(),
         filter_reason=reasons,doublet_status=np.where(prior,'predicted_doublet','predicted_singlet'),
         prior_doublet_method='published_author_flag' if 'doublet' in a.obs else 'not_previously_called'))
    if len(set(original))!=len(original):
        raise ValueError('Native original cell IDs not unique within study')
    a.obs['original_native_row_position']=np.arange(a.n_obs)
    a.obs_names=[args.dataset+'::'+str(i)+'::'+str(c) for i,c in enumerate(original)]
    a.obs['dataset']=args.dataset;a.obs['total_counts']=total;a.obs['n_genes']=ng;a.obs['pct_counts_mt']=pct
    a.obs['doublet_score']=a.obs.get('source_obs_doublet_score',a.obs.get('doublet_score',np.nan))
    a.obs['doublet_status']=ledger.doublet_status.to_numpy();a.obs['predicted_doublet']=prior
    a.layers['counts']=x
    a.uns['raw_counts_source']=str(args.native_object)
    a.uns['native_QC_author_PASS_reused']=bool(author_pass)
    return a[keep].copy(),x,ledger,thresholds


def dataset_job(args):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse
    sc.settings.n_jobs=1
    source = Path(args.collection); dest=Path(args.output)/args.dataset
    dest.mkdir(parents=True,exist_ok=False)
    started = now(); event(args.dataset,'load_QC')
    if args.native_object:
        a,x,ledger,thresholds=native_contract(args)
        save(dest/'native_QC_input_contract.json',dict(status='sample_QC_applied_pending_independent_audit',
             native_source_ref=args.native_object,n_raw=len(ledger),n_QC=len(a),thresholds=thresholds,
             author_PASS_reused=bool(a.uns['native_QC_author_PASS_reused']),source_final_AC=False))
        ledger.to_csv(dest/'native_QC_ledger_private.tsv.gz',sep='\t',index=False)
    else:
        a = ad.read_zarr(source/'objects'/('qc_'+args.dataset+'.zarr'))
        proof=json.loads((source/'manifests'/('accepted_'+args.dataset+'.json')).read_text())
        if proof['n_qc_nonduplicate']!=a.n_obs:
            raise ValueError('QC shard and registered cell count mismatch')
    if args.dataset=='GSE286827':
        path=source/'artifacts/QC_patient_mapping_repair_20261011_v2/GSE286827_QC_cell_metadata_overlay_v2.tsv.gz'
        if sha(path)!='4d291759b53727cf8d946651f66095eea994a558c0de3a394774cedce82cb576':
            raise ValueError('Patient overlay SHA mismatch')
        overlay=pd.read_csv(path,sep='\t',keep_default_na=False).set_index('stable_cell_id')
        if not overlay.index.is_unique or not overlay.index.isin(a.obs_names).all():
            raise ValueError('Patient overlay ID mapping not exact')
        for field in ['patient_uid','timepoint','patient_identity_status']:
            a.obs[field]=a.obs[field].astype(str)
            a.obs.loc[overlay.index,field]=overlay[field]
    if not args.native_object:
        native_path=Path(a.uns['raw_counts_source']);event(args.dataset,'load_pre_Scrublet_counts')
        native=ad.read_h5ad(native_path)
        x=sparse.csr_matrix(native.layers['counts'] if 'counts' in native.layers else native.X)
        ledger=pd.read_csv(source/'manifests'/(args.dataset+'_cell_qc_and_filter_reason.tsv.gz'),sep='\t',keep_default_na=False)
        if len(ledger)!=native.n_obs or not np.array_equal(ledger.original_cell_id.astype(str),native.obs_names.astype(str)):
            raise ValueError('Original ledger order does not match raw object')
    good=ledger.filter_reason.isin(['','sample_doublet','duplicate_original_counts_and_barcode']).to_numpy()
    groups=a.obs.groupby(a.obs.original_library_id.astype(str),sort=False).indices
    results=[]
    if args.library_limit:
        groups=dict(list(groups.items())[:args.library_limit])
    global CONTEXT
    CONTEXT=dict(a=a,x=x,ledger=ledger,good=good,args=args,dest=dest)
    tasks=[(i,lib,pos) for i,(lib,pos) in enumerate(groups.items())]
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('fork')) as pool:
      pending=[pool.submit(context_job,item) for item in tasks]
      for i,future in enumerate(as_completed(pending)):
        results.append(future.result())
        save(dest/('progress_%03d.json'%i),dict(dataset=args.dataset,utc=now(),completed=i+1,
            total=len(groups),succeeded=sum(r['status']=='success' for r in results),
            failed=sum(r['status']=='failed' for r in results),output_cells=sum(r.get('n_output',0) for r in results)))
    final=dict(schema='library_preprocess_dataset_v1',dataset=args.dataset,start_time=started,end_time=now(),
               status='success' if all(r['status']=='success' for r in results) else 'partial_failure',
               n_libraries=len(results),n_input_QC=a.n_obs,n_output=sum(r.get('n_output',0) for r in results),
               n_additional_doublets_removed=sum(r.get('n_additional_removed',0) for r in results),
               pilot_only=bool(args.library_limit),source_final_AC=False,cross_study_integrated=False,libraries=results)
    save(dest/'dataset_summary.json',final)
    return 0 if final['status']=='success' else 1


def main():
    p=argparse.ArgumentParser();p.add_argument('--collection',required=True);p.add_argument('--output',required=True)
    p.add_argument('--dataset',required=True);p.add_argument('--rscript',required=True);p.add_argument('--rworker',required=True)
    p.add_argument('--library-limit',type=int,default=0,help='Real library pilot only; never claims full coverage')
    p.add_argument('--workers',type=int,default=1)
    p.add_argument('--native-object',help='Audited native-count Zarr needing initial per-library QC')
    args=p.parse_args()
    event(args.dataset,'environment',python=sys.version.split()[0],versions={name:importlib.metadata.version(name)
          for name in ['scanpy','anndata','numpy','scipy','doubletdetection']},script_sha=sha(__file__))
    try:
        return dataset_job(args)
    except Exception as exc:
        output=Path(args.output); output.mkdir(parents=True,exist_ok=True)
        save(output/('fatal_'+args.dataset+'.json'),dict(utc=now(),error=repr(exc),traceback=traceback.format_exc()))
        event(args.dataset,'dataset_fatal',error=repr(exc));return 1


if __name__=='__main__':
    sys.exit(main())
