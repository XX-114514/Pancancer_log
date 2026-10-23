"""Frozen feature factor and graph-engine comparison on existing pilot cells.

Counts remain in the native collection. Stored graphs/models are private bounded
sidecars. No clinical outcomes, final labels, or engine winner are generated.
"""
import argparse
import datetime
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score

from annotation_representation_pilot_v1 import select_features,pca_fit,harmonize,save
from annotation_graph_diagnostics_v1 import (checked_graph,unit_out_graph_state,
    identity_neighbor_profiles,rare_neighbor_profiles)
from paired_granularity_review_v1 import paired_scores
from annotation_granularity_metrics_v1 import granularity
from full_context_atlas_v1 import FullContextAtlas,sha256_file
from load_annotation_pilot_native_v1 import load_pilot


VARIANTS=['shared_Harmony','lineage_uncorrected','lineage_Harmony','lineage_BBKNN']
CONTRASTS=[('shared_Harmony','lineage_Harmony'),('lineage_uncorrected','lineage_Harmony'),
           ('lineage_Harmony','lineage_BBKNN')]


def qualify_expression_units(source, values):
    # NumPy 1.x does not implement scalar string + Unicode ndarray via add.
    # Explicit string construction preserves leading zeros and source namespace.
    return np.asarray([str(source)+'::'+str(value) for value in values],dtype=str)


def build_graph(embedding,source,engine,seed,parameters):
    import anndata as ad
    import scanpy as sc
    import numba
    sc.settings.n_jobs=parameters['threads'];numba.set_num_threads(parameters['threads'])
    a=ad.AnnData(sp.csr_matrix((len(embedding),0)));a.obsm['X_pilot']=embedding
    if engine=='BBKNN':
        import bbknn.matrix
        _,counts=np.unique(source,return_counts=True)
        if len(counts)<2 or counts.min()<10:raise ValueError('BBKNN source coverage unqualified')
        distance,graph,params=bbknn.matrix.bbknn(embedding,np.asarray(source),neighbors_within_batch=5,
            n_pcs=embedding.shape[1],trim=50*len(counts),computation='pynndescent',
            pynndescent_n_neighbors=30,pynndescent_random_state=seed,metric='euclidean',
            set_op_mix_ratio=1.,local_connectivity=1.)
        a.obsp['distances']=distance;a.obsp['connectivities']=graph
        a.uns['neighbors']={'distances_key':'distances','connectivities_key':'connectivities','params':params}
    else:
        sc.pp.neighbors(a,use_rep='X_pilot',n_neighbors=min(parameters['neighbors'],len(embedding)-1),
            random_state=seed,transformer='pynndescent',metric='euclidean')
    sc.tl.leiden(a,resolution=1.,random_state=seed,flavor='igraph',n_iterations=2,directed=False)
    _,labels=np.unique(a.obs['leiden'].to_numpy().astype(str),return_inverse=True)
    w=checked_graph(a.obsp['connectivities'])
    difference=w-w.T
    if difference.nnz and np.max(np.abs(difference.data))>1e-6:raise ValueError('Expected undirected graph is asymmetric')
    degree=np.diff(w.indptr)
    info={'engine':engine,'n_cells':len(embedding),'graph_nnz':w.nnz,
          'actual_degree_min':int(degree.min()),'actual_degree_median':float(np.median(degree)),
          'actual_degree_max':int(degree.max()),'n_isolated':int((degree==0).sum()),
          'actual_csr_uncompressed_bytes':w.data.nbytes+w.indices.nbytes+w.indptr.nbytes,
          'symmetric_within_1e6':True,'self_edges_removed_for_diagnostics':True}
    return w,labels.astype(np.int32),info


def actual_engine_smoke(seed,parameters):
    rng=np.random.default_rng(seed);x=rng.normal(size=(160,8)).astype(np.float32)
    source=np.repeat(['synthetic_a','synthetic_b'],80)
    records=[]
    for engine in ['Scanpy','BBKNN']:
        graph,labels,meta=build_graph(x,source,engine,seed,parameters)
        if graph.shape!=(160,160) or len(labels)!=160 or graph.nnz==0:raise ValueError('Actual graph API smoke failed')
        records.append(meta)
    return records


def evaluate_graphs(graphs,source,units,depth,states,direct,gamma,cfg,checkpoint=None):
    records=[];paired=[];identity=[];rare=[];groups=[]
    for variant,graph in graphs.items():
        identity.extend(dict(variant=variant,**r) for r in identity_neighbor_profiles(graph,source,direct,cfg['minimum_shared_type_cells']))
        rare.extend(dict(variant=variant,**r) for r in rare_neighbor_profiles(graph,source,gamma))
    for ds in np.unique(source):
        scopes=[('all_route_cells_not_pure_type',source==ds)]
        scopes.extend((str(label),(source==ds)&(direct==label)) for label in np.unique(direct[source==ds]) if label!='unresolved')
        for identity_scope,mask in scopes:
            groups.append({'source':str(ds),'identity_scope':identity_scope,'n_cells':int(mask.sum()),
                'status':'evaluated' if mask.sum()>=cfg['minimum_stratum_cells'] else 'below_minimum_cells'})
            if mask.sum()<cfg['minimum_stratum_cells']:continue
            indices=np.flatnonzero(mask);u=units[indices];d=depth[indices]
            cache={name:{} for name in graphs}
            for variant,graph in graphs.items():
                sub=graph[indices][:,indices]
                for state,values in states.items():
                    result=unit_out_graph_state(sub,u,d,values[indices],cfg['minimum_other_unit_neighbors'],
                        cfg['minimum_train_cells'],return_predictions=True)
                    cache[variant][state]=result
                    records.append({'source':str(ds),'identity_scope':identity_scope,'variant':variant,'state':state,
                        **{k:v for k,v in result.items() if k not in ['prediction','baseline']}})
            for first,second in CONTRASTS:
                if first not in cache or second not in cache:continue
                for state,values in states.items():
                    paired.append({'source':str(ds),'identity_scope':identity_scope,'state':state,
                        'first_variant':first,'second_variant':second,
                        **paired_scores(values[indices],u,cache[first][state],cache[second][state])})
            if checkpoint:checkpoint('graph_source_identity_diagnostics_complete',source=str(ds),identity_scope=identity_scope)
    return {'state_profiles':records,'paired_state_profiles':paired,'identity_neighbor_profiles':identity,
            'rare_RNA_neighbor_profiles':rare,'stratum_coverage':groups}


def main(root):
    start=time.monotonic();cfg=json.loads((root/'config/graph_comparison_spec_v2.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Graph comparison input changed')
    for path in cfg['accepted_dependencies']:
        sup=json.loads((Path(path)/'provenance/supervisor_terminal_v1.json').read_text())
        if sup['actual_child_exit_code']!=0 or sup['wall_timeout']:raise ValueError('Graph comparison dependency incomplete')
    def checkpoint(stage,**extra):
        save(root/'checkpoints'/(str(time.time_ns())+'.json'),{'observed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'stage':stage,'wall_seconds':time.monotonic()-start,**extra})
    t=subprocess.run([sys.executable,'-B',str(root/'tests/test_annotation_graph_diagnostics_v2.py')],capture_output=True,text=True,timeout=180)
    (root/'provenance/test_stdout.txt').write_text(t.stdout);(root/'provenance/test_stderr.txt').write_text(t.stderr)
    if t.returncode!=0 or 'Ran 8 tests' not in t.stderr or 'skipped' in t.stderr.lower():raise ValueError('Graph diagnostic counterexamples failed')
    params=json.loads((Path(cfg['discovery'])/'config/model_spec_v1.json').read_text())
    smoke=actual_engine_smoke(cfg['seed'],params)
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':8,'no_skips':True,'actual_graph_API_smoke':smoke})
    checkpoint('numerical_and_actual_API_tests_passed')
    pilot=Path(cfg['pilot']);sl=Path(cfg['source_lineage']);confidence=Path(cfg['confidence'])
    data=load_pilot(pilot,cfg['manifest'],checkpoint=checkpoint)
    with np.load(sl/'outputs/derived_axis_and_model_mask.npz',allow_pickle=False) as f:
        source,route,positive=f['source'],f['provisional_native_route'],f['modeled']
        np.testing.assert_array_equal(source,data['source']);np.testing.assert_array_equal(positive,data['native_total']>0)
    with np.load(sl/'outputs/candidate_partitions_2026102301.npz',allow_pickle=False) as f:oldc=f['C']
    inventory=json.loads((pilot/'outputs/pilot_input_inventory_v1.json').read_text())
    confidence_index=json.loads((confidence/'outputs/cell_identity_confidence_index_v2.json').read_text())
    reader=FullContextAtlas(cfg['manifest']);units=[];direct=[];axis_records=[];offset=0
    for inv,entry in zip(inventory,confidence_index):
        ds=inv['dataset'];n=inv['n_selected_cells'];positions=np.load(pilot/inv['selector_file'],allow_pickle=False)
        if entry['source']!=ds or entry['n_pilot_cells']!=n or not np.all(source[offset:offset+n]==ds):raise ValueError('Confidence source axes differ')
        side=confidence/'outputs'/entry['sidecar'];meta=confidence/'outputs'/entry['metadata']
        if sha256_file(side)!=entry['sidecar_sha256'] or sha256_file(meta)!=entry['metadata_sha256']:raise ValueError('Accepted confidence sidecar changed')
        with gzip.open(meta,'rt') as f:md=json.load(f)
        with np.load(side,allow_pickle=False) as f:
            np.testing.assert_array_equal(f['native_row'],reader.rows(ds)[positions])
            labels=np.array(md['categorical_dictionaries']['C__direct_fine_candidate'])[f['C__direct_fine_candidate_code']]
            for method in ['A','B']:
                other=np.array(md['categorical_dictionaries'][method+'__direct_fine_candidate'])[f[method+'__direct_fine_candidate_code']]
                np.testing.assert_array_equal(labels,other)
        units.extend(qualify_expression_units(ds,reader.frozen_expression_units(ds)[positions]).tolist());direct.extend(labels.tolist())
        axis_records.append({'source':ds,'n_cells':n,'selector_sha256':inv['selector_sha256'],'cell_axis_sha256':entry['cell_axis_sha256']});offset+=n
    if offset!=len(source):raise ValueError('Confidence records incomplete')
    del reader
    units=np.asarray(units);direct=np.asarray(direct)
    shared=json.loads((Path(cfg['discovery'])/'outputs/discovery_feature_selection_v1.json').read_text())
    shared_idx=np.asarray(shared['selected_feature_indices'],np.int64)
    if [data['genes'][i] for i in shared_idx]!=shared['selected_feature_names']:raise ValueError('Shared feature axis changed')
    save(root/'outputs/pilot_axis_binding_v1.json',axis_records)
    output_index=[];features=[];fit_records=[];all_profiles=[];rare_coverage=[]
    for ds in np.unique(source):
        measured=(source==ds)&np.isfinite(data['gamma_RNA_candidate']);rare=measured&(data['gamma_RNA_candidate']>0)
        rare_coverage.append({'source':str(ds),'n_measured_all_pilot':int(measured.sum()),'n_RNA_candidates_all_pilot':int(rare.sum()),
            'n_RNA_candidates_in_requested_routes':int((rare&np.isin(route,cfg['routes'])&positive).sum()),
            'n_RNA_candidates_outside_requested_routes':int((rare&~(np.isin(route,cfg['routes'])&positive)).sum())})
    for parent in cfg['routes']:
        idx=np.flatnonzero(positive&(route==parent));s=source[idx];u=units[idx];df=direct[idx]
        if len(idx)<params['minimum_fit_group']:
            fit_records.append({'route':parent,'status':'below_minimum_fit_cells'});continue
        lineage_idx=select_features(data['x'][idx],s,cfg['number_features'])
        raw,pca=pca_fit(data['x'][idx][:,lineage_idx],cfg['seed'],params['components'],params['pca_iterations'])
        feature_file=root/'models'/(parent+'_lineage_PCA_v1.npz')
        with feature_file.open('xb') as f:np.savez_compressed(f,mean=pca.mean_,components=pca.components_,features=lineage_idx)
        features.append({'route':parent,'lineage_feature_indices':lineage_idx.tolist(),'lineage_feature_names':[data['genes'][i] for i in lineage_idx],
            'shared_feature_overlap':len(set(lineage_idx)&set(shared_idx)),'feature_estimator':'mean_within_source_log_dispersion_rank_not_Seurat_HVG',
            'PCA_file':feature_file.name,'PCA_sha256':sha256_file(feature_file),'all_pilot_sources_used':True})
        raw_sha=hashlib.sha256(raw.tobytes()).hexdigest()
        graphs={};partitions={}
        for variant in VARIANTS:
            began=time.monotonic()
            if variant=='shared_Harmony':
                shared_raw,_=pca_fit(data['x'][idx][:,shared_idx],cfg['seed'],params['components'],params['pca_iterations'])
                embedding,info=harmonize(shared_raw,s,cfg['seed'],params);del shared_raw
                engine='Scanpy'
            elif variant=='lineage_Harmony':embedding,info=harmonize(raw.copy(),s,cfg['seed'],params);engine='Scanpy'
            elif variant=='lineage_uncorrected':embedding=raw.copy();info={'status':'uncorrected_PCA_diagnostic'};engine='Scanpy'
            else:
                _,count=np.unique(s,return_counts=True)
                if len(count)<2 or count.min()<10:
                    fit_records.append({'route':parent,'variant':variant,'status':'BBKNN_source_coverage_unqualified_no_source_dropped'});continue
                embedding=raw.copy();info={'status':'BBKNN_graph_not_corrected_embedding_no_out_of_sample_transform_claim'};engine='BBKNN'
            graph,labels,meta=build_graph(embedding,s,engine,cfg['seed'],params)
            if hashlib.sha256(raw.tobytes()).hexdigest()!=raw_sha:raise ValueError('Shared lineage PCA mutated across engines')
            if variant=='shared_Harmony':
                ari=float(adjusted_rand_score(oldc[idx],labels))
                if ari<1-1e-12:raise ValueError('Shared Harmony C partition not reproduced')
                meta['ARI_with_frozen_C_resolution1']=ari
            path=root/'graphs'/(parent+'__'+variant+'.npz')
            with path.open('xb') as f:sp.save_npz(f,graph,compressed=True)
            back=sp.load_npz(path)
            for key in ['indptr','indices','data']:np.testing.assert_array_equal(getattr(back,key),getattr(graph,key))
            graphs[variant]=graph;partitions[variant]=labels
            record={'route':parent,'variant':variant,'status':'fit_complete','wall_seconds':time.monotonic()-began,
                'graph_file':path.name,'graph_sha256':sha256_file(path),'graph_compressed_bytes':path.stat().st_size,
                'embedding_processing':info,'granularity':granularity(labels),**meta}
            fit_records.append(record);checkpoint('graph_variant_complete',route=parent,variant=variant)
            del embedding
        axes=root/'outputs'/(parent+'_axis_and_partitions_v1.npz')
        with axes.open('xb') as f:np.savez_compressed(f,pilot_row=idx,source=s,expression_unit=u,direct_candidate=df,
            **{'partition__'+k:v for k,v in partitions.items()})
        profile=evaluate_graphs(graphs,s,u,np.log1p(data['native_total'][idx]),{k:v[idx] for k,v in data['state'].items()},
            df,data['gamma_RNA_candidate'][idx],cfg,checkpoint=checkpoint)
        profile.update(route=parent,state_gene_sets=data['state_genes'],native_RNA_unchanged=True,
            conditional_reference_and_routing_may_encode_states=True,new_authoritative_labels=False)
        path=root/'outputs'/(parent+'_graph_profiles_v1.json.gz')
        with gzip.open(path,'xt') as f:json.dump(profile,f,allow_nan=False)
        output_index.append({'route':parent,'n_pilot_cells':len(idx),'axes':axes.name,'axes_sha256':sha256_file(axes),
            'profile':path.name,'profile_sha256':sha256_file(path),'variants':list(graphs)})
        del graphs,partitions,raw,pca
        size=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
        if size>cfg['output_cap_bytes']:raise ValueError('Graph output cap exceeded')
        checkpoint('route_graph_comparison_complete',route=parent,written_bytes=size)
    save(root/'outputs/feature_factor_selection_v1.json',features)
    save(root/'outputs/graph_fit_resource_records_v1.json',fit_records)
    save(root/'outputs/gamma_RNA_route_coverage_v1.json',rare_coverage)
    index=root/'outputs/graph_comparison_index_v1.json';save(index,output_index)
    size=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    if size>cfg['output_cap_bytes']:raise ValueError('Graph output cap exceeded')
    save(root/'provenance/terminal_v1.json',{'status':'LINEAGE_FEATURE_AND_GRAPH_FACTOR_COMPARISON_COMPLETE_INDEPENDENT_AUDIT_PENDING',
        'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_tests_actual_no_skips':8,
        'n_route_branches':len(output_index),'index_sha256':sha256_file(index),'wall_seconds':time.monotonic()-start,
        'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'written_bytes_before_terminal':size,
        'all_shared_Harmony_C_partitions_reproduced':True,'native_counts_copied':False,
        'engine_winner_selected':False,'new_authoritative_labels':False,'source_holdout_replaced':False,
        'versions':{n:importlib.metadata.version(n) for n in ['bbknn','scanpy','numpy','scipy','pynndescent','scikit-learn']}})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)})
        raise
