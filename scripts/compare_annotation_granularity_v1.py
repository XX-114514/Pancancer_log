"""Frozen-resolution curves on unchanged A/B/C pilot graph definitions.

All-source discovery is transductive. Unit-out RNA diagnostics exclude query
units from depth/group prediction, but do not make discovery source-heldout.
Native counts and graphs remain in RAM; only small labels/profiles are saved.
"""
import argparse
import datetime
import gzip
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score

from annotation_representation_pilot_v1 import pca_fit,harmonize,save
from annotation_granularity_metrics_v1 import granularity,unit_out_state,gamma_group_retention
from full_context_atlas_v1 import FullContextAtlas,sha256_file
from load_annotation_pilot_native_v1 import load_pilot


def graph_grid(embedding,seed,parameters,resolutions):
    import anndata as ad
    import scanpy as sc
    import numba
    sc.settings.n_jobs=parameters['threads'];numba.set_num_threads(parameters['threads'])
    a=ad.AnnData(sp.csr_matrix((len(embedding),0)));a.obsm['X_pilot']=embedding
    sc.pp.neighbors(a,use_rep='X_pilot',n_neighbors=min(parameters['neighbors'],len(embedding)-1),
                    random_state=seed,transformer='pynndescent',metric='euclidean')
    out={}
    for resolution in resolutions:
        sc.tl.leiden(a,resolution=resolution,random_state=seed,flavor='igraph',n_iterations=2,directed=False)
        _,codes=np.unique(a.obs['leiden'].to_numpy().astype(str),return_inverse=True)
        out[str(resolution)]=codes.astype(np.int32)
    return out


def main(root):
    start=time.monotonic();cfg=json.loads((root/'config/granularity_spec_v1.json').read_text())
    for item in cfg['bindings'].values():
        if sha256_file(item['path'])!=item['sha256']:raise ValueError('Frozen granularity input changed')
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_annotation_granularity_metrics_v1.py')],capture_output=True,text=True,timeout=300)
    (root/'provenance/test_stdout.txt').write_text(test.stdout);(root/'provenance/test_stderr.txt').write_text(test.stderr)
    ok=test.returncode==0 and 'Ran 6 tests' in test.stderr and 'skipped' not in test.stderr.lower()
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':test.returncode,'tests':6,'no_skips':ok})
    if not ok:raise ValueError('Actual granularity counterexamples did not pass')
    def checkpoint(stage,**extra):
        save(root/'checkpoints'/(str(time.time_ns())+'.json'),{'observed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'stage':stage,'wall_seconds':time.monotonic()-start,**extra})
    data=load_pilot(cfg['pilot'],cfg['manifest'],checkpoint=checkpoint)
    discovery=Path(cfg['discovery']);source,route=data['source'],data['route'];positive=data['native_total']>0
    oldcfg=json.loads((discovery/'config/model_spec_v1.json').read_text())
    selected=json.loads((discovery/'outputs/discovery_feature_selection_v1.json').read_text())
    hvg=np.asarray(selected['selected_feature_indices'],np.int64)
    if [data['genes'][i] for i in hvg]!=selected['selected_feature_names']:raise ValueError('Model feature axis changed')
    with np.load(discovery/'outputs/derived_axis_and_model_mask.npz',allow_pickle=False) as f:
        if not np.array_equal(source,f['source']) or not np.array_equal(route,f['inherited_route']) or not np.array_equal(positive,f['modeled']):
            raise ValueError('Frozen graph input axes changed')
    x=data.pop('x')[:,hvg];seed=cfg['seed'];resolutions=cfg['resolutions']
    reference=FullContextAtlas(cfg['manifest'])
    inv=json.loads((Path(cfg['pilot'])/'outputs/pilot_input_inventory_v1.json').read_text());units=[]
    for row in inv:
        ds=row['dataset'];positions=np.load(Path(cfg['pilot'])/row['selector_file'],allow_pickle=False)
        units.extend(reference.frozen_expression_units(ds)[positions].astype(str))
    units=np.asarray(units)
    if len(units)!=len(source):raise ValueError('Expression-unit axis differs')
    # Evaluation is within source, so source-specific unit namespaces never mix.
    del reference
    partitions={m:{str(r):np.full(len(source),-1,np.int32) for r in resolutions} for m in ['A','B','C']}
    global_raw,pca=pca_fit(x[positive],seed,oldcfg['components'],oldcfg['pca_iterations'])
    with np.load(discovery/'models'/('global_discovery_PCA_'+str(seed)+'.npz'),allow_pickle=False) as f:
        np.testing.assert_allclose(pca.mean_,f['mean'],atol=1e-6,rtol=1e-6)
        np.testing.assert_allclose(pca.components_,f['components'],atol=1e-5,rtol=1e-5)
    corrected,_=harmonize(global_raw,source[positive],seed,oldcfg)
    for resolution,labels in graph_grid(corrected,seed,oldcfg,resolutions).items():partitions['B'][resolution][positive]=labels
    checkpoint('B_same_graph_resolution_grid_complete')
    del global_raw,corrected,pca
    for method,grouping,correct in [('A',source,False),('C',route,True)]:
        offsets={str(r):0 for r in resolutions}
        for group in np.unique(grouping):
            member=np.flatnonzero(positive & (grouping==group))
            if len(member)<oldcfg['minimum_fit_group']:continue
            embedding,_=pca_fit(x[member],seed,oldcfg['components'],oldcfg['pca_iterations'])
            if correct:embedding,_=harmonize(embedding,source[member],seed,oldcfg)
            for resolution,labels in graph_grid(embedding,seed,oldcfg,resolutions).items():
                partitions[method][resolution][member]=labels+offsets[resolution]
                offsets[resolution]+=len(np.unique(labels))
            checkpoint('group_same_graph_resolution_grid_complete',method=method,group=str(group),n_cells=len(member))
    del x
    reproductions={}
    with np.load(discovery/'outputs'/('candidate_partitions_'+str(seed)+'.npz'),allow_pickle=False) as f:
        for method in ['A','B','C']:
            a,b=partitions[method]['1.0'],f[method]
            if not np.array_equal(a>=0,b>=0):raise ValueError('Original partition eligibility not reproduced')
            ari=float(adjusted_rand_score(a[a>=0],b[b>=0]))
            reproductions[method]={'adjusted_Rand_with_original_resolution1':ari}
            if ari<1-1e-12:raise ValueError('Original resolution1 graph partition not reproduced')
    save(root/'outputs/original_partition_reproduction_v1.json',reproductions)
    side={m+'__'+r:a for m,v in partitions.items() for r,a in v.items()}
    with (root/'outputs/granularity_partitions_v1.npz').open('xb') as f:np.savez_compressed(f,**side)
    profiles=[];overlaps=[]
    for ds in np.unique(source):
        for parent in np.unique(route[source==ds]):
            member=(source==ds)&(route==parent)&positive
            if member.sum()<cfg['minimum_stratum_cells']:continue
            yy={k:v[member] for k,v in data['state'].items()}
            stratum_profiles={}
            for method,variants in partitions.items():
                for resolution,labels in variants.items():
                    lab=labels[member];profile=granularity(lab)
                    record={'source':str(ds),'inherited_parent_not_truth':str(parent),'method':method,'resolution':float(resolution),
                            'granularity':profile,'state_prediction':{},'n_expression_units_not_patients':len(np.unique(units[member])),
                            'gamma_RNA_candidate':gamma_group_retention(lab,data['gamma_RNA_candidate'][member])}
                    for name,values in yy.items():
                        record['state_prediction'][name]=unit_out_state(lab,units[member],np.log1p(data['native_total'][member]),values,
                                                                       minimum_other_cells=cfg['minimum_other_unit_group_cells'])
                    profiles.append(record);stratum_profiles[(method,resolution)]=profile
            for a,b in [('A','B'),('A','C'),('B','C')]:
                matches=[]
                for ra in resolutions:
                    for rb in resolutions:
                        pa,pb=stratum_profiles[(a,str(ra))],stratum_profiles[(b,str(rb))]
                        if pa['status']!='estimated' or pb['status']!='estimated':continue
                        ratio=max(pa['effective_group_count'],pb['effective_group_count'])/min(pa['effective_group_count'],pb['effective_group_count'])
                        if ratio<=cfg['maximum_effective_group_ratio']:
                            matches.append({'resolution_'+a:ra,'resolution_'+b:rb,'effective_group_count_ratio':ratio})
                overlaps.append({'source':str(ds),'inherited_parent_not_truth':str(parent),'methods':[a,b],
                                 'status':'granularity_overlap_only_not_biological_equivalence' if matches else 'no_overlap_in_frozen_grid',
                                 'overlapping_grid_pairs':matches,'outcome_used_for_matching':False})
            checkpoint('source_parent_profiles_complete',source=str(ds),inherited_parent=str(parent))
    payload={'profiles':profiles,'granularity_overlap':overlaps,
             'state_gene_sets':data['state_genes'],'discovery_all_sources':True,'new_method_winner_selected':False,
             'source_unit_independence_not_assumed':True,'C_inherited_route_not_new_identity':True,
             'marker_readouts_excluded_from_model_but_prior_labels_may_encode_them':True,
             'first_seed_only_grid_three_seed_original_discovery_preserved':True}
    path=root/'outputs/source_parent_granularity_profiles_v1.json.gz'
    with gzip.open(path,'xt') as f:json.dump(payload,f,allow_nan=False)
    size=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    if size>cfg['output_cap_bytes']:raise ValueError('Granularity output cap exceeded')
    save(root/'provenance/terminal_v1.json',{'status':'SAME_INPUT_GRAPH_RESOLUTION_CURVES_COMPLETE_NOT_METHOD_SELECTION',
         'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'original_resolution1_reproductions':reproductions,
         'output_sha256':sha256_file(path),'wall_seconds':time.monotonic()-start,
         'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
         'written_bytes_before_terminal':size,'n_profile_records':len(profiles),'n_tests_actual_no_skips':6,
         'native_counts_copied':False,'new_authoritative_labels':False,'source_holdout_replaced':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)})
        raise
