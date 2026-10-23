"""Audit graph factors and replay conditional RNA diagnostics without refitting.

All aggregates are reproduced; the smallest complete pilot source is additionally
checked with independent centered OLS and direct row-neighbor calculations.
"""
import argparse
import datetime
import gzip
import json
from pathlib import Path
import resource
import time
import numpy as np
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score

from annotation_representation_pilot_v1 import save
from annotation_graph_diagnostics_v1 import checked_graph,unit_out_graph_state,direct_unit_out_graph_state
from compare_lineage_graph_engines_v2 import evaluate_graphs,qualify_expression_units
from paired_granularity_review_v1 import compare_json_numeric
from full_context_atlas_v1 import FullContextAtlas,sha256_file
from load_annotation_pilot_native_v1 import load_pilot


def main(root,producer):
    started=time.monotonic();cfg=json.loads((producer/'config/graph_comparison_spec_v2.json').read_text())
    sup=json.loads((producer/'provenance/supervisor_terminal_v1.json').read_text())
    if sup['actual_child_exit_code']!=0 or sup['wall_timeout']:raise ValueError('Graph producer not complete')
    term=json.loads((producer/'provenance/terminal_v1.json').read_text())
    index_path=producer/'outputs/graph_comparison_index_v1.json'
    if sha256_file(index_path)!=term['index_sha256']:raise ValueError('Graph output index changed')
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Graph input changed before audit')
    def checkpoint(stage,**extra):
        save(root/'checkpoints'/(str(time.time_ns())+'.json'),{'observed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'stage':stage,**extra})
    data=load_pilot(cfg['pilot'],cfg['manifest'],checkpoint=checkpoint)
    sl=Path(cfg['source_lineage']);confidence=Path(cfg['confidence']);pilot=Path(cfg['pilot'])
    with np.load(sl/'outputs/derived_axis_and_model_mask.npz',allow_pickle=False) as f:
        route,positive=f['provisional_native_route'],f['modeled'];np.testing.assert_array_equal(data['source'],f['source'])
    with np.load(sl/'outputs/candidate_partitions_2026102301.npz',allow_pickle=False) as f:oldc=f['C']
    inv=json.loads((pilot/'outputs/pilot_input_inventory_v1.json').read_text())
    smallest=min(inv,key=lambda r:(r['n_selected_cells'],r['dataset']))['dataset']
    ci=json.loads((confidence/'outputs/cell_identity_confidence_index_v2.json').read_text())
    reader=FullContextAtlas(cfg['manifest']);units=[];direct=[];expected_axes=[]
    for v,e in zip(inv,ci):
        ds=v['dataset'];positions=np.load(pilot/v['selector_file'],allow_pickle=False)
        units.extend(qualify_expression_units(ds,reader.frozen_expression_units(ds)[positions]).tolist())
        meta=confidence/'outputs'/e['metadata'];side=confidence/'outputs'/e['sidecar']
        if sha256_file(meta)!=e['metadata_sha256'] or sha256_file(side)!=e['sidecar_sha256']:raise ValueError('Confidence sidecar changed')
        with gzip.open(meta,'rt') as f:md=json.load(f)
        with np.load(side,allow_pickle=False) as f:
            np.testing.assert_array_equal(f['native_row'],reader.rows(ds)[positions])
            direct.extend(np.array(md['categorical_dictionaries']['C__direct_fine_candidate'])[f['C__direct_fine_candidate_code']].tolist())
        expected_axes.append({'source':ds,'n_cells':len(positions),'selector_sha256':v['selector_sha256'],'cell_axis_sha256':e['cell_axis_sha256']})
    del reader
    if expected_axes!=json.loads((producer/'outputs/pilot_axis_binding_v1.json').read_text()):raise ValueError('Graph source axis binding differs')
    units=np.asarray(units);direct=np.asarray(direct);source=data['source']
    features=json.loads((producer/'outputs/feature_factor_selection_v1.json').read_text())
    feature_by={r['route']:r for r in features}
    fits=json.loads((producer/'outputs/graph_fit_resource_records_v1.json').read_text())
    fit_by={(r['route'],r['variant']):r for r in fits if r['status']=='fit_complete'}
    index=json.loads(index_path.read_text());reports=[];n_cases=0;max_error=0.;replayed=0
    for entry in index:
        parent=entry['route'];ap=producer/'outputs'/entry['axes'];pp=producer/'outputs'/entry['profile']
        if sha256_file(ap)!=entry['axes_sha256'] or sha256_file(pp)!=entry['profile_sha256']:raise ValueError('Graph profile or axis changed')
        with np.load(ap,allow_pickle=False) as f:axis={k:f[k] for k in f.files}
        expected=np.flatnonzero(positive&(route==parent));idx=axis['pilot_row']
        np.testing.assert_array_equal(idx,expected)
        for name,whole in [('source',source),('expression_unit',units),('direct_candidate',direct)]:np.testing.assert_array_equal(axis[name],whole[idx])
        fr=feature_by[parent];fi=np.asarray(fr['lineage_feature_indices'])
        if len(fi)!=cfg['number_features'] or len(set(fi))!=len(fi) or [data['genes'][k] for k in fi]!=fr['lineage_feature_names']:
            raise ValueError('Lineage feature axis differs')
        # load_pilot validates measured genes and physical heldout aliases in every source.
        mp=producer/'models'/fr['PCA_file']
        if sha256_file(mp)!=fr['PCA_sha256']:raise ValueError('Saved lineage PCA changed')
        with np.load(mp,allow_pickle=False) as f:
            np.testing.assert_array_equal(f['features'],fi)
            if f['mean'].shape!=(len(fi),) or f['components'].shape[1]!=len(fi):raise ValueError('PCA feature axis malformed')
        graphs={}
        for variant in entry['variants']:
            fit=fit_by[parent,variant];gp=producer/'graphs'/fit['graph_file']
            if sha256_file(gp)!=fit['graph_sha256']:raise ValueError('Accepted graph changed')
            raw=sp.load_npz(gp);w=checked_graph(raw)
            if raw.shape!=(len(idx),len(idx)) or raw.nnz!=fit['graph_nnz']:raise ValueError('Graph dimensions differ')
            for field in ['data','indices','indptr']:np.testing.assert_array_equal(getattr(raw,field),getattr(w,field))
            delta=w-w.T
            if delta.nnz and np.max(np.abs(delta.data))>1e-6:raise ValueError('Stored graph not symmetric')
            if variant=='shared_Harmony' and adjusted_rand_score(axis['partition__'+variant],oldc[idx])<1-1e-12:
                raise ValueError('Frozen C identity partition reproduction differs')
            graphs[variant]=w
        values={k:v[idx] for k,v in data['state'].items()};s=source[idx];u=units[idx];d=np.log1p(data['native_total'][idx]);labels=direct[idx]
        replay=evaluate_graphs(graphs,s,u,d,values,labels,data['gamma_RNA_candidate'][idx],cfg,checkpoint=checkpoint)
        with gzip.open(pp,'rt') as f:original=json.load(f)
        for key,value in replay.items():
            if isinstance(value,list):
                if len(value)!=len(original[key]):raise ValueError('Diagnostic record count differs')
                for before,after in zip(original[key],value):compare_json_numeric(before,after)
                replayed+=len(value)
            else:compare_json_numeric(original[key],value)
        scopes=[('all_route_cells_not_pure_type',s==smallest)]
        scopes.extend((str(label),(s==smallest)&(labels==label)) for label in np.unique(labels[s==smallest]) if label!='unresolved')
        local_cases=0
        for scope,mask in scopes:
            if mask.sum()<cfg['minimum_stratum_cells']:continue
            q=np.flatnonzero(mask)
            for variant,graph in graphs.items():
                sub=graph[q][:,q]
                for state,y in values.items():
                    first=unit_out_graph_state(sub,u[q],d[q],y[q],cfg['minimum_other_unit_neighbors'],cfg['minimum_train_cells'],True)
                    p,b=direct_unit_out_graph_state(sub,u[q],d[q],y[q],cfg['minimum_other_unit_neighbors'],cfg['minimum_train_cells'])
                    for a,z in [(first['prediction'],p),(first['baseline'],b)]:
                        np.testing.assert_array_equal(np.isfinite(a),np.isfinite(z))
                        ok=np.isfinite(a)
                        np.testing.assert_allclose(a[ok],z[ok],atol=1e-9,rtol=1e-8)
                        if ok.any():max_error=max(max_error,float(np.max(np.abs(a[ok]-z[ok]))))
                    n_cases+=1;local_cases+=1
        reports.append({'route':parent,'all_graph_axes_features_readbacks_and_profile_records_passed':True,
            'independent_prediction_cases':local_cases,'independent_source':smallest})
        checkpoint('route_graph_independent_audit_complete',route=parent)
        del graphs
    expected_rare=[]
    for ds in np.unique(source):
        measured=(source==ds)&np.isfinite(data['gamma_RNA_candidate']);rare=measured&(data['gamma_RNA_candidate']>0)
        expected_rare.append({'source':str(ds),'n_measured_all_pilot':int(measured.sum()),'n_RNA_candidates_all_pilot':int(rare.sum()),
            'n_RNA_candidates_in_requested_routes':int((rare&np.isin(route,cfg['routes'])&positive).sum()),
            'n_RNA_candidates_outside_requested_routes':int((rare&~(np.isin(route,cfg['routes'])&positive)).sum())})
    if expected_rare!=json.loads((producer/'outputs/gamma_RNA_route_coverage_v1.json').read_text()):raise ValueError('Rare candidate full denominator differs')
    out=root/'outputs/graph_comparison_independent_audit_v1.json'
    save(out,{'route_reports':reports,'replayed_aggregate_records':replayed,'independent_prediction_cases':n_cases,
        'independent_prediction_max_absolute_error':max_error,'numeric_source_selection':'smallest_frozen_pilot_source',
        'all_rare_RNA_route_denominators_checked':True,'native_features_and_heldout_physical_aliases_validated':True,
        'graph_rebuild_not_claimed':True,'new_source_holdout_not_claimed':True,'biological_accuracy_not_claimed':True})
    save(root/'provenance/terminal_v1.json',{'status':'LINEAGE_GRAPH_FACTOR_AXES_FULL_AGGREGATES_AND_SELECTED_DIRECT_NEIGHBOR_REPLAY_ACCEPTED',
        'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'output_sha256':sha256_file(out),
        'replayed_aggregate_records':replayed,'independent_prediction_cases':n_cases,'independent_prediction_max_absolute_error':max_error,
        'wall_seconds':time.monotonic()-started,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        'no_graph_or_model_refit':True,'native_counts_copied':False,'new_authoritative_labels':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--producer',required=True);args=p.parse_args()
    root=Path(args.output)
    try:main(root,Path(args.producer))
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)})
        raise
