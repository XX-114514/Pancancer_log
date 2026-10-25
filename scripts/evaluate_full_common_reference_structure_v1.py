"""Paired, source-conditioned native RNA diagnostics for common model candidates.

All-source reference discovery is transductive; guarded types are conditioning
strata, not accuracy targets. These diagnostics do not assign final identities.
"""
import argparse,datetime,gc,json,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
from annotation_pilot_sampling_v1 import balanced_cell_selection,stable_seed
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from common_reference_structure_metrics_v1 import paired_native_states,exact_neighbor_graph
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import row_axis_sha256
from native_panel_aggregate_v1 import STATE_PANELS
from summarize_full_native_annotation_evidence_v1 import read_source


def gamma_summary(graph,units,values):
    graph=graph.copy();row=np.repeat(np.arange(len(units)),np.diff(graph.indptr));keep=units[row]!=units[graph.indices];graph.data[~keep]=0;graph.eliminate_zeros()
    denominator=np.asarray(graph.sum(axis=1)).ravel();numerator=np.asarray(graph@values.astype(float)).ravel();q=values&(denominator>0)
    return {'n_RNA_candidates':int(values.sum()),'n_RNA_candidates_with_other_unit_neighbors':int(q.sum()),'candidate_prevalence':float(values.mean()),'candidate_neighbor_fraction':float(np.mean(numerator[q]/denominator[q])) if q.any() else None,'orthogonal_gamma_identity_validation':False,'same_unit_and_self_excluded':True}


def main(root):
    start=time.monotonic();cfg=json.loads((root/'config/common_structure_spec_v1.json').read_text())
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Frozen structure diagnostic binding changed')
    producer,auditor=(Path(cfg[k]) for k in ['producer','model_audit'])
    for stage in [producer,auditor]:
        s=json.loads((stage/'provenance/supervisor_terminal_v1.json').read_text())
        if s['actual_child_exit_code'] or s['wall_timeout']:raise ValueError('Successful model and computational audit required')
    audit=json.loads((auditor/'provenance/terminal_v1.json').read_text())
    if audit['producer_terminal_sha256']!=sha256_file(producer/'provenance/terminal_v1.json'):raise ValueError('Model audit belongs to another producer')
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_common_reference_structure_metrics_v1.py')],capture_output=True,text=True,timeout=180)
    (root/'provenance/tests_stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 3 tests' not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Paired structure numeric counterexamples failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':3,'no_skips':True})
    spec=json.loads((producer/'config/common_reference_model_spec_v1.json').read_text());inputs=Path(spec['inputs']);inputspec=json.loads((inputs/'config/common_reference_input_spec_v1.json').read_text());candidate=Path(inputspec['candidate'])
    index=index_at(candidate,'full_native_cell_identity_candidate_index_v1.json');metadata=json.loads((producer/'outputs/reference_union_axis_metadata_v1.json').read_text());sources=metadata['sources']
    with np.load(producer/'outputs/reference_union_axes_and_native_readouts_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
    n=len(axes['native_row']);labels=np.full(n,'unresolved',dtype='U80');parent=np.full(n,-1,np.int32);gamma=np.zeros(n,bool)
    for sc,ds in enumerate(sources):
        m,a=read_source(candidate,index[ds]);u=np.flatnonzero(axes['source_code']==sc);positions=axes['native_source_position'][u]
        np.testing.assert_array_equal(a['native_row'][positions],axes['native_row'][u]);labels[u]=np.asarray(m['dictionaries']['consensus_celltype'])[a['consensus_celltype_code'][positions]];parent[u]=a['parent_source_local_group_code'][positions]
        gamma[u]=(a['own_fine_bitmask'][positions]&(1<<m['fine_candidate_order'].index('T_gamma_delta_RNA_candidate')))!=0
        del a;gc.collect()
    records=json.loads((producer/'outputs/common_reference_candidate_model_index_v1.json').read_text());routes=sorted({r['route'] for r in records});states=metadata['state_order'];coverage=[];profiles=[];paired=[];rare=[];sampling=[];local_group_profiles=[];independent_replays=0;max_error=0.
    for route in routes:
        storage_gate(root,cfg);models=[r for r in records if r['route']==route and r['status']=='candidate_model_fitted_not_selected_for_release'];coordinates={};rows=None
        if not models:
            coverage.append({'route':route,'status':'no_fitted_model_small_reference'});continue
        for r in models:
            p=producer/'outputs'/r['coordinates']
            if sha256_file(p)!=r['coordinates_sha256']:raise ValueError('Diagnostic model coordinate hash differs')
            with np.load(p,allow_pickle=False) as f:
                current=f['reference_union_row'];raw=f['raw_PCA_fit_embedding'];corrected=f['reference_embedding']
            if rows is None:rows=current
            else:np.testing.assert_array_equal(rows,current)
            coordinates[r['model_id']+'__raw_PCA']=raw
            if r['reference_kind']=='source_Harmony_Symphony_candidate':coordinates[r['model_id']+'__source_Harmony']=corrected
        source=axes['source_code'][rows];fine=labels[rows]
        for sc in np.unique(source):
            ds=sources[int(sc)]
            for label in np.unique(fine[source==sc]):
                eligible=(source==sc)&(fine==label);size=int(eligible.sum());record={'route':route,'source':ds,'conditioning_guarded_type_not_truth':str(label),'n_reference_stratum_rows':size}
                if size<cfg['minimum_stratum_cells']:
                    coverage.append({**record,'status':'insufficient_reference_stratum_rows'});continue
                sample=balanced_cell_selection(axes['expression_unit_code'][rows],parent[rows].astype(str),eligible,cfg['maximum_stratum_cells'],stable_seed(cfg['seed'],ds,route+'::'+str(label)),rare_reserve=2)
                q=rows[sample];units=axes['expression_unit_code'][q];depth=np.log1p(axes['native_total_RNA'][q]);values=axes['state_readouts'][q].astype(float)
                # Exact common panel is primary; reduced but measured panel remains a separately named sensitivity.
                full=[]
                for j,state in enumerate(states):
                    measured=metadata['source_state_measurement_coverage'][ds][state]['measured_genes'];available=set(measured)==set(STATE_PANELS[state]);full.append(available)
                primary=values.copy();primary[:,~np.asarray(full)]=np.nan
                z={k:v[sample] for k,v in coordinates.items()};audit_now=independent_replays<cfg['n_independent_stratum_replays'];graphs={k:exact_neighbor_graph(x) for k,x in z.items()}
                result,contrast,error=paired_native_states(z,units,depth,primary,states,cfg['minimum_other_unit_neighbors'],cfg['minimum_train_cells'],audit_now,graphs=graphs)
                if audit_now and any(r['n_evaluable']>0 for r in result):independent_replays+=1;max_error=max(max_error,error)
                for r in result:profiles.append({**record,'n_sampled_reference_rows':len(sample),'panel_scope':'complete_frozen_panel_primary',**r})
                for r in contrast:paired.append({**record,'n_sampled_reference_rows':len(sample),'panel_scope':'complete_frozen_panel_primary',**r})
                if not all(full):
                    result2,contrast2,_=paired_native_states(z,units,depth,values,states,cfg['minimum_other_unit_neighbors'],cfg['minimum_train_cells'],False,graphs=graphs)
                    for r in result2:profiles.append({**record,'n_sampled_reference_rows':len(sample),'panel_scope':'measured_at_least_three_genes_sensitivity',**r})
                    for r in contrast2:paired.append({**record,'n_sampled_reference_rows':len(sample),'panel_scope':'measured_at_least_three_genes_sensitivity',**r})
                for variant,x in z.items():
                    graph=graphs[variant].copy();row=np.repeat(np.arange(len(q)),np.diff(graph.indptr));keep=units[row]!=units[graph.indices];graph.data[~keep]=0;graph.eliminate_zeros();rr=np.repeat(np.arange(len(q)),np.diff(graph.indptr));den=np.asarray(graph.sum(axis=1)).ravel();same=parent[q][rr]==parent[q][graph.indices];weighted=graph.copy();weighted.data*=same;num=np.asarray(weighted.sum(axis=1)).ravel();valid=den>0
                    local_group_profiles.append({**record,'variant':variant,'n_sampled_reference_rows':len(q),'n_original_local_parent_groups':len(np.unique(parent[q])),'n_other_unit_neighbor_supported':int(valid.sum()),'same_original_local_group_neighbor_fraction':float(np.mean(num[valid]/den[valid])) if valid.any() else None,'local_group_preservation_not_identity_accuracy':True})
                if route in ['T_NK','global_overview']:
                    for variant,x in z.items():rare.append({**record,'variant':variant,**gamma_summary(graphs[variant],units,gamma[q])})
                sampling.append({**record,'union_reference_rows':q.tolist(),'row_axis_sha256':row_axis_sha256(q),'n_expression_units_not_patients':len(np.unique(units))})
                coverage.append({**record,'status':'diagnostic_computed_not_biological_validation','n_sampled_reference_rows':len(sample),'n_expression_units_not_patients':len(np.unique(units)),'complete_state_panels':dict(zip(states,map(bool,full)))})
            save(root/'checkpoints'/(route+'__'+ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'route':route,'source':ds,'stage':'conditional_native_structure_diagnostics_complete'})
        del coordinates;gc.collect()
    output=root/'outputs/common_reference_native_structure_diagnostics_v1.json';save(output,{'coverage':coverage,'state_profiles':profiles,'paired_state_profiles':paired,'gamma_RNA_profiles':rare,'local_group_profiles':local_group_profiles,'model_selected':False,'source_conditioned_transductive_reference_diagnostics':True,'guarded_conditioning_types_not_independent_truth':True,'clinical_or_CM_model_selection':False,'independent_replays':independent_replays,'independent_prediction_max_absolute_error':max_error})
    save(root/'outputs/private_fixed_diagnostic_sampling_v1.json',sampling)
    save(root/'provenance/terminal_v1.json',{'status':'PAIRED_COMMON_REFERENCE_NATIVE_STRUCTURE_DIAGNOSTICS_COMPUTED_SELECTION_AND_FULL_MAPPING_PENDING','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'output_sha256':sha256_file(output),'n_profiles':len(profiles),'n_paired_profiles':len(paired),'n_coverage_records':len(coverage),'n_independent_stratum_replays':independent_replays,'independent_prediction_max_absolute_error':max_error,'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'new_final_identity':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
