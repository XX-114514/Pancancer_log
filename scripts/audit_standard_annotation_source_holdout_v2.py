"""All-fold training-HVG/reference axes and full query transfer-score readback.

Independent distances/weights replay every query's five prediction probes.
Reference mapping reuses pinned Symphony; C prototype construction membership
is not saved by the producer and is explicitly not independently replayed.
"""
import argparse,collections,datetime,json,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
from full_context_atlas_v1 import sha256_file
from annotation_pilot_sampling_v1 import balanced_cell_selection,stable_seed
from standard_annotation_order_primitives_v1 import select_features,coarse_routes_from_local
from standard_annotation_order_readback_v1 import compare_record
from standard_study_holdout_readback_v1 import restore_reference,prototype_prediction,independent_scores,check_score_bounds
from frozen_symphony_reference_standard_v2 import load_symphonypy,map_query,reference_fingerprint
from load_standard_annotation_pilot_native_v2 import load_pilot

def load(p):return json.loads(Path(p).read_text())
def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def wait_file(p,producer):
    p=Path(p)
    while not p.exists():
        t=producer/'provenance/supervisor_terminal_v1.json'
        if t.exists():
            assert load(t)['actual_child_exit_code']==0,'Study holdout producer actually failed'
            raise ValueError('Successful producer omitted artifact')
        time.sleep(5)
def main(R):
    start=time.monotonic();acfg=load(R/'config/spec_v1.json');P=Path(acfg['producer']);cfg=load(P/'config/source_holdout_spec_v1.json');req=load(P/'config/supervisor_request_v1.json')
    for name,digest in req['script_SHA256'].items():assert sha256_file(P/name)==digest
    for b in cfg['bindings'].values():assert sha256_file(b['path'])==b['sha256']
    launch=load(P/'provenance/launch_receipt_v1.json');assert sha256_file(P/'config/source_holdout_spec_v1.json')==launch['config_SHA256']
    assert cfg['no_clinical_outcome_inputs'] and cfg['HVG_PCA_Harmony_and_reference_partitions_training_only']
    sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,acfg)
    test=subprocess.run([sys.executable,'-B',str(R/'tests/test_standard_study_holdout_readback_v1.py')],capture_output=True,text=True,timeout=180)
    (R/'provenance/test_stdout.txt').write_text(test.stdout);(R/'provenance/test_stderr.txt').write_text(test.stderr)
    assert test.returncode==0 and 'Ran 2 tests' in test.stderr and 'skipped' not in test.stderr
    save(R/'provenance/test_receipt_v1.json',dict(actual_exit_code=0,n_tests=2,no_skips=True))
    def check(stage,**v):save(R/'checkpoints'/(str(time.time_ns())+'.json'),dict(observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage=stage,**v))
    data=load_pilot(Path(cfg['bindings']['pilot_spec']['path']).parents[1],cfg['bindings']['native_manifest']['path'],check,cfg['read_block_rows'],discovery=cfg['discovery'],native_auditor=cfg['native_auditor'])
    source=data['source'];total=data['native_total'];programs=sorted(data['state']);target=np.column_stack([data['state'][p] for p in programs]);folds=sorted(set(data['fold']));assert folds==cfg['expected_fold_groups'];symphony=load_symphonypy(cfg['wheel_path'],cfg['wheel_SHA256']);results=[];fold_receipts=[]
    for number,fold in enumerate(folds):
        train=np.flatnonzero((total>0)&(data['fold']!=fold));query=np.flatnonzero((total>0)&(data['fold']==fold));assert not set(source[train])&set(source[query])
        selected,_=select_features(data['counts'][train],source[train],cfg['n_features'],data['genes'],cfg['candidate_HVGs']);genes=[data['genes'][j] for j in selected]
        F=P/'models'/('fold_'+str(number));wait_file(F/'reference_axes_and_parameters_v1.json',P);wait_file(F/'private_fold_axes_and_partitions_v2.npz',P)
        meta=load(F/'reference_axes_and_parameters_v1.json');np.testing.assert_array_equal(selected,meta['selected_model_feature_indices']);assert meta['training_sources']==sorted(set(source[train])) and meta['query_sources']==sorted(set(source[query]))
        with np.load(F/'private_fold_axes_and_partitions_v2.npz',allow_pickle=False) as f:a={k:f[k] for k in f.files}
        np.testing.assert_array_equal(a['train_rows'],train);np.testing.assert_array_equal(a['query_rows'],query)
        troute,tr=coarse_routes_from_local(source[train],a['A_train'],data['native_flags'][train],data['source_assessed']);qroute,qr=coarse_routes_from_local(source[query],a['A_query'],data['native_flags'][query],data['source_assessed']);np.testing.assert_array_equal(troute,a['training_coarse_RNA_route']);np.testing.assert_array_equal(qroute,a['query_coarse_RNA_route'])
        routes=load(F/'source_local_native_coarse_routes_v2.json');assert routes['training']==tr and routes['query']==qr
        ref,rm=restore_reference(F/'A_B_shared_minimal_reference.npz',genes,selected);before=reference_fingerprint(ref);assert before==meta['reference_fingerprint_before_query']==meta['reference_fingerprint_after_query']
        np.testing.assert_allclose(ref.var['mean'],np.asarray(data['x'][train][:,selected].mean(axis=0)).ravel(),rtol=3e-4,atol=3e-4)
        _,mapped,supported=map_query(ref,data['x'][query][:,selected],genes,source[query],symphony);assert reference_fingerprint(ref)==before
        fold_records=[];replayed_queries=0;route_receipts=[]
        for wanted in cfg['reference_routes']:
            lt=np.flatnonzero(troute==wanted);lq=np.flatnonzero(qroute==wanted);nt,nq=len(lt),len(lq)
            if nt<cfg['minimum_training_group'] or nq==0:
                route_receipts.append(dict(route=wanted,status='reference_or_query_scope_unavailable',n_training_cells=nt,n_query_cells=nq));continue
            wait_file(F/(wanted+'_private_scoring_readback_v2.npz'),P)
            cref,cm=restore_reference(F/(wanted+'_minimal_reference.npz'),genes,selected);cbefore=reference_fingerprint(cref)
            np.testing.assert_allclose(cref.var['mean'],np.asarray(data['x'][train[lt]][:,selected].mean(axis=0)).ravel(),rtol=3e-4,atol=3e-4)
            _,cmap,csupport=map_query(cref,data['x'][query[lq]][:,selected],genes,source[query[lq]],symphony);assert reference_fingerprint(cref)==cbefore;common=supported[lq]&csupport
            with np.load(F/(wanted+'_private_scoring_readback_v2.npz'),allow_pickle=False) as f:side={k:f[k] for k in f.files}
            anchors=balanced_cell_selection(source[train[lt]],data['fine'][train[lt]],np.ones(nt,bool),cfg['anchor_cap_per_route'],stable_seed(cfg['seed'],fold,wanted),rare_reserve=2);anchor_rows=train[lt[anchors]];qrows=query[lq]
            np.testing.assert_array_equal(side['query_pilot_rows'],qrows);np.testing.assert_array_equal(side['anchor_pilot_rows'],anchor_rows);np.testing.assert_array_equal(side['query_source'],source[qrows]);np.testing.assert_array_equal(side['common_supported'],common);np.testing.assert_array_equal(side['query_native_total'],total[qrows]);np.testing.assert_allclose(side['query_RNA'],target[qrows],equal_nan=True,atol=0,rtol=0);np.testing.assert_allclose(side['anchor_RNA'],target[anchor_rows],equal_nan=True,atol=0,rtol=0)
            np.testing.assert_allclose(side['B_query'],mapped[lq],equal_nan=True,rtol=2e-5,atol=2e-5);np.testing.assert_allclose(side['C_query'],cmap,equal_nan=True,rtol=2e-5,atol=2e-5)
            predictions={}
            for name,zref,zq in [('A_B_shared_reference_neighbors',side['B_anchor'],side['B_query']),('C_reference_neighbors',side['C_anchor'],side['C_query'])]:
                pred=np.full(target[qrows].shape,np.nan);pred[common]=prototype_prediction(zref,target[anchor_rows],zq[common],cfg['state_neighbors']);predictions[name]=pred
            for name,prefix,zq in [('A_conditional_partition_prototypes','A_conditional_partition_prototypes',mapped[lq]),('B_conditional_partition_prototypes','B_conditional_partition_prototypes',mapped[lq]),('C_conditional_partition_prototypes','C_partition_prototypes',cmap)]:
                with np.load(F/(prefix+'_'+wanted+'.npz'),allow_pickle=False) as f:centers=f['centers'];outcomes=f['native_RNA_program_means']
                pred=np.full(target[qrows].shape,np.nan);pred[common]=prototype_prediction(centers,outcomes,zq[common],cfg['prototype_neighbors']);predictions[name]=pred
            for name,pred in predictions.items():
                np.testing.assert_allclose(pred,side[name],rtol=3e-5,atol=3e-5,equal_nan=True)
                for ds in np.unique(source[qrows]):
                    member=(source[qrows]==ds)&common
                    for j,program in enumerate(programs):
                        record=dict(heldout_fold=fold,query_source=str(ds),source_local_native_coarse_route=wanted,probe=name,program=program,n_query_route_cells=int((source[qrows]==ds).sum()),same_evaluable_query_subset_across_compared_probes=True,**independent_scores(target[qrows[member],j],pred[member,j],total[qrows[member]],cfg['minimum_query_evaluation']));check_score_bounds(record);fold_records.append(record)
            replayed_queries+=int(common.sum());route_receipts.append(dict(route=wanted,status='all_query_neighbor_and_prototype_predictions_recomputed',n_query_cells=nq,n_common_supported=int(common.sum()),n_training_anchors=len(anchor_rows),C_prototype_construction_membership_not_independently_replayed=True))
            check('full_route_query_predictions_independently_replayed',fold=fold,route=wanted,n_common_queries=int(common.sum()));storage_gate(R,acfg)
        wait_file(P/'outputs'/('fold_'+str(number)+'_qualification.json'),P);coverage=load(P/'outputs'/('fold_'+str(number)+'_qualification.json'));assert coverage['training_sources']==sorted(set(source[train])) and coverage['query_sources']==sorted(set(source[query]))
        byroute={r['route']:r for r in coverage['routes']}
        for r in route_receipts:
            observed=byroute[r['route']];assert observed['n_training_cells']==int((troute==r['route']).sum()) and observed['n_query_cells']==r['n_query_cells']
            if r['status']=='reference_or_query_scope_unavailable':assert observed['status']==r['status']
            else:assert observed['n_common_numerically_supported_query_cells']==r['n_common_supported']
        results.extend(fold_records);fold_receipts.append(dict(fold=fold,training_only_HVG_recomputed=True,n_queries_independently_predicted=replayed_queries,routes=route_receipts));check('study_fold_train_features_and_all_query_transfer_replayed',fold=fold)
    wait_file(P/'provenance/supervisor_terminal_v1.json',P);assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    actual=load(P/'outputs/source_heldout_RNA_structure_v1.json');key=lambda r:tuple(r[k] for k in ['heldout_fold','query_source','source_local_native_coarse_route','probe','program']);expected={key(r):r for r in results};assert len(expected)==len(actual)==len(results)
    for record in actual:compare_record(record,expected[key(record)],atol=3e-5)
    byte_count=sum(p.stat().st_size for p in P.rglob('*') if p.is_file());assert byte_count<=cfg['output_cap_bytes']
    save(R/'outputs/standard_study_holdout_readback_acceptance_v1.json',dict(folds=fold_receipts,n_result_records=len(results),all_training_only_feature_selections_recomputed=True,all_query_predictions_and_scores_recomputed=True,reference_mapping_backend_shared_with_producer=True,all_native_readouts_previously_independently_replayed=True,C_prototype_construction_membership_not_independently_replayed=True,primary_cell_neighbor_reference_geometry_and_secondary_prototypes_reported_separately=True,producer_output_bytes=byte_count,producer_terminal_SHA256=sha256_file(P/'provenance/terminal_v1.json'),scientific_summary_SHA256=sha256_file(P/'outputs/source_heldout_RNA_structure_v1.json'),method_winner_or_final_identity_selected=False))
    save(R/'provenance/terminal_v1.json',dict(status='STANDARD_STUDY_HOLDOUT_TRAIN_FEATURES_REFERENCE_AXES_ALL_QUERY_PREDICTIONS_AND_SCORES_ACCEPTED_SECONDARY_C_PROTOTYPE_CONSTRUCTION_NOT_REPLAYED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_folds=len(fold_receipts),n_scientific_records=len(results),acceptance_SHA256=sha256_file(R/'outputs/standard_study_holdout_readback_acceptance_v1.json'),wall_seconds=time.monotonic()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
