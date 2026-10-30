"""Native all-row evaluation replay and all-seed saved-representation checks.

Uses independent sparse-window reconstruction and explicit Moran/correlation
formulae. Frozen biological RNA policy is shared and not independent truth.
No counts are copied to disk and no production identity is changed.
"""
import argparse,collections,datetime,hashlib,json,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
from full_context_atlas_v1 import FullContextAtlas,sha256_file
from native_archive_collection_v1 import NativeArchiveCollection
from native_annotation_review_v1 import native_measured_genes
from TME_all_lineage_cluster_policy_v1 import GATES,GATE_ORDER,coherent_flags,annotate_cluster
from TME_assay_aware_T_anchor_policy_v1 import ACCESSORY,SIGNAL,CHAIN,assay_aware_flags
from standard_annotation_order_primitives_v1 import COARSE_ROUTING_GENES,technical
from standard_annotation_order_readback_v1 import gather_windows,independent_conditional_metrics,compare_record

def load(p):return json.loads(Path(p).read_text())
def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def wait_file(path,producer):
    path=Path(path)
    while not path.exists():
        t=producer/'provenance/supervisor_terminal_v1.json'
        if t.exists():
            assert load(t)['actual_child_exit_code']==0,'Producer failed; no accepted substitute'
            raise ValueError('Successful producer omitted required artifact')
        time.sleep(5)

def main(R):
    started=time.monotonic();cfg=load(R/'config/spec_v1.json');P=Path(cfg['producer']);pc=load(P/'config/model_spec_v1.json')
    req=load(P/'config/supervisor_request_v1.json')
    for name,digest in req['script_SHA256'].items():assert sha256_file(P/name)==digest
    for item in pc['bindings'].values():assert sha256_file(item['path'])==item['sha256']
    assert pc['clinical_outcomes_read'] is False and pc['native_count_disk_copies'] is False
    pilot=Path(pc['bindings']['pilot_spec']['path']).parents[1];inventory=load(pc['bindings']['pilot_inventory']['path']);plan=load(pc['bindings']['features_and_readouts']['path'])
    def check(stage,**details):save(R/'checkpoints'/(str(time.time_ns())+'.json'),dict(observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage=stage,**details))
    sys.path.insert(0,str(Path(pc['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg)
    test=subprocess.run([sys.executable,'-B',str(R/'tests/test_standard_annotation_order_readback_v1.py')],capture_output=True,text=True,timeout=180)
    (R/'provenance/test_stdout.txt').write_text(test.stdout);(R/'provenance/test_stderr.txt').write_text(test.stderr)
    assert test.returncode==0 and 'Ran 3 tests' in test.stderr and 'skipped' not in test.stderr
    save(R/'provenance/test_receipt_v1.json',dict(actual_exit_code=0,n_tests=3,no_skips=True))
    for name in ['derived_axis_and_model_mask.npz','discovery_feature_selection_v1.json','evaluation_RNA_contract_v1.json','private_shared_evaluation_vectors_v1.npz']:
        wait_file(P/'outputs'/name,P)
    native=NativeArchiveCollection(pc['native_archive_manifest']);reader=FullContextAtlas(pc['bindings']['native_manifest']['path'])
    with np.load(P/'outputs/derived_axis_and_model_mask.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
    with np.load(P/'outputs/private_shared_evaluation_vectors_v1.npz',allow_pickle=False) as f:vectors={k:f[k] for k in f.files}
    contract=load(P/'outputs/evaluation_RNA_contract_v1.json');selection=load(P/'outputs/discovery_feature_selection_v1.json');selected=selection['selected_feature_names']
    assert len(selected)==len(set(selected))==pc['n_features'] and not any(technical(g) for g in selected)
    tab=pd.read_csv(P/'outputs/standard_batch_HVG_v1.tsv',sep='\t',index_col=0);assert tab.selected_HVG.sum()==len(selected)
    expected=tab.loc[~tab.technical_gene&tab.highly_variable_rank.notna()].sort_values(['highly_variable_rank','highly_variable_nbatches'],ascending=[True,False],kind='mergesort').index[:pc['n_features']].tolist()
    assert expected==selected and set(selected)<=set(selection['candidate_gene_universe'])
    held={g for panels in plan['heldout_RNA_readouts'].values() for genes in panels.values() for g in genes}
    assert not set(selected)&held
    all_names={}
    for entry in inventory:
        ds=entry['dataset'];names,cols=native_measured_genes(native,ds);all_names[ds]=dict(zip(names,map(int,cols)));native._fields.clear()
    states={name:[] for name in contract['state']};totals=[];flags=[];gamma=[];sources=[];probe_rows=[];probe_x=[];cursor=0
    for entry in inventory:
        ds=entry['dataset'];selector=pilot/entry['selector_file'];assert sha256_file(selector)==entry['selector_sha256'];positions=np.load(selector,allow_pickle=False)
        physical=reader.rows(ds)[positions];assert len(physical)==len(np.unique(physical)) and np.isin(physical,native.rows(ds)).all()
        names=all_names[ds];lookup=reader.feature_lookup(ds);modelcols=[names[g] for g in selected]
        assert len(set(modelcols))==len(selected)
        reserved={next(iter(lookup[g])) for g in held if len(lookup.get(g,set()))==1};assert not set(modelcols)&reserved
        wanted=set(ACCESSORY+SIGNAL+CHAIN+['CD3D','CD3E'])|{g for modules in GATES.values() for clause in modules for g in clause}
        typing=[g for g in names if g in wanted];typingcols=[names[g] for g in typing]
        for name,r in contract['state'].items():
            expected=[g for g in plan['heldout_RNA_readouts']['state'][name] if g not in COARSE_ROUTING_GENES and all(g in by for by in all_names.values())]
            assert expected==r['used_gene_names'] and not set(expected)&set(selected)
        z=reader.store(ds);key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        probes=set(np.unique(np.linspace(0,len(positions)-1,min(12,len(positions))).astype(int)).tolist())
        for low in range(0,len(positions),1024):
            high=min(len(positions),low+1024);raw=gather_windows(z[key],physical[low:high]);assert np.equal(raw.data,np.floor(raw.data)).all()
            total=np.asarray(raw.sum(axis=1)).ravel();totals.append(total)
            analysis=raw
            if 'ambient_removed_counts' in z['layers']:analysis=(raw-gather_windows(z['layers/ambient_removed_counts'],physical[low:high])).tocsr()
            analysis.eliminate_zeros();assert analysis.data.min(initial=0)>=-1e-8;analysis.data[analysis.data<0]=0;analysis.eliminate_zeros()
            RNA=np.asarray(analysis.sum(axis=1)).ravel();panel=analysis[:,typingcols];observed,assessed=coherent_flags(panel,typing);_,_,mode=assay_aware_flags(panel,typing)
            assert assessed==contract['source_assessed'][ds] and mode==contract['source_T_modes'][ds];flags.append(observed)
            g=np.full(high-low,np.nan)
            if assessed['gamma_delta_T_RNA_candidate']:g[:]=observed[:,GATE_ORDER.index('gamma_delta_T_RNA_candidate')]
            gamma.append(g)
            for name,r in contract['state'].items():
                used=r['used_gene_names'];value=np.full(high-low,np.nan)
                if len(used)>=pc['minimum_state_genes_measured']:
                    scale=np.zeros(len(RNA));ok=RNA>0;scale[ok]=10000/RNA[ok]
                    value=np.log1p(analysis[:,[names[g] for g in used]].toarray()*scale[:,None]).mean(axis=1);value[~ok]=np.nan
                states[name].append(value)
            local=np.asarray(sorted(probes&set(range(low,high))),np.int64)-low
            if len(local):
                values=raw[local][:,modelcols].toarray();scale=np.zeros(len(local));ok=total[local]>0;scale[ok]=10000/total[local][ok]
                probe_x.append(np.log1p(values*scale[:,None]));probe_rows.extend((cursor+low+local).tolist())
        sources.extend([ds]*len(positions));cursor+=len(positions);check('full_native_evaluation_source_replayed',source=ds,n_rows=len(positions));native._fields.clear();reader._fields.clear()
    source=np.asarray(sources);total=np.concatenate(totals);flag=np.vstack(flags);rare=np.concatenate(gamma);state={k:np.concatenate(v) for k,v in states.items()};probe_x=np.vstack(probe_x);probe_rows=np.asarray(probe_rows)
    np.testing.assert_array_equal(source,axes['source']);np.testing.assert_array_equal(total>0,axes['modeled']);np.testing.assert_array_equal(total,vectors['native_total']);np.testing.assert_array_equal(flag,vectors['native_coherent_flags']);np.testing.assert_allclose(rare,vectors['gamma_RNA_candidate'],rtol=0,atol=0,equal_nan=True)
    for name,v in state.items():np.testing.assert_allclose(v,vectors['state_'+name],rtol=2e-5,atol=2e-5,equal_nan=True)
    save(R/'outputs/native_evaluation_full_replay_acceptance_v1.json',dict(n_records=len(source),n_sources=len(inventory),all_raw_total_and_native_RNA_readouts_replayed=True,measurement_masks_and_unmeasured_gamma_replayed=True,frozen_RNA_policy_shared_not_identity_truth=True,model_feature_domain_and_reserved_alias_checks_passed=True,official_HVG_rank_not_recomputed_in_this_audit=True,n_native_PCA_probes=len(probe_rows)))
    rng=np.random.default_rng(pc['diagnostic_sampling_seed']);anchors={}
    for ds in np.unique(source):
        for route in np.unique(axes['inherited_route']):
            wh=np.flatnonzero((source==ds)&(axes['inherited_route']==route)&axes['modeled'])
            if len(wh)>=pc['minimum_diagnostic_stratum']:anchors[(ds,route)]=np.sort(rng.choice(wh,min(len(wh),pc['diagnostic_cells_per_stratum']),replace=False))
    seed_receipts=[]
    for seed in pc['model_seeds']:
        for name in ['candidate_partitions_'+str(seed)+'.npz','representation_diagnostics_'+str(seed)+'.json','source_local_coarse_RNA_routes_'+str(seed)+'.json','private_embeddings_and_C_spaces_'+str(seed)+'.npz']:
            wait_file(P/'outputs'/name,P)
        with np.load(P/'outputs'/('candidate_partitions_'+str(seed)+'.npz'),allow_pickle=False) as f:part={k:f[k] for k in f.files}
        with np.load(P/'outputs'/('private_embeddings_and_C_spaces_'+str(seed)+'.npz'),allow_pickle=False) as f:emb={k:f[k] for k in f.files}
        cspace=emb.pop('C_space_id');assert len(cspace)==len(source)
        for method,groups in [('A',source),('B',np.repeat('joint',len(source))),('C',cspace)]:
            assert len(part[method])==len(source) and (part[method]>=-1).all()
            for k in np.unique(part[method][part[method]>=0]):assert len(np.unique(groups[part[method]==k]))==1
        routes=np.full(len(source),'unresolved_context',dtype='<U40');expected_records=[]
        for ds in np.unique(source):
            for k in np.unique(part['A'][source==ds]):
                if k<0:continue
                wh=(source==ds)&(part['A']==k);count=dict(zip(GATE_ORDER,flag[wh].sum(axis=0).astype(int).tolist()));proposal=annotate_cluster(int(wh.sum()),count,contract['source_assessed'][ds]);major=proposal['major'];route='T_NK' if major in ['T','NK','T_NK'] else 'unresolved_context' if major=='unresolved' else major;routes[wh]=route
                expected_records.append(dict(source=str(ds),source_local_partition=int(k),n_cells=int(wh.sum()),coarse_RNA_route=route,status=proposal['status'],not_independent_identity_truth=True))
        np.testing.assert_array_equal(routes,cspace);assert expected_records==load(P/'outputs'/('source_local_coarse_RNA_routes_'+str(seed)+'.json'))
        with np.load(P/'models'/('global_discovery_PCA_'+str(seed)+'.npz'),allow_pickle=False) as f:components=f['components'];mean=f['mean'];variance=f['explained_variance']
        np.testing.assert_allclose(components@components.T,np.eye(len(components)),rtol=3e-4,atol=3e-4);assert np.isfinite(variance).all() and (variance>=0).all()
        np.testing.assert_allclose((probe_x-mean)@components.T,emb['raw_global_control'][probe_rows],rtol=3e-4,atol=3e-3)
        actual=load(P/'outputs'/('representation_diagnostics_'+str(seed)+'.json'))['conditional_RNA_diagnostics'];expected=[]
        for (ds,route),fixed in anchors.items():
            for cr in np.unique(cspace[fixed]):
                rows=fixed[cspace[fixed]==cr]
                for name,z in emb.items():
                    assert z.shape==(len(source),pc['components'])
                    if not np.isfinite(z[rows]).all():metrics=dict(status='representation_unavailable_small_group')
                    else:metrics=independent_conditional_metrics(z,state,total,rare,rows,pc['diagnostic_neighbors'])
                    expected.append(dict(source=str(ds),inherited_route=str(route),source_local_C_RNA_route=str(cr),n_frozen_base_anchors=len(fixed),every_method_uses_same_anchor_subset=True,representation=name,**metrics))
        assert len(actual)==len(expected)
        for a,e in zip(actual,expected):compare_record(a,e)
        seed_receipts.append(dict(seed=seed,n_conditional_records_replayed=len(expected),all_partition_and_route_axes_checked=True,native_PCA_probes=len(probe_rows)))
        check('seed_saved_axes_native_PCA_and_conditional_metrics_replayed',**seed_receipts[-1]);storage_gate(R,cfg)
    wait_file(P/'provenance/supervisor_terminal_v1.json',P);assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    total_bytes=sum(p.stat().st_size for p in P.rglob('*') if p.is_file());assert total_bytes<=pc['output_cap_bytes']
    files=[dict(path=str(p.relative_to(P)),bytes=p.stat().st_size,SHA256=sha256_file(p)) for p in sorted(P.rglob('*')) if p.is_file() and p.parent.name in ['outputs','models']]
    save(R/'outputs/standard_ABC_readback_acceptance_v1.json',dict(n_records=len(source),n_sources=len(inventory),seeds=seed_receipts,producer_output_bytes=total_bytes,producer_file_SHA256=files,producer_terminal_SHA256=sha256_file(P/'provenance/terminal_v1.json'),native_evaluation_all_rows_replayed=True,saved_coordinate_conditional_metrics_independently_recomputed=True,official_HVG_ranking_not_independently_recomputed=True,biological_identity_accuracy_and_study_holdout_not_proven=True,method_winner_selected=False))
    save(R/'provenance/terminal_v1.json',dict(status='STANDARD_ABC_NATIVE_READOUTS_ALL_SAVED_AXES_AND_CONDITIONAL_STATISTICS_ACCEPTED_REQUIRES_BIOLOGICAL_AND_STUDY_HOLDOUT_REVIEW',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_records=len(source),n_seeds=len(seed_receipts),acceptance_SHA256=sha256_file(R/'outputs/standard_ABC_readback_acceptance_v1.json'),wall_seconds=time.monotonic()-started,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
