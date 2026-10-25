"""Bounded common-coordinate candidate models from accepted native references.

Native local identities and full expression remain authoritative. These models
are transductive representation candidates, not final taxonomy or clinical
validation. Context-sensitive epithelial/neural routes retain raw PCA only.
"""
import argparse
import datetime
import gc
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
import psutil
import scipy.sparse as sp
from annotation_representation_pilot_v1 import normalized_log_counts,select_features,pca_fit,save
from audit_full_source_local_graphs_v1 import index_at
from extract_full_source_local_marker_evidence_v1 import storage_gate
from frozen_symphony_reference_v1 import load_symphonypy,build_reference,save_minimal_reference,reference_fingerprint
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,row_axis_sha256,text_axis_sha256
from native_lineage_routing_v1 import FAMILIES
from native_panel_aggregate_v1 import gather_csr_rows
from summarize_full_native_annotation_evidence_v1 import read_source


def save_arrays(path,**arrays):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('xb') as f:np.savez_compressed(f,**arrays)
    with np.load(path,allow_pickle=False) as f:
        if set(f.files)!=set(arrays):raise ValueError('Reference model serialization keys differ')
        for k,v in arrays.items():np.testing.assert_array_equal(v,f[k])


def main(root):
    began=time.monotonic();cfg=json.loads((root/'config/common_reference_model_spec_v1.json').read_text())
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Reference-model input binding changed')
    inputs,audit,reference=(Path(cfg[k]) for k in ['inputs','input_audit','native_evidence'])
    for stage in [inputs,audit]:
        sup=json.loads((stage/'provenance/supervisor_terminal_v1.json').read_text())
        if sup['actual_child_exit_code'] or sup['wall_timeout']:raise ValueError('Accepted full reference input/audit required')
    at=json.loads((audit/'provenance/terminal_v1.json').read_text())
    if at['producer_terminal_sha256']!=sha256_file(inputs/'provenance/terminal_v1.json'):raise ValueError('Reference input audit binds another producer')
    params=cfg['parameters'];env=os.environ.copy();env.update(ATLAS_SYMPHONY_WHEEL=cfg['wheel_path'],ATLAS_SYMPHONY_SHA256=cfg['wheel_sha256'])
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_frozen_symphony_reference_v1.py')],env=env,capture_output=True,text=True,timeout=240)
    (root/'provenance/reference_tests.stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 4 tests' not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Real cached reference API tests failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':4,'no_skips':True})
    symphony=load_symphonypy(cfg['wheel_path'],cfg['wheel_sha256']);index=index_at(inputs,'full_common_reference_input_index_v1.json');ri=index_at(reference,'full_archive_reference_index_v1.json')
    reader=NativeArchiveCollection(cfg['archive_manifest']);sources=sorted(index)
    if set(sources)!=set(reader.archive_partitions) or set(ri)!=set(sources):raise ValueError('Reference-model source scope differs')
    featurefile=inputs/'outputs/common_measured_model_feature_candidates_v1.json';genes=json.loads(featurefile.read_text())['genes'];p=len(genes)
    reference_union=0
    for ds in sources:
        meta,side=read_source(inputs,index[ds]);reference_union+=int((side['global_reference_fit']|side['lineage_reference_fit']).sum())
    estimated=reference_union*p*16+cfg['RAM_workspace_reserve_bytes']
    if psutil.virtual_memory().available<estimated:raise MemoryError('Common reference sparse/dense-upper-bound RAM gate failed')
    save(root/'provenance/reference_RAM_preflight_v1.json',{'observed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_union_reference_rows':reference_union,
      'n_common_candidates':p,'conservative_peak_bound_bytes':estimated,'available_RAM_bytes':int(psutil.virtual_memory().available),'all_archive_dense_matrix_requested':False})
    chunks=[];axis={name:[] for name in ['source_code','native_row','native_source_position','global_reference_fit','lineage_reference_fit','route_code','unit_name','native_total_RNA','state_readouts']}
    readout_order=None;state_coverage={};nnz=0
    for source_code,ds in enumerate(sources):
        storage_gate(root,cfg);meta,side=read_source(inputs,index[ds]);rm,rv=read_source(reference,ri[ds]);rows=reader.rows(ds);n=len(rows)
        np.testing.assert_array_equal(rows,side['native_row']);np.testing.assert_array_equal(rows,rv['native_row'])
        if meta['native_ids_sha256']!=text_axis_sha256(reader.native_ids(ds)) or rm['native_ids_sha256']!=meta['native_ids_sha256']:raise ValueError('Reference-model native identity changed')
        f=inputs/'outputs'/meta['native_feature_columns']
        if sha256_file(f)!=meta['native_feature_columns_sha256'] or sha256_file(featurefile)!=meta['common_feature_list_sha256']:raise ValueError('Reference-model feature axes changed')
        with np.load(f,allow_pickle=False) as a:columns=a['native_columns']
        lookup=reader.feature_lookup(ds)
        if len(columns)!=p or len(set(columns))!=p or any(lookup.get(g,set())!={int(c)} for g,c in zip(genes,columns)):raise ValueError('Reference-model measured native columns differ')
        chosen=np.flatnonzero(side['global_reference_fit']|side['lineage_reference_fit'])
        if np.any(~side['fit_pool_candidate'][chosen]):raise ValueError('Archive exclusion entered model reference fit')
        native=reader.store(ds);key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        for low in range(0,len(chosen),cfg['read_block_size']):
            positions=chosen[low:low+cfg['read_block_size']];raw=gather_csr_rows(native[key],rows[positions]);totals=np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel()
            if not np.isfinite(raw.data).all() or np.any(raw.data<0) or not np.equal(raw.data,np.floor(raw.data)).all():raise ValueError('Reference-model native integer count domain differs')
            np.testing.assert_array_equal(totals,rv['native_total_RNA'][positions]);log,positive=normalized_log_counts(raw[:,columns],totals)
            if not positive.all():raise ValueError('Zero native RNA entered selected model reference')
            chunks.append(log);nnz+=raw.nnz
        if readout_order is None:readout_order=rm['state_order']
        if readout_order!=rm['state_order']:raise ValueError('Native readout dictionary differs')
        state_coverage[ds]=rm['state_measurement_coverage'];units=reader.frozen_expression_units(ds).astype(str)
        axis['source_code'].append(np.full(len(chosen),source_code,np.int16));axis['native_row'].append(rows[chosen]);axis['native_source_position'].append(chosen)
        axis['global_reference_fit'].append(side['global_reference_fit'][chosen]);axis['lineage_reference_fit'].append(side['lineage_reference_fit'][chosen]);axis['route_code'].append(side['guarded_route_code'][chosen])
        axis['unit_name'].append(np.char.add(ds+'::',units[chosen]));axis['native_total_RNA'].append(rv['native_total_RNA'][chosen]);axis['state_readouts'].append(rv['state_native_logCP10K'][chosen])
        save(root/'checkpoints'/('input_'+ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'stage':'union_reference_native_rows_loaded','n_reference_rows':len(chosen)})
        reader._fields.clear();del side,rv,units;gc.collect()
        if len(chosen):del raw,log
    x=sp.vstack(chunks,format='csr');del chunks;axis={k:np.concatenate(v) for k,v in axis.items()}
    if x.shape!=(reference_union,p) or any(len(v)!=reference_union for v in axis.values()):raise ValueError('Reference expression/row concatenation differs')
    unit_names,unit_code=np.unique(axis.pop('unit_name'),return_inverse=True);axis['expression_unit_code']=unit_code.astype(np.int32)
    axes=root/'outputs/reference_union_axes_and_native_readouts_v1.npz';save_arrays(axes,**axis)
    save(root/'outputs/reference_union_axis_metadata_v1.json',{'sources':sources,'expression_unit_dictionary':unit_names.tolist(),
      'unit_dictionary_sha256':text_axis_sha256(unit_names),'state_order':readout_order,'source_state_measurement_coverage':state_coverage,
      'expression_units_not_assumed_independent_patients':True,'source_reference_input_index_sha256':sha256_file(inputs/'outputs/full_common_reference_input_index_v1.json'),
      'reference_native_nnz_read':int(nnz),'native_expression_counts_saved':False,'all_source_transductive_discovery':True})
    global_rows=np.flatnonzero(axis['global_reference_fit']);shared=select_features(x[global_rows],axis['source_code'][global_rows],params['n_features'])
    jobs=[('global_overview','shared',global_rows,shared,False)]
    for family in FAMILIES:
        rows=np.flatnonzero(axis['lineage_reference_fit']&(axis['route_code']==FAMILIES.index(family)))
        if len(rows)<cfg['minimum_reference_rows']:
            jobs.append((family,'not_fitted_small_reference',rows,None,True));continue
        lineage=select_features(x[rows],axis['source_code'][rows],params['n_features'])
        if family in cfg['shared_feature_comparison_routes']:jobs.append((family,'shared',rows,shared,family in cfg['raw_only_context_routes']))
        jobs.append((family,'lineage',rows,lineage,family in cfg['raw_only_context_routes']))
    model_index=[]
    for family,variant,rows,selected,raw_only in jobs:
        storage_gate(root,cfg);identifier=family+'__'+variant;started=time.monotonic()
        if selected is None:
            model_index.append({'model_id':identifier,'route':family,'variant':variant,'status':'reference_too_small_retained_for_global_overview','n_reference_rows':len(rows)});continue
        selected=np.asarray(selected,np.int64);selected_genes=[genes[j] for j in selected];source=axis['source_code'][rows]
        raw_only=raw_only or len(np.unique(source))<2;kind='raw_PCA_context_preserving' if raw_only else 'source_Harmony_Symphony_candidate'
        model_path=root/'models'/(identifier+'_reference.npz');model_path.parent.mkdir(parents=True,exist_ok=True)
        if raw_only:
            raw,pca=pca_fit(x[rows][:,selected],cfg['seed'],params['components'],params['pca_iterations']);corrected=raw.copy();fingerprint=None;converged=None
            save_arrays(model_path,genes=np.asarray(selected_genes),mean=pca.mean_,std=np.ones(len(selected),np.float32),loadings=pca.components_.T,
              feature_candidate_indices=selected,explained_variance=pca.explained_variance_)
        else:
            reference_model,raw,corrected,pca=build_reference(x[rows][:,selected],source,selected_genes,cfg['seed'],params,symphony)
            fingerprint=reference_fingerprint(reference_model);h=reference_model.uns['harmony'];converged=None if h['converged'] is None else bool(h['converged'])
            save_minimal_reference(model_path,reference_model,selected,{'model_id':identifier,'scope':'transductive_common_coordinate_candidate_not_identity_truth',
              'reference_row_axis_sha256':row_axis_sha256(rows),'source_domain_transfer_not_yet_qualified':True,'reference_fingerprint':fingerprint})
            with np.load(model_path,allow_pickle=False) as saved:
                np.testing.assert_array_equal(saved['genes'],selected_genes);np.testing.assert_array_equal(saved['mean'],reference_model.var['mean'].to_numpy());np.testing.assert_array_equal(saved['loadings'],reference_model.varm['PCs'])
            if reference_fingerprint(reference_model)!=fingerprint:raise ValueError('Saving common model mutated frozen reference')
            del reference_model
        coordinates=root/'outputs'/(identifier+'_reference_coordinates_v1.npz');save_arrays(coordinates,reference_union_row=rows,raw_PCA_fit_embedding=raw,reference_embedding=corrected)
        record={'model_id':identifier,'route':family,'variant':variant,'status':'candidate_model_fitted_not_selected_for_release','reference_kind':kind,
          'n_reference_rows':len(rows),'n_reference_sources':len(np.unique(source)),'n_selected_features':len(selected),'selected_feature_names':selected_genes,
          'feature_selection_scope':'global_reference_rows' if variant=='shared' else 'lineage_reference_rows','model':str(model_path.relative_to(root)),
          'model_sha256':sha256_file(model_path),'coordinates':coordinates.name,'coordinates_sha256':sha256_file(coordinates),'reference_row_axis_sha256':row_axis_sha256(rows),
          'reference_fingerprint':fingerprint,'harmony_converged':converged,'raw_PCA_retained':True,'raw_fit_transform_not_assumed_equal_to_later_transform':True,
          'clinical_or_CM_hyperparameter_selection':False,'source_cancer_modality_or_patient_effects_not_claimed_separated':True,
          'malignant_clone_context_not_forced_into_normal_cell_identity':True,'local_fine_identity_replacement':False,
          'query_batch_adaptation_scope_must_be_explicit_before_full_mapping':True,'reference_diagnostics_and_full_mapping_pending':True,
          'wall_seconds':time.monotonic()-started}
        mp=root/'models'/(identifier+'_metadata.json');save(mp,record);model_index.append({**record,'metadata':str(mp.relative_to(root)),'metadata_sha256':sha256_file(mp)})
        save(root/'checkpoints'/('model_'+identifier+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'model_id':identifier,'stage':'common_reference_candidate_fitted','n_models_completed':sum(r['status']=='candidate_model_fitted_not_selected_for_release' for r in model_index)})
        del raw,corrected,pca;gc.collect()
    output=root/'outputs/common_reference_candidate_model_index_v1.json';save(output,model_index)
    save(root/'provenance/terminal_v1.json',{'status':'COMMON_REFERENCE_CANDIDATE_MODELS_FITTED_AUDIT_TRANSFER_AND_SELECTION_PENDING','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'n_sources':len(sources),'n_union_reference_rows':reference_union,'n_model_entries':len(model_index),'n_fitted_models':sum(r['status']=='candidate_model_fitted_not_selected_for_release' for r in model_index),
      'reference_native_nnz_read':int(nnz),'index_sha256':sha256_file(output),'reference_axes_sha256':sha256_file(axes),'input_auditor_terminal_sha256':sha256_file(audit/'provenance/terminal_v1.json'),
      'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'new_final_identity':False,'full_archive_mapping_complete':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
