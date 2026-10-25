"""Full reference/model axis audit plus native query projection counterexamples.

Computational acceptance does not qualify cross-source biological transfer.
Query tests use a declared bounded batch, never pretend to map the full atlas.
"""
import argparse,datetime,gc,json,resource,time
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from annotation_representation_pilot_v1 import save,normalized_log_counts
from audit_full_source_local_graphs_v1 import index_at
from extract_full_source_local_marker_evidence_v1 import storage_gate
from frozen_symphony_reference_v1 import load_symphonypy,map_query,reference_fingerprint
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,row_axis_sha256,text_axis_sha256
from native_lineage_routing_v1 import FAMILIES
from native_panel_aggregate_v1 import gather_csr_rows
from summarize_full_native_annotation_evidence_v1 import read_source


def restore(v):
    import anndata as ad
    ref=ad.AnnData(sp.csr_matrix((0,len(v['genes'])),dtype=np.float32),var=pd.DataFrame({'mean':v['mean'],'std':v['std']},index=v['genes']))
    ref.varm['PCs']=v['loadings'];ref.uns['harmony']={'C':v['C'],'Nr':v['Nr'],'K':len(v['Nr']),'sigma':v['sigma'],'ref_basis_loadings':'PCs','ref_basis_adjusted':'X_pca_harmony','converged':None}
    return ref


def main(root):
    started=time.monotonic();cfg=json.loads((root/'config/common_reference_model_audit_spec_v1.json').read_text());producer=Path(cfg['producer'])
    sup=json.loads((producer/'provenance/supervisor_terminal_v1.json').read_text());terminal=json.loads((producer/'provenance/terminal_v1.json').read_text())
    if sup['actual_child_exit_code'] or sup['wall_timeout']:raise ValueError('Successful full model producer required')
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Model audit input binding changed')
    spec=json.loads((producer/'config/common_reference_model_spec_v1.json').read_text());inputs=Path(spec['inputs']);evidence=Path(spec['native_evidence'])
    indexfile=producer/'outputs/common_reference_candidate_model_index_v1.json';axesfile=producer/'outputs/reference_union_axes_and_native_readouts_v1.npz'
    if sha256_file(indexfile)!=terminal['index_sha256'] or sha256_file(axesfile)!=terminal['reference_axes_sha256']:raise ValueError('Full model output hashes differ')
    with np.load(axesfile,allow_pickle=False) as f:axes={k:f[k] for k in f.files}
    am=json.loads((producer/'outputs/reference_union_axis_metadata_v1.json').read_text());sources=am['sources'];units=np.asarray(am['expression_unit_dictionary'])
    if text_axis_sha256(units)!=am['unit_dictionary_sha256'] or any(len(v)!=terminal['n_union_reference_rows'] for v in axes.values()):raise ValueError('Union reference axes differ')
    ci=index_at(inputs,'full_common_reference_input_index_v1.json');ri=index_at(evidence,'full_archive_reference_index_v1.json');reader=NativeArchiveCollection(spec['archive_manifest'])
    if sources!=sorted(ci) or set(ri)!=set(sources) or set(reader.archive_partitions)!=set(sources):raise ValueError('Audit source scope differs')
    genes=np.asarray(json.loads((inputs/'outputs/common_measured_model_feature_candidates_v1.json').read_text())['genes']);cursor=0
    for sc,ds in enumerate(sources):
        storage_gate(root,cfg);meta,side=read_source(inputs,ci[ds]);em,ev=read_source(evidence,ri[ds]);chosen=np.flatnonzero(side['global_reference_fit']|side['lineage_reference_fit']);sl=slice(cursor,cursor+len(chosen))
        if not np.all(axes['source_code'][sl]==sc):raise ValueError('Source order differs')
        np.testing.assert_array_equal(axes['native_source_position'][sl],chosen);np.testing.assert_array_equal(axes['native_row'][sl],reader.rows(ds)[chosen])
        np.testing.assert_array_equal(axes['global_reference_fit'][sl],side['global_reference_fit'][chosen]);np.testing.assert_array_equal(axes['lineage_reference_fit'][sl],side['lineage_reference_fit'][chosen]);np.testing.assert_array_equal(axes['route_code'][sl],side['guarded_route_code'][chosen])
        np.testing.assert_array_equal(axes['native_total_RNA'][sl],ev['native_total_RNA'][chosen]);np.testing.assert_array_equal(axes['state_readouts'][sl],ev['state_native_logCP10K'][chosen])
        np.testing.assert_array_equal(units[axes['expression_unit_code'][sl]],np.char.add(ds+'::',reader.frozen_expression_units(ds).astype(str)[chosen]))
        if meta['native_ids_sha256']!=text_axis_sha256(reader.native_ids(ds)) or em['state_order']!=am['state_order'] or em['state_measurement_coverage']!=am['source_state_measurement_coverage'][ds]:raise ValueError('Native readout/identity metadata differs')
        cursor+=len(chosen);reader._fields.clear();del side,ev;gc.collect()
    if cursor!=len(axes['native_row']):raise ValueError('Unvisited union reference rows')
    records=json.loads(indexfile.read_text());identifiers=[r['model_id'] for r in records]
    expected=['global_overview__shared']
    for family in FAMILIES:
        rows=np.flatnonzero(axes['lineage_reference_fit']&(axes['route_code']==FAMILIES.index(family)))
        if len(rows)<spec['minimum_reference_rows']:expected.append(family+'__not_fitted_small_reference')
        else:
            if family in spec['shared_feature_comparison_routes']:expected.append(family+'__shared')
            expected.append(family+'__lineage')
    if identifiers!=expected:raise ValueError('Missing, duplicated or reordered model candidates')
    symphony=load_symphonypy(spec['wheel_path'],spec['wheel_sha256']);audit=[]
    for record in records:
        storage_gate(root,cfg);family=record['route'];rows=np.flatnonzero(axes['global_reference_fit']) if family=='global_overview' else np.flatnonzero(axes['lineage_reference_fit']&(axes['route_code']==FAMILIES.index(family)))
        if len(rows)!=record['n_reference_rows']:raise ValueError('Model reference denominator differs')
        if record['status']=='reference_too_small_retained_for_global_overview':
            if len(rows)>=spec['minimum_reference_rows']:raise ValueError('Improper small-pool fallback')
            audit.append({'model_id':record['model_id'],'status':'small_pool_fallback_verified'});continue
        model=producer/record['model'];coords=producer/'outputs'/record['coordinates'];metadata=producer/record['metadata']
        if any(sha256_file(p)!=record[k] for p,k in [(model,'model_sha256'),(coords,'coordinates_sha256'),(metadata,'metadata_sha256')]):raise ValueError('Saved model or coordinate hash differs')
        if row_axis_sha256(rows)!=record['reference_row_axis_sha256']:raise ValueError('Model row hash differs')
        with np.load(model,allow_pickle=False) as f:v={k:f[k].copy() for k in f.files}
        selected=v['feature_candidate_indices'];np.testing.assert_array_equal(v['genes'],genes[selected]);np.testing.assert_array_equal(v['genes'],record['selected_feature_names'])
        if len(set(selected))!=len(selected) or v['loadings'].shape[0]!=len(selected) or v['mean'].shape!=(len(selected),) or not np.all(v['std']==1):raise ValueError('Model feature/scaling axes invalid')
        for k in ['mean','std','loadings']:
            if not np.isfinite(v[k]).all():raise ValueError('Nonfinite model parameters')
        np.testing.assert_allclose(v['loadings'].T@v['loadings'],np.eye(v['loadings'].shape[1]),rtol=2e-4,atol=2e-4)
        with np.load(coords,allow_pickle=False) as f:
            np.testing.assert_array_equal(f['reference_union_row'],rows);raw=f['raw_PCA_fit_embedding'];corrected=f['reference_embedding']
            if raw.shape!=corrected.shape or raw.shape!=(len(rows),v['loadings'].shape[1]) or not np.isfinite(raw).all() or not np.isfinite(corrected).all():raise ValueError('Coordinate axes or finite support differ')
            if record['reference_kind']=='raw_PCA_context_preserving':np.testing.assert_array_equal(raw,corrected)
        # One predetermined native source/batch per model; order fixed before reading expression.
        probe=rows[:min(32,len(rows))];native=[]
        for sc in np.unique(axes['source_code'][probe]):
            ds=sources[int(sc)];positions=np.flatnonzero(axes['source_code'][probe]==sc);target=probe[positions];lookup=reader.feature_lookup(ds)
            columns=[next(iter(lookup[str(g)])) for g in v['genes']]
            if any(lookup[str(g)]!={int(c)} for g,c in zip(v['genes'],columns)):raise ValueError('Probe native feature ambiguity')
            store=reader.store(ds);key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key'];counts=gather_csr_rows(store[key],axes['native_row'][target]);total=np.asarray(counts.sum(axis=1,dtype=np.float64)).ravel()
            np.testing.assert_array_equal(total,axes['native_total_RNA'][target]);x,positive=normalized_log_counts(counts[:,columns],total)
            if not positive.all():raise ValueError('Nonpositive native test reference')
            native.append((positions,x))
        x=sp.vstack([next(x[j] for pos,x in native for j,k in enumerate(pos) if k==i) for i in range(len(probe))],format='csr')
        projected=np.asarray(x@v['loadings'])-v['mean']@v['loadings'];error=None
        if record['reference_kind']=='source_Harmony_Symphony_candidate':
            ref=restore(v);before=reference_fingerprint(ref)
            if before!=record['reference_fingerprint']:raise ValueError('Restored frozen reference differs')
            batch=axes['source_code'][probe].astype(str);direct,mapped,valid=map_query(ref,x,v['genes'].tolist(),batch,symphony);order=np.arange(len(probe))[::-1]
            reverse,rmapped,rvalid=map_query(ref,x[order],v['genes'].tolist(),batch[order],symphony)
            np.testing.assert_array_equal(valid,rvalid[order]);np.testing.assert_allclose(direct,projected,rtol=2e-5,atol=2e-4);np.testing.assert_allclose(mapped[valid],rmapped[order][valid],rtol=2e-5,atol=2e-4)
            if not valid.any() or reference_fingerprint(ref)!=before:raise ValueError('Probe numerical mapping unsupported or reference mutated')
            error=float(np.max(np.abs(mapped[valid]-rmapped[order][valid])));del ref
        audit.append({'model_id':record['model_id'],'status':'model_axes_and_native_projection_verified','n_native_probe_rows':len(probe),'query_order_max_absolute_difference':error,'probe_is_not_full_source_transfer_qualification':True})
        save(root/'checkpoints'/(record['model_id']+'.json'),audit[-1]);del v,raw,corrected,x;reader._fields.clear();gc.collect()
    output=root/'outputs/common_reference_model_audit_v1.json';save(output,{'all_union_reference_axes_verified':True,'model_checks':audit,'full_archive_mapping_or_biological_transfer_qualified':False})
    save(root/'provenance/terminal_v1.json',{'status':'COMMON_REFERENCE_COMPUTATIONAL_AUDIT_ACCEPTED_TRANSFER_AND_SELECTION_PENDING','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(sources),'n_union_reference_rows':cursor,'n_models_checked':len(audit),'output_sha256':sha256_file(output),'producer_terminal_sha256':sha256_file(producer/'provenance/terminal_v1.json'),'wall_seconds':time.monotonic()-started,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
