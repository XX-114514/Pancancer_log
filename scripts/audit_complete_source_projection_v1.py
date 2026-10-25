"""Full source masks and independent complete-query ridge replay.

Only bounded native rows and geometric probes are reread. All included query
rows contribute to one sufficient-statistic fit, with no per-chunk refitting.
"""
import argparse,datetime,gc,json,resource,subprocess,sys,time,os
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from annotation_representation_pilot_v1 import save,normalized_log_counts
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
from native_lineage_routing_v1 import FAMILIES
from native_panel_aggregate_v1 import gather_csr_rows
from single_source_symphony_replay_v1 import fit_single_source_offsets,apply_single_source_offsets
from summarize_full_native_annotation_evidence_v1 import read_source


def main(root):
    started=time.monotonic();cfg=json.loads((root/'config/complete_projection_audit_spec_v1.json').read_text());producer=Path(cfg['producer']);t=json.loads((producer/'provenance/terminal_v1.json').read_text());s=json.loads((producer/'provenance/supervisor_terminal_v1.json').read_text())
    if s['actual_child_exit_code'] or s['wall_timeout']:raise ValueError('Successful complete-source projection required')
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Projection audit input binding changed')
    pspec=json.loads((producer/'config/complete_projection_spec_v1.json').read_text());scope=Path(pspec['qualification_inputs']);plan=json.loads((scope/'outputs/complete_source_projection_qualification_inputs_v1.json').read_text());inputs=Path(plan['reference_input_stage']);models=Path(plan['model_stage']);mspec=json.loads((models/'config/common_reference_model_spec_v1.json').read_text());reader=NativeArchiveCollection(mspec['archive_manifest'])
    env=os.environ.copy();env.update(ATLAS_SYMPHONY_WHEEL=mspec['wheel_path'],ATLAS_SYMPHONY_SHA256=mspec['wheel_sha256'])
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_single_source_symphony_replay_v1.py')],env=env,capture_output=True,text=True,timeout=180);(root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 2 tests' not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Independent real-library ridge replay tests failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':2,'no_skips':True})
    indexfile=producer/'outputs/complete_source_projection_qualification_index_v1.json'
    if sha256_file(indexfile)!=t['index_sha256']:raise ValueError('Projection index changed')
    index=json.loads(indexfile.read_text());source_plan={r['source']:r for r in plan['sources']}
    if {(r['source'],r['model_id']) for r in index}!={(ds,m) for ds in source_plan for m in plan['selected_models']} or len(index)!=len(source_plan)*len(plan['selected_models']):raise ValueError('Complete source/model scope differs')
    am=json.loads((models/'outputs/reference_union_axis_metadata_v1.json').read_text())
    with np.load(models/'outputs/reference_union_axes_and_native_readouts_v1.npz',allow_pickle=False) as f:ref_source=f['source_code'];ref_position=f['native_source_position']
    reports=[]
    for record in index:
        storage_gate(root,cfg);ds=record['source'];identifier=record['model_id'];p=producer/'outputs'/record['metadata'];a=producer/'outputs'/record['arrays']
        if sha256_file(p)!=record['metadata_sha256'] or sha256_file(a)!=record['arrays_sha256']:raise ValueError('Projection output bytes changed')
        meta=json.loads(p.read_text());source=source_plan[ds];cm,side=read_source(inputs,source['reference_row_axis_record']);model=plan['selected_models'][identifier]
        mp=models/model['model']
        if sha256_file(mp)!=model['model_sha256']:raise ValueError('Projection reference parameters changed')
        with np.load(mp,allow_pickle=False) as f:v={k:f[k].copy() for k in ['genes','mean','loadings','C','Nr','sigma']}
        with np.load(a,allow_pickle=False) as f:values={k:f[k] for k in f.files}
        n=len(reader.rows(ds));np.testing.assert_array_equal(values['native_row'],reader.rows(ds));np.testing.assert_array_equal(values['native_row'],side['native_row'])
        if meta['native_ids_sha256']!=text_axis_sha256(reader.native_ids(ds)) or n!=source['n_archive_rows']:raise ValueError('Projection source identity or denominator changed')
        domain=np.ones(n,bool) if model['route']=='global_overview' else side['guarded_route_code']==FAMILIES.index(model['route']);positive=values['native_total_RNA']>0
        np.testing.assert_array_equal(values['domain_candidate'],domain);np.testing.assert_array_equal(values['fit_pool_candidate'],side['fit_pool_candidate']);np.testing.assert_array_equal(values['model_gene_signal'],values['n_detected_model_genes']>0)
        eligible=domain&positive&side['fit_pool_candidate']&values['model_gene_signal'];raw=values['raw_PCA'];corrected=values['candidate_Harmony'];raw_rows=domain&positive
        if raw.shape!=(n,v['loadings'].shape[1]) or corrected.shape!=raw.shape or not np.isfinite(raw[raw_rows]).all() or not np.isnan(raw[~raw_rows]).all():raise ValueError('Raw projection axes/finite scope invalid')
        numerical=np.zeros(n,bool);q=np.flatnonzero(raw_rows);numerical[q]=(np.linalg.norm(raw[q],axis=1)>1e-10)&(raw[q].max(axis=1)>0);included=eligible&numerical
        np.testing.assert_array_equal(values['adaptation_eligible'],eligible);np.testing.assert_array_equal(values['projection_numerical_supported'],numerical);np.testing.assert_array_equal(values['adaptation_included'],included)
        exposure=side['global_reference_fit'] if model['route']=='global_overview' else side['lineage_reference_fit']&domain;np.testing.assert_array_equal(values['reference_fit_exposure'],exposure)
        if not np.isfinite(corrected[included]).all() or not np.isnan(corrected[~included]).all() or np.any(included&~side['fit_pool_candidate']):raise ValueError('Excluded/unsupported query entered corrected coordinates')
        # Fit once from the complete included query domain. Closed-form ridge
        # algebra and chunked sums differ from the library's per-cluster inverse.
        maximum=0.;rms=0.
        if included.any():
            qr=raw[included];offsets=fit_single_source_offsets(qr,v['C'],v['Nr'],v['sigma'],1.,cfg['replay_block_rows']);expected=apply_single_source_offsets(qr,v['C'],v['sigma'],offsets,cfg['replay_block_rows']);observed=corrected[included];maximum=float(np.max(np.abs(expected-observed)));rms=float(np.sqrt(np.mean((expected-observed)**2)))
            np.testing.assert_allclose(expected,observed,atol=cfg['corrected_replay_atol'],rtol=cfg['corrected_replay_rtol'],err_msg='Complete-domain ridge replay differs; do not refit arbitrary sub-batches');del expected,observed,qr,offsets
        # Direct native count/normalization/projection checks include archive
        # exclusions, reference exposures, no-signal and both source endpoints.
        selected=np.unique(np.r_[np.linspace(0,n-1,min(64,n),dtype=np.int64),np.flatnonzero(domain&~side['fit_pool_candidate'])[:8],np.flatnonzero(exposure)[:8],np.flatnonzero(domain&~values['model_gene_signal'])[:8]])
        native=reader.store(ds);key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key'];counts=gather_csr_rows(native[key],reader.rows(ds)[selected]);total=np.asarray(counts.sum(axis=1,dtype=np.float64)).ravel();np.testing.assert_array_equal(total,values['native_total_RNA'][selected]);lookup=reader.feature_lookup(ds)
        if any(len(lookup.get(str(g),set()))!=1 for g in v['genes']):raise ValueError('Native projection genes missing or ambiguous')
        columns=[next(iter(lookup[str(g)])) for g in v['genes']];mc=counts[:,columns];det=np.asarray((mc>0).sum(axis=1)).ravel();np.testing.assert_array_equal(det,values['n_detected_model_genes'][selected]);fraction=np.divide(np.asarray(mc.sum(axis=1,dtype=np.float64)).ravel(),total,out=np.zeros(len(total)),where=total>0);np.testing.assert_allclose(fraction,values['native_RNA_fraction_in_model_genes'][selected],rtol=1e-6,atol=1e-7)
        log,_=normalized_log_counts(mc,total);direct=(log.toarray().astype(np.float64)-v['mean'].astype(np.float64))@v['loadings'].astype(np.float64);valid=raw_rows[selected];np.testing.assert_allclose(direct[valid],raw[selected[valid]],atol=2e-4,rtol=2e-5)
        geometry=producer/'outputs'/meta['geometry_probes']
        if sha256_file(geometry)!=meta['geometry_probes_sha256']:raise ValueError('Projection geometric probes changed')
        with np.load(geometry,allow_pickle=False) as f:g={k:f[k] for k in f.files}
        probe=g['source_native_position']
        if len(set(probe))!=len(probe) or np.any(probe<0) or np.any(probe>=n):raise ValueError('Geometric probe axes invalid')
        cp=models/'outputs'/model['coordinates']
        with np.load(cp,allow_pickle=False) as f:refrows=f['reference_union_row'];rz=f['raw_PCA_fit_embedding'];cz=f['reference_embedding']
        excluded=np.full(n,-1,np.int64);match=ref_source[refrows]==am['sources'].index(ds);excluded[ref_position[refrows[match]]]=np.flatnonzero(match);distance_error=0.
        for name,z,reference,good in [('raw',raw,rz,raw_rows&values['model_gene_signal']),('candidate_corrected',corrected,cz,included)]:
            vector=g['nearest_other_reference_raw_distance' if name=='raw' else 'nearest_other_reference_candidate_corrected_distance'];check=good[probe]
            np.testing.assert_array_equal(np.isfinite(vector),check);summary=meta['raw_probe_geometry' if name=='raw' else 'corrected_probe_geometry']
            if int(check.sum())!=summary['n_numeric_probe_rows'] or int((vector[check]>summary['reference_leave_self_out_radius_q99']).sum())!=summary['n_probe_rows_above_reference_radius_q99']:raise ValueError('Geometric probe qualification denominator differs')
            take=np.flatnonzero(check);take=take[np.unique(np.linspace(0,len(take)-1,min(16,len(take)),dtype=np.int64))] if len(take) else take
            if len(take):
                distances=cdist(z[probe[take]].astype(float),reference.astype(float))
                for i,qp in enumerate(probe[take]):
                    if excluded[qp]>=0:distances[i,excluded[qp]]=np.inf
                nearest=distances.min(axis=1);distance_error=max(distance_error,float(np.max(np.abs(nearest-vector[take]))));np.testing.assert_allclose(nearest,vector[take],rtol=2e-5,atol=2e-4)
        report={'source':ds,'model_id':identifier,'n_archive_rows':n,'n_included_complete_query_rows_replayed':int(included.sum()),'n_archive_excluded_domain_rows_raw_only':int((domain&~side['fit_pool_candidate']).sum()),'whole_domain_corrected_replay_max_absolute_error':maximum,'whole_domain_corrected_replay_RMS_error':rms,'n_direct_native_rows':len(selected),'geometric_probe_max_absolute_error':distance_error,'no_per_chunk_batch_refit':True,'identity_or_biological_transfer_not_accepted_by_numerical_audit':True}
        reports.append(report);save(root/'checkpoints'/(ds+'__'+identifier+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),**report});reader._fields.clear();del side,values,raw,corrected,counts,mc,log,rz,cz,v;gc.collect()
    output=root/'outputs/complete_source_projection_audit_v1.json';save(output,{'records':reports,'n_complete_sources':len(source_plan),'full_archive_qualification_complete':False,'source_comparability_and_identity_release_pending':True})
    save(root/'provenance/terminal_v1.json',{'status':'COMPLETE_SOURCE_PROJECTION_MASKS_NATIVE_PROBES_AND_WHOLE_QUERY_RIDGE_REPLAY_AUDITED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_source_model_records':len(reports),'output_sha256':sha256_file(output),'producer_terminal_sha256':sha256_file(producer/'provenance/terminal_v1.json'),'wall_seconds':time.monotonic()-started,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
