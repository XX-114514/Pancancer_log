"""Audit explicit lineage axes, native probes and complete-source query fits.

The producer is required to have actually exited successfully. Numerical
acceptance is separate from neighbor precision, identity and source eligibility.
No expression matrix, individual record or learned parameter is published.
"""
import argparse, datetime, gc, json, os, resource, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from annotation_representation_pilot_v1 import save, normalized_log_counts
from audit_full_source_local_graphs_v1 import index_at
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, text_axis_sha256
from native_lineage_routing_v1 import FAMILIES
from native_panel_aggregate_v1 import gather_csr_rows
from single_source_symphony_replay_v1 import fit_single_source_offsets, apply_single_source_offsets
from summarize_full_native_annotation_evidence_v1 import read_source


def exact_distances(query, reference, excluded):
    out = np.empty(len(query), float)
    for low in range(0, len(query), 16):
        distance = cdist(query[low:low+16].astype(float), reference.astype(float))
        for j, index in enumerate(excluded[low:low+16]):
            if index >= 0:
                distance[j, index] = np.inf
        out[low:low+len(distance)] = distance.min(axis=1)
    return out


def main(root):
    began = time.monotonic()
    cfg = json.loads((root/'config/full_projection_audit_spec_v1.json').read_text())
    producer = Path(cfg['producer'])
    supervisor = json.loads((producer/'provenance/supervisor_terminal_v1.json').read_text())
    terminal = json.loads((producer/'provenance/terminal_v1.json').read_text())
    if supervisor['actual_child_exit_code'] or supervisor['wall_timeout']:
        raise ValueError('Actual successful complete producer required')
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']:
            raise ValueError('Full projection audit binding changed')
    spec = json.loads((producer/'config/full_archive_projection_spec_v1.json').read_text())
    scope = Path(spec['projection_inputs'])
    plan = json.loads((scope/'outputs/full_archive_projection_inputs_v1.json').read_text())
    inputs, models = Path(plan['reference_input_stage']), Path(plan['model_stage'])
    evidence = Path(plan['native_evidence'])
    evidence_index = index_at(evidence, 'full_archive_reference_index_v1.json')
    mspec = json.loads((models/'config/common_reference_model_spec_v1.json').read_text())
    reader = NativeArchiveCollection(plan['archive_manifest'])
    env = os.environ.copy()
    env.update(ATLAS_SYMPHONY_WHEEL=mspec['wheel_path'], ATLAS_SYMPHONY_SHA256=mspec['wheel_sha256'])
    test = subprocess.run([sys.executable, '-B', str(root/'tests/test_single_source_symphony_replay_v1.py')], env=env, capture_output=True, text=True, timeout=180)
    (root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 2 tests' not in test.stderr or 'skipped' in test.stderr.lower():
        raise ValueError('Independent actual-library ridge replay tests failed')
    save(root/'provenance/numerical_tests_v1.json', {'actual_exit_code':0, 'n_tests':2, 'no_skips':True})
    indexfile = producer/'outputs/full_archive_projection_index_v1.json'
    if sha256_file(indexfile) != terminal['index_sha256']:
        raise ValueError('Complete projection index changed')
    index = json.loads(indexfile.read_text())
    sources = {r['source']:r for r in plan['sources']}
    expected = {(ds, m) for ds in sources for m in plan['selected_models']}
    if {(r['source'], r['model_id']) for r in index} != expected or len(index) != len(expected):
        raise ValueError('Full source/model Cartesian coverage differs')
    if terminal['n_archive_rows'] != sum(r['n_archive_rows'] for r in sources.values()):
        raise ValueError('Full archive denominator differs')
    am = json.loads((models/'outputs/reference_union_axis_metadata_v1.json').read_text())
    with np.load(models/'outputs/reference_union_axes_and_native_readouts_v1.npz', allow_pickle=False) as f:
        refsource, refposition = f['source_code'], f['native_source_position']
    reports = []
    for ds in sorted(sources):
        source = sources[ds]
        _, side = read_source(inputs, source['reference_row_axis_record'])
        _, native_evidence = read_source(evidence, evidence_index[ds])
        native_rows = reader.rows(ds)
        np.testing.assert_array_equal(native_rows, side['native_row'])
        np.testing.assert_array_equal(native_rows, native_evidence['native_row'])
        if len(native_rows) != source['n_archive_rows'] or text_axis_sha256(reader.native_ids(ds)) != source['native_ids_sha256']:
            raise ValueError('Full archive native identity changed')
        for record in (r for r in index if r['source'] == ds):
            storage_gate(root, cfg)
            identifier = record['model_id']
            model = plan['selected_models'][identifier]
            meta_path = producer/'outputs'/record['metadata']
            if sha256_file(meta_path) != record['metadata_sha256']:
                raise ValueError('Projection metadata changed')
            meta = json.loads(meta_path.read_text())
            q = np.arange(len(native_rows)) if model['route'] == 'global_overview' else np.flatnonzero(side['guarded_route_code'] == FAMILIES.index(model['route']))
            if len(q) != record['n_model_domain_rows'] or meta['n_model_domain_rows'] != len(q) or meta['source'] != ds or meta['model_id'] != identifier:
                raise ValueError('Model domain denominator/key differs')
            if not len(q):
                if 'arrays' in record or meta['status'] != 'no_guarded_rows_in_this_model_route_global_archive_retained':
                    raise ValueError('Empty lineage fabricated arrays')
                reports.append({'source':ds, 'model_id':identifier, 'n_domain_rows':0, 'status':'explicit_empty_domain_checked', 'neighbor_precision_passed':True})
                continue
            array_path = producer/'outputs'/record['arrays']
            if sha256_file(array_path) != record['arrays_sha256'] or meta['arrays_sha256'] != record['arrays_sha256']:
                raise ValueError('Full projection arrays changed')
            with np.load(array_path, allow_pickle=False) as f:
                v = {k:f[k] for k in f.files}
            if any(len(a) != len(q) for a in v.values()):
                raise ValueError('Explicit domain array lengths differ')
            np.testing.assert_array_equal(v['source_native_position'], q)
            np.testing.assert_array_equal(v['native_row'], native_rows[q])
            np.testing.assert_array_equal(v['native_total_RNA'], native_evidence['native_total_RNA'][q])
            np.testing.assert_array_equal(v['fit_pool_candidate'], side['fit_pool_candidate'][q])
            exposure = side['global_reference_fit'][q] if model['route'] == 'global_overview' else side['lineage_reference_fit'][q]
            np.testing.assert_array_equal(v['reference_fit_exposure'], exposure)
            if meta['native_ids_sha256'] != source['native_ids_sha256'] or meta['array_axis'] != 'explicit_source_native_position_model_domain_rows':
                raise ValueError('Projection source axis identity changed')
            model_path = models/model['model']
            if sha256_file(model_path) != model['model_sha256']:
                raise ValueError('Frozen model bytes changed')
            with np.load(model_path, allow_pickle=False) as f:
                params = {k:f[k].copy() for k in f.files}
            raw, corrected = v['raw_PCA'], v['candidate_Harmony']
            positive = v['native_total_RNA'] > 0
            np.testing.assert_array_equal(v['model_gene_signal'], v['n_detected_model_genes'] > 0)
            if raw.shape != (len(q), params['loadings'].shape[1]) or corrected.shape != raw.shape or not np.isfinite(raw[positive]).all() or not np.isnan(raw[~positive]).all():
                raise ValueError('Raw coordinate shape or measurement mask differs')
            raw_only = model['reference_kind'] == 'raw_PCA_context_preserving'
            if meta['raw_context_only'] != raw_only:
                raise ValueError('Context model policy changed')
            numerical = np.zeros(len(q), bool)
            eligible = np.zeros(len(q), bool)
            if not raw_only:
                numerical[positive] = (np.linalg.norm(raw[positive], axis=1) > 1e-10) & (raw[positive].max(axis=1) > 0)
                eligible = positive & v['fit_pool_candidate'] & v['model_gene_signal']
            included = eligible & numerical
            if int(eligible.sum()) != meta['n_adaptation_eligible'] or int(included.sum()) != meta['n_adaptation_included'] or int((~v['fit_pool_candidate']).sum()) != meta['n_archive_excluded_domain_rows_raw_only']:
                raise ValueError('Adaptation denominator metadata differs')
            for key, value in [('projection_numerical_supported', numerical), ('adaptation_eligible', eligible), ('adaptation_included', included)]:
                np.testing.assert_array_equal(v[key], value)
            if not np.isfinite(corrected[included]).all() or not np.isnan(corrected[~included]).all():
                raise ValueError('Raw-context, excluded or unsupported row acquired correction')
            maximum, rms = 0., 0.
            if included.any():
                qr = raw[included]
                offsets = fit_single_source_offsets(qr, params['C'], params['Nr'], params['sigma'], 1., cfg['replay_block_rows'])
                replay = apply_single_source_offsets(qr, params['C'], params['sigma'], offsets, cfg['replay_block_rows'])
                maximum = float(np.max(np.abs(replay-corrected[included])))
                rms = float(np.sqrt(np.mean((replay-corrected[included])**2)))
                np.testing.assert_allclose(replay, corrected[included], atol=cfg['corrected_replay_atol'], rtol=cfg['corrected_replay_rtol'])
                del qr, offsets, replay
            selected = np.unique(np.r_[np.linspace(0,len(q)-1,min(64,len(q)),dtype=np.int64), np.flatnonzero(~v['fit_pool_candidate'])[:8], np.flatnonzero(exposure)[:8], np.flatnonzero(~v['model_gene_signal'])[:8], np.flatnonzero(~positive)[:8]])
            store = reader.store(ds)
            key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
            count = gather_csr_rows(store[key], native_rows[q[selected]])
            total = np.asarray(count.sum(axis=1,dtype=np.float64)).ravel()
            np.testing.assert_array_equal(total, v['native_total_RNA'][selected])
            lookup = reader.feature_lookup(ds)
            if any(len(lookup.get(str(g),set())) != 1 for g in params['genes']):
                raise ValueError('Native model feature missing/ambiguous')
            columns = [next(iter(lookup[str(g)])) for g in params['genes']]
            mc = count[:,columns]
            np.testing.assert_array_equal(np.asarray((mc>0).sum(axis=1)).ravel(),v['n_detected_model_genes'][selected])
            fraction = np.divide(np.asarray(mc.sum(axis=1,dtype=np.float64)).ravel(),total,out=np.zeros(len(total)),where=total>0)
            np.testing.assert_allclose(fraction,v['native_RNA_fraction_in_model_genes'][selected],atol=1e-7,rtol=1e-6)
            log,_ = normalized_log_counts(mc,total)
            direct = (log.toarray().astype(float)-params['mean'].astype(float)) @ params['loadings'].astype(float)
            np.testing.assert_allclose(direct[positive[selected]],raw[selected[positive[selected]]],atol=2e-4,rtol=2e-5)
            cp = models/'outputs'/model['coordinates']
            if sha256_file(cp) != model['coordinates_sha256']:
                raise ValueError('Frozen reference coordinates changed')
            with np.load(cp,allow_pickle=False) as f:
                refrows, rz, cz = f['reference_union_row'], f['raw_PCA_fit_embedding'], f['reference_embedding']
            same = refsource[refrows] == am['sources'].index(ds)
            exclude = np.full(len(native_rows),-1,np.int64)
            exclude[refposition[refrows[same]]] = np.flatnonzero(same)
            geometry_reports = []
            for basis,z,reference,good in [('raw',raw,rz,positive&v['model_gene_signal']),('candidate_corrected',corrected,cz,included)]:
                profile = next(p for p in meta['geometry_profiles'] if p['basis'] == basis)
                distance = v[basis+'_reference_distance']
                np.testing.assert_array_equal(v[basis+'_reference_distance_available'],good)
                np.testing.assert_array_equal(np.isfinite(distance),good)
                radius = profile['reference_other_neighbor_radius_q99']
                outside = good & (distance>radius) if radius is not None else np.zeros(len(q),bool)
                np.testing.assert_array_equal(v[basis+'_outside_reference_radius_q99'],outside)
                if int(good.sum()) != profile['n_queried_model_domain_rows'] or int(outside.sum()) != profile['n_rows_above_reference_radius_q99']:
                    raise ValueError('Geometry denominator/flag differs')
                if not good.any():
                    if radius is not None or profile['n_exact_probes']:
                        raise ValueError('Unavailable geometry fabricated support')
                    geometry_reports.append(True)
                    continue
                local = np.flatnonzero(good)
                count_probes = min(spec['n_independent_exact_probes'],len(local))
                base = np.unique(np.linspace(0,len(local)-1,count_probes,dtype=np.int64))
                midpoint = np.unique(np.floor((np.arange(count_probes)+.5)*len(local)/count_probes).astype(np.int64))
                passed = True
                for positions, summary in [(base,profile),(np.setdiff1d(midpoint,base),profile.get('additional_midpoint_exact_probes'))]:
                    if not len(positions):
                        continue
                    if summary is None or summary['n_exact_probes'] != len(positions):
                        raise ValueError('Exact probe domain/count differs')
                    take = local[positions]
                    exact = exact_distances(z[take],reference,exclude[q[take]])
                    error = distance[take]-exact
                    if np.any(error < -(2e-4+2e-5*exact)):
                        raise ValueError('Geometry distance under independent exact minimum')
                    matched = float(np.isclose(distance[take],exact,atol=2e-4,rtol=2e-5).mean())
                    if matched != summary['distance_exact_match_fraction']:
                        raise ValueError('Producer exact-distance fraction not independently reproduced')
                    passed = passed and matched >= spec['neighbor_parameters']['minimum_exact_probe_fraction']
                if passed != profile['exact_probe_gate_passed']:
                    raise ValueError('Declared neighbor gate differs')
                geometry_reports.append(passed)
            precision = all(geometry_reports)
            if precision != record['neighbor_precision_all_queried_bases_passed'] or precision != meta['neighbor_precision_all_queried_bases_passed']:
                raise ValueError('Record neighbor gate differs')
            report = {'source':ds,'model_id':identifier,'n_domain_rows':len(q),'n_complete_query_rows_replayed':int(included.sum()),'raw_context_only':raw_only,'max_corrected_replay_error':maximum,'RMS_corrected_replay_error':rms,'n_native_probe_rows':len(selected),'neighbor_precision_passed':precision,'status':'full_domain_numerical_checks_passed_identity_release_pending'}
            reports.append(report)
            save(root/'checkpoints'/(ds+'__'+identifier+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),**report})
            del v,raw,corrected,params,count,mc,log,direct,rz,cz
            gc.collect()
        reader._fields.clear()
        del side, native_evidence
        gc.collect()
    precision = all(r['neighbor_precision_passed'] for r in reports)
    if precision != terminal['all_queried_neighbor_precision_gates_passed']:
        raise ValueError('Terminal precision acceptance differs')
    out = root/'outputs/full_archive_projection_audit_v1.json'
    save(out,{'records':reports,'n_sources':len(sources),'all_neighbor_precision_gates_passed':precision,'new_final_identity':False,'clinical_source_views_or_standalone_Zarr_accepted':False})
    save(root/'provenance/terminal_v1.json',{'status':'FULL_ARCHIVE_DOMAIN_AXES_NATIVE_PROBES_AND_COMPLETE_QUERY_FITS_AUDITED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(sources),'n_source_model_records':len(reports),'all_neighbor_precision_gates_passed':precision,'output_sha256':sha256_file(out),'producer_terminal_sha256':sha256_file(producer/'provenance/terminal_v1.json'),'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True)
    root=Path(p.parse_args().output)
    try:
        main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)})
        raise
