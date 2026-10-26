"""Freeze a readonly, source-backed atlas delivery and exercise every source.

The delivery is explicit about candidate identities, existing measurement
coverage and held scientific views. It does not materialize a counts Zarr.
"""
import argparse,datetime,gc,json,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
from annotation_representation_pilot_v1 import save
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from joined_native_atlas_delivery_v1 import JoinedNativeAtlas,completed
from audit_full_source_local_graphs_v1 import direct_native_rows


def main(root):
    start=time.monotonic();cfg=json.loads((root/'config/joined_delivery_spec_v1.json').read_text())
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Delivery input changed')
    for stage in cfg['successful_prerequisite_stages']:completed(stage)
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_joined_native_atlas_delivery_v1.py')],capture_output=True,text=True,timeout=120)
    (root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 4 tests' not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Actual delivery join counterexamples failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':4,'no_skips':True})
    manifest={k:cfg[k] for k in ['archive_manifest','catalog','unit_views','native_evidence','local_review','qualified_projection','projection_qualifier','bindings']}
    manifest.update(schema='joined_native_atlas_delivery_v1',native_counts_copied=False,standalone_Zarr=False,
        new_authoritative_taxonomy=False,new_clinical_or_composition_eligibility=False,
        source_local_structural_IDs_are_not_common_biological_types=True)
    mp=root/'outputs/joined_native_atlas_delivery_v1.json';save(mp,manifest)
    reader=JoinedNativeAtlas(mp);records=[];genes=['CD3D','CD3E','TRAC','LST1','LYZ','EPCAM','COL1A1','PECAM1','__UNMEASURED_SENTINEL__']
    for ds in reader.sources:
        storage_gate(root,cfg);s=reader.source(ds);n=len(s['rows']);p=np.asarray([n-1,0,min(1,n-1),n-1],np.int64)
        m=reader.metadata(ds,p);cols=m['columns'];np.testing.assert_array_equal(cols['native_row'],s['rows'][p])
        if any(len(v)!=len(p) for v in cols.values()):raise ValueError('Joined metadata row axis differs')
        for f in ['author_label','local_major','local_subcluster_id','consensus_major','consensus_celltype','global_subcluster_id','confidence','label_version','marker','marker_evidence_type','marker_origin']:
            if f not in cols:raise ValueError('Required delivery field absent')
        counts=reader.read_counts(ds,p,genes);native=reader.native;z=native.store(ds);key=counts['count_origin'];direct=direct_native_rows(z[key],s['rows'][p])
        np.testing.assert_array_equal(counts['native_total_counts'],np.asarray(direct.sum(axis=1,dtype=np.float64)).ravel())
        lookup=native.feature_lookup(ds)
        for j,gene in enumerate(genes):
            physical=lookup.get(gene,set());measured=len(physical)==1
            if counts['measured'][j]!=measured:raise ValueError('Native measurement mask changed')
            expected=direct[:,next(iter(physical))].toarray().ravel() if measured else np.zeros(len(p))
            np.testing.assert_array_equal(counts['counts'][:,j].toarray().ravel(),expected)
        projection=reader.projection(ds,'global_overview__shared',p)
        if not projection['model_domain_available'].all():raise ValueError('Global model lost archive cells')
        np.testing.assert_array_equal(projection['arrays']['native_row'],s['rows'][p])
        if counts['measured'][-1]:raise ValueError('Unknown gene was fabricated as measured')
        row={'source':ds,'n_archive_rows':n,'n_expression_units_not_patients':s['um']['n_expression_units_not_patients'],
             'n_required_metadata_columns':len(cols),'n_native_count_probe_rows':len(p),
             'native_count_and_mask_probe_exact':True,'metadata_duplicate_order_preserved':True,
             'native_fields_present':m['provenance']['original_carried_fields'],
             'native_fields_absent':m['provenance']['absent_carried_fields'],
             'RNA_state_order':s['em']['state_order'],'candidate_dictionaries':s['cm']['candidate_dictionaries'],
             'clinical_or_composition_eligibility_newly_accepted':False}
        records.append(row);save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'n_complete_sources':len(records)})
        del m,cols,counts,direct,projection;reader._cache.clear();reader.native._fields.clear();reader.native._native_arrays.clear();gc.collect()
    # Explicitly exercise every exact patch, including its complete native
    # subaxis, instead of relying on a global-coordinate sample to hit it.
    patched=0
    for entry in reader.projections['records']:
        for patch in entry['exact_neighbor_overlays']:
            with np.load(patch['repair_arrays'],allow_pickle=False) as f:v={k:f[k] for k in f.files}
            output=reader.projection(entry['source'],entry['model_id'],v['source_native_position']);b=patch['basis']
            np.testing.assert_array_equal(output['exact_neighbor_overlay_applied'],True)
            np.testing.assert_allclose(output['arrays'][b+'_reference_distance'],v['exact_nearest_other_reference_distance'],atol=0,rtol=0)
            np.testing.assert_array_equal(output['arrays'][b+'_outside_reference_radius_q99'],v['outside_original_reference_radius_q99']);patched+=1
    coverage=root/'outputs/joined_delivery_source_coverage_v1.json';save(coverage,{'sources':records,'n_sources':len(records),'n_archive_rows':sum(r['n_archive_rows'] for r in records),'n_units_not_patients':sum(r['n_expression_units_not_patients'] for r in records),'n_exact_overlay_domains_exercised':patched,'count_readback_scope':'bounded_direct_native_probes_every_source_not_allcount_recomputation','full_candidate_catalog_and_state_evidence_audits_bound':True})
    save(root/'provenance/terminal_v1.json',{'status':'ALL_SOURCE_READONLY_JOINED_NATIVE_ATLAS_DELIVERY_AND_BOUNDED_READBACK_COMPLETE_AUDIT_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(records),'n_archive_rows':sum(r['n_archive_rows'] for r in records),'manifest_sha256':sha256_file(mp),'coverage_sha256':sha256_file(coverage),'n_exact_overlay_domains_exercised':patched,'native_counts_copied':False,'standalone_Zarr':False,'new_authoritative_taxonomy':False,'new_clinical_or_composition_eligibility':False,'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as e:save(root/'provenance/failure_v1.json',{'error_class':type(e).__name__,'error':str(e)});raise
