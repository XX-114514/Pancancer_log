"""Native whole-source projection qualification with frozen models and scope.

No expression matrices or authoritative identities are written. Query adaptation
is performed once per complete eligible source/model domain, never per block.
"""
import argparse,datetime,gc,json,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
import psutil
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors
from annotation_representation_pilot_v1 import normalized_log_counts,save
from audit_annotation_source_holdout_v1 import restore_reference
from audit_full_source_local_graphs_v1 import index_at
from extract_full_source_local_marker_evidence_v1 import storage_gate
from frozen_symphony_reference_v1 import load_symphonypy,reference_fingerprint
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
from native_lineage_routing_v1 import FAMILIES
from native_panel_aggregate_v1 import gather_csr_rows
from project_frozen_query_domain_v1 import project_domain
from summarize_full_native_annotation_evidence_v1 import read_source


def save_arrays(path,arrays):
    with path.open('xb') as f:np.savez_compressed(f,**arrays)
    with np.load(path,allow_pickle=False) as f:
        if set(f.files)!=set(arrays):raise ValueError('Projection serialization keys differ')
        for key,value in arrays.items():np.testing.assert_array_equal(value,f[key])


def nearest_excluding(neighbors,z,excluded):
    d,j=neighbors.kneighbors(z);out=np.empty(len(z),float)
    for i in range(len(z)):
        allowed=j[i]!=excluded[i]
        if not allowed.any():raise ValueError('No other reference row in nearest candidates')
        out[i]=d[i,allowed].min()
    return out


def probe_geometry(reference_coords,query_coords,probe,eligible,excluded_reference,cfg):
    count=min(cfg['reference_radius_calibration_rows'],len(reference_coords));cal=np.unique(np.linspace(0,len(reference_coords)-1,count,dtype=np.int64))
    near=NearestNeighbors(n_neighbors=2,n_jobs=4,algorithm='brute').fit(reference_coords)
    radius=nearest_excluding(near,reference_coords[cal],cal);threshold=float(np.quantile(radius,.99));out=np.full(len(probe),np.nan,np.float64)
    valid=eligible[probe]&np.isfinite(query_coords[probe]).all(axis=1)
    if valid.any():out[valid]=nearest_excluding(near,query_coords[probe[valid]],excluded_reference[probe[valid]])
    return out,{'n_reference_calibration_rows':len(cal),'reference_leave_self_out_radius_q99':threshold,'n_numeric_probe_rows':int(valid.sum()),'n_probe_rows_above_reference_radius_q99':int((out[valid]>threshold).sum()),'radius_threshold_is_relative_geometry_not_biological_or_identity_confidence':True,'reference_and_query_may_share_training_exposure_not_independent_validation':True}


def main(root):
    start=time.monotonic();cfg=json.loads((root/'config/complete_projection_spec_v1.json').read_text())
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Whole-source projection input changed')
    scope=Path(cfg['qualification_inputs']);sf=scope/'outputs/complete_source_projection_qualification_inputs_v1.json';st=json.loads((scope/'provenance/terminal_v1.json').read_text())
    if sha256_file(sf)!=st['output_sha256']:raise ValueError('Projection qualification scope hash differs')
    plan=json.loads(sf.read_text());models=Path(plan['model_stage']);inputs=Path(plan['reference_input_stage']);mspec=json.loads((models/'config/common_reference_model_spec_v1.json').read_text());evidence=Path(mspec['native_evidence'])
    for b in plan['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Bound model/input auditor terminal changed')
    env=dict(__import__('os').environ);env.update(ATLAS_SYMPHONY_WHEEL=mspec['wheel_path'],ATLAS_SYMPHONY_SHA256=mspec['wheel_sha256'])
    for test,count in [('test_project_frozen_query_domain_v1.py',2),('test_frozen_symphony_reference_v1.py',4)]:
        result=subprocess.run([sys.executable,'-B',str(root/'tests'/test)],env=env,capture_output=True,text=True,timeout=240);(root/'provenance'/(test+'.stderr.txt')).write_text(result.stderr)
        if result.returncode or ('Ran %d tests'%count) not in result.stderr or 'skipped' in result.stderr.lower():raise ValueError('Actual projection-domain/cached-library tests failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':6,'no_skips':True})
    symphony=load_symphonypy(mspec['wheel_path'],mspec['wheel_sha256']);reader=NativeArchiveCollection(mspec['archive_manifest']);ri=index_at(evidence,'full_archive_reference_index_v1.json')
    genes=json.loads((inputs/'outputs/common_measured_model_feature_candidates_v1.json').read_text())['genes'];am=json.loads((models/'outputs/reference_union_axis_metadata_v1.json').read_text())
    with np.load(models/'outputs/reference_union_axes_and_native_readouts_v1.npz',allow_pickle=False) as f:ref_source=f['source_code'];ref_position=f['native_source_position']
    selected_models=plan['selected_models'];model_values={}
    for identifier,record in selected_models.items():
        path=models/record['model']
        if sha256_file(path)!=record['model_sha256']:raise ValueError('Selected projection model changed')
        with np.load(path,allow_pickle=False) as f:model_values[identifier]={k:f[k].copy() for k in ['genes','feature_candidate_indices']}
        v=model_values[identifier];np.testing.assert_array_equal(v['genes'],np.asarray(genes)[v['feature_candidate_indices']])
    records=[];native_nnz=0
    for source_entry in plan['sources']:
        began=time.monotonic();ds=source_entry['source'];storage_gate(root,cfg);meta,side=read_source(inputs,source_entry['reference_row_axis_record']);em,ev=read_source(evidence,ri[ds]);rows=reader.rows(ds);n=len(rows)
        np.testing.assert_array_equal(rows,side['native_row']);np.testing.assert_array_equal(rows,ev['native_row'])
        if n!=source_entry['n_archive_rows'] or text_axis_sha256(reader.native_ids(ds))!=source_entry['native_ids_sha256'] or em['native_ids_sha256']!=source_entry['native_ids_sha256']:raise ValueError('Complete source native identity/denominator changed')
        required=n*len(genes)*16+cfg['RAM_workspace_reserve_bytes']
        if psutil.virtual_memory().available<required:raise MemoryError('Whole-source sparse/model workspace RAM gate failed')
        cf=inputs/'outputs'/source_entry['native_feature_columns']
        if sha256_file(cf)!=source_entry['native_feature_columns_sha256']:raise ValueError('Source feature columns changed')
        with np.load(cf,allow_pickle=False) as f:columns=f['native_columns']
        lookup=reader.feature_lookup(ds)
        if len(columns)!=len(genes) or any(lookup.get(g,set())!={int(c)} for g,c in zip(genes,columns)):raise ValueError('Measured common native gene axes differ')
        native=reader.store(ds);key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key'];chunks=[];detected={k:np.zeros(n,np.int16) for k in selected_models};fraction={k:np.zeros(n,np.float32) for k in selected_models};totals=np.zeros(n,np.float64)
        for low in range(0,n,cfg['read_block_rows']):
            high=min(n,low+cfg['read_block_rows']);raw=gather_csr_rows(native[key],rows[low:high]);native_nnz+=raw.nnz
            if not np.isfinite(raw.data).all() or np.any(raw.data<0) or not np.equal(raw.data,np.floor(raw.data)).all():raise ValueError('Whole-source native count domain changed')
            total=np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel();np.testing.assert_array_equal(total,ev['native_total_RNA'][low:high]);totals[low:high]=total
            common=raw[:,columns];log,positive=normalized_log_counts(common,total);chunks.append(log)
            for identifier,v in model_values.items():
                model_counts=common[:,v['feature_candidate_indices']];detected[identifier][low:high]=np.asarray((model_counts>0).sum(axis=1)).ravel();mass=np.asarray(model_counts.sum(axis=1,dtype=np.float64)).ravel();fraction[identifier][low:high]=np.divide(mass,total,out=np.zeros(len(total)),where=total>0)
            if low==0 or high==n or (low//cfg['read_block_rows'])%32==0:save(root/'checkpoints'/(ds+'_read_'+str(high)+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'stage':'whole_source_native_reference_feature_read','completed_native_rows':high,'n_archive_rows':n})
        x=sp.vstack(chunks,format='csr');del chunks,raw,common,log,model_counts;gc.collect()
        if x.shape!=(n,len(genes)):raise ValueError('Whole source sparse query matrix axis differs')
        for identifier,record in selected_models.items():
            storage_gate(root,cfg);v=model_values[identifier];ref,rm=restore_reference(models/record['model'],v['genes'].tolist(),v['feature_candidate_indices'])
            if reference_fingerprint(ref)!=record['reference_fingerprint']:raise ValueError('Restored projection reference fingerprint differs')
            family=record['route'];domain=np.ones(n,bool) if family=='global_overview' else side['guarded_route_code']==FAMILIES.index(family);positive=totals>0;model_x=x[:,v['feature_candidate_indices']]
            result=project_domain(ref,model_x,v['genes'].tolist(),ds,domain,side['fit_pool_candidate'],positive,symphony=symphony)
            np.testing.assert_array_equal(result['model_gene_signal'],detected[identifier]>0)
            if not np.isfinite(fraction[identifier]).all() or np.any(fraction[identifier]<0) or np.any(fraction[identifier]>1+1e-6):raise ValueError('Native model RNA mass fractions invalid')
            fit_exposure=side['global_reference_fit'] if family=='global_overview' else side['lineage_reference_fit']&domain
            candidates=np.flatnonzero(domain&positive);probe=candidates[np.unique(np.linspace(0,len(candidates)-1,min(cfg['query_geometry_probe_rows'],len(candidates)),dtype=np.int64))] if len(candidates) else np.array([],np.int64)
            stress=np.unique(np.r_[np.flatnonzero(domain&~side['fit_pool_candidate'])[:16],np.flatnonzero(fit_exposure)[:16],np.flatnonzero(domain&~result['model_gene_signal'])[:16]])
            supplemental=np.setdiff1d(stress,probe);probe=np.r_[probe,supplemental];probe_role=np.r_[np.zeros(len(probe)-len(supplemental),np.int8),np.ones(len(supplemental),np.int8)]
            cp=models/'outputs'/record['coordinates']
            if sha256_file(cp)!=record['coordinates_sha256']:raise ValueError('Frozen reference coordinates changed')
            with np.load(cp,allow_pickle=False) as f:ref_rows=f['reference_union_row'];raw_reference=f['raw_PCA_fit_embedding'];corrected_reference=f['reference_embedding']
            excluded=np.full(n,-1,np.int64);match=ref_source[ref_rows]==am['sources'].index(ds);excluded[ref_position[ref_rows[match]]]=np.flatnonzero(match)
            raw_supported=domain&positive&result['model_gene_signal']&np.isfinite(result['raw_PCA']).all(axis=1)
            raw_distance,raw_geometry=probe_geometry(raw_reference,result['raw_PCA'],probe,raw_supported,excluded,cfg);corrected_distance,corrected_geometry=probe_geometry(corrected_reference,result['candidate_Harmony'],probe,result['adaptation_included'],excluded,cfg)
            arrays={**result,'native_row':rows,'domain_candidate':domain,'fit_pool_candidate':side['fit_pool_candidate'],'reference_fit_exposure':fit_exposure,'native_total_RNA':totals,'n_detected_model_genes':detected[identifier],'native_RNA_fraction_in_model_genes':fraction[identifier]}
            out=root/'outputs'/(ds+'__'+identifier+'_projection_v1.npz');save_arrays(out,arrays)
            probes=root/'outputs'/(ds+'__'+identifier+'_geometry_probes_v1.npz');save_arrays(probes,{'source_native_position':probe,'probe_role':probe_role,'nearest_other_reference_raw_distance':raw_distance,'nearest_other_reference_candidate_corrected_distance':corrected_distance})
            info={'source':ds,'model_id':identifier,'n_archive_rows':n,'n_domain_rows':int(domain.sum()),'n_adaptation_eligible':int(result['adaptation_eligible'].sum()),'n_adaptation_included':int(result['adaptation_included'].sum()),'n_domain_archive_exclusions_raw_only':int((domain&~side['fit_pool_candidate']).sum()),'native_ids_sha256':source_entry['native_ids_sha256'],'arrays':out.name,'arrays_sha256':sha256_file(out),'geometry_probes':probes.name,'geometry_probes_sha256':sha256_file(probes),'raw_probe_geometry':raw_geometry,'corrected_probe_geometry':corrected_geometry,'query_scope':'complete_source_eligible_positive_RNA_and_model_signal_domain_one_adaptation_call','raw_reference_fit_transform_not_assumed_identical_to_query_transform':True,'source_modality_capture_patient_and_clinical_transfer_not_accepted_by_numeric_projection':True,'source_context':source_entry['source_modality_context'],'reference_fingerprint_unchanged':reference_fingerprint(ref)==record['reference_fingerprint'],'all_reference_sources_were_in_training_transductive_qualification':True,'no_new_cell_identity':True}
            mp=out.with_suffix('.json');save(mp,info);records.append({'source':ds,'model_id':identifier,'n_archive_rows':n,'metadata':mp.name,'metadata_sha256':sha256_file(mp),'arrays':out.name,'arrays_sha256':sha256_file(out)})
            save(root/'checkpoints'/(ds+'__'+identifier+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'model_id':identifier,'stage':'whole_source_model_projection_and_probe_geometry_complete','n_adaptation_included':int(result['adaptation_included'].sum())})
            del ref,model_x,result,arrays,raw_reference,corrected_reference;gc.collect()
        save(root/'checkpoints'/(ds+'_complete.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'stage':'complete_source_all_selected_models_projected','wall_seconds':time.monotonic()-began});reader._fields.clear();del x,ev,side;gc.collect()
    output=root/'outputs/complete_source_projection_qualification_index_v1.json';save(output,records)
    save(root/'provenance/terminal_v1.json',{'status':'COMPLETE_SOURCE_PROJECTION_AND_BOUNDED_GEOMETRY_COMPUTED_INDEPENDENT_AUDIT_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_complete_sources':len(plan['sources']),'n_source_model_records':len(records),'n_archive_source_rows':sum(r['n_archive_rows'] for r in plan['sources']),'native_nnz_read':int(native_nnz),'index_sha256':sha256_file(output),'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'new_final_identity':False,'full_archive_mapping_complete':False,'native_expression_matrices_saved':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
