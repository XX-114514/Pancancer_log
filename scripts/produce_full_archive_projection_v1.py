"""Selected full-archive projections with per-cell reference-distance flags.

Counts are read once per source into a temporary sparse common-feature cache.
Lineage coordinates have explicit source-position axes; the global model retains
all archive rows. No authoritative identities or clinical eligibility change.
"""
import argparse,datetime,gc,json,os,resource,subprocess,sys,time,importlib.metadata
from pathlib import Path
import numpy as np
import psutil
import scipy.sparse as sp
from annotation_representation_pilot_v1 import save,normalized_log_counts
from audit_annotation_source_holdout_v1 import restore_reference
from audit_full_source_local_graphs_v1 import index_at
from extract_full_source_local_marker_evidence_v1 import storage_gate
from frozen_reference_neighbors_v1 import FrozenReferenceNeighbors,actual_runtime_smoke
from frozen_symphony_reference_v1 import load_symphonypy,reference_fingerprint
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
from native_lineage_routing_v1 import FAMILIES
from native_panel_aggregate_v1 import gather_csr_rows
from project_frozen_query_domain_v1 import project_domain,project_raw_context
from summarize_full_native_annotation_evidence_v1 import read_source


def save_arrays(path,arrays):
    with path.open('xb') as f:np.savez_compressed(f,**arrays)
    with np.load(path,allow_pickle=False) as f:
        if set(f.files)!=set(arrays):raise ValueError('Full projection serialized keys differ')
        for k,v in arrays.items():np.testing.assert_array_equal(v,f[k])


def main(root):
    began=time.monotonic();cfg=json.loads((root/'config/full_archive_projection_spec_v1.json').read_text())
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Full archive projection binding changed')
    ir=Path(cfg['projection_inputs']);ip=ir/'outputs/full_archive_projection_inputs_v1.json';it=json.loads((ir/'provenance/terminal_v1.json').read_text())
    if sha256_file(ip)!=it['output_sha256']:raise ValueError('Full archive input scope changed')
    plan=json.loads(ip.read_text());neighbor_qualification=Path(plan['neighbor_qualification_stage']);ns=json.loads((neighbor_qualification/'provenance/supervisor_terminal_v1.json').read_text());nt=json.loads((neighbor_qualification/'provenance/terminal_v1.json').read_text())
    if ns['actual_child_exit_code'] or ns['wall_timeout'] or not nt['all_exact_probe_gates_passed']:raise ValueError('Actual complete-source neighbor precision qualification required')
    for b in plan['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Selected model or qualification request changed')
    if importlib.metadata.version('pynndescent')!=cfg['pynndescent_version']:raise ValueError('Qualified neighbor runtime version changed')
    models=Path(plan['model_stage']);inputs=Path(plan['reference_input_stage']);mspec=json.loads((models/'config/common_reference_model_spec_v1.json').read_text());evidence=Path(plan['native_evidence']);reader=NativeArchiveCollection(plan['archive_manifest']);ri=index_at(evidence,'full_archive_reference_index_v1.json')
    env=os.environ.copy();env.update(ATLAS_SYMPHONY_WHEEL=mspec['wheel_path'],ATLAS_SYMPHONY_SHA256=mspec['wheel_sha256'])
    for test,n in [('test_project_frozen_query_domain_v1.py',3),('test_frozen_symphony_reference_v1.py',4)]:
        r=subprocess.run([sys.executable,'-B',str(root/'tests'/test)],env=env,capture_output=True,text=True,timeout=240);(root/'provenance'/(test+'.stderr.txt')).write_text(r.stderr)
        if r.returncode or ('Ran %d tests'%n) not in r.stderr or 'skipped' in r.stderr.lower():raise ValueError('Actual full-domain/context/reference tests failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':7,'no_skips':True});save(root/'provenance/actual_neighbor_runtime_smoke_v1.json',actual_runtime_smoke(cfg['seed'],cfg['neighbor_parameters']))
    symphony=load_symphonypy(mspec['wheel_path'],mspec['wheel_sha256']);genes=json.loads((inputs/'outputs/common_measured_model_feature_candidates_v1.json').read_text())['genes'];am=json.loads((models/'outputs/reference_union_axis_metadata_v1.json').read_text())
    if {r['source'] for r in plan['sources']}!=set(reader.archive_partitions) or set(ri)!=set(reader.archive_partitions):raise ValueError('Full archive source coverage changed')
    with np.load(models/'outputs/reference_union_axes_and_native_readouts_v1.npz',allow_pickle=False) as f:refsource=f['source_code'];refposition=f['native_source_position']
    model_values={};search_cache={};index=[];native_nnz=0;source_records=[]
    for identifier,record in plan['selected_models'].items():
        p=models/record['model']
        if sha256_file(p)!=record['model_sha256']:raise ValueError('Selected model parameter bytes changed')
        with np.load(p,allow_pickle=False) as f:model_values[identifier]={k:f[k].copy() for k in ['genes','feature_candidate_indices','mean','loadings']}
        np.testing.assert_array_equal(model_values[identifier]['genes'],np.asarray(genes)[model_values[identifier]['feature_candidate_indices']])
    for source_entry in sorted(plan['sources'],key=lambda r:(r['n_archive_rows'],r['source'])):
        source_began=time.monotonic();storage_gate(root,cfg);ds=source_entry['source'];m,side=read_source(inputs,source_entry['reference_row_axis_record']);em,e=read_source(evidence,ri[ds]);rows=reader.rows(ds);n=len(rows)
        np.testing.assert_array_equal(rows,side['native_row']);np.testing.assert_array_equal(rows,e['native_row'])
        if n!=source_entry['n_archive_rows'] or text_axis_sha256(reader.native_ids(ds))!=source_entry['native_ids_sha256'] or em['native_ids_sha256']!=source_entry['native_ids_sha256']:raise ValueError('Full native source identity changed')
        if psutil.virtual_memory().available<n*len(genes)*16+cfg['RAM_workspace_reserve_bytes']:raise MemoryError('Full archive single-source RAM gate failed')
        feature=inputs/'outputs'/source_entry['native_feature_columns']
        if sha256_file(feature)!=source_entry['native_feature_columns_sha256']:raise ValueError('Full-source model measured gene axes changed')
        with np.load(feature,allow_pickle=False) as f:columns=f['native_columns']
        lookup=reader.feature_lookup(ds)
        if len(columns)!=len(genes) or len(set(columns))!=len(columns) or any(lookup.get(g,set())!={int(c)} for g,c in zip(genes,columns)):raise ValueError('Missing, ambiguous or aliased native model gene')
        domains={k:np.arange(n,dtype=np.int64) if r['route']=='global_overview' else np.flatnonzero(side['guarded_route_code']==FAMILIES.index(r['route'])) for k,r in plan['selected_models'].items()};detected={k:np.zeros(len(q),np.int16) for k,q in domains.items()};fraction={k:np.zeros(len(q),np.float32) for k,q in domains.items()};totals=np.zeros(n,np.float64);chunks=[];native=reader.store(ds);countkey=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key'];source_nnz=0
        for low in range(0,n,cfg['read_block_rows']):
            high=min(n,low+cfg['read_block_rows']);raw=gather_csr_rows(native[countkey],rows[low:high]);source_nnz+=raw.nnz
            if not np.isfinite(raw.data).all() or np.any(raw.data<0) or not np.equal(raw.data,np.floor(raw.data)).all():raise ValueError('Native full archive count domain changed')
            total=np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel();np.testing.assert_array_equal(total,e['native_total_RNA'][low:high]);totals[low:high]=total;common=raw[:,columns];log,_=normalized_log_counts(common,total);chunks.append(log)
            for identifier,q in domains.items():
                lo,hi=np.searchsorted(q,[low,high]);local=q[lo:hi]-low
                if not len(local):continue
                c=common[local][:,model_values[identifier]['feature_candidate_indices']];detected[identifier][lo:hi]=np.asarray((c>0).sum(axis=1)).ravel();mass=np.asarray(c.sum(axis=1,dtype=np.float64)).ravel();fraction[identifier][lo:hi]=np.divide(mass,total[local],out=np.zeros(len(local)),where=total[local]>0)
            if low==0 or high==n or (low//cfg['read_block_rows'])%64==0:save(root/'checkpoints'/(ds+'_native_read_'+str(high)+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'stage':'full_archive_source_native_model_gene_read','completed_source_rows':high,'n_archive_rows':n})
        native_nnz+=source_nnz;x=sp.vstack(chunks,format='csr');del chunks,raw,common,log,c;gc.collect()
        for identifier,model in plan['selected_models'].items():
            storage_gate(root,cfg);q=domains[identifier];v=model_values[identifier];kind=model['reference_kind'];metadata={'source':ds,'model_id':identifier,'n_source_archive_rows':n,'n_model_domain_rows':len(q),'native_ids_sha256':source_entry['native_ids_sha256'],'array_axis':'explicit_source_native_position_model_domain_rows','raw_context_only':kind=='raw_PCA_context_preserving','source_context':source_entry['source_modality_context'],'source_alias_patient_capture_modality_or_clinical_qualification_not_promoted':True,'new_authoritative_identity':False}
            if not len(q):
                metadata.update(status='no_guarded_rows_in_this_model_route_global_archive_retained',neighbor_precision_all_queried_bases_passed=True);mp=root/'outputs'/(ds+'__'+identifier+'_empty_domain_v1.json');save(mp,metadata);index.append({'source':ds,'model_id':identifier,'n_model_domain_rows':0,'metadata':mp.name,'metadata_sha256':sha256_file(mp),'neighbor_precision_all_queried_bases_passed':True});continue
            model_x=x[q][:,v['feature_candidate_indices']];positive=totals[q]>0
            if kind=='raw_PCA_context_preserving':result=project_raw_context(model_x,v['mean'],v['loadings'],positive)
            else:
                ref,_=restore_reference(models/model['model'],v['genes'].tolist(),v['feature_candidate_indices'])
                if reference_fingerprint(ref)!=model['reference_fingerprint']:raise ValueError('Frozen whole-archive reference fingerprint differs')
                result=project_domain(ref,model_x,v['genes'].tolist(),ds,np.ones(len(q),bool),side['fit_pool_candidate'][q],positive,symphony=symphony)
                if reference_fingerprint(ref)!=model['reference_fingerprint']:raise ValueError('Whole source adaptation mutated reference')
                del ref
            np.testing.assert_array_equal(result['model_gene_signal'],detected[identifier]>0)
            if np.any(fraction[identifier]<0) or np.any(fraction[identifier]>1+1e-6):raise ValueError('Native model RNA fraction invalid')
            arrays={**result,'source_native_position':q,'native_row':rows[q],'fit_pool_candidate':side['fit_pool_candidate'][q],'reference_fit_exposure':side['global_reference_fit'][q] if model['route']=='global_overview' else side['lineage_reference_fit'][q],'native_total_RNA':totals[q],'n_detected_model_genes':detected[identifier],'native_RNA_fraction_in_model_genes':fraction[identifier]};geometry=[];cp=models/'outputs'/model['coordinates']
            if sha256_file(cp)!=model['coordinates_sha256']:raise ValueError('Frozen reference coordinates changed')
            with np.load(cp,allow_pickle=False) as f:refrows=f['reference_union_row'];rawref=f['raw_PCA_fit_embedding'];correctedref=f['reference_embedding']
            source_excluded=np.full(n,-1,np.int64);same=refsource[refrows]==am['sources'].index(ds);source_excluded[refposition[refrows[same]]]=np.flatnonzero(same);excluded=source_excluded[q]
            for basis,key,reference,good in [('raw','raw_PCA',rawref,positive&result['model_gene_signal']),('candidate_corrected','candidate_Harmony',correctedref,result['adaptation_included'])]:
                distance=np.full(len(q),np.nan,np.float32);available=good.copy();qualified=True;threshold=None;exact={'n_exact_probes':0,'distance_exact_match_fraction':None}
                if good.any():
                    ck=(identifier,basis)
                    if ck not in search_cache:
                        search=FrozenReferenceNeighbors(reference,cfg['seed'],cfg['neighbor_parameters']);cal=np.unique(np.linspace(0,len(reference)-1,min(cfg['reference_calibration_rows'],len(reference)),dtype=np.int64));caldistance=search.distances(reference[cal],cal,cfg['neighbor_query_block_rows']);search_cache[ck]=(search,float(np.quantile(caldistance,.99)))
                    search,threshold=search_cache[ck];local=np.flatnonzero(good);distance[local]=search.distances(result[key][local],excluded[local],cfg['neighbor_query_block_rows']);probe=local[np.unique(np.linspace(0,len(local)-1,min(cfg['n_independent_exact_probes'],len(local)),dtype=np.int64))];exact=search.exact_probe_check(result[key][probe],excluded[probe],distance[probe]);qualified=exact['distance_exact_match_fraction']>=cfg['neighbor_parameters']['minimum_exact_probe_fraction']
                    count=min(cfg['n_independent_exact_probes'],len(local));base=np.unique(np.linspace(0,len(local)-1,count,dtype=np.int64));mid=np.unique(np.floor((np.arange(count)+.5)*len(local)/count).astype(np.int64));extra=local[np.setdiff1d(mid,base)]
                    if len(extra):
                        confirmation=search.exact_probe_check(result[key][extra],excluded[extra],distance[extra]);exact['additional_midpoint_exact_probes']=confirmation;qualified=qualified and confirmation['distance_exact_match_fraction']>=cfg['neighbor_parameters']['minimum_exact_probe_fraction']
                outside=good&(distance>threshold) if threshold is not None else np.zeros(len(q),bool);arrays[basis+'_reference_distance']=distance;arrays[basis+'_reference_distance_available']=available;arrays[basis+'_outside_reference_radius_q99']=outside
                geometry.append({'basis':basis,'n_queried_model_domain_rows':int(good.sum()),'reference_other_neighbor_radius_q99':threshold,'n_rows_above_reference_radius_q99':int(outside.sum()),'exact_probe_gate_passed':bool(qualified),'status':'estimated_approximate_reference_geometry' if good.any() else 'not_applicable_or_no_supported_rows',**exact})
            output=root/'outputs'/(ds+'__'+identifier+'_full_projection_v1.npz');save_arrays(output,arrays);metadata.update(status='full_source_model_domain_projected_not_identity_release',arrays=output.name,arrays_sha256=sha256_file(output),n_adaptation_eligible=int(result['adaptation_eligible'].sum()),n_adaptation_included=int(result['adaptation_included'].sum()),n_archive_excluded_domain_rows_raw_only=int((~side['fit_pool_candidate'][q]).sum()),geometry_profiles=geometry,neighbor_precision_all_queried_bases_passed=all(r['exact_probe_gate_passed'] for r in geometry),radius_flags_are_not_biological_or_identity_confidence=True,query_adaptation_scope='complete_eligible_source_model_domain_single_call_or_context_raw_only',raw_fit_transform_not_assumed_equal_to_later_transform=True)
            mp=output.with_suffix('.json');save(mp,metadata);index.append({'source':ds,'model_id':identifier,'n_model_domain_rows':len(q),'metadata':mp.name,'metadata_sha256':sha256_file(mp),'arrays':output.name,'arrays_sha256':sha256_file(output),'neighbor_precision_all_queried_bases_passed':metadata['neighbor_precision_all_queried_bases_passed']});save(root/'checkpoints'/(ds+'__'+identifier+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'model_id':identifier,'stage':'full_source_model_domain_projection_and_neighbor_flags_complete','n_model_domain_rows':len(q),'neighbor_precision_all_queried_bases_passed':metadata['neighbor_precision_all_queried_bases_passed']});del result,arrays,model_x,rawref,correctedref;gc.collect()
        source_records.append({'source':ds,'n_archive_rows':n,'native_nnz_read':int(source_nnz)});save(root/'checkpoints'/(ds+'_complete.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'stage':'full_source_all_selected_models_projected','n_archive_rows':n,'wall_seconds':time.monotonic()-source_began,'n_complete_sources':len(source_records)});reader._fields.clear();del x,side,e;gc.collect()
    output=root/'outputs/full_archive_projection_index_v1.json';save(output,index);save(root/'outputs/source_native_count_read_coverage_v1.json',source_records)
    save(root/'provenance/terminal_v1.json',{'status':'FULL_ARCHIVE_SELECTED_PROJECTIONS_AND_GEOMETRIC_FLAGS_COMPUTED_INDEPENDENT_AUDIT_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(source_records),'n_archive_rows':sum(r['n_archive_rows'] for r in source_records),'n_source_model_records':len(index),'all_queried_neighbor_precision_gates_passed':all(r['neighbor_precision_all_queried_bases_passed'] for r in index),'native_nnz_read':int(native_nnz),'index_sha256':sha256_file(output),'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'native_count_matrices_saved':False,'new_final_identity':False,'standalone_Zarr':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
