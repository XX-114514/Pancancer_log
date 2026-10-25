"""Full-cell reference distance flags on accepted whole-source projections.

The reference-radius flags are geometric diagnostics, never identity accuracy,
patient validation, or a biological rejection rule. No counts are read.
"""
import argparse,datetime,gc,json,resource,time
from pathlib import Path
import numpy as np
from annotation_representation_pilot_v1 import save
from extract_full_source_local_marker_evidence_v1 import storage_gate
from frozen_reference_neighbors_v1 import FrozenReferenceNeighbors,actual_runtime_smoke
from full_context_atlas_v1 import sha256_file


def main(root):
    start=time.monotonic();cfg=json.loads((root/'config/neighbor_qualification_spec_v1.json').read_text());producer=Path(cfg['projection']);auditor=Path(cfg['projection_audit'])
    for stage in [producer,auditor]:
        sup=json.loads((stage/'provenance/supervisor_terminal_v1.json').read_text())
        if sup['actual_child_exit_code'] or sup['wall_timeout']:raise ValueError('Successful projection and independent full-query audit required')
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Neighbor qualification input changed')
    at=json.loads((auditor/'provenance/terminal_v1.json').read_text())
    if at['producer_terminal_sha256']!=sha256_file(producer/'provenance/terminal_v1.json'):raise ValueError('Projection auditor belongs to another producer')
    smoke=actual_runtime_smoke(cfg['seed'],cfg['neighbor_parameters']);save(root/'provenance/actual_neighbor_runtime_smoke_v1.json',smoke)
    pspec=json.loads((producer/'config/complete_projection_spec_v1.json').read_text());plan=json.loads((Path(pspec['qualification_inputs'])/'outputs/complete_source_projection_qualification_inputs_v1.json').read_text());models=Path(plan['model_stage']);am=json.loads((models/'outputs/reference_union_axis_metadata_v1.json').read_text())
    with np.load(models/'outputs/reference_union_axes_and_native_readouts_v1.npz',allow_pickle=False) as f:refsource=f['source_code'];refposition=f['native_source_position']
    t=json.loads((producer/'provenance/terminal_v1.json').read_text());ip=producer/'outputs/complete_source_projection_qualification_index_v1.json'
    if sha256_file(ip)!=t['index_sha256']:raise ValueError('Projection index changed')
    records=json.loads(ip.read_text());cache={};results=[]
    for record in records:
        storage_gate(root,cfg);identifier=record['model_id'];ds=record['source'];model=plan['selected_models'][identifier];cp=models/'outputs'/model['coordinates']
        if sha256_file(cp)!=model['coordinates_sha256']:raise ValueError('Reference coordinate bytes changed')
        meta=json.loads((producer/'outputs'/record['metadata']).read_text());p=producer/'outputs'/record['arrays']
        if sha256_file(p)!=record['arrays_sha256']:raise ValueError('Complete-source coordinates changed')
        with np.load(p,allow_pickle=False) as f:a={k:f[k] for k in ['native_row','raw_PCA','candidate_Harmony','domain_candidate','native_total_RNA','model_gene_signal','adaptation_included','reference_fit_exposure']}
        with np.load(cp,allow_pickle=False) as f:refrows=f['reference_union_row'];rawref=f['raw_PCA_fit_embedding'];corrref=f['reference_embedding']
        n=len(a['native_row']);exclude=np.full(n,-1,np.int64);same=refsource[refrows]==am['sources'].index(ds);exclude[refposition[refrows[same]]]=np.flatnonzero(same);arrays={'native_row':a['native_row'],'reference_fit_exposure':a['reference_fit_exposure']};report={'source':ds,'model_id':identifier,'n_archive_rows':n,'basis_profiles':[],'radius_flags_are_not_biological_or_identity_confidence':True,'source_alias_patient_capture_or_modality_qualification_not_promoted':True}
        for basis,key,reference,good in [('raw','raw_PCA',rawref,a['domain_candidate']&(a['native_total_RNA']>0)&a['model_gene_signal']),('candidate_corrected','candidate_Harmony',corrref,a['adaptation_included'])]:
            cache_key=(identifier,basis)
            if cache_key not in cache:
                search=FrozenReferenceNeighbors(reference,cfg['seed'],cfg['neighbor_parameters']);cal=np.unique(np.linspace(0,len(reference)-1,min(cfg['reference_calibration_rows'],len(reference)),dtype=np.int64));radius=search.distances(reference[cal],cal,cfg['query_block_rows']);threshold=float(np.quantile(radius,.99));cache[cache_key]=(search,threshold)
            search,threshold=cache[cache_key];q=np.flatnonzero(good);distance=np.full(n,np.nan,np.float32)
            if len(q):distance[q]=search.distances(a[key][q],exclude[q],cfg['query_block_rows'])
            probes=q[np.unique(np.linspace(0,len(q)-1,min(cfg['n_independent_exact_probes'],len(q)),dtype=np.int64))] if len(q) else np.array([],np.int64)
            exact=search.exact_probe_check(a[key][probes],exclude[probes],distance[probes]);qualified=exact['distance_exact_match_fraction'] is not None and exact['distance_exact_match_fraction']>=cfg['neighbor_parameters']['minimum_exact_probe_fraction']
            confirmation=None
            if cfg.get('additional_midpoint_probes',False) and len(q):
                count=min(cfg['n_independent_exact_probes'],len(q));base=np.unique(np.linspace(0,len(q)-1,count,dtype=np.int64));mid=np.unique(np.floor((np.arange(count)+.5)*len(q)/count).astype(np.int64));extra=q[np.setdiff1d(mid,base)]
                if len(extra):
                    confirmation=search.exact_probe_check(a[key][extra],exclude[extra],distance[extra]);qualified=qualified and confirmation['distance_exact_match_fraction']>=cfg['neighbor_parameters']['minimum_exact_probe_fraction']
            # Preserve failed exact-search qualifications; a later all-archive
            # launch must require these gates and cannot silently accept them.
            arrays[basis+'_reference_distance']=distance;arrays[basis+'_distance_available']=good;arrays[basis+'_outside_reference_radius_q99']=good&(distance>threshold)
            report['basis_profiles'].append({'basis':basis,'n_full_source_numeric_rows_queried':int(good.sum()),'reference_other_neighbor_radius_q99':threshold,'n_rows_above_reference_radius_q99':int((good&(distance>threshold)).sum()),'exact_probe_gate_passed':qualified,'additional_midpoint_exact_probes':confirmation,**exact})
            save(root/'checkpoints'/(ds+'__'+identifier+'__'+basis+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'model_id':identifier,'basis':basis,'stage':'all_valid_source_rows_reference_neighbors_and_exact_probes_complete','exact_probe_gate_passed':qualified})
        output=root/'outputs'/(ds+'__'+identifier+'_reference_neighbor_flags_v1.npz')
        with output.open('xb') as f:np.savez_compressed(f,**arrays)
        with np.load(output,allow_pickle=False) as f:
            for k,v in arrays.items():np.testing.assert_array_equal(v,f[k])
        report['arrays']=output.name;report['arrays_sha256']=sha256_file(output);mf=output.with_suffix('.json');save(mf,report);results.append({'source':ds,'model_id':identifier,'metadata':mf.name,'metadata_sha256':sha256_file(mf),'all_basis_exact_probe_gates_passed':all(r['exact_probe_gate_passed'] for r in report['basis_profiles'])});del a,arrays,rawref,corrref;gc.collect()
    output=root/'outputs/complete_source_neighbor_qualification_index_v1.json';save(output,results);passed=all(r['all_basis_exact_probe_gates_passed'] for r in results)
    save(root/'provenance/terminal_v1.json',{'status':'COMPLETE_SOURCE_NEIGHBOR_FLAGS_AND_EXACT_PROBE_GATES_COMPUTED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_source_model_records':len(results),'all_exact_probe_gates_passed':passed,'index_sha256':sha256_file(output),'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'native_counts_read':False,'new_final_identity':False,'full_archive_mapping_complete':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
