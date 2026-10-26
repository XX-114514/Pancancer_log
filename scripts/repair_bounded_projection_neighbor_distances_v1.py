"""Exact coordinate-only repair for bounded failed reference-neighbor domains.

Preserves original arrays and failed precision records. Full-query projection
arithmetic and native marker audits remain separate acceptance requirements.
"""
import argparse,datetime,gc,json,resource,time
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from annotation_representation_pilot_v1 import save
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file


def exact_two_implementations(query,reference,excluded):
    """SciPy distance followed by independent blocked sum-of-squares replay."""
    first=cdist(query.astype(float),reference.astype(float))
    for j,index in enumerate(excluded):
        if index>=0:first[j,index]=np.inf
    exact=first.min(axis=1);del first
    replay=np.full(len(query),np.inf)
    for low in range(0,len(reference),2048):
        difference=query.astype(float)[:,None,:]-reference[None,low:low+2048,:].astype(float)
        squared=np.einsum('ijk,ijk->ij',difference,difference)
        for j,index in enumerate(excluded):
            if low<=index<low+squared.shape[1]:squared[j,int(index)-low]=np.inf
        replay=np.minimum(replay,np.sqrt(squared.min(axis=1)))
    np.testing.assert_allclose(exact,replay,rtol=1e-12,atol=1e-12)
    if not np.isfinite(exact).all():raise ValueError('No finite other-reference exact neighbor')
    return exact,float(np.max(np.abs(exact-replay)))


def main(root):
    began=time.monotonic();cfg=json.loads((root/'config/bounded_neighbor_repair_spec_v1.json').read_text());producer=Path(cfg['producer'])
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Neighbor repair input binding changed')
    s=json.loads((producer/'provenance/supervisor_terminal_v1.json').read_text());t=json.loads((producer/'provenance/terminal_v1.json').read_text());ip=producer/'outputs/full_archive_projection_index_v1.json'
    if s['actual_child_exit_code'] or s['wall_timeout'] or sha256_file(ip)!=t['index_sha256']:raise ValueError('Completed frozen full projection required')
    index=json.loads(ip.read_text());spec=json.loads((producer/'config/full_archive_projection_spec_v1.json').read_text());plan=json.loads((Path(spec['projection_inputs'])/'outputs/full_archive_projection_inputs_v1.json').read_text());models=Path(plan['model_stage']);axis=json.loads((models/'outputs/reference_union_axis_metadata_v1.json').read_text())
    with np.load(models/'outputs/reference_union_axes_and_native_readouts_v1.npz',allow_pickle=False) as f:refsource,refposition=f['source_code'],f['native_source_position']
    failed=[r for r in index if not r['neighbor_precision_all_queried_bases_passed']]
    expected={(r['source'],r['model_id']) for r in cfg['failed_records']}
    if len(failed)!=len(expected) or {(r['source'],r['model_id']) for r in failed}!=expected:raise ValueError('Frozen failed-record scope changed')
    reports=[]
    for record in failed:
        storage_gate(root,cfg);ds,identifier=record['source'],record['model_id'];mp=producer/'outputs'/record['metadata'];ap=producer/'outputs'/record['arrays']
        if sha256_file(mp)!=record['metadata_sha256'] or sha256_file(ap)!=record['arrays_sha256']:raise ValueError('Original failed-record bytes changed')
        meta=json.loads(mp.read_text());model=plan['selected_models'][identifier];cp=models/'outputs'/model['coordinates']
        if sha256_file(cp)!=model['coordinates_sha256']:raise ValueError('Frozen reference geometry changed')
        with np.load(cp,allow_pickle=False) as f:refrows,rawref,correctedref=f['reference_union_row'],f['raw_PCA_fit_embedding'],f['reference_embedding']
        with np.load(ap,allow_pickle=False) as f:arrays={k:f[k] for k in f.files}
        position=arrays['source_native_position'];same=refsource[refrows]==axis['sources'].index(ds);mapping={int(p):int(i) for p,i in zip(refposition[refrows[same]],np.flatnonzero(same))};exclude=np.asarray([mapping.get(int(p),-1) for p in position],np.int64)
        for profile in meta['geometry_profiles']:
            if profile['exact_probe_gate_passed']:continue
            basis=profile['basis'];good=arrays[basis+'_reference_distance_available'];selected=np.flatnonzero(good)
            if not 0<len(selected)<=cfg['maximum_exact_repair_rows']:raise ValueError('Failed coordinate domain exceeds frozen bounded exact repair')
            reference=rawref if basis=='raw' else correctedref;query=arrays['raw_PCA' if basis=='raw' else 'candidate_Harmony'][selected]
            exact,error=exact_two_implementations(query,reference,exclude[selected]);radius=profile['reference_other_neighbor_radius_q99'];outside=exact>radius
            output=root/'outputs'/(ds+'__'+identifier+'__'+basis+'_exact_neighbor_patch_v1.npz');values={'source_native_position':position[selected],'native_row':arrays['native_row'][selected],'exact_nearest_other_reference_distance':exact,'distance_available':np.ones(len(selected),bool),'outside_original_reference_radius_q99':outside}
            with output.open('xb') as handle:np.savez_compressed(handle,**values)
            with np.load(output,allow_pickle=False) as f:
                for key,value in values.items():np.testing.assert_array_equal(f[key],value)
            reports.append({'source':ds,'model_id':identifier,'basis':basis,'n_complete_failed_basis_rows':len(selected),'arrays':output.name,'arrays_sha256':sha256_file(output),'original_record':record,'original_precision_profile':profile,'reference_coordinates_sha256':model['coordinates_sha256'],'original_reference_radius_q99_preserved':radius,'two_exact_implementations_maximum_absolute_difference':error,'effective_exact_match_fraction':1.0,'precision_gate_passed':True,'threshold_not_lowered':True,'no_native_counts_read_or_model_refit':True})
        del arrays,rawref,correctedref;gc.collect()
    out=root/'outputs/bounded_neighbor_exact_repair_index_v1.json';save(out,reports)
    save(root/'provenance/terminal_v1.json',{'status':'BOUNDED_FAILED_NEIGHBOR_DOMAINS_EXACTLY_RECOMPUTED_FULL_BASE_AUDIT_ACCEPTANCE_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_failed_source_model_records':len(failed),'n_repaired_bases':len(reports),'n_exact_rows':sum(r['n_complete_failed_basis_rows'] for r in reports),'index_sha256':sha256_file(out),'producer_terminal_sha256':sha256_file(producer/'provenance/terminal_v1.json'),'all_repaired_precision_gates_passed':all(r['precision_gate_passed'] for r in reports),'native_counts_read':False,'new_final_identity':False,'base_full_projection_audit_still_required':True,'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
