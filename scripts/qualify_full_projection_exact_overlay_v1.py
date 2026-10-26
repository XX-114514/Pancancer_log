"""Accept a full computational projection after its audited exact-distance patch.

No expression, coordinate or original failure file is rewritten. The manifest
binds original per-record projections and explicit, complete failed-basis patches.
It is not a final identity, patient eligibility or biological transfer release.
"""
import argparse,datetime,json
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import sha256_file


def completed(stage):
    s=json.loads((stage/'provenance/supervisor_terminal_v1.json').read_text())
    if s['actual_child_exit_code'] or s['wall_timeout']:raise ValueError('Actual successful completed stage required')
    return json.loads((stage/'provenance/terminal_v1.json').read_text())


def main(root):
    cfg=json.loads((root/'config/exact_overlay_qualification_spec_v1.json').read_text())
    producer,auditor,repair=map(Path,[cfg['producer'],cfg['auditor'],cfg['repair']])
    t,a,r=map(completed,[producer,auditor,repair])
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Qualification bound input changed')
    producer_sha=sha256_file(producer/'provenance/terminal_v1.json')
    if a['producer_terminal_sha256']!=producer_sha or r['producer_terminal_sha256']!=producer_sha or not r['all_repaired_precision_gates_passed']:raise ValueError('Independent audit/repair belongs to another projection')
    ip=producer/'outputs/full_archive_projection_index_v1.json';ap=auditor/'outputs/full_archive_projection_audit_v1.json';rp=repair/'outputs/bounded_neighbor_exact_repair_index_v1.json'
    for p,digest in [(ip,t['index_sha256']),(ap,a['output_sha256']),(rp,r['index_sha256'])]:
        if sha256_file(p)!=digest:raise ValueError('Full audit or repair output changed')
    original=json.loads(ip.read_text());audit=json.loads(ap.read_text());patches=json.loads(rp.read_text());records={(v['source'],v['model_id']):v for v in original};ar={(v['source'],v['model_id']):v for v in audit['records']}
    if len(records)!=len(original) or set(ar)!=set(records) or len(ar)!=len(audit['records']):raise ValueError('Complete projection audit scope differs')
    mapping={}
    for v in patches:
        key=(v['source'],v['model_id'],v['basis'])
        if key in mapping:raise ValueError('Duplicate exact patch basis')
        mapping[key]=v
    spec=json.loads((producer/'config/full_archive_projection_spec_v1.json').read_text());plan=json.loads((Path(spec['projection_inputs'])/'outputs/full_archive_projection_inputs_v1.json').read_text());models=Path(plan['model_stage']);axis=json.loads((models/'outputs/reference_union_axis_metadata_v1.json').read_text())
    with np.load(models/'outputs/reference_union_axes_and_native_readouts_v1.npz',allow_pickle=False) as f:refsource,refposition=f['source_code'],f['native_source_position']
    effective=[];consumed=set()
    for key,record in records.items():
        if ar[key]['neighbor_precision_passed']!=record['neighbor_precision_all_queried_bases_passed']:raise ValueError('Base independent precision audit differs')
        selected=[]
        if not record['neighbor_precision_all_queried_bases_passed']:
            mp=producer/'outputs'/record['metadata'];fp=producer/'outputs'/record['arrays']
            if sha256_file(mp)!=record['metadata_sha256'] or sha256_file(fp)!=record['arrays_sha256']:raise ValueError('Original failed projection changed')
            meta=json.loads(mp.read_text());model=plan['selected_models'][key[1]];cp=models/'outputs'/model['coordinates']
            if sha256_file(cp)!=model['coordinates_sha256']:raise ValueError('Frozen reference coordinate bytes changed')
            with np.load(cp,allow_pickle=False) as f:refrows=f['reference_union_row'];references={'raw':f['raw_PCA_fit_embedding'],'candidate_corrected':f['reference_embedding']}
            with np.load(fp,allow_pickle=False) as f:base={k:f[k] for k in f.files}
            same=refsource[refrows]==axis['sources'].index(key[0]);excluded={int(p):int(i) for p,i in zip(refposition[refrows[same]],np.flatnonzero(same))}
            for profile in meta['geometry_profiles']:
                if profile['exact_probe_gate_passed']:continue
                pkey=(*key,profile['basis'])
                if pkey not in mapping:raise ValueError('Missing complete failed-basis exact repair')
                patch=mapping[pkey];pp=repair/'outputs'/patch['arrays']
                if patch['original_record']!=record or patch['original_precision_profile']!=profile or sha256_file(pp)!=patch['arrays_sha256']:raise ValueError('Exact patch provenance differs')
                with np.load(pp,allow_pickle=False) as f:values={k:f[k] for k in f.files}
                local=np.flatnonzero(base[profile['basis']+'_reference_distance_available'])
                np.testing.assert_array_equal(values['source_native_position'],base['source_native_position'][local]);np.testing.assert_array_equal(values['native_row'],base['native_row'][local]);np.testing.assert_array_equal(values['distance_available'],np.ones(len(local),bool))
                query=base['raw_PCA' if profile['basis']=='raw' else 'candidate_Harmony'][local];distance=cdist(query.astype(float),references[profile['basis']].astype(float))
                for j,position in enumerate(values['source_native_position']):
                    if int(position) in excluded:distance[j,excluded[int(position)]]=np.inf
                exact=distance.min(axis=1);np.testing.assert_allclose(values['exact_nearest_other_reference_distance'],exact,rtol=1e-12,atol=1e-12);np.testing.assert_array_equal(values['outside_original_reference_radius_q99'],exact>profile['reference_other_neighbor_radius_q99'])
                if patch['original_reference_radius_q99_preserved']!=profile['reference_other_neighbor_radius_q99'] or patch['effective_exact_match_fraction']!=1.0 or not patch['precision_gate_passed']:raise ValueError('Patch changed radius or failed precision')
                selected.append({'basis':profile['basis'],'repair_arrays':str(pp),'repair_arrays_sha256':patch['arrays_sha256'],'complete_basis_domain_replaced':True});consumed.add(pkey)
        effective.append({'source':key[0],'model_id':key[1],'base_projection_record':record,'exact_neighbor_overlays':selected,'effective_reference_neighbor_precision_qualified':True})
    if consumed!=set(mapping):raise ValueError('Extra or unconsumed exact-distance patch')
    out=root/'outputs/qualified_full_projection_overlay_index_v1.json';save(out,{'schema':'qualified_full_projection_overlay_index_v1','base_projection':str(producer),'base_full_auditor':str(auditor),'exact_repair':str(repair),'records':effective,'coordinates_and_native_measurement_not_changed':True,'original_failures_preserved':True,'new_final_identity':False})
    save(root/'provenance/terminal_v1.json',{'status':'FULL_COMPUTATIONAL_PROJECTION_WITH_EXACT_NEIGHBOR_OVERLAY_QUALIFIED_NOT_IDENTITY_RELEASE','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_source_model_records':len(effective),'n_exactly_repaired_bases':len(consumed),'output_sha256':sha256_file(out),'producer_terminal_sha256':producer_sha,'base_auditor_terminal_sha256':sha256_file(auditor/'provenance/terminal_v1.json'),'repair_terminal_sha256':sha256_file(repair/'provenance/terminal_v1.json'),'new_final_identity':False,'source_patient_and_biological_transfer_qualification_not_promoted':True})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
