"""Freeze one common broad compartment domain before comparing fine schemes.

The same native marker/reference scope is used for every comparison arm.
Compartments are explicit denominators, not new claims that T and NK or
fibroblast and mural cells have identical biological identity.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from TME_shared_compartment_scope_v1 import shared_compartment,counterexamples

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def main(R):
    cfg=load(R/'config/spec_v1.json');save(R/'provenance/policy_counterexamples_v1.json',counterexamples());sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg);S=Path(cfg['semantics']);C=Path(cfg['coexpression']);CA=Path(cfg['coexpression_auditor']);H=Path(cfg['model']);B=Path(cfg['naming']);policy=load(R/'config/frozen_compartment_policy_v1.json');assert sha(R/'config/frozen_compartment_policy_v1.json')==cfg['policy_SHA256'];deadline=time.monotonic()+cfg['dependency_wait_seconds']
    for root in [S,C,CA,B]:
        p=root/'provenance/supervisor_terminal_v1.json'
        while not p.exists():
            assert time.monotonic()<deadline;time.sleep(5)
        assert load(p)['actual_child_exit_code']==0
    sp=S/'outputs/accepted_reference_semantics_overlay_v1.json';assert sha(sp)==load(S/'provenance/terminal_v1.json')['acceptance_SHA256'];semantic=load(sp);cat=load(B/'outputs/marker_first_cluster_catalog_v1.json');groups=np.load(B/'outputs/private_full_fit_cluster_code_v1.npy');N=cat['n_fit'];assert N==semantic['n_fit']==len(groups)
    cm=C/'outputs/joint_cluster_coexpression_manifest_v1.json';manifest=load(cm);audit=load(CA/'outputs/accepted_joint_cluster_coexpression_v1.json');assert audit['producer_terminal_SHA256']==sha(C/'provenance/terminal_v1.json');assert audit['all_cell_coexpression_bits_independently_reproduced'];assert manifest['input_axes_SHA256']==sha(H/'outputs/private_cell_axes_v1.npz') and manifest['n_fit']==N
    bitpath=C/'outputs/private_all_cell_coherent_lineage_bitmask_v1.npy';assert sha(bitpath)==next(r['SHA256'] for r in manifest['outputs'] if r['path']==bitpath.name);bits=np.load(bitpath,mmap_mode='r');assert bits.shape==(N,) and bits.dtype==np.uint32;gate={name:i for i,name in enumerate(manifest['gate_order'])};names=list(policy['compartments']);by={name:i for i,name in enumerate(names)}
    family_of={parent:family for family,item in policy['compartments'].items() for parent in item['identity_parents']};assert sum(len(v['identity_parents']) for v in policy['compartments'].values())==len(family_of)
    masks=np.asarray([sum(1<<gate[g] for g in policy['compartments'][name]['native_gates']) for name in names],np.uint32)
    cluster_map=np.full(len(cat['clusters']),-1,np.int16)
    for row in cat['clusters']:
        k=row['cluster'];assert k<len(cluster_map)
        if k not in policy['ambiguous_marker_cluster_ids']:cluster_map[k]=by.get(family_of.get(row['major_marker_first']),-1)
    marker=cluster_map[groups];semparents=semantic['parents'];reference_map=np.asarray([by.get(family_of.get(p),-1) for p in semparents],np.int16);references={};raw_bindings=[]
    for key in ['immune_low','immune_high']:
        row=next(r for r in semantic['models'] if r['model']==key);binding=row['parent_proposal'];assert sha(binding['path'])==binding['SHA256'];codes=np.load(binding['path'],mmap_mode='r');assert codes.shape==(N,);v=np.full(N,-1,np.int16);good=codes>=0;v[good]=reference_map[codes[good]];references[key]=v;raw_bindings.append(binding)
    final,why=shared_compartment(marker,references['immune_low'],references['immune_high'],bits,masks);fp=R/'outputs/private_shared_compartment_code_v1.npy';rp=R/'outputs/private_shared_compartment_reason_v1.npy';np.save(fp,final);np.save(rp,why);read=np.load(fp,mmap_mode='r');np.testing.assert_array_equal(read,final)
    # One scalar case for every observed rule outcome, with explicit marker
    # ambiguity and native/ref support visible in the private technical proof.
    probes=[]
    for value in np.unique(why):
        i=int(np.flatnonzero(why==value)[0]);lo=int(references['immune_low'][i]);hi=int(references['immune_high'][i]);ma=int(marker[i]);has=lo>=0 and hi==lo and bool(int(bits[i])&int(masks[lo]));expected=lo if has else ma;assert int(read[i])==expected;probes.append(dict(rule_code=int(value),scalar_case_replayed=True))
    counts=np.bincount(final[final>=0],minlength=len(names));reasons=np.bincount(why,minlength=5);assert int(counts.sum())+int((final<0).sum())==N;assert (final>=-1).all() and (final<len(names)).all()
    accepted=dict(schema='accepted_TME_shared_compartment_scope_v1',created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),accepted_purpose='qualified_workunit_annotation_comparison_exploration',n_fit=N,n_compartment_assessed=int((final>=0).sum()),n_compartment_unresolved=int((final<0).sum()),compartments=names,compartment_counts=counts.tolist(),reason_counts=reasons.tolist(),reason_definitions=policy['reason_definitions'],codes=dict(path=str(fp),SHA256=sha(fp)),reasons=dict(path=str(rp),SHA256=sha(rp)),policy_SHA256=sha(R/'config/frozen_compartment_policy_v1.json'),reference_semantics_acceptance=dict(path=str(sp),SHA256=sha(sp)),reference_bindings=raw_bindings,coexpression_acceptance_SHA256=sha(CA/'outputs/accepted_joint_cluster_coexpression_v1.json'),marker_catalog_SHA256=sha(B/'outputs/marker_first_cluster_catalog_v1.json'),input_axis_SHA256=sha(H/'outputs/private_cell_axes_v1.npz'),scalar_rule_probes=probes,compartments_are_declared_denominators_not_identical_cell_types=True,marker_fallback_is_working_cluster_identity_not_per_cell_truth=True,cell_state_doublet_or_malignancy_not_reestimated=True,all_archive_and_raw_model_labels_retained=True,independent_patients_not_certified=True,final_identity_or_clinical_module_eligibility_not_granted=True)
    save(R/'outputs/accepted_shared_compartment_scope_v1.json',accepted);storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='SHARED_MARKER_REFERENCE_NATIVE_COMPARTMENT_SCOPE_FROZEN_FOR_EXPLORATORY_COMPARISON',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(R/'outputs/accepted_shared_compartment_scope_v1.json')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
