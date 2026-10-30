"""Exact reference-label semantics, separate from final biological identity.

Raw inference codes and learned numerical objects are retained unchanged.
The overlay corrects lexical routing and records phenotype/state caveats; it
does not grant tissue, patient, biological or CoVarNet input qualification.
"""
import argparse,collections,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def main(R):
    cfg=load(R/'config/spec_v1.json');C=Path(cfg['automatic_controller']);deadline=time.monotonic()+cfg['dependency_wait_seconds']
    sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg);assert load(C/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;pointer=load(C/'outputs/automatic_execution_pointer_v1.json');A=Path(pointer['automatic']);D=Path(pointer['auditor'])
    for root in [A,D]:
        p=root/'provenance/supervisor_terminal_v1.json'
        while not p.exists():
            assert time.monotonic()<deadline,'Official inference/audit dependency did not complete';time.sleep(5)
        assert load(p)['actual_child_exit_code']==0,'Inference failure is not a semantic correction'
    ap=D/'outputs/accepted_TME_CellTypist_cluster_crosscheck_v1.json';assert sha(ap)==load(D/'provenance/terminal_v1.json')['acceptance_SHA256'];accepted=load(ap);autop=Path(accepted['automatic_crosscheck']['path']);assert sha(autop)==accepted['automatic_crosscheck']['SHA256'];auto=load(autop);N=auto['n_fit'];completed_models={};completion_binding=None
    if cfg.get('core_coverage_completion'):
        completion=Path(cfg['core_coverage_completion']);p=completion/'provenance/supervisor_terminal_v1.json'
        while not p.exists():
            assert time.monotonic()<deadline;time.sleep(5)
        assert load(p)['actual_child_exit_code']==0
        cp=completion/'outputs/accepted_core_coverage_completion_v1.json';assert sha(cp)==load(completion/'provenance/terminal_v1.json')['acceptance_SHA256'];record=load(cp);assert record['schema']=='accepted_TME_core_reference_coverage_completion_v1' and record['n_fit']==N;assert record['baseline_execution_acceptance']['SHA256']==sha(ap);assert record['baseline_automatic']==str(A) and record['all_existing_predictions_and_model_numbers_unchanged'];completed_models={r['model']:r for r in record['models']};completion_binding=dict(path=str(cp),SHA256=sha(cp))
    op=R/'config/frozen_reference_label_ontology_v2.json';assert sha(op)==cfg['ontology_SHA256'];ontology=load(op);assert ontology['schema']=='frozen_exact_TME_reference_label_ontology_v2';assert ontology['raw_prediction_codes_and_model_values_unchanged'];models=ontology['models'];assert set(models)==set(auto['models'])
    parents=sorted({r['coarse_reference_parent'] for m in models.values() for r in m['rows']} - {'Reference_unresolved'});parent_index={name:i for i,name in enumerate(parents)};identities=sorted({(r['coarse_reference_parent'],r['canonical_identity']) for m in models.values() for r in m['rows'] if r['coarse_reference_parent']!='Reference_unresolved'});identity_index={key:i for i,key in enumerate(identities)};summary=[];probes=[]
    for key,m in models.items():
        actual=auto['models'][key];assert actual['labels']==m['labels'] and actual['model_SHA256']==m['model_SHA256'];assert m['label_axis_SHA256']==hashlib.sha256('\0'.join(m['labels']).encode()).hexdigest();assert [r['label'] for r in m['rows']]==m['labels'];assert len(set(m['labels']))==len(m['labels'])
        rawpath=A/'outputs'/('private_'+key+'_best_code_v1.npy')
        if key in completed_models:
            c=completed_models[key];assert c['labels']==m['labels'] and c['model_SHA256']==m['model_SHA256'];rawpath=Path(c['fields']['best_code']['path']);assert sha(rawpath)==c['fields']['best_code']['SHA256']
        before=sha(rawpath);raw=np.load(rawpath,mmap_mode='r');assert raw.shape==(N,) and ((raw>=-1)&(raw<len(m['rows']))).all();good=raw>=0;pmap=np.asarray([parent_index.get(r['coarse_reference_parent'],-1) for r in m['rows']],np.int16);fmap=np.asarray([identity_index.get((r['coarse_reference_parent'],r['canonical_identity']),-1) for r in m['rows']],np.int16);p=np.full(N,-1,np.int16);f=np.full(N,-1,np.int16);p[good]=pmap[raw[good]];f[good]=fmap[raw[good]];assert np.array_equal(p>=0,f>=0)
        pp=R/'outputs'/('private_'+key+'_reference_parent_proposal_v1.npy');fp=R/'outputs'/('private_'+key+'_reference_identity_proposal_v1.npy');np.save(pp,p);np.save(fp,f);read_parent=np.load(pp,mmap_mode='r');read_fine=np.load(fp,mmap_mode='r')
        # Independent scalar string-key lookup after disk readback, including
        # every actually observed reference class and missing predictions.
        by_label={r['label']:r for r in m['rows']};indices=[]
        for code in np.unique(raw):
            at=int(np.flatnonzero(raw==code)[0]);indices.append(at)
            if code<0:assert read_parent[at]==read_fine[at]==-1;continue
            row=by_label[actual['labels'][int(code)]];assert read_parent[at]==parent_index.get(row['coarse_reference_parent'],-1);assert read_fine[at]==identity_index.get((row['coarse_reference_parent'],row['canonical_identity']),-1)
        for label in ['Mast cells','Mast']:
            if label in m['labels']:
                assert by_label[label]['coarse_reference_parent']=='Mast';mask=raw==m['labels'].index(label);assert (read_parent[mask]==parent_index['Mast']).all()
        assert sha(rawpath)==before;counts=np.bincount(p[p>=0],minlength=len(parents));assert int(counts.sum())+int((p<0).sum())==N
        summary.append(dict(model=key,n_raw_predicted=int(good.sum()),n_reference_parent_resolved=int((p>=0).sum()),n_no_prediction=int((~good).sum()),n_unresolved_reference_parent_among_predicted=int((good&(p<0)).sum()),reference_parent_proposal_counts=dict(zip(parents,counts.tolist())),raw_prediction_codes=dict(path=str(rawpath),SHA256=before),parent_proposal=dict(path=str(pp),SHA256=sha(pp)),identity_proposal=dict(path=str(fp),SHA256=sha(fp)),reference_model_SHA256=m['model_SHA256'],proposals_not_final_cell_identities=True));probes.append(dict(model=key,n_scalar_raw_label_cases=len(indices),every_observed_class_and_unpredicted_case_replayed=True))
    # Ontology-only counterexamples do not depend on whether these classes
    # occur in this atlas. Progenitors and state words cannot acquire identity.
    low={r['label']:r for r in models['immune_low']['rows']};assert low['Mast cells']['coarse_reference_parent']=='Mast';assert low['Megakaryocyte-erythroid-mast cell progenitor']['coarse_reference_parent']=='Reference_unresolved';assert low['Cycling DCs']['coarse_reference_parent']=='Myeloid';assert low['Tem/Trm cytotoxic T cells']['canonical_identity_is_fallback']
    lung={r['label']:r for r in models['lung']['rows']};assert not lung['Hematopoietic stem cells']['reference_state_hints'];assert not lung['EC venous systemic']['reference_state_hints'];assert lung['Mast cells']['coarse_reference_parent']=='Mast'
    manifest=dict(schema='accepted_TME_reference_semantics_overlay_v1',created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_fit=N,parents=parents,identities=[dict(parent=p,canonical_identity=f) for p,f in identities],models=summary,raw_execution_acceptance=dict(path=str(ap),SHA256=sha(ap)),core_coverage_completion_acceptance=completion_binding,exact_ontology=dict(path=str(op),SHA256=sha(op)),full_vocabulary_rows=sum(len(m['rows']) for m in models.values()),independent_scalar_label_replays=probes,raw_prediction_codes_unchanged=True,substring_router_not_used=True,reference_state_hints_not_new_state_scores=True,biological_identity_and_CM_input_eligibility_not_granted=True)
    save(R/'outputs/accepted_reference_semantics_overlay_v1.json',manifest);storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='EXACT_REFERENCE_LABEL_SEMANTICS_FULL_AXIS_OVERLAY_ACCEPTED_IDENTITY_GATE_SEPARATE',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(R/'outputs/accepted_reference_semantics_overlay_v1.json')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
