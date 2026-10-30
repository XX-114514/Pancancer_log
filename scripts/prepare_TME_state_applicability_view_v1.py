"""Bind reference identity contexts to exact resource-program applicability.

This sidecar never rewrites cell identities or raw programme scores. Primary
scope follows explicit resource parents; legacy-use intersection is retained.
Every use also requires its separate measured-gene mask.
"""
import argparse,datetime,hashlib,json,sys
from pathlib import Path
import numpy as np
from TME_state_program_scope_policy_v1 import applicable

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def binding(p):return dict(path=str(p),SHA256=sha(p))
def bound(b):
    assert sha(b['path'])==b['SHA256'];return Path(b['path'])

def main(R):
    cfg=load(R/'config/spec_v1.json');sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg);S=Path(cfg['scope']);M=Path(cfg['semantics']);P=Path(cfg['preparation']);C=Path(cfg['coverage'])
    for root in [S,M,P,C]:assert load(root/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    sp=S/'outputs/accepted_shared_compartment_scope_v1.json';scope=load(sp)
    assert sha(sp)==load(S/'provenance/terminal_v1.json')['acceptance_SHA256']
    sem=load(bound(scope['reference_semantics_acceptance']))
    assert Path(scope['reference_semantics_acceptance']['path']).parent.parent==M
    low=next(x for x in sem['models'] if x['model']=='immune_low')
    identity=np.load(bound(low['identity_proposal']));p=np.load(bound(scope['codes']))
    selected=np.load(P/'outputs/private_common_assessed_cell_mask_v1.npy')
    assert selected.shape==p.shape==identity.shape and selected.dtype.kind=='b';assert (p[selected]>=0).all()
    policy=load(S/'config/frozen_compartment_policy_v1.json');assert sha(S/'config/frozen_compartment_policy_v1.json')==scope['policy_SHA256']
    family={x:k for k,r in policy['compartments'].items() for x in r['identity_parents']}
    rows=sem['identities'];imap=np.array([scope['compartments'].index(family[x['parent']]) if x['parent'] in family else -1 for x in rows],np.int16)
    valid=(identity>=0)&selected;good=valid.copy();good[valid]&=imap[identity[valid]]==p[valid]
    reviewed=np.where(good,identity,-1)
    signatures=np.column_stack([p[selected],reviewed[selected]])
    pairs,inverse=np.unique(signatures,axis=0,return_inverse=True)
    context=np.full(len(p),-1,np.int16);context[selected]=inverse
    counts=np.bincount(inverse,minlength=len(pairs));contexts=[]
    for k,(parent,i) in enumerate(pairs):
        parent=int(parent);i=int(i);r=rows[i] if i>=0 else {}
        contexts.append(dict(code=k,compartment=scope['compartments'][parent],
            reference_parent=r.get('parent'),canonical_identity=r.get('canonical_identity'),
            identity_scope='in_scope_reference_proposal' if i>=0 else 'working_compartment_only',
            n_cells=int(counts[k]),reference_identity_not_independent_biological_truth=True))
    cp=R/'outputs/private_reference_state_context_code_v1.npy';np.save(cp,context)
    catalog=load(bound(cfg['catalog']));curation=load(bound(cfg['curation']))
    cmap={r['program_id']:r for r in curation['program_metadata']}
    covp=C/'outputs/state_native_coverage_preparation_v1.json';cov=load(covp)
    assert sha(covp)==load(C/'provenance/terminal_v1.json')['manifest_SHA256']
    assert cfg['catalog']['SHA256']==cov['catalog_SHA256']==curation['scored_catalog_SHA256']
    assert cfg['curation']['SHA256']==cov['curation_SHA256']
    mode=cfg.get('scope_policy_mode','legacy_default_intersection')
    programs=catalog['programs'];matrix=np.zeros((len(contexts),len(programs)),bool);legacy=np.zeros_like(matrix);ledger=[]
    for j,program in enumerate(programs):
        meta=cmap[program['program_id']];reasons=[]
        for i,c in enumerate(contexts):
            ok,reason=applicable(program['reference_cell_type'],meta,c,mode=mode);matrix[i,j]=ok;reasons.append(reason)
            legacy[i,j]=applicable(program['reference_cell_type'],meta,c,mode='legacy_default_intersection')[0]
        ledger.append(dict(program_id=program['program_id'],default_resource_scope=program['reference_cell_type'],
            Cytopus_parents=meta['reference_parents'],exact_resource_gene_match=meta.get('exact_Cytopus_gene_match',False),
            n_eligible_contexts=int(matrix[:,j].sum()),n_legacy_eligible_contexts=int(legacy[:,j].sum()),n_eligible_cells_before_gene_measurement=int(counts[matrix[:,j]].sum()),
            scope_reasons=sorted(set(reasons)),original_article_citation_status=meta['citation_status'],
            function_validated=False,module_admitted=False))
    mp=R/'outputs/private_program_by_reference_context_applicability_v1.npy';np.save(mp,matrix)
    lp=R/'outputs/private_legacy_default_program_applicability_v1.npy';np.save(lp,legacy)
    assert not np.any(legacy&~matrix)
    # Check every realized context after readback. Outside scope is a masked
    # observation even when the package produced a large raw numerical score.
    read=np.load(mp)
    for i,c in enumerate(contexts):
        for j,program in enumerate(programs):
            assert bool(read[i,j])==applicable(program['reference_cell_type'],cmap[program['program_id']],c,mode=mode)[0]
    np.testing.assert_array_equal(np.load(cp)[selected],inverse)
    out=R/'outputs/accepted_state_program_applicability_view_v1.json'
    save(out,dict(schema='accepted_TME_state_program_reference_applicability_v1',
        accepted_purpose='diagnostic_expression_program_scope_before_readout_and_module_admission',
        n_fit=len(p),n_assessed=int(selected.sum()),n_contexts=len(contexts),n_programs=len(programs),
        scope_policy_mode=mode,legacy_default_sensitivity=binding(lp),
        context_codes=binding(cp),context_program_applicability=binding(mp),contexts=contexts,programs=ledger,
        catalog=cfg['catalog'],curation=cfg['curation'],source_gene_coverage_manifest=binding(covp),
        scope_acceptance=binding(sp),reference_semantics_acceptance=scope['reference_semantics_acceptance'],
        common_cell_mask_SHA256=sha(P/'outputs/private_common_assessed_cell_mask_v1.npy'),
        measurement_and_context_mask_must_both_pass=True,outside_scope_is_NaN_not_negative_state=True,
        no_automatic_state_names_used_as_new_score=True,no_identity_or_raw_score_overwrite=True,
        readout_reliability_and_functional_and_module_admission_not_granted=True))
    storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='EXACT_RESOURCE_AND_REFERENCE_IDENTITY_STATE_SCOPE_INTERSECTION_ACCEPTED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(out)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
