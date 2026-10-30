"""Build frozen complete-denominator panels from accepted annotation memberships.

This entry does not infer identities, select a clinical endpoint, impute absent
parents, or choose a parent panel from its downstream CM performance. Input
membership acceptance and the requested parent panels must already be frozen.
"""
import argparse,datetime,hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def bound(binding):
    assert sha(binding['path'])==binding['SHA256'],'Input binding changed'
    return Path(binding['path'])
def main(R):
    cfg=load(R/'config/spec_v1.json');sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    from TME_CM_comparison_contract_v1 import frequency_with_measurement_mask
    storage_gate(R,cfg);design=Path(cfg['design']);unit_closure=Path(cfg['unit_closure'])
    for D in [design,unit_closure]:assert load(D/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    acceptance=load(unit_closure/'outputs/accepted_fit_qualification_unit_bridge_v3.json');unit=np.load(bound(acceptance['unit_codes']));N=len(unit)
    table=pd.read_csv(design/'outputs/private_analysis_unit_context_and_coverage_v1.tsv',sep='\t');assert table.qualification_unit_code.tolist()==list(range(len(table)));nu=len(table)
    technical=np.load(design/'outputs/private_technical_material_candidate_mask_v1.npy');assert technical.shape==(N,) and technical.dtype.kind=='b'
    training=np.load(bound(cfg['training_cell_mask']));assert training.shape==(N,) and training.dtype.kind=='b';assert cfg['training_cell_mask']['frozen_before_module_fitting']
    # A whole material is one analysis observation. Never split cells of one
    # material across discovery/held-out roles by a random cell mask.
    exposed=np.bincount(unit[training],minlength=nu);all_n=np.bincount(unit,minlength=nu);assert np.all((exposed==0)|(exposed==all_n))
    shared=technical.copy();schemes=[];axis=None
    for binding in cfg['schemes']:
        ap=bound(binding);a=load(ap);assert a['schema']=='accepted_TME_identity_scheme_membership_v1';assert a['accepted_purpose']==cfg['accepted_purpose'];assert a['final_membership_codes_not_display_names'] and a['unresolved_fine_identity_has_explicit_fallback_node']
        p=np.load(bound(a['parent_codes']));n=np.load(bound(a['node_codes']));assert p.shape==n.shape==(N,) and p.dtype.kind in 'iu' and n.dtype.kind in 'iu';assert (p>=-1).all() and (n>=-1).all()
        parents=a['parents'];nodes=a['nodes'];assert [r['code'] for r in parents]==list(range(len(parents))) and [r['code'] for r in nodes]==list(range(len(nodes)))
        assert len({r['name'] for r in parents})==len(parents) and len({r['name'] for r in nodes})==len(nodes);assert (p<len(parents)).all() and (n<len(nodes)).all();assert np.array_equal(p>=0,n>=0)
        mapping=np.asarray([r['parent_code'] for r in nodes],np.int64);assert ((mapping>=0)&(mapping<len(parents))).all();assigned=n>=0;np.testing.assert_array_equal(mapping[n[assigned]],p[assigned])
        if axis is None:axis=a['input_axis_SHA256']
        assert a['input_axis_SHA256']==axis
        if a.get('node_measurement_mask'):
            opportunity=np.load(bound(a['node_measurement_mask']));assert opportunity.shape==(len(nodes),nu) and opportunity.dtype.kind=='b'
        else:
            assert cfg['accepted_purpose']=='synthetic_validation','Production schemes require explicit per-node measurement opportunity'
            opportunity=np.ones((len(nodes),nu),bool)
        shared&=assigned;schemes.append(dict(manifest=a,manifest_binding=binding,parent=p,node=n,mapping=mapping,opportunity=opportunity))
    # Shared assessed-cell domain is identical across schemes. Each scheme's
    # complete parent denominator is retained before omitting discovery nodes.
    np.save(R/'outputs/private_common_assessed_cell_mask_v1.npy',shared);common_total=np.bincount(unit[shared],minlength=nu);report=[];candidate_bindings=[]
    for s in schemes:
        a=s['manifest'];nn=len(a['nodes']);np_=len(a['parents']);count=np.bincount(s['node'][shared].astype(np.int64)*nu+unit[shared],minlength=nn*nu).reshape(nn,nu);parent=np.bincount(s['parent'][shared].astype(np.int64)*nu+unit[shared],minlength=np_*nu).reshape(np_,nu)
        for k in range(np_):np.testing.assert_array_equal(count[s['mapping']==k].sum(axis=0),parent[k])
        assert (count[~s['opportunity']]==0).all(),'A class called in an unmeasurable node/unit needs its definition corrected'
        s.update(count=count,parent_counts=parent);np.save(R/'outputs'/('private_'+a['scheme_id']+'_complete_parent_counts_v1.npy'),parent)
    for panel in cfg['panels']:
        base=np.load(design/'outputs'/('private_'+panel['material_panel']+'_unit_mask_v1.npy'));assert base.shape==(nu,) and base.dtype.kind=='b';use=base&(common_total>=cfg['minimum_total_cells'])&(exposed>0);requested=panel['parents'];assert requested and len(set(requested))==len(requested)
        reason=None
        for s in schemes:
            by={r['name']:r['code'] for r in s['manifest']['parents']}
            if not set(requested)<=set(by):reason='requested_parent_absent_from_scheme_ontology';break
            use&=(s['parent_counts'][[by[n] for n in requested]]>=cfg['minimum_parent_cells']).all(axis=0)
        # Freeze a pre-coverage node pool using only genuinely measured units.
        # Unapplied/unsupported classifiers must not create apparent variation
        # by contributing artificial zeros outside their observation scope.
        pools=[]
        if reason is None:
            preliminary=use.copy()
            for s in schemes:
                a=s['manifest'];f,measured=frequency_with_measurement_mask(s['count'],s['mapping'],s['parent_counts'],cfg['minimum_parent_cells']);pool=[]
                for j,node in enumerate(a['nodes']):
                    if a['parents'][node['parent_code']]['name'] not in requested or not node['discovery_eligible']:continue
                    observed=preliminary&s['opportunity'][j]&measured[j]
                    if s['count'][j,observed].sum()>=cfg['minimum_node_cells'] and (s['count'][j,observed]>0).sum()>=cfg['minimum_node_observed_units'] and observed.sum() and np.ptp(f[j,observed])>0:pool.append(j)
                pools.append(pool)
                if pool:use&=s['opportunity'][pool].all(axis=0)
        else:pools=[[] for s in schemes]
        nblocks=int(table.loc[use,'patient_candidate_block'].nunique())
        if reason is None and (use.sum()<cfg['minimum_units'] or nblocks<cfg['minimum_blocks']):reason='insufficient_complete_comparable_units_or_candidate_blocks'
        rows=[]
        for s,pool in zip(schemes,pools):
            a=s['manifest'];f,measured=frequency_with_measurement_mask(s['count'],s['mapping'],s['parent_counts'],cfg['minimum_parent_cells']);eligible=[]
            for j,node in enumerate(a['nodes']):
                if j not in pool:continue
                if a['parents'][node['parent_code']]['name'] not in requested or not node['discovery_eligible']:continue
                assert s['opportunity'][j,use].all()
                if s['count'][j,use].sum()<cfg['minimum_node_cells'] or (s['count'][j,use]>0).sum()<cfg['minimum_node_observed_units']:continue
                if use.sum() and np.isfinite(f[j,use]).all() and np.ptp(f[j,use])>0:eligible.append(j)
            status=reason or ('insufficient_variable_discovery_nodes' if len(eligible)<cfg['minimum_nodes'] else 'prepared_waiting_independent_input_audit')
            row=dict(panel=panel['id'],scheme=a['scheme_id'],status=status,n_units=int(use.sum()),n_candidate_blocks=nblocks,n_sources=int(table.loc[use,'source'].nunique()),n_nodes=len(eligible),n_precoverage_discovery_nodes=len(pool),node_opportunity_applied_without_missing_to_zero=True,parents_requested=requested,n_common_cells=int(common_total[use].sum()),candidate_blocks_are_not_certified_independent_patients=True)
            report.append(row);rows.append((s,eligible,row))
        # No scheme-specific sample fallback: a head-to-head panel is usable
        # only when every scheme passes on exactly these same observations.
        if any(row['status']!='prepared_waiting_independent_input_audit' for s,e,row in rows):
            for s,e,row in rows:row['comparison_ready']=False
            continue
        for s,eligible,row in rows:
            a=s['manifest'];prefix=panel['id']+'__'+a['scheme_id'];fq,mask=frequency_with_measurement_mask(s['count'],s['mapping'],s['parent_counts'],cfg['minimum_parent_cells']);assert mask[np.ix_(eligible,use)].all();out=fq[np.ix_(eligible,use)];assert np.isfinite(out).all()
            ids=['UNIT_'+str(i) for i in np.flatnonzero(use)];frame=pd.DataFrame(out,index=[a['nodes'][j]['name'] for j in eligible],columns=ids);frame.insert(0,'majorCluster',[a['parents'][s['mapping'][j]]['name'] for j in eligible]);fp=R/'outputs'/('private_'+prefix+'_frequency_v1.tsv');frame.to_csv(fp,sep='\t',index_label='subCluster')
            meta=table.loc[use].copy();meta.insert(0,'sampleID',ids);mp=R/'outputs'/('private_'+prefix+'_metadata_v1.tsv');meta.to_csv(mp,sep='\t',index=False)
            row['comparison_ready']=True;candidate=dict(schema='candidate_TME_CM_scheme_input_v1',accepted_purpose=cfg['accepted_purpose'],independent_input_audit_pending=True,complete_measured_denominators_no_missing_to_zero=True,annotation_and_unit_masks_bound=True,node_measurement_mask_binding=a.get('node_measurement_mask'),frequency=dict(path=str(fp),SHA256=sha(fp)),metadata=dict(path=str(mp),SHA256=sha(mp)),annotation_acceptance=s['manifest_binding'],common_assessed_cell_mask_SHA256=sha(R/'outputs/private_common_assessed_cell_mask_v1.npy'),training_mask_binding=cfg['training_cell_mask'],requested_parent_panel=panel,coverage=row,rank_upper_bound_by_candidate_blocks=min(cfg['maximum_rank'],len(eligible)-1,nblocks//4))
            cp=R/'outputs'/('private_'+prefix+'_candidate_input_v1.json');save(cp,candidate);candidate_bindings.append(dict(path=str(cp),SHA256=sha(cp)))
    storage_gate(R,cfg);save(R/'outputs/comparable_CM_panel_preparation_v1.json',dict(schema='TME_CM_comparable_panel_preparation_v1',n_fit=N,n_common_assessed_cells=int(shared.sum()),panels=report,candidate_inputs=candidate_bindings,accepted_analysis_inputs_not_yet_granted=True,input_axis_SHA256=axis,no_missing_to_zero=True,shared_cell_domain_and_unit_panels=True))
    save(R/'provenance/terminal_v1.json',dict(status='COMPARABLE_CM_PANELS_PREPARED_PENDING_INDEPENDENT_INPUT_AUDIT',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_candidate_inputs=len(candidate_bindings),manifest_SHA256=sha(R/'outputs/comparable_CM_panel_preparation_v1.json')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
