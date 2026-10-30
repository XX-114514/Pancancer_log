"""Recount private CM inputs independently from accepted cell memberships.

Preparation alone cannot grant execution eligibility. This auditor retains
all denominator children, checks actual measurement opportunity, and never
turns unavailable observations into negative cell counts.
"""
import argparse,datetime,hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def binding(p):return dict(path=str(p),SHA256=sha(p))
def bound(b):
    assert sha(b['path'])==b['SHA256'];return Path(b['path'])

def main(R):
    cfg=load(R/'config/spec_v1.json');P=Path(cfg['preparation'])
    sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg)
    assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    prep=load(P/'outputs/comparable_CM_panel_preparation_v1.json')
    spec=load(P/'config/spec_v1.json')
    # Bind the group-reference table to its independently replayed execution,
    # rather than accepting a successful directory terminal by itself.
    sr=Path(cfg['SingleR_auditor'])
    assert load(sr/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    sap=sr/'outputs/accepted_SingleR_group_execution_v1.json'
    assert sha(sap)==load(sr/'provenance/terminal_v1.json')['acceptance_SHA256']
    sr_acceptance=load(sap)
    closure=Path(spec['unit_closure']);design=Path(spec['design'])
    rec=load(closure/'outputs/accepted_fit_qualification_unit_bridge_v3.json')
    unit=np.load(bound(rec['unit_codes']));N=len(unit)
    table=pd.read_csv(design/'outputs/private_analysis_unit_context_and_coverage_v1.tsv',sep='\t')
    nu=len(table);technical=np.load(design/'outputs/private_technical_material_candidate_mask_v1.npy')
    train=np.load(bound(spec['training_cell_mask']));common=technical.copy();schemes={}
    for b in spec['schemes']:
        a=load(bound(b));p=np.load(bound(a['parent_codes']));n=np.load(bound(a['node_codes']))
        assert a['SingleR_proposals']['SHA256']==sr_acceptance['table_SHA256']
        bound(a['SingleR_proposals'])
        opportunity=np.load(bound(a['node_measurement_mask']))
        assert len(p)==len(n)==N and np.array_equal(p>=0,n>=0)
        common&=p>=0;schemes[a['scheme_id']]=dict(a=a,p=p,n=n,opportunity=opportunity)
    np.testing.assert_array_equal(common,np.load(P/'outputs/private_common_assessed_cell_mask_v1.npy'))
    # Alternative aggregation implementation: indexed addition rather than
    # the preparation driver's packed node-by-unit bincount.
    rows=np.flatnonzero(common)
    for s in schemes.values():
        a=s['a'];s['counts']=np.zeros((len(a['nodes']),nu),np.int64)
        s['parent_counts']=np.zeros((len(a['parents']),nu),np.int64)
        np.add.at(s['counts'],(s['n'][rows],unit[rows]),1)
        np.add.at(s['parent_counts'],(s['p'][rows],unit[rows]),1)
        assert not s['counts'][~s['opportunity']].any()
    accepted=[];panel_samples={};probes=[]
    for b in prep['candidate_inputs']:
        cp=bound(b);a=load(cp);assert a['schema']=='candidate_TME_CM_scheme_input_v1'
        mb=bound(a['annotation_acceptance']);membership=load(mb)
        sid=membership['scheme_id'];s=schemes[sid]
        fq=pd.read_csv(bound(a['frequency']),sep='\t',index_col=0)
        meta=pd.read_csv(bound(a['metadata']),sep='\t')
        assert fq.columns[1:].tolist()==meta.sampleID.tolist()
        assert not meta.sampleID.duplicated().any() and fq.index.is_unique
        indices=np.array([int(x.removeprefix('UNIT_')) for x in meta.sampleID],np.int64)
        np.testing.assert_array_equal(indices,meta.qualification_unit_code.to_numpy())
        assert (indices>=0).all() and (indices<nu).all()
        for col in table.columns:
            left=meta[col].fillna('').astype(str).tolist()
            right=table.iloc[indices][col].fillna('').astype(str).tolist()
            assert left==right,'Source/material/block metadata changed'
        key=a['requested_parent_panel']['id']
        if key not in panel_samples:panel_samples[key]=indices.tolist()
        assert panel_samples[key]==indices.tolist(),'Scheme-specific sample fallback'
        allowed=np.load(design/'outputs'/('private_'+a['requested_parent_panel']['material_panel']+'_unit_mask_v1.npy'))
        assert allowed[indices].all()
        for u in indices:
            assert train[unit==u].all(),'Partial work-unit training leakage'
            assert np.count_nonzero(common&(unit==u))>=spec['minimum_total_cells']
        nodemap={n['name']:n for n in membership['nodes']}
        parents=membership['parents'];observed=fq.iloc[:,1:].to_numpy(float)
        assert np.isfinite(observed).all() and ((observed>=0)&(observed<=1)).all()
        nzeros=0
        for j,name in enumerate(fq.index):
            node=nodemap[name];k=node['code'];parent=node['parent_code']
            assert node['discovery_eligible']
            assert parents[parent]['name']==fq.iloc[j,0]
            assert parents[parent]['name'] in a['requested_parent_panel']['parents']
            assert s['opportunity'][k,indices].all()
            denominator=s['parent_counts'][parent,indices]
            numerator=s['counts'][k,indices]
            assert (denominator>=spec['minimum_parent_cells']).all()
            expected=numerator/denominator
            np.testing.assert_allclose(observed[j],expected,rtol=1e-12,atol=1e-14)
            assert numerator.sum()>=spec['minimum_node_cells']
            assert np.count_nonzero(numerator)>=(spec['minimum_node_observed_units'])
            assert np.ptp(expected)>0
            # Direct scalar checks bypass both full-array counting methods.
            for mask,case in [(numerator>0,'observed_positive'),(numerator==0,'observed_zero')]:
                at=np.flatnonzero(mask)
                if at.size and not any(x['case']==case and x['panel']==key and x['scheme']==sid for x in probes):
                    u=indices[at[0]];c=common&(unit==u)
                    num=np.count_nonzero(c&(s['n']==k));den=np.count_nonzero(c&(s['p']==parent))
                    assert np.isclose(observed[j,at[0]],num/den)
                    probes.append(dict(panel=key,scheme=sid,case=case,direct_cell_recount_pass=True))
            nzeros+=int(np.count_nonzero(numerator==0))
        blocks=meta.patient_candidate_block.nunique()
        assert len(indices)>=spec['minimum_units'] and blocks>=spec['minimum_blocks']
        assert len(fq)>=spec['minimum_nodes']
        assert a['rank_upper_bound_by_candidate_blocks']==min(spec['maximum_rank'],len(fq)-1,blocks//4)
        result=dict(a);result.update(schema='accepted_TME_CM_scheme_input_v1',
            independent_input_audit_pending=False,independent_all_cell_recount_pass=True,
            input_candidate=binding(cp),n_explicit_observed_zero_frequencies=nzeros,
            uncertainty_nodes_kept_in_complete_denominator=True,
            work_units_and_candidate_blocks_not_certified_independent_patients=True,
            audited_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
        op=R/'outputs'/cp.name.replace('_candidate_input_','_accepted_input_');save(op,result)
        accepted.append(dict(panel=key,scheme=sid,acceptance=binding(op),
            n_units=len(indices),n_nodes=len(fq),maximum_rank=a['rank_upper_bound_by_candidate_blocks']))
    out=R/'outputs/accepted_comparable_CM_inputs_v1.json'
    save(out,dict(schema='accepted_TME_comparable_CM_inputs_v1',inputs=accepted,
        n_fit=N,n_common_assessed_cells=int(common.sum()),independent_count_implementation='numpy_add_at',
        direct_scalar_positive_and_observed_zero_recounts=probes,
        all_comparison_panels_have_identical_observations=True,
        input_preparation=binding(P/'outputs/comparable_CM_panel_preparation_v1.json'),
        no_clinical_or_independent_patient_validation_claim=True))
    storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(
        status='COMPARABLE_CM_INPUTS_INDEPENDENTLY_RECOUNTED_ACCEPTED',
        completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        n_accepted_inputs=len(accepted),acceptance_SHA256=sha(out)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
