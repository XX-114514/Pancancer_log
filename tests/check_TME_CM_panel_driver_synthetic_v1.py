"""Hand-counted end-to-end panel I/O counterexamples; no human inputs."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v))
def binding(p):return dict(path=str(p),SHA256=sha(p))
def array(p,v):p.parent.mkdir(parents=True,exist_ok=True);np.save(p,v);return binding(p)
def main(R):
    sys.path.insert(0,str(R/'scripts'))
    from prepare_TME_CM_comparable_panels_v1 import main as prepare
    cfg=json.loads((R/'config/spec_v1.json').read_text());base=R/'outputs/synthetic_only';base.mkdir()
    units=[];parent=[];nodes0=[];nodes1=[];nu=16
    for i in range(nu):
        a0=10+5*(i%3)
        for j in range(60):
            units.append(i);parent.append(0);nodes0.append(0 if j<a0 else 1);nodes1.append((0 if j<3 else 1) if j<a0 else 2)
        if i!=2:
            b0=0 if i==0 else 5+5*(i%4)
            for j in range(30):units.append(i);parent.append(1);nodes0.append(2 if j<b0 else 3);nodes1.append(3 if j<b0 else 4)
    units=np.asarray(units,np.int32);parent=np.asarray(parent,np.int16);N=len(units);D=base/'design';U=base/'unit_bridge'
    for d in [D,U]:put(d/'provenance/supervisor_terminal_v1.json',dict(actual_child_exit_code=0))
    ub=array(U/'outputs/private_units.npy',units);put(U/'outputs/accepted_fit_qualification_unit_bridge_v3.json',dict(unit_codes=ub,n_fit=N))
    mask=np.ones(N,bool);mask[:2]=False;array(D/'outputs/private_technical_material_candidate_mask_v1.npy',mask);array(D/'outputs/private_all_qualified_tumour_unit_mask_v1.npy',np.ones(nu,bool))
    pd.DataFrame(dict(qualification_unit_code=np.arange(nu),source=['ARTIFICIAL_'+str(i//4) for i in range(nu)],patient_candidate_block=['ARTIFICIAL_BLOCK_'+str(i) for i in range(nu)])).to_csv(D/'outputs/private_analysis_unit_context_and_coverage_v1.tsv',sep='\t',index=False)
    training=array(base/'training_mask.npy',units<15);training['frozen_before_module_fitting']=True;schemes=[]
    for k,(node,mapping) in enumerate([(nodes0,[0,0,1,1]),(nodes1,[0,0,0,1,1])]):
        node=np.asarray(node,np.int16);ap=base/('scheme_'+str(k)+'.json');p=array(base/('parent_'+str(k)+'.npy'),parent);n=array(base/('node_'+str(k)+'.npy'),node)
        put(ap,dict(schema='accepted_TME_identity_scheme_membership_v1',accepted_purpose='synthetic_validation',scheme_id='ARTIFICIAL_SCHEME_'+str(k),final_membership_codes_not_display_names=True,unresolved_fine_identity_has_explicit_fallback_node=True,parent_codes=p,node_codes=n,parents=[dict(code=j,name='ARTIFICIAL_PARENT_'+str(j)) for j in range(2)],nodes=[dict(code=j,name='ARTIFICIAL_NODE_'+str(j),parent_code=v,discovery_eligible=True) for j,v in enumerate(mapping)],input_axis_SHA256='ARTIFICIAL_AXIS_NOT_HUMAN'))
        schemes.append(binding(ap))
    spec=dict(cfg,design=str(D),unit_closure=str(U),training_cell_mask=training,schemes=schemes,accepted_purpose='synthetic_validation',panels=[dict(id='ARTIFICIAL_COMPLETE',material_panel='all_qualified_tumour',parents=['ARTIFICIAL_PARENT_0','ARTIFICIAL_PARENT_1'])],minimum_total_cells=20,minimum_parent_cells=20,minimum_units=12,minimum_blocks=8,minimum_node_cells=10,minimum_node_observed_units=4,minimum_nodes=4,maximum_rank=12)
    target=base/'complete';(target/'config').mkdir(parents=True);(target/'outputs').mkdir();(target/'provenance').mkdir();put(target/'config/spec_v1.json',spec);prepare(target)
    index=json.loads((target/'outputs/comparable_CM_panel_preparation_v1.json').read_text());assert len(index['candidate_inputs'])==2
    for j,b in enumerate(index['candidate_inputs']):
        a=json.loads(Path(b['path']).read_text());assert a['schema']=='candidate_TME_CM_scheme_input_v1' and a['independent_input_audit_pending'];f=pd.read_csv(a['frequency']['path'],sep='\t',index_col=0)
        assert 'UNIT_2' not in f and 'UNIT_15' not in f;assert f.columns[1:].tolist()==['UNIT_'+str(i) for i in range(15) if i!=2]
        assert np.isclose(float(f.loc['ARTIFICIAL_NODE_0','UNIT_0']),8/58 if j==0 else 1/58)
        assert float(f.loc['ARTIFICIAL_NODE_'+str(2 if j==0 else 3),'UNIT_0'])==0
        for major,g in f.groupby('majorCluster'):np.testing.assert_allclose(g.iloc[:,1:].sum(axis=0),1,atol=1e-12)
        assert a['rank_upper_bound_by_candidate_blocks']==3
    # A single cell of a held-out material cannot leak into discovery input.
    bad=base/'partial_material';(bad/'config').mkdir(parents=True);(bad/'outputs').mkdir();(bad/'provenance').mkdir();partial=units<15;partial[np.flatnonzero(units==15)[0]]=True;b=array(base/'partial_training_mask.npy',partial);b['frozen_before_module_fitting']=True;put(bad/'config/spec_v1.json',dict(spec,training_cell_mask=b))
    try:prepare(bad)
    except AssertionError:pass
    else:raise AssertionError('Partial material train/test leakage was accepted')
    assert not (bad/'provenance/terminal_v1.json').exists()
    # The same observed zero is retained when measured above, but cannot be
    # counted as a negative class if that node was not measurable in the unit.
    measurement_schemes=[]
    for k,original in enumerate(schemes):
        a=json.loads(Path(original['path']).read_text());mask=np.ones((len(a['nodes']),nu),bool)
        if k==0:mask[2,0]=False
        a['node_measurement_mask']=array(base/('opportunity_'+str(k)+'.npy'),mask);ap=base/('measured_scheme_'+str(k)+'.json');put(ap,a);measurement_schemes.append(binding(ap))
    measured=base/'measurement_opportunity';(measured/'config').mkdir(parents=True);(measured/'outputs').mkdir();(measured/'provenance').mkdir();put(measured/'config/spec_v1.json',dict(spec,schemes=measurement_schemes));prepare(measured)
    check=json.loads((measured/'outputs/comparable_CM_panel_preparation_v1.json').read_text());assert len(check['candidate_inputs'])==2;axes=[]
    for b in check['candidate_inputs']:
        a=json.loads(Path(b['path']).read_text());fq=pd.read_csv(a['frequency']['path'],sep='\t',index_col=0);axes.append(fq.columns.tolist());assert 'UNIT_0' not in fq;assert a['coverage']['node_opportunity_applied_without_missing_to_zero']
    assert axes[0]==axes[1]
    put(R/'outputs/hand_counted_panel_driver_acceptance_v1.json',dict(synthetic_only=True,whole_material_split_checked=True,absent_parent_not_subgroup_zero=True,observed_subgroup_zero_retained=True,unmeasured_node_zero_excludes_unit_from_both_schemes=True,technical_hold_changes_complete_denominator=True,split_node_comparison_same_units=True,candidate_not_accepted_input=True,rank_limited_by_candidate_blocks=True,n_hand_counted_schemes=2))
    put(R/'provenance/terminal_v1.json',dict(status='HAND_COUNTED_COMPARABLE_PANEL_DRIVER_COUNTEREXAMPLES_PASS',synthetic_only=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);main(Path(p.parse_args().output))
