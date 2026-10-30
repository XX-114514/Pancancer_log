"""Compare completed composition modules on a training-only common refinement.

The transported quantities are normalized feature-loading masses, not cell
counts. Membership transport is empirical and does not validate identities.
The comparison preserves unmatched factors and separates nsNMF W from W S.
"""
import argparse,datetime,hashlib,itertools,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from TME_CM_comparison_contract_v1 import partition_transport,map_module_mass,match_modules,weighted_jaccard

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def bound(b):
    assert sha(b['path'])==b['SHA256'];return Path(b['path'])

def main(R):
    cfg=load(R/'config/spec_v1.json');sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg);deadline=time.monotonic()+cfg['dependency_wait_seconds']
    items=[];units=None;common=None;unit_codes=None
    for root_name in cfg['fits']:
        root=Path(root_name);terminal=root/'provenance/supervisor_terminal_v1.json'
        while not terminal.exists():
            assert time.monotonic()<deadline,'Fit dependency remains incomplete';time.sleep(5)
        assert load(terminal)['actual_child_exit_code']==0,'A failed fit is not an unmatched module'
        run_cfg=load(root/'config/spec_v1.json');assert run_cfg['analysis_role']==cfg['accepted_purpose']
        acc=load(bound(run_cfg['input_acceptance']));membership=load(bound(acc['annotation_acceptance']))
        frequency=pd.read_csv(bound(acc['frequency']),sep='\t',index_col=0)
        metadata=pd.read_csv(bound(acc['metadata']),sep='\t')
        assert frequency.columns[1:].tolist()==metadata.sampleID.tolist()
        if units is None:units=metadata.sampleID.tolist()
        assert units==metadata.sampleID.tolist()
        nodemap={n['name']:n['code'] for n in membership['nodes']}
        rows=np.array([nodemap[n] for n in frequency.index],np.int64)
        rawnode=np.load(bound(membership['node_codes']));parents=np.load(bound(membership['parent_codes']))
        if common is None:
            prep=Path(cfg['preparation']);common=np.load(prep/'outputs/private_common_assessed_cell_mask_v1.npy')
            closure=Path(load(prep/'config/spec_v1.json')['unit_closure'])
            unit_codes=np.load(bound(load(closure/'outputs/accepted_fit_qualification_unit_bridge_v3.json')['unit_codes']))
        train=np.load(bound(acc['training_mask_binding']))
        selected_units=metadata.qualification_unit_code.to_numpy()
        selected_parents=[p['code'] for p in membership['parents'] if p['name'] in acc['requested_parent_panel']['parents']]
        mask=common&train&np.isin(unit_codes,selected_units)&np.isin(parents,selected_parents)
        if items:np.testing.assert_array_equal(mask,items[0]['mask'])
        assert (rawnode[mask]>=0).all()
        report=load(root/'outputs/CoVarNet_execution_report_v1.json')
        items.append(dict(root=root,scheme=membership['scheme_id'],node=rawnode,
            membership=membership,rows=rows,mask=mask,report=report,
            seeds=run_cfg['R_parameters']['final_seeds'],input_acceptance=run_cfg['input_acceptance']))
    assert len({x['scheme'] for x in items})==len(items)
    use=items[0]['mask'];signatures=np.column_stack([x['node'][use] for x in items])
    unique,atom=np.unique(signatures,axis=0,return_inverse=True)
    n_atoms=len(unique);atom_full=np.full(len(use),-1,np.int32);atom_full[use]=atom
    atom_path=R/'outputs/private_training_common_refinement_atom_codes_v1.npy';np.save(atom_path,atom_full)
    models=[];crosswalks=[]
    for item in items:
        T,measured,diag=partition_transport(atom_full,item['node'],n_atoms,len(item['membership']['nodes']),use)
        T=T[:,item['rows']];assert measured[item['rows']].all()
        # Direct member counting for one atom in every fitted node supplies
        # a second check on column-normalized transport and its orientation.
        for j,node in enumerate(item['rows']):
            at=np.flatnonzero(T[:,j]>0)[0];node_mask=use&(item['node']==node)
            scalar=np.count_nonzero(node_mask&(atom_full==at))/np.count_nonzero(node_mask)
            assert np.isclose(T[at,j],scalar)
        np.testing.assert_allclose(T.sum(axis=0),1,atol=1e-12)
        crosswalks.append(dict(scheme=item['scheme'],n_training_cells=diag['n_jointly_mapped'],
            n_atoms=n_atoms,n_fitted_nodes=len(item['rows']),one_scalar_membership_replay_per_node=True))
        for seed in item['seeds']:
            D=item['root']/'outputs'/('seed_'+str(seed));part={}
            for kind in ['W','A']:
                f=pd.read_csv(D/('private_'+kind+'_v1.tsv'),sep='\t',index_col=0)
                assert f.index.tolist()==[item['membership']['nodes'][j]['name'] for j in item['rows']]
                part[kind]=map_module_mass(f.to_numpy(float),T)
            h=pd.read_csv(D/'private_H_v1.tsv',sep='\t',index_col=0)
            assert h.columns.tolist()==units
            models.append(dict(scheme=item['scheme'],seed=seed,masses=part,H=h.to_numpy(float),rank=len(h)))
    rows=[];unmatched=[]
    for left,right in itertools.combinations(models,2):
        if left['scheme']!=right['scheme'] and left['seed']!=right['seed']:continue
        kind_of='same_scheme_seed_repeat' if left['scheme']==right['scheme'] else 'cross_scheme_same_seed'
        for basis in ['A','W']:
            match=match_modules(left['masses'][basis],right['masses'][basis],cfg['minimum_cosine_similarity'])
            for i,j,cosine in match['matches']:
                correlation=spearmanr(left['H'][i],right['H'][j]).statistic
                rows.append(dict(comparison=kind_of,left_scheme=left['scheme'],right_scheme=right['scheme'],
                    left_seed=left['seed'],right_seed=right['seed'],basis=basis,left_module=i+1,right_module=j+1,
                    cosine=cosine,weighted_jaccard=weighted_jaccard(left['masses'][basis][:,i],right['masses'][basis][:,j]),
                    activity_Spearman=float(correlation) if np.isfinite(correlation) else None,
                    same_samples=len(units),same_rank=left['rank']==right['rank']))
            unmatched.append(dict(comparison=kind_of,left_scheme=left['scheme'],right_scheme=right['scheme'],
                left_seed=left['seed'],right_seed=right['seed'],basis=basis,
                unmatched_left=[i+1 for i in match['unmatched_A']],unmatched_right=[i+1 for i in match['unmatched_B']],
                n_matched=len(match['matches']),left_rank=left['rank'],right_rank=right['rank'],
                full_pairwise_cosine=match['cosine_matrix'].tolist()))
    pd.DataFrame(rows).to_csv(R/'outputs/private_matched_module_comparison_v1.tsv',sep='\t',index=False)
    output=R/'outputs/accepted_membership_module_comparison_v1.json'
    save(output,dict(schema='accepted_TME_membership_module_comparison_v1',
        comparison_role=cfg['comparison_role'],n_schemes=len(items),n_models=len(models),
        n_common_training_cells=int(use.sum()),n_common_training_atoms=n_atoms,n_same_analysis_units=len(units),
        minimum_cosine_similarity=cfg['minimum_cosine_similarity'],crosswalk_checks=crosswalks,
        matches=rows,unmatched=unmatched,matching='one_to_one_Hungarian_with_unmatched_dummies',
        fitted_input_acceptances=[x['input_acceptance'] for x in items],
        mass_is_normalized_feature_loading_not_cell_count=True,
        common_refinement_is_empirical_training_only_not_external_validation=True,
        source_within_stratum_and_common_rank_checks_not_replaced=True,
        no_clinical_state_or_independent_patient_validation_claim=True))
    storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(
        status='TRAINING_MEMBERSHIP_CM_COMPARISON_COMPLETE_SOURCE_AND_STATE_GATES_SEPARATE',
        completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(output)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
