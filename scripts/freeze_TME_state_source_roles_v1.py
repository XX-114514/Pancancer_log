"""Inherit only frozen source roles, never legacy lineage/patient qualification."""
import argparse,datetime,hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def binding(p):return dict(path=str(p),SHA256=sha(p))

def main(R):
    cfg=load(R/'config/spec_v1.json');sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg);H=Path(cfg['model']);P=Path(cfg['preparation']);b=cfg['existing_source_split']
    assert sha(b['path'])==b['SHA256'];old=pd.read_csv(b['path'],sep='\t')
    assert not old.source.duplicated().any()
    sources=load(H/'outputs/scope_and_common_genes_v1.json')['sources']
    roles=old.set_index('source').split.to_dict();assert set(sources)<=set(roles)
    assert set(roles[s] for s in sources)<={'development','lockbox'}
    assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    selected=np.load(P/'outputs/private_common_assessed_cell_mask_v1.npy')
    with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:source=f['source_code']
    assert source.shape==selected.shape and selected.dtype.kind=='b'
    dev=selected&np.isin(source,[i for i,s in enumerate(sources) if roles[s]=='development'])
    locked=selected&np.isin(source,[i for i,s in enumerate(sources) if roles[s]=='lockbox'])
    np.testing.assert_array_equal(dev|locked,selected);assert not np.any(dev&locked)
    dp=R/'outputs/private_state_development_cell_mask_v1.npy';lp=R/'outputs/private_state_lockbox_cell_mask_v1.npy'
    np.save(dp,dev);np.save(lp,locked)
    out=R/'outputs/accepted_state_source_roles_v1.json'
    save(out,dict(schema='accepted_TME_state_learning_source_roles_v1',
        inherited_source_split=b,n_fit=len(source),n_common_assessed=int(selected.sum()),
        n_development_cells=int(dev.sum()),n_lockbox_cells=int(locked.sum()),
        n_development_sources=sum(roles[s]=='development' for s in sources),
        n_lockbox_sources=sum(roles[s]=='lockbox' for s in sources),
        source_roles=[dict(source=s,role=roles[s]) for s in sources],
        development_mask=binding(dp),lockbox_mask=binding(lp),
        input_axis_SHA256=sha(H/'outputs/private_cell_axes_v1.npz'),
        common_cell_mask_SHA256=sha(P/'outputs/private_common_assessed_cell_mask_v1.npy'),
        only_source_role_inherited_not_legacy_TNK_scope_or_patient_counts=True,
        fixed_program_quantification_may_cover_both_roles=True,
        state_feature_selection_and_pattern_training_development_only=True,
        state_feature_and_pattern_freeze_required_before_lockbox_biology_review=True,
        all_source_numeric_audit_not_biological_feature_selection=True,
        heldout_is_conditional_on_already_frozen_all_source_atlas=True,
        not_end_to_end_external_or_independent_patient_validation=True,
        clinical_outcomes_used=False))
    storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='EXISTING_STATE_SOURCE_LEARNING_ROLES_FROZEN_FOR_NEW_TME',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(out)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
