"""Finalize measurement candidates from development-only readout diagnostics.

This stage excludes lucky partial splits and unsupported one-material signals.
It does not admit function, clinical mechanism or a final state feature model.
"""
import argparse,datetime,hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from TME_state_readout_candidate_policy_v1 import candidate_status

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def binding(p):return dict(path=str(p),SHA256=sha(p))
def bound(b):
    assert sha(b['path'])==b['SHA256'];return Path(b['path'])

def main(R):
    cfg=load(R/'config/spec_v1.json');P=Path(cfg['readouts'])
    sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg);assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    ip=P/'outputs/accepted_development_state_readouts_v1.json';assert sha(ip)==load(P/'provenance/terminal_v1.json')['acceptance_SHA256']
    m=load(ip);spec=load(P/'config/spec_v1.json');assert m['n_lockbox_score_arrays_read']==0
    app=load(bound(m['applicability']));roles=load(bound(m['source_roles']));assert app['scope_policy_mode']=='resource_hierarchy_primary'
    allowed={r['source'] for r in roles['source_roles'] if r['role']=='development'}
    d=pd.read_csv(bound(m['diagnostics']),sep='\t');names=m['program_ids'];by_name={s:i for i,s in enumerate(names)}
    assert set(d.source)<=allowed;source_index={x['source']:x for x in m['readout_sources']}
    extra=[]
    for source,group in d.groupby('source',sort=True):
        record=source_index[source]
        with np.load(bound(record['readouts']),allow_pickle=False) as f:
            contexts=f['context_code'];counts=f['n_cells'];unit_ok=f['unit_eligible_for_module_analysis'];means=f['means']
            for row in group.to_dict('records'):
                use=(contexts==row['context_code'])&(counts>=spec['minimum_split_cells'])&unit_ok
                j=by_name[row['program_id']];positive=int(np.count_nonzero(means[use,j]>0))
                assert int(use.sum())==int(row['n_units'])
                passed,reason=candidate_status(row,spec['n_splits'],positive,
                    minimum_units=spec['minimum_reliability_units'],minimum_nonzero_units=cfg['minimum_nonzero_units'],
                    minimum_finite_fraction=cfg['minimum_finite_split_fraction'],minimum_median=spec['minimum_split_median_Spearman'],
                    minimum_positive_fraction=spec['minimum_positive_split_fraction'])
                extra.append(dict(row,n_units_observed_positive_score=positive,measurement_candidate=passed,candidate_reason=reason))
    table=R/'outputs/private_state_measurement_candidate_diagnostics_v1.tsv';out=pd.DataFrame(extra);out.to_csv(table,sep='\t',index=False)
    candidates=[]
    for (context,programme),group in out[out.measurement_candidate].groupby(['context_code','program_id'],sort=True):
        sources=sorted(group.source.tolist());assert len(sources)==len(set(sources))
        candidates.append(dict(context_code=int(context),program_id=programme,n_supporting_development_sources=len(sources),
            supporting_sources=sources,scope='multiple_source_measurement_support' if len(sources)>=2 else 'source_limited_measurement_support',
            biological_function_not_validated=True,final_state_feature_model_not_frozen=True))
    op=R/'outputs/accepted_development_state_measurement_candidates_v1.json'
    save(op,dict(schema='accepted_TME_development_state_measurement_candidates_v1',
        readout_acceptance=binding(ip),applicability=m['applicability'],source_roles=m['source_roles'],
        diagnostics=binding(table),n_source_context_programme_diagnostics=len(out),
        n_measurement_candidate_records=int(out.measurement_candidate.sum()),
        context_programme_candidates=candidates,
        minimum_finite_split_fraction=cfg['minimum_finite_split_fraction'],minimum_nonzero_units=cfg['minimum_nonzero_units'],
        development_only=True,score_arrays_or_heldout_outcomes_not_opened=True,
        no_candidate_selected_by_clinical_or_CM_effect=True,
        source_limited_candidates_retained_not_automatically_rejected=True,
        depth_thinning_signature_overlap_and_final_feature_freeze_not_replaced=True))
    storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='DEVELOPMENT_STATE_MEASUREMENT_CANDIDATES_REQUIRE_INFORMATIVE_REPEATED_SPLITS',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(op)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
