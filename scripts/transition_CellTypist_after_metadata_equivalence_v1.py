"""Bounded, explicit owned transition after numerical-equivalence acceptance.

Healthy execution is interrupted only for this recorded performance change,
after ordinary/Unicode official inference agrees and a complete-source
checkpoint is available. Other failures are never treated as this transition.
"""
import argparse,datetime,hashlib,json,os,runpy,signal,statistics,time
from pathlib import Path

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2)
def exact_owned(pid,ticks):
    p=Path('/proc')/str(pid)/'stat'
    return p.exists() and p.read_text().rsplit(')',1)[1].split()[19]==ticks
def main(R):
    cfg=load(R/'config/spec_v1.json');E=cfg['equivalence_acceptance'];ep=Path(E['path']);assert sha(ep)==E['SHA256'];e=load(ep);assert load(ep.parent.parent/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    assert e['expression_model_numbers_and_official_functions_unchanged'] and e['no_global_numpy_or_package_monkeypatch'];probes=e['probes'];assert len(probes)>=12
    assert all(p['decision_and_probability_bitwise_identical'] and p['labels_identical'] and p['original_model_unchanged'] for p in probes)
    speedup=statistics.median(p['ordinary_seconds']/p['unicode_metadata_seconds'] for p in probes);assert speedup>=cfg['minimum_median_probe_speedup'],'Performance transition not justified'
    P=Path(cfg['previous']);terminal=P/'provenance/supervisor_terminal_v1.json';child=load(P/'provenance/child_launch_v1.json');request=load(P/'config/supervisor_request_v1.json');entry='scripts/'+request['entrypoint'];assert sha(P/entry)==request['script_SHA256'][entry]==cfg['previous_frozen_script_SHA256'];deadline=time.monotonic()+cfg['wait_seconds']
    while not terminal.exists():
        assert exact_owned(child['PID'],child['startticks']),'Previous owned child disappeared without a terminal'
        completed=[]
        for p in sorted((P/'outputs').glob('private_source_coverage_*_v1.json')):
            try:
                value=load(p)
                assert isinstance(value['source'],str) and value['n_fit']>0 and 'reference_predictions' in value
            except (json.JSONDecodeError,KeyError,AssertionError):continue
            completed.append(p)
        if len(completed)>=cfg['minimum_completed_sources']:break
        assert time.monotonic()<deadline,'Complete-source checkpoint wait expired';time.sleep(5)
    stopped=False
    if not terminal.exists():
        transition=dict(reason='metadata_only_gene_match_optimization_after_bitwise_API_equivalence',requested_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),controller_stage=str(R),old_child_PID=child['PID'],old_child_startticks=child['startticks'],equivalence_acceptance_SHA256=sha(ep),n_completed_source_checkpoints=len(completed),complete_checkpoint_bindings=[dict(path=str(p),SHA256=sha(p)) for p in completed],partial_source_will_not_be_reused=True,median_probe_speedup=speedup,scientific_models_input_normalization_scope_and_batch_size_unchanged=True)
        save(P/'provenance/requested_unicode_metadata_transition_v1.json',transition);assert exact_owned(child['PID'],child['startticks']);os.kill(child['PID'],signal.SIGTERM);stopped=True
        stop_deadline=time.monotonic()+60
        while not terminal.exists():assert time.monotonic()<stop_deadline,'Owned child terminal not observed';time.sleep(1)
    t=load(terminal)
    if not stopped and t['actual_child_exit_code']==0:
        save(R/'outputs/automatic_execution_pointer_v1.json',dict(automatic=str(P),auditor=cfg['previous_auditor'],resume_needed=False));status='PREVIOUS_CLASSIFIER_COMPLETED_TRANSITION_NOT_NEEDED'
    else:
        assert stopped and t['actual_child_exit_code']==-15 and not t['wall_timeout'];assert not exact_owned(child['PID'],child['startticks'])
        for key in ['resume_launcher','auditor_launcher']:
            b=cfg[key];p=R/'scripts'/b['file'];assert sha(p)==b['SHA256'];runpy.run_path(str(p),run_name='__main__')
        save(R/'outputs/automatic_execution_pointer_v1.json',dict(automatic=cfg['new_automatic'],auditor=cfg['new_auditor'],resume_needed=True,reason='explicit_bitwise_validated_gene_metadata_performance_transition',old_checkpoint_outputs_retained=True,partial_source_not_reused=True,CM_input_gates_still_separate=True));status='VERIFIED_METADATA_TRANSITION_COMPLETE_SOURCE_RESUME_AND_AUDIT_LAUNCHED'
    save(R/'provenance/terminal_v1.json',dict(status=status,completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),owned_child_intentionally_stopped=stopped,median_probe_speedup=speedup,old_supervisor_terminal_SHA256=sha(terminal)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
