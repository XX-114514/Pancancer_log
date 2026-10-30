"""Bounded queue orchestration fixture; no expression, R or NMF is run."""
import argparse,datetime,hashlib,importlib.util,json,tempfile
from pathlib import Path

def put(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v))
def main():
    assert Path('/proc/self/stat').exists(),'Linux /proc exact process identity is required for this orchestration test'
    module=Path(__file__).resolve().parents[1]/'scripts/run_TME_CM_frozen_fit_queue_v1.py'
    spec=importlib.util.spec_from_file_location('queue_under_test',module)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    base=Path(tempfile.mkdtemp(prefix='TME_synthetic_finite_queue_'))
    joined=base/'ARTIFICIAL_joined';(joined/'scripts').mkdir(parents=True)
    (joined/'scripts/extract_full_source_local_marker_evidence_v1.py').write_text('def storage_gate(root,cfg): pass\n')
    initial=base/'ARTIFICIAL_initial';put(initial/'provenance/supervisor_terminal_v1.json',dict(actual_child_exit_code=0))
    auditor=base/'ARTIFICIAL_auditor';ap=auditor/'outputs/accepted_comparable_CM_inputs_v1.json'
    put(ap,dict(inputs=[dict(panel='ARTIFICIAL_PANEL',scheme='ARTIFICIAL_SCHEME',maximum_rank=4,acceptance=dict(path='ARTIFICIAL_UNUSED',SHA256='ARTIFICIAL_UNUSED'))]))
    put(auditor/'provenance/supervisor_terminal_v1.json',dict(actual_child_exit_code=0))
    put(auditor/'provenance/terminal_v1.json',dict(acceptance_SHA256=hashlib.sha256(ap.read_bytes()).hexdigest()))
    stub='''import argparse,json,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--output');r=Path(p.parse_args().output)
time.sleep(.15)
code=9 if 'FAIL' in r.name else 0
(r/'outputs/numeric_readback_acceptance_v1.json').write_text('{}')
(r/'provenance/supervisor_terminal_v1.json').write_text(json.dumps({'actual_child_exit_code':code}))
'''
    for failure in [False,True]:
        root=base/('ARTIFICIAL_FAILURE' if failure else 'ARTIFICIAL_SUCCESS')
        for n in ['config','scripts','outputs','provenance','checkpoints']:(root/n).mkdir(parents=True)
        for n in ['execute_TME_CM_scheme_v1.py','run_TME_CoVarNet_scheme_v1.R','frozen_author_CoVarNet_v1.R']:(root/'scripts'/n).write_text('ARTIFICIAL_UNUSED')
        (root/'scripts/supervise_atlas_stage_v1.py').write_text(stub)
        targets=base/(root.name+'_children');targets.mkdir()
        requests=[dict(stage_name=('FAIL_' if failure and i==0 else 'ARTIFICIAL_')+str(i),panel='ARTIFICIAL_PANEL',scheme='ARTIFICIAL_SCHEME',rank_mode='fixed_common_rank' if i==0 else 'own_survey_rank') for i in range(4)]
        cfg=dict(joined=str(joined),initial_fits=[str(initial)],input_auditor=str(auditor),requests=requests,runs_parent=str(targets),base_fit_config=dict(threads=1,R_parameters={}),maximum_concurrent_fits=2,fixed_common_rank=4,queue_wall_seconds=40,per_fit_supervisor_wall_seconds=5,frozen_before_module_result_review=True,clinical_outcomes_used=False)
        put(root/'config/spec_v1.json',cfg)
        if failure:
            try:m.main(root)
            except AssertionError as e:assert 'queued fit failed' in str(e)
            else:raise AssertionError('Scientific failure was silently retried/accepted')
            assert len(list(targets.iterdir()))==2
            assert not (root/'provenance/terminal_v1.json').exists()
        else:
            m.main(root);out=json.loads((root/'outputs/completed_frozen_fit_queue_v1.json').read_text());assert len(out['finished'])==4
            events=json.loads((root/'checkpoints/queue_progress_v1.json').read_text())['events'];active=0
            for e in events:
                active+=1 if e['event']=='launch' else -1
                assert 0<=active<=2
            assert active==0
            first=json.loads((targets/'ARTIFICIAL_0/config/spec_v1.json').read_text())
            assert first['R_parameters']['rank_override']==4
    result=dict(status='FINITE_QUEUE_CONCURRENCY_FAILURE_AND_FIXED_RANK_COUNTEREXAMPLES_PASS',synthetic_only=True,fixture=str(base))
    print(json.dumps(result));return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output');args=p.parse_args()
    result=main()
    if args.output:
        root=Path(args.output);put(root/'outputs/synthetic_queue_acceptance_v1.json',result)
        put(root/'provenance/terminal_v1.json',dict(status=result['status'],synthetic_only=True,completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat()))
