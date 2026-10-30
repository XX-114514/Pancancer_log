"""Finite, bounded scientific fit queue for prespecified annotation comparisons.

No model, rank or panel is selected from a clinical result. Existing attempts
are observed but never overwritten or retried automatically after failure.
Each child has its own immutable request, exact process record and supervisor.
"""
import argparse,datetime,hashlib,json,os,subprocess,sys,time
from pathlib import Path

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()

def main(R):
    cfg=load(R/'config/spec_v1.json');sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg);deadline=time.monotonic()+cfg['queue_wall_seconds']
    assert cfg['frozen_before_module_result_review'] and not cfg['clinical_outcomes_used']
    # Deliberately wait for the first panel's actual numeric acceptance before
    # expanding the compute queue; mere process disappearance is insufficient.
    for root in cfg['initial_fits']:
        p=Path(root)/'provenance/supervisor_terminal_v1.json'
        while not p.exists():
            assert time.monotonic()<deadline;time.sleep(5)
        assert load(p)['actual_child_exit_code']==0,'Initial fit failed; investigate before expansion'
    A=Path(cfg['input_auditor']);ap=A/'outputs/accepted_comparable_CM_inputs_v1.json'
    assert load(A/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    assert sha(ap)==load(A/'provenance/terminal_v1.json')['acceptance_SHA256']
    inputs=load(ap)['inputs'];pending=list(cfg['requests']);running=[];finished=[];events=[]
    while pending or running:
        assert time.monotonic()<deadline,'Finite queue wall limit reached; child supervisors remain separately recorded'
        for item in list(running):
            p=Path(item['root'])/'provenance/supervisor_terminal_v1.json'
            if p.exists():
                result=load(p);events.append(dict(UTC=now(),event='terminal',root=item['root'],actual_child_exit_code=result['actual_child_exit_code']))
                assert result['actual_child_exit_code']==0,'A queued fit failed; no automatic scientific retry'
                assert (Path(item['root'])/'outputs/numeric_readback_acceptance_v1.json').exists()
                finished.append(item);running.remove(item)
        while pending and len(running)<cfg['maximum_concurrent_fits']:
            request=pending.pop(0);target=Path(cfg['runs_parent'])/request['stage_name']
            assert not target.exists(),'Existing stage must be inspected, not overwritten'
            entry=next(x for x in inputs if x['panel']==request['panel'] and x['scheme']==request['scheme'])
            c=json.loads(json.dumps(cfg['base_fit_config']));c['input_acceptance']=entry['acceptance']
            if request['rank_mode']=='fixed_common_rank':
                rank=cfg['fixed_common_rank'];assert 2<=rank<=entry['maximum_rank']
                c['R_parameters'].update(minimum_rank=rank,maximum_rank=rank,rank_override=rank)
            else:
                assert request['rank_mode']=='own_survey_rank'
                c['R_parameters']['maximum_rank']=entry['maximum_rank'];c['R_parameters'].pop('rank_override',None)
            os.umask(0o077);target.mkdir()
            for n in ['scripts','config','outputs','checkpoints','provenance']:(target/n).mkdir()
            hashes={}
            for name in ['execute_TME_CM_scheme_v1.py','run_TME_CoVarNet_scheme_v1.R','frozen_author_CoVarNet_v1.R','supervise_atlas_stage_v1.py']:
                original=R/'scripts'/name;destination=target/'scripts'/name
                destination.write_bytes(original.read_bytes());hashes['scripts/'+name]=sha(destination)
            save(target/'config/spec_v1.json',c)
            save(target/'config/supervisor_request_v1.json',dict(entrypoint='execute_TME_CM_scheme_v1.py',arguments=['--output',str(target)],wall_seconds=cfg['per_fit_supervisor_wall_seconds'],script_SHA256=hashes))
            env=os.environ.copy();env.update(OMP_NUM_THREADS=str(c['threads']),OPENBLAS_NUM_THREADS=str(c['threads']),MKL_NUM_THREADS=str(c['threads']),PYTHONDONTWRITEBYTECODE='1')
            with (target/'provenance/supervisor_stdout.txt').open('x') as out,(target/'provenance/supervisor_stderr.txt').open('x') as err:
                proc=subprocess.Popen([sys.executable,'-B',str(target/'scripts/supervise_atlas_stage_v1.py'),'--output',str(target)],stdout=out,stderr=err,env=env,start_new_session=True)
            ticks=(Path('/proc')/str(proc.pid)/'stat').read_text().rsplit(')',1)[1].split()[19]
            item=dict(root=str(target),PID=proc.pid,startticks=ticks,request=request)
            save(target/'provenance/queue_launch_v1.json',dict(item,UTC=now(),queue=str(R)))
            running.append(item);events.append(dict(UTC=now(),event='launch',**item))
            # Event records preserve the exact sequence. Checkpoint is only a
            # live convenience view and does not replace immutable terminals.
        (R/'checkpoints/queue_progress_v1.json').write_text(json.dumps(dict(UTC=now(),n_pending=len(pending),running=running,finished=finished,events=events),indent=2))
        storage_gate(R,cfg)
        if pending or running:time.sleep(5)
    out=R/'outputs/completed_frozen_fit_queue_v1.json'
    save(out,dict(schema='completed_TME_frozen_fit_queue_v1',finished=finished,
        requests=cfg['requests'],clinical_outcomes_used=False,all_actual_terminals_checked=True,
        maximum_concurrent_fits=cfg['maximum_concurrent_fits']))
    save(R/'provenance/terminal_v1.json',dict(status='ALL_PRESPECIFIED_QUEUED_CM_FITS_NUMERICALLY_ACCEPTED',completed_UTC=now(),n_completed_fits=len(finished),acceptance_SHA256=sha(out)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
