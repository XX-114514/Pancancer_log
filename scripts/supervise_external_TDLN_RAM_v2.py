"""Finite detached source supervisor; only validated aggregate frames on disk.

The worker alone holds native data in RAM. It has a separate process group so
timeouts or explicit termination clean up the worker and R descendants.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import resource
import selectors
import signal
import subprocess
import sys
import traceback


KINDS = {'SOURCE_FIXED_WORKFLOW_STARTED','PUBLIC_INPUT_RAM_BYTES_PROGRESS',
         'COMPLETE_INPUT_IDENTITY_ACCEPTED','COMPLETE_EXTERNAL_FIXED_WORKFLOW_RAM_RELEASED',
         'EXTERNAL_FIXED_WORKFLOW_NOT_ACCEPTED','EXTERNAL_FIXED_WORKFLOW_NATIVE_GATE_FAILURE'}


def now():
    return datetime.now(timezone.utc).isoformat()


def ticks(pid):
    try:
        return int((Path('/proc')/str(pid)/'stat').read_text().split()[21])
    except FileNotFoundError:
        return None


def write_new(path,value):
    with path.open('x') as stream:
        json.dump(value,stream,indent=2)
        stream.write('\n')


def supervise(root):
    os.umask(0o077)
    plan=json.loads((root/'source_plan.json').read_text())
    child=None;actual=None;failure=None;complete=0;frames=0;received=0;phase='VERIFY_GENERIC_CODE'
    output=root/'source_aggregate_frames.jsonl'
    rejected_shape=None
    known_codes={'UNEXPECTED_INDIVIDUAL_OR_CLINICAL_PAYLOAD','UNEXPECTED_IDENTIFIER_VECTOR','UNEXPECTED_LONG_ARRAY','SOURCE_AGGREGATE_LOG_LIMIT','SOURCE_FRAME_LIMIT','UNRECOGNIZED_SOURCE_FRAME','INCOMPLETE_SOURCE_FRAME'}
    def expire(_sig,_frame):
        raise TimeoutError('FINITE_SOURCE_WORKFLOW_BUDGET')
    for sig in (signal.SIGALRM,signal.SIGTERM,signal.SIGHUP):
        signal.signal(sig,expire)
    signal.alarm(plan['maximum_download_seconds']+plan['maximum_R_seconds']+plan.get('supervisor_grace_seconds',120))
    try:
        for name,sha in plan['source_files_SHA256'].items():
            assert Path(name).name==name and name not in ('.','..')
            assert hashlib.sha256((root/name).read_bytes()).hexdigest()==sha
        resource.setrlimit(resource.RLIMIT_AS,(plan['supervisor_address_cap_bytes'],plan['inherited_address_hard_cap_bytes']))
        from accept_external_TDLN_fixed_spatial_v1 import reject_individual_payload
        write_new(root/'source_supervisor_identity.json',dict(started_UTC=now(),
            source_pid=os.getpid(),source_start_ticks=ticks(os.getpid()),
            source_plan_SHA256=hashlib.sha256((root/'source_plan.json').read_bytes()).hexdigest()))
        phase='RUN_FINITE_SOURCE_WORKER'
        child=subprocess.Popen([sys.executable,'-B',str(root/plan['source_worker_file']),str(root)],
            stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,start_new_session=True)
        write_new(root/'source_worker_identity.json',dict(started_UTC=now(),
            source_pid=child.pid,source_start_ticks=ticks(child.pid),source_worker_file=plan['source_worker_file']))
        with output.open('x') as stream:
            pending=bytearray();group_cleaned=False
            with selectors.DefaultSelector() as reader:
                reader.register(child.stdout,selectors.EVENT_READ)
                while True:
                    if child.poll() is not None and not group_cleaned:
                        # Descendants can inherit stdout after the worker exits.
                        # Kill the dedicated group before waiting for pipe EOF.
                        try:
                            os.killpg(child.pid,signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                        group_cleaned=True
                    if not reader.select(timeout=.25):
                        continue
                    chunk=os.read(child.stdout.fileno(),65536)
                    if not chunk:
                        assert not pending, 'INCOMPLETE_SOURCE_FRAME'
                        break
                    received+=len(chunk);assert received<=2*1024**2, 'SOURCE_AGGREGATE_LOG_LIMIT'
                    pending.extend(chunk)
                    while b'\n' in pending:
                        position=pending.index(b'\n')+1
                        assert position<=1024**2, 'SOURCE_FRAME_LIMIT'
                        line=bytes(pending[:position]);del pending[:position]
                        frame=json.loads(line);assert frame['kind'] in KINDS, 'UNRECOGNIZED_SOURCE_FRAME'
                        try:
                            reject_individual_payload(frame)
                        except AssertionError:
                            # Only constant field roles and types: never values,
                            # arbitrary dictionary keys, rejected text or payload.
                            issues=[]
                            def shape(value):
                                if isinstance(value,dict):
                                    forbidden={'internal','case_values','case_matrix','patient_ids','sample_ids','cell_ids','Patient.ID','ROI.ID','Image.name','sample_id','ObjectNumber','pg_clusters_corrected_40','coordinates','neighbors','bootstrap_indices','Group','Group.location','RFS.time'}
                                    for key in sorted(forbidden.intersection(value)):
                                        issues.append(dict(rule='forbidden_key',field_role=key))
                                    for key in ('fine','image','case','cell_rows'):
                                        if key in value and value[key] is not None and not (type(value[key]) is int and value[key]>=0):
                                            issues.append(dict(rule='noninteger_aggregate_count',field_role=key,value_type=type(value[key]).__name__))
                                    for child in value.values():shape(child)
                                elif isinstance(value,list):
                                    if len(value)>128:issues.append(dict(rule='long_array',length=len(value)))
                                    for child in value:shape(child)
                            shape(frame)
                            rejected_shape=issues[:32]
                            raise
                        if frame['kind']=='COMPLETE_EXTERNAL_FIXED_WORKFLOW_RAM_RELEASED':
                            complete+=1;assert complete==1
                        stream.write(json.dumps(frame)+'\n');stream.flush();frames+=1
                    assert len(pending)<=1024**2, 'SOURCE_FRAME_LIMIT'
        actual=child.wait(timeout=30)
    except Exception as error:
        locations=[dict(module=Path(v.filename).name,line=v.lineno) for v in traceback.extract_tb(error.__traceback__)]
        failure=dict(phase=phase,error_class=type(error).__name__,
            error_code=str(error) if str(error) in known_codes else 'REDACTED_OR_NO_LITERAL_CODE',
            generic_code_locations=locations[-6:],rejected_aggregate_shape=rejected_shape)
    finally:
        signal.alarm(0)
        if child is not None:
            if child.poll() is None:
                try:
                    os.killpg(child.pid,signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    child.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(child.pid,signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    child.wait(timeout=10)
            actual=child.returncode
            # A worker may exit before a descendant. Its dedicated group is
            # always terminated, including on apparent normal completion.
            try:
                os.killpg(child.pid,signal.SIGKILL)
            except ProcessLookupError:
                pass
    effective=actual if failure is None and actual not in (None,0) else (0 if failure is None and actual==0 and complete==1 else 2)
    result=dict(status='SOURCE_SUPERVISOR_ACTUAL_WORKER_EXIT_COLLECTED',observed_UTC=now(),
        actual_worker_exit_code=actual,effective_exit_code=effective,complete_source_frames=complete,
        persisted_validated_aggregate_frames=frames,failure=failure,
        aggregate_frames_SHA256=hashlib.sha256(output.read_bytes()).hexdigest() if output.exists() else None,
        source_plan_SHA256=hashlib.sha256((root/'source_plan.json').read_bytes()).hexdigest(),
        native_individual_coordinate_state_graph_and_resample_files_not_created=True,
        worker_and_R_descendants_waited_or_terminated=True)
    write_new(root/'source_actual_exit.json',result)
    return effective


if __name__=='__main__':
    raise SystemExit(supervise(Path(sys.argv[1])))
