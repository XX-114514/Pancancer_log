"""One dependency-aware analysis step, not a recurring scheduler.

Keeps the existing input job. It starts the outcome-blind coverage child only
after a real complete terminal and matching data/producer hashes. No rows or
arrays are written, and no failed/missing dependency is silently restarted.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os,subprocess,sys,time,traceback

PRODUCER='eb142bcdeb172f068f9d182a6ba89866d3d07e852c75f0b976ec4c5b43cd436d'
COMPLETE='COMPLETE_SELECTED_PROCESSED_CELL_MEMBER_MEASUREMENT_QUALIFICATION'


def sha(b):return hashlib.sha256(b).hexdigest()


def ready(accepted):
    terminal=accepted.parent/'private_terminal.json'
    if not accepted.exists() or not terminal.exists():return False
    try:
        raw=accepted.read_bytes();d=json.loads(raw);t=json.loads(terminal.read_bytes())
    except json.JSONDecodeError:return False  # A concurrent small terminal write.
    assert d['status']==t['status']==COMPLETE
    assert d['producer_SHA256']==t['producer_SHA256']==PRODUCER
    assert t['files'][accepted.name]['SHA256']==sha(raw)
    assert d['complete_member_CRC_checked'] and d['complete_independent_CSV_numeric_ROI_and_label_counter_replay']
    assert d['published_metadata_and_directory_replayed'] and d['all_literal_ROI_keys_exact']
    assert d['author_label_counts']['celltype20'] is not None
    return True


def parent_matches(pid,ticks):
    p=Path('/proc')/str(pid)/'stat'
    try:fields=p.read_text().split()
    except FileNotFoundError:return False
    return fields[21]==ticks and fields[2]!='Z'


def main():
    root=Path(__file__).parent;record=root/'private_waiter_status.json'
    assert not record.exists()
    def status(state,**extra):
        value=dict(status=state,observed_UTC=datetime.now(timezone.utc).isoformat(),
            source_SHA256=sha(Path(__file__).read_bytes()),no_individual_rows_or_arrays_written=True,**extra)
        record.write_text(json.dumps(value,indent=2)+'\n');return value
    accepted=Path(os.environ['PANCANCER_IMC_ACCEPTED_MEASUREMENT'])
    child=Path(os.environ['PANCANCER_IMC_COVERAGE_CODE'])
    assert sha(child.read_bytes())==os.environ['PANCANCER_IMC_COVERAGE_SHA256']
    pid=int(os.environ['PANCANCER_IMC_PARENT_PID']);ticks=os.environ['PANCANCER_IMC_PARENT_STARTTICKS']
    status('WAITING_FOR_SAME_INPUT_PRODUCER',parent_pid=pid,parent_startticks=ticks,
        coverage_code_SHA256=sha(child.read_bytes()),no_second_full_member_read_started=True)
    started=time.monotonic()
    while not ready(accepted):
        assert time.monotonic()-started<=90*60,'bounded_dependency_wait_expired'
        assert parent_matches(pid,ticks),'input_producer_missing_or_terminal_without_accepted_measurement'
        assert not (accepted.parent/'private_failed_terminal.json').exists(),'input_dependency_failed'
        time.sleep(5)
    # Child repeats all acceptance gates before its first source request.
    with (root/'private_coverage_stdout.txt').open('wb') as out,(root/'private_coverage_stderr.txt').open('wb') as err:
        process=subprocess.Popen([sys.executable,'-B',str(child)],env=os.environ.copy(),stdout=out,stderr=err)
    status('OUTCOME_BLIND_COVERAGE_CHILD_RUNNING',child_pid=process.pid,
        child_source_SHA256=sha(child.read_bytes()),input_measurement_SHA256=sha(accepted.read_bytes()),
        no_clinical_or_spatial_test_started=True)
    exit_code=process.wait()
    child_terminal=Path(os.environ['PANCANCER_RAM_OUTPUT'])/'private_terminal.json'
    assert exit_code==0 and child_terminal.exists(),'coverage_child_not_complete'
    terminal=json.loads(child_terminal.read_text())
    assert terminal['status']=='COMPLETE_OUTCOME_BLIND_PATIENT_TYPE_MARKER_COVERAGE'
    result=status('COMPLETE_OUTCOME_BLIND_COVERAGE_DEPENDENCY_STEP',child_exit_code=exit_code,
        child_terminal_SHA256=sha(child_terminal.read_bytes()),no_clinical_or_spatial_effects=True)
    print(json.dumps(dict(status=result['status'],receipt_SHA256=sha(record.read_bytes()))))


if __name__=='__main__':
    try:main()
    except Exception as error:
        root=Path(__file__).parent
        value=dict(status='DEPENDENCY_OR_COVERAGE_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],
            no_source_values_in_error=True,no_failed_dependency_restarted=True)
        (root/'private_waiter_failed_terminal.json').write_text(json.dumps(value,indent=2)+'\n')
        print(json.dumps(dict(status=value['status'],error_class=value['error_class'])))
        raise SystemExit(1)
