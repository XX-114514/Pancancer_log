"""One bounded existing-audit dependency, then fixed inference and full replay."""
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def now():return datetime.now(timezone.utc).isoformat()


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,16*1024**2)
    code=Path(__file__).parent
    stages=[('test_immune_RNA_MCP_directional_associations_v2.py','55ddfe50f7d4ceb90f15f2e6b42439bff1940bb6dd3a5defda0a15031a2beb2a','immune_RNA_MCP_directional_associations_20261018_v2'),
        ('audit_immune_RNA_MCP_directional_associations_v2.py','6d5af3485a7c383e8096ee88f783c4b22d7c5e0983c8ea3a292798d646cb09aa','immune_RNA_MCP_directional_associations_audit_20261018_v2')]
    for name,digest,_ in stages:assert sha(code/name)==digest
    dependency=root/'immune_RNA_MCP_disjoint_receiver_preflight_audit_20261018_v1/private_terminal.json'
    launch=root/'HNSCC_receiver_followup_job_control_20261016_v2/audit_immune_RNA_MCP_disjoint_receiver_preflight_v1/launch_receipt.json'
    receipt=json.loads(launch.read_text());assert receipt['pid']==288375 and receipt['proc_start_ticks']=='314635735'
    assert receipt['source_sha256']=='8dae2aa90f1f22893622abb75bbc27fd0c68a55000468144d2e2219749620d2a'
    out=root/'immune_RNA_MCP_directional_dependency_pipeline_20261018_v2';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_single_dependency_pipeline.json').write_text(json.dumps(dict(
        created_UTC=now(),source_sha256=sha(__file__),budget=budget,existing_audit_launch_sha256=sha(launch),
        expected_existing_pid=288375,expected_startticks='314635735',dependency_wait_cap_seconds=1800,
        sequence=stages,no_original_audit_or_model_restart=True,no_partial_scientific_effect_review=True),indent=2)+'\n')
    completed=[];deadline=time.monotonic()+1800
    try:
        while not dependency.exists():
            p=Path('/proc')/str(receipt['pid'])/'stat'
            assert p.exists(),'EXISTING_AUDIT_HANDLE_GONE_WITHOUT_TERMINAL'
            fields=p.read_text().rsplit(')',1)[1].split()
            assert fields[19]==receipt['proc_start_ticks'] and fields[0]!='Z','EXISTING_AUDIT_HANDLE_NOT_LIVE'
            assert time.monotonic()<deadline,'DEPENDENCY_WAIT_TIMEOUT_NO_RESTART'
            (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=now(),status='WAITING_SAME_LIVE_INDEPENDENT_PREFLIGHT_AUDIT',
                pid=receipt['pid'],startticks=fields[19],state=fields[0]))+'\n')
            time.sleep(15)
        proof=json.loads(dependency.read_text())
        assert proof['status']=='ACCEPTED_ALL_DISJOINT_READOUTS_NODE_BASELINES_ENDPOINTS_AND_COVERAGE'
        assert proof['source_sha256']==receipt['source_sha256']
        assert proof['producer_terminal_sha256']=='469df1863a6770168a7c67fd78078455b586f654ef7be3e350ede9a7f3592f2b'
        for index,(name,digest,resultdir) in enumerate(stages):
            assert sha(code/name)==digest and not (root/resultdir).exists()
            started=now()
            with (out/('private_stage_%d_stdout_stderr.log'%index)).open('xb',buffering=0) as log:
                child=subprocess.Popen([sys.executable,str(code/name)],env=os.environ.copy(),stdout=log,stderr=subprocess.STDOUT)
                p=Path('/proc')/str(child.pid)/'stat';ticks=p.read_text().rsplit(')',1)[1].split()[19]
                stage=dict(stage=index,source=name,source_sha256=digest,pid=child.pid,startticks=ticks,started_UTC=started,status='RUNNING')
                (out/'private_active_stage.json').write_text(json.dumps(stage,indent=2)+'\n')
                exit_code=child.wait()
            stage.update(completed_UTC=now(),exit_code=exit_code,status='TERMINAL')
            assert exit_code==0,'FIXED_CHILD_FAILED_NO_AUTOMATIC_RETRY'
            terminal=root/resultdir/'private_terminal.json';assert terminal.exists()
            stage['terminal_sha256']=sha(terminal);completed.append(stage)
            (out/('private_stage_%d_receipt.json'%index)).write_text(json.dumps(stage,indent=2)+'\n')
        result=dict(status='DIRECTIONAL_FAMILY_AND_COMPLETE_RANDOM_REPLAY_ACCEPTED_SCIENTIFIC_REVIEW_PENDING',
            completed_UTC=now(),source_sha256=sha(__file__),preflight_acceptance_sha256=sha(dependency),stages=completed,
            original_audit_or_model_restarted=False,clinical_or_causal_validation=False)
    except Exception as e:
        result=dict(status='DIRECTIONAL_DEPENDENCY_PIPELINE_FAILED_NO_AUTOMATIC_RETRY',completed_UTC=now(),
            source_sha256=sha(__file__),stages=completed,error_class=type(e).__name__,private_error_message=str(e))
        (out/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');raise
    (out/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n')
    print('DIRECTIONAL_FIXED_PIPELINE_COMPLETE_SCIENTIFIC_REVIEW_PENDING')


if __name__=='__main__':main()
