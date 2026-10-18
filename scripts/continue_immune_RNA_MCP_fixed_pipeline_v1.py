"""Single bounded dependency pipeline for the already-running fixed family."""
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def utc():
    return datetime.now(timezone.utc).isoformat()


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,96*1024**2)
    models=root/'immune_RNA_MCP_official_models_20261018_v1'
    original_driver=Path(__file__).with_name('run_immune_RNA_MCP_official_family_v1.py')
    assert sha(original_driver)=='24ed7c348a96100a6b00d63a2809d1e7f39ba10f842253320280d515698090bc'
    original_receipt=root/'HNSCC_receiver_followup_job_control_20261016_v2/run_immune_RNA_MCP_official_family_v1_launchfix2/launch_receipt.json'
    receipt=json.loads(original_receipt.read_text())
    assert receipt['pid']==392166 and str(receipt['proc_start_ticks'])=='314029788'
    stages=[
        ('project_immune_RNA_MCP_official_family_v1.py','fcb603ababa242460b0a780718807a5232d907f13a1c257d015a6c85483014e9','immune_RNA_MCP_frozen_projection_20261018_v1',7200),
        ('audit_immune_RNA_MCP_frozen_projection_v1.py','297fe7bb6b3fc136bf9d453ad5f91f303f5c6c2d27b31539803d3cae903688b1','immune_RNA_MCP_frozen_projection_audit_20261018_v1',7200),
        ('run_immune_RNA_MCP_accepted_postfit_v1.py','336d44c01c49b363ae671186a1aa63c554a06e12caecfbfe99e6ae5276017d2c','immune_RNA_MCP_accepted_postfit_execution_20261018_v1',18000),
        ('audit_immune_RNA_MCP_stability_and_coordination_v1.py','705f2ffb34748aefebce46f337f6ec4b566320e8e05455f638f0a928cc409c1c','immune_RNA_MCP_stability_coordination_audit_20261018_v1',14400)]
    for name,digest,destination,_ in stages:
        assert sha(Path(__file__).with_name(name))==digest and not (root/destination).exists()
    rscript=Path(os.environ.get('PANCANCER_DIALOGUE_RSCRIPT',
        str(Path(sys.executable).parents[2]/'scrna_r/bin/Rscript')))
    assert rscript.is_file()
    os.environ['PANCANCER_DIALOGUE_RSCRIPT']=str(rscript)
    output=root/'immune_RNA_MCP_fixed_dependency_pipeline_20261018_v1'
    assert not output.exists();output.mkdir(mode=0o700)
    plan=dict(created_UTC=utc(),source_sha256=sha(__file__),budget=budget,
        upstream_driver_pid=receipt['pid'],upstream_start_ticks=receipt['proc_start_ticks'],
        original_launch_receipt_sha256=sha(original_receipt),stages=stages,
        wait='Single existing upstream family, at most31hours (five original6h limits+1h margin). Verify PID startticks while waiting; missing/reused/nonlive upstream without terminal stops. Never restart any model.',
        gates='All5 original fit attempts and manifest hashes, frozen contract and numeric gates are required. Execute projection, independent score replay, full fixed stability/source tests and independent full random-stream replay in order. Errors preserved, no automatic scientific retry.',
        not_recurring_or_scheduled=True,no_new_download_or_clinical_endpoint=True)
    (output/'plan_before_dependency_pipeline.json').write_text(json.dumps(plan,indent=2)+'\n')
    status=dict(status='WAITING_FOR_CURRENT_UPSTREAM_FAMILY',started_UTC=utc(),finished_stages=[])
    status_path=output/'private_current_stage.json'
    try:
        deadline=time.monotonic()+31*3600
        while not (models/'private_terminal.json').exists():
            proc=Path('/proc')/str(receipt['pid'])/'stat'
            try:
                fields=proc.read_text().rsplit(')',1)[1].split()
            except FileNotFoundError:
                if (models/'private_terminal.json').exists():break
                raise AssertionError('UPSTREAM_MISSING_WITHOUT_TERMINAL_NO_RESTART')
            if fields[19]!=str(receipt['proc_start_ticks']) or fields[0] in ['Z','X']:
                if (models/'private_terminal.json').exists():break
                raise AssertionError('UPSTREAM_NONLIVE_OR_REUSED_NO_RESTART')
            assert time.monotonic()<deadline, 'BOUNDED_DEPENDENCY_WAIT_EXPIRED_NO_RESTART'
            status.update(observed_UTC=utc(),upstream_state=fields[0])
            status_path.write_text(json.dumps(status,indent=2)+'\n')
            time.sleep(15)
        terminal=None
        for attempt in range(3):
            try:
                terminal=json.loads((models/'private_terminal.json').read_text());break
            except json.JSONDecodeError:
                time.sleep(1)
        assert terminal is not None, 'UPSTREAM_TERMINAL_INCOMPLETE_NO_RESTART'
        assert terminal['status']=='ALL_FIVE_PREDECLARED_OFFICIAL_RNA_MCP_FITS_ATTEMPTED_VALIDATION_PENDING'
        assert terminal['source_sha256']==sha(original_driver)
        fits=json.loads((models/'private_all_planned_fit_status.json').read_text())
        assert len(fits)==5 and all('exit_code' in s for s in fits)
        for name,item in terminal['files'].items():assert sha(models/name)==item['sha256']
        for name,digest,destination,timeout in stages:
            resource(root,80*1024**2)
            assert sha(Path(__file__).with_name(name))==digest and not (root/destination).exists()
            stage=dict(module=name,source_sha256=digest,started_UTC=utc())
            status.update(status='EXECUTING_FIXED_STAGE',current_module=name,observed_UTC=utc())
            status_path.write_text(json.dumps(status,indent=2)+'\n')
            with (output/(name+'.log')).open('xb') as log:
                result=subprocess.run([sys.executable,str(Path(__file__).with_name(name))],
                    stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
            proof=root/destination/'private_terminal.json'
            stage.update(exit_code=result.returncode,finished_UTC=utc(),
                terminal_sha256=sha(proof) if proof.exists() else None)
            status['finished_stages'].append(stage)
            status_path.write_text(json.dumps(status,indent=2)+'\n')
            assert result.returncode==0 and proof.exists(), 'DOWNSTREAM_TECHNICAL_FAILURE_PRESERVED_NO_AUTOMATIC_RETRY'
        status.update(status='FIXED_PIPELINE_AND_INDEPENDENT_NUMERIC_AUDIT_COMPLETE_SCIENTIFIC_REVIEW_PENDING',finished_UTC=utc())
    except Exception as error:
        status.update(status='STOPPED_REVIEW_REQUIRED_ALL_ORIGINAL_OUTPUTS_PRESERVED',
            error_type=type(error).__name__,message=str(error),finished_UTC=utc())
        raise
    finally:
        status_path.write_text(json.dumps(status,indent=2)+'\n')
        if 'finished_UTC' in status:
            record=dict(status,source_sha256=sha(__file__),biological_clinical_or_direct_interaction_validation=False)
            (output/'private_terminal.json').write_text(json.dumps(record,indent=2)+'\n')
    print('FIXED_PIPELINE_NUMERICALLY_COMPLETE_SCIENTIFIC_REVIEW_PENDING')


if __name__=='__main__':main()
