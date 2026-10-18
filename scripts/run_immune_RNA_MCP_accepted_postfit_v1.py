"""Validate every frozen input and acceptance before stability/source inference."""
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 80 * 1024**2)
    contract = root / 'immune_RNA_MCP_postfit_contract_20261018_v1/plan_before_learned_outputs_or_projection.json'
    assert sha(contract) == '1172d887cffb7c3b0c28a21583154fbe3718a1f13d79edb85500343ead5e2b09'
    scope = root / 'immune_RNA_MCP_discovery_scope_20261018_v1'
    assert sha(scope / 'private_terminal.json') == '0f4833a07bd5cdf5e8f6806d06805777dddf78adad60a6151eca991dfdda6196'
    for name, item in json.loads((scope / 'private_terminal.json').read_text())['files'].items():
        assert sha(scope / name) == item['sha256']
    inputs = root / 'immune_RNA_MCP_native_expression_inputs_20261018_v1'
    assert sha(inputs / 'private_terminal.json') == '2474355a003d96232b752ad2c77ba92a506f4c68829c7e0adf434045cf03d656'
    for name, item in json.loads((inputs / 'private_terminal.json').read_text())['files'].items():
        assert sha(inputs / name) == item['sha256']
    projection = root / 'immune_RNA_MCP_frozen_projection_20261018_v1'
    audit = root / 'immune_RNA_MCP_frozen_projection_audit_20261018_v1/private_terminal.json'
    proof = json.loads(audit.read_text())
    assert proof['status'] == 'ALL_PLANNED_SLOTS_AND_ALL_SUPPORTED_PROJECTION_VALUES_ACCEPTED'
    assert sha(projection / 'private_terminal.json') == proof['producer_terminal_sha256']
    for name, digest in json.loads((projection / 'private_terminal.json').read_text())['files'].items():
        assert sha(projection / name) == digest
    plan = json.loads((projection / 'plan_before_official_output_read.json').read_text())
    models = root / 'immune_RNA_MCP_official_models_20261018_v1'
    assert sha(models / 'private_terminal.json') == plan['original_model_terminal_sha256']
    for name, digest in plan['original_model_files_sha256'].items():
        assert sha(models / name) == digest
    math_path = root / 'immune_RNA_MCP_inference_math_check_20261018_v1/private_terminal.json'
    math_proof = json.loads(math_path.read_text())
    assert math_proof['status'] == 'PASS_SYNTHETIC_SIGNED_UNAVAILABLE_ALL_RANDOM_REPLAY'
    sources = {'review_immune_RNA_MCP_seed_stability_v1.py': '141a8ddffe10c8eb0808b8359f0a95cb11ea62dab11c5adfcd818f713e3b60d3',
        'test_immune_RNA_MCP_source_coordination_v1.py': '4bff1ebe975983237035de8b401dc9d9e13ff38c321276643fa73c2caaaf37ff'}
    for name, digest in sources.items():
        assert sha(Path(__file__).with_name(name)) == digest == math_proof['production_sources'][name]
    output = root / 'immune_RNA_MCP_accepted_postfit_execution_20261018_v1'
    assert not output.exists()
    destinations = ['immune_RNA_MCP_seed_stability_20261018_v1', 'immune_RNA_MCP_source_coordination_20261018_v1']
    assert not any((root / n).exists() for n in destinations), 'Preserve old outputs; no implicit rerun.'
    output.mkdir(mode=0o700)
    (output / 'plan_before_accepted_postfit_execution.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),budget=budget,
        source_sha256_by_module=sources,projection_acceptance_sha256=sha(audit),
        synthetic_math_acceptance_sha256=sha(math_path),contract_sha256=sha(contract),
        gates='Current exact file hashes for all frozen case/source fields, native inputs, projection outputs and original RDS. Execute every fixed stability and source slot; no choosing by result. Independent real-data random-stream audit remains required.'),indent=2)+'\n')
    statuses=[]
    for (name,digest),destination in zip(sources.items(),destinations):
        resource(root,72*1024**2)
        status=dict(module=name,source_sha256=digest,started_UTC=datetime.now(timezone.utc).isoformat())
        with (output/(name+'.log')).open('xb') as log:
            code=subprocess.run([sys.executable,str(Path(__file__).with_name(name))],stdout=log,
                stderr=subprocess.STDOUT,timeout=7200).returncode
        status.update(exit_code=code,finished_UTC=datetime.now(timezone.utc).isoformat())
        terminal=root/destination/'private_terminal.json'
        status['terminal_sha256']=sha(terminal) if terminal.exists() else None
        statuses.append(status)
        (output/'private_stage_status.json').write_text(json.dumps(statuses,indent=2)+'\n')
        assert code==0 and terminal.exists(), 'Technical failure preserved; do not replace parameters or restart model.'
    (output/'private_terminal.json').write_text(json.dumps(dict(
        status='FIXED_STABILITY_AND_SOURCE_INFERENCE_DONE_REAL_NUMERIC_ACCEPTANCE_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),stages=statuses,
        independent_biological_clinical_or_direct_interaction_validation=False),indent=2)+'\n')
    print('FIXED_POSTFIT_STAGES_DONE_INDEPENDENT_REAL_DATA_AUDIT_PENDING')


if __name__=='__main__':
    main()
