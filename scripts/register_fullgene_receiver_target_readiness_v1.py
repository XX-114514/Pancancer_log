"""Accept audited preparation only; never promotes a communication mechanism.

Roots are supplied by PANCANCER_PREFLIGHT_ROOT/PANCANCER_ANNOTATION_ROOT.
All checks precede state mutation; an existing milestone is never overwritten.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import subprocess

import pandas as pd


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8388608), b''):
            digest.update(block)
    return digest.hexdigest()


def write(path, value):
    temporary = Path(str(path) + '.atomic_tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    allow_nan=False) + '\n')
    os.replace(temporary, path)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    annotation = Path(os.environ['PANCANCER_ANNOTATION_ROOT'])
    output = root / 'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1'
    milestone = root / 'immune_same_name_family_source_review_v1/fullgene_receiver_target_milestone_v1'
    register = root / 'CM_program_screen_evidence_register_v1.json'
    checkpoint = root / 'cross_source_cell_identity_programme_localization_v1/current_execution_checkpoint.json'
    control = annotation / 'control/handoffs/20261011_GSE243013_preflight_v1.json'
    secondary = annotation / 'control/handoffs/20261011_fullgene_receiver_target_milestone_v1.json'
    old = read(register)
    assert old['schema'] == 'native_LR_program_evidence_register_v16'
    assert old['active_jobs'] == [] and not milestone.exists() and not secondary.exists()
    plan = read(output / 'predeclared_plan.json')
    assert sha(old['latest_verified_handoff']) == plan['parent_handoff_sha256']
    inventory = read(output / 'inventory.json')
    audit = read(output / 'independent_audit.json')
    freeze = read(root / 'fullgene_receiver_target_code_v1/fullgene_code_freeze_v3.json')
    assert inventory['status'] == 'PRODUCER_COMPLETE_PENDING_FULLGENE_AND_PRIOR_AUDIT'
    assert audit['status'] == 'PASS_ALL336_15_1008_PRIOR_AND_FULLGENE_ORIGINAL_NODE_UNIT_AGGREGATES'
    assert sha(output / 'inventory.json') == audit['producer_inventory_sha256']
    assert inventory['plan_sha256'] == freeze['plan_sha256'] == sha(output / 'predeclared_plan.json')
    assert inventory['producer_sha256'] == freeze['code_shas']['fullgene_receiver_target_readiness_v3.py']
    assert audit['auditor_sha256'] == freeze['code_shas']['audit_fullgene_receiver_target_readiness_v3.py']
    for name, digest in inventory['outputs'].items():
        assert sha(output / name) == digest, name
    for path, value in old['verified_evidence'].items():
        assert sha(path) == value['sha256'], path
    executions = {}
    expected = {
        'assets': ('TERMINAL_SUCCESS_RESOURCE_ONLY', None),
        'production': ('TERMINAL_SUCCESS_PENDING_FULLGENE_AUDIT', inventory['producer_sha256']),
        'audit': ('TERMINAL_SUCCESS_INDEPENDENT_FULLGENE_AUDIT', audit['auditor_sha256']),
    }
    for phase, (status, script_sha) in expected.items():
        suffix = 'v3' if phase == 'assets' else 'v1'
        path = root / ('family_LR_fullgene_' + phase + '_' + suffix + '.execution.json')
        execution = read(path)
        assert execution['status'] == status and execution['exit_code'] == 0
        if script_sha:
            assert execution['script_sha256'] == script_sha
        for key in ['pid', 'child_pid']:
            assert key in execution
            # Fail closed on reused PIDs; /proc is not affected by ps header output.
            assert not Path('/proc/' + str(execution[key])).exists(), (phase, key)
        executions[str(path)] = dict(manifest=execution, worker_and_child_absent=True)
    receipts = read(output / 'verified_prior_inputs_v3/download_receipts.json')
    assert receipts['status'] == 'ALL3_OFFICIAL_MD5_VERIFIED'
    for item in plan['prior']['files']:
        for folder in ['verified_prior_inputs_v2', 'verified_prior_inputs_v3']:
            path = output / folder / item['key']
            assert path.stat().st_size == item['size']
            assert hashlib.md5(path.read_bytes()).hexdigest() == item['checksum'].split(':')[1]
            assert sha(path) == inventory['prior_files'][item['key']]['sha256']
    candidates = pd.read_csv(output / 'all336_prior_coverage.tsv', sep='\t')
    followup = pd.read_csv(output / 'all15_followup_prior_coverage.tsv', sep='\t')
    coverage = pd.read_csv(output / 'all1008_source_target_coverage.tsv', sep='\t')
    sources = pd.read_csv(output / 'fullgene_source_summary.tsv', sep='\t')
    checks = pd.read_csv(output / 'independent_all_source_checks.tsv', sep='\t')
    assert len(candidates) == audit['candidate_records'] == 336
    assert len(followup) == audit['followup_records'] == 15
    assert candidates.model_key.is_unique and followup.model_key.is_unique
    assert set(followup.model_key).issubset(set(candidates.model_key))
    assert len(coverage) == audit['source_comparisons'] == 1008
    assert not coverage[['model_key', 'native_gse_id']].duplicated().any()
    assert len(sources) == len(checks) == 3 and sources.gse_id.is_unique
    assert set(sources.gse_id) == set(plan['scope']['full_native_sources']) == set(checks.gse_id)
    assert int(sources.selected_cells.sum()) == audit['selected_cells'] == 285360
    assert int(sources.libraries.sum()) == 164 and int(sources.source_cells.sum()) == 809784
    assert int(sources.unit_rows.sum()) == int(checks.unit_rows.sum()) == audit['unit_rows']
    assert int(checks.aggregate_gene_values_checked.sum()) == audit['aggregate_gene_values_checked']
    assert inventory['scope']['original_catalogued_nodes'] == 13
    assert inventory['scope']['observed_original_node_union'] == 9
    for value in [inventory, audit]:
        for field in ['statistical_models_fitted', 'formal_states', 'scientific_release', 'axis_validation_complete']:
            assert value[field] is False
    assert inventory['ligand_tf_prior_is_not_a_curated_TF_inventory'] is True
    summary = dict(
        source_summary=sources.to_dict('records'), records=336, unique_LR=209,
        followup_records=15, source_comparisons=1008, selected_cells=285360,
        libraries=164, source_cells=809784, catalogued_original_nodes=13,
        observed_original_node_union=9, unit_rows=audit['unit_rows'],
        aggregate_gene_values_checked=audit['aggregate_gene_values_checked'],
        LR_status_counts=candidates.LR_status.value_counts().to_dict(),
        target_coverage_status_counts=coverage.target_coverage_status.value_counts().to_dict(),
        exact_target_prior_records=int(candidates.exact_ligand_target_column.sum()),
        prior_dimensions=inventory['prior_dimensions'],
        statistical_models_fitted=False, formal_states=False,
        scientific_release=False, axis_validation_complete=False)
    evidence = dict(old['verified_evidence'])
    for name in ['inventory.json', 'independent_audit.json', 'independent_all_source_checks.tsv',
                 'all336_prior_coverage.tsv', 'all15_followup_prior_coverage.tsv',
                 'all1008_source_target_coverage.tsv', 'fullgene_source_summary.tsv']:
        evidence[str(output / name)] = dict(sha256=sha(output / name),
                                           status='ACCEPTED_FULLGENE_PREPARATION_NOT_AXIS_VALIDATION')
    decision = dict(old['scientific_decision'])
    decision['fullgene_receiver_target_readiness_zh'] = (
        '固定336/15/1008及285360细胞全基因原编号单位、官方先验经独立复算验收；'
        '当前3来源观察9/目录13节点，缺捕获和缺先验保持NA。仅输入准备，'
        '未拟合新统计模型、未赋正式状态、未验证轴或原9CM driver。')
    next_task = dict(
        question='Across reviewed cancer/source contexts, which fixed cell-pair hypotheses have reproducible ligand-associated receiver responses with receptor and patient-unit evidence?',
        scope='Predeclare a separate full45 reviewed-source readiness and response scope after actual gene/identity/clinical-unit coverage; retain all336 record mapping and all original13 node identities',
        steps=[
            'Audit native inputs and exact original node/library alignment across all45 reviewed sources; file presence alone is not data acceptance',
            'Retain numbered nodes separately, source/cancer-limited ecological findings and measured zero versus unavailable distinctions',
            'Freeze target-response estimands, full testing family, covariates and patient/sample blocks before new significance tests',
            'Use separate receptor/path-specific priors or signed curated regulons for receptor/TF claims; ligand-target potential alone is not receptor-specific',
            'Identify truly unused baseline/nonoverlapping clinical or spatial validation cohorts; current3 sources were discovery inputs',
            'Retain previous negative fixed spatial results and0direct originalCM edges; do not promote the family extension to aCMdriver'],
        acceptance='Actual independently audited response and validation evidence; coverage or code passing does not complete the scientific goal')
    stamp = datetime.now(timezone.utc).isoformat()
    handoff = dict(
        run_id='20261011_fullgene_receiver_target_milestone_v1', verified_utc=stamp,
        biological_question=plan['question'], scope=inventory['scope'],
        execution_status='COMPLETE_FULLGENE_PRIOR_READINESS_AUDITED_NOT_AXIS_VALIDATION',
        outcome='accepted_input_and_unsigned_prior_coverage_not_new_statistical_or_mechanistic_evidence',
        accepted_baseline=old['baseline'], parent_version=old['latest_verified_handoff'],
        verified_evidence=evidence, verified_execution=executions, summary=summary,
        scientific_decision=decision, identity_decision='Original numbered identities retained',
        malignancy_decision='No changes', state_decision='No formal states assigned',
        scientific_limits=plan['limits'], active_jobs=[], scientific_release=False,
        axis_validation_complete=False, next_bounded_task=next_task,
        files_to_avoid_modifying='CanonicalZarr172nodes9CM and all original/historical/failed inputs and outputs')
    milestone.mkdir()
    for name, path in [('previous_evidence_register_v16.json', register),
                       ('previous_execution_checkpoint.json', checkpoint),
                       ('previous_root_control_handoff.json', control)]:
        (milestone / name).write_bytes(path.read_bytes())
    write(milestone / 'milestone_handoff.json', handoff)
    write(register, dict(schema='native_LR_program_evidence_register_v17', baseline=old['baseline'],
                         latest_verified_handoff=str(milestone / 'milestone_handoff.json'),
                         parent_verified_handoff=old['latest_verified_handoff'], verified_evidence=evidence,
                         active_jobs=[], scientific_decision=decision, scientific_release=False,
                         axis_validation_complete=False, next_bounded_task=next_task))
    write(checkpoint, dict(verified_utc=stamp, latest_complete_handoff=str(milestone / 'milestone_handoff.json'),
                           verified_execution=executions, active_jobs=[], latest_bounded_stage_complete=True,
                           fullgene_preparation_complete=True, scientific_release=False, axis_validation_complete=False))
    root_control = read(control)
    root_control.update(execution_status=handoff['execution_status'], active_jobs=[], last_verified_utc=stamp,
                        scientific_release=False, latest_verified_handoff=str(milestone / 'milestone_handoff.json'),
                        latest_verified_evidence=evidence, next_task=next_task)
    root_control.pop('current_execution_handoff', None)
    write(control, root_control)
    write(secondary, handoff)
    subprocess.run(['python3', str(annotation / 'control/refresh_state.py')], check=True)
    conflicts = read(annotation / 'control/STATUS_CONFLICTS.json')
    write(milestone / 'state_refresh_result.json', conflicts)
    assert not conflicts['conflicts']
    print(json.dumps(dict(status=handoff['execution_status'], register='v17_evidence_not_annotation',
                          registered_UTC=stamp, conflicts=0,
                          handoff_sha256=sha(milestone / 'milestone_handoff.json'), summary=summary)), flush=True)


if __name__ == '__main__':
    main()
