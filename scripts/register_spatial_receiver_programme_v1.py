from pathlib import Path
from datetime import datetime, timezone
import os, json, hashlib, subprocess
import pandas as pd

R = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); P = R / 'immune_same_name_family_source_review_v1'
O = P / 'spatial_receiver_programme_exposure_v1'; M = P / 'spatial_receiver_programme_milestone_v1'
AR = Path(os.environ['PANCANCER_ANNOTATION_ROOT']); ROOT = R / 'cross_source_cell_identity_programme_localization_v1'
def read(p): return json.loads(Path(p).read_text())
def write(p,x): Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()

reg = R / 'CM_program_screen_evidence_register_v1.json'; old = read(reg)
assert old['schema'] == 'native_LR_program_evidence_register_v15' and old['active_jobs'] == [] and not M.exists()
v = read(O / 'inventory.json'); a = read(O / 'independent_audit.json'); c = read(O / 'accepted_primary_consistency_checks.json')
assert a['status'] == 'PASS_ALL70_NATIVE_PROGRAMME_GEOMETRY_EXPOSURE_QR_AND_FULL10000_CONDITIONAL_NULL'
assert c['status'] == 'PASS_ALL70_MAPPING_DENOMINATORS_DIAGNOSTIC_STATUSES_AND_LOFO'
assert a['models'] == c['models'] == 70 and a['cells'] == 254458 and a['mapped_rows'] == c['mapped_rows'] == 75
assert a['native_extra_values_checked'] == 1272290 and a['unique_ligand_draws_replayed'] == 150000
assert sha(O / 'inventory.json') == a['producer_inventory_sha256']
for name,h in v['outputs'].items(): assert sha(O/name) == h
for path,x in old['verified_evidence'].items(): assert sha(path) == x['sha256']
d = pd.read_csv(O / 'all70_primary_programme_exposure.tsv',sep='\t',float_precision='round_trip')
assert len(d) == 70 and not d[['sample_id','programme_interface_id']].duplicated().any()
assert d.programme_interface_id.nunique() == 14 and d.sample_id.nunique() == 5
assert d.loc[d.status != 'ESTIMABLE_CONDITIONAL',['p_two_sided','BH_full70']].isna().all().all()
execution = {}
for phase in ['production','audit']:
    path = R / ('family_LR_spatial_receiver_'+phase+'_v1.execution.json'); j = read(path)
    assert j['status'].startswith('TERMINAL_SUCCESS') and j['exit_code'] == 0
    assert j['script_sha256'] == (v['producer_sha256'] if phase == 'production' else a['auditor_sha256'])
    for role in ['worker_pid','child_pid']:
        check = subprocess.run(['/bin/ps','-p',str(j[role]),'-o','pid=','-o','stat=','-o','etime=','-o','args='],capture_output=True,text=True)
        assert not check.stdout.strip(), 'Inspect any live or reused PID before acceptance'
    execution[str(path)] = dict(manifest=j,all_roles_absent=True)
evidence = dict(old['verified_evidence'])
statuses = {
    'inventory.json':v['status'],'independent_audit.json':a['status'],
    'accepted_primary_consistency_checks.json':c['status'],
    'all70_primary_programme_exposure.tsv':'ACCEPTED_ALL70_FIXED_CONDITIONAL_PROGRAMME_RESULTS',
    'all75_record_mapping.tsv':'ACCEPTED_ALL75_FIXED_RECORD_MAPPING_NO_DUPLICATE_EVIDENCE',
    'all140_radius_sensitivity.tsv':'ACCEPTED_ALL140_PRESPECIFIED_OBSERVED_RADIUS_DIAGNOSTICS',
    'all70_programme_name_sensitivity.tsv':'ACCEPTED_ALL70_PRESPECIFIED_PROGRAMME_MARKER_DIAGNOSTICS',
    'all70_receptor_phenotype.tsv':'ACCEPTED_ALL70_DESCRIPTIVE_RECEPTOR_PHENOTYPES',
    'native_gene_coverage.tsv':'ACCEPTED_ALL35_NEW_GENE_COVERAGE_ROWS',
    'accepted_all70_consistency_checks.tsv':c['status'],
}
for name,status in statuses.items(): evidence[str(O/name)] = dict(sha256=sha(O/name),status=status)
summary = dict(tests=70,mapped_rows=75,programme_interfaces=14,files=5,estimable=int(d.status.eq('ESTIMABLE_CONDITIONAL').sum()),status_counts={str(k):int(n) for k,n in d.status.value_counts().items()},nominal_p05=int((d.p_two_sided<.05).sum()),primary_BH_q05=int((d.BH_full70<.05).sum()),minimal_nominal_p=float(d.p_two_sided.min()),minimal_primary_BH=float(d.BH_full70.min()),all_scope_retained=True)
policy = dict(old['scientific_decision'])
policy['spatial_receiver_programme_exposure_zh'] = '全15固定记录去重复为14定义×5=70主模型，22可估计/48不可估计，10000次L-only条件置换和独立原生/几何/QR复算通过；0未校正p<0.05/0完整主BH显著。固定通用程序未支持配体邻居效应，不能据计算通过宣布轴验证。缺失不是阴性；来源/癌种限定生态价值仍保留，未赋正式状态或改原CM。'
next_task = dict(question='Which existing cell-pair hypotheses can be tested through receptor-appropriate target/TF responses and independent full-gene sources without reselecting current generic programmes to significance?',scope='Start with all15 frozen records and existing originalCM/cellpair evidence; audit prior mechanism and immune-treatment cohort work before choosing a new prespecified validation scope',steps=['Review existing negative mechanism, cell-identity and source-specific branches so failed analyses are not silently repeated','Inventory full-gene native and immune-treatment cohorts, true clinical/patient blocks, baseline overlap and endpoint identity; subset matrices are not full-gene response evidence','Retrieve primary ligand-receptor-target/signalling/TF priors and record exact prior version, orientation, complexes and coverage; generic programme association is not receptor-specific downstream response','Before new testing freeze all eligible targets/TFs, exclusions, within-identity and source/cancer estimands, complete multiple-testing family and independent validation units','Retain cancer/source-limited cell-pair ecology and original9CM membership; current familyextension has0direct originalCMedges and is not aCMdriver','Do not change current radius/programme/endpoints after negative results; report all missing and negative evidence'],acceptance='Independent full-gene/identity/block/prior coverage checks and actual validation results; programme, binding, ecological association and clinical effects distinguished')
t = datetime.now(timezone.utc).isoformat()
h = dict(run_id='20261011_spatial_receiver_programme_milestone_v1',verified_utc=t,biological_question='Within receptor-positive reviewed family identities, do fixed receiver programmes associate with ligand-positive sender neighbourhood exposure after identity/RNA/density/family-count adjustment?',scope=dict(cells=254458,files=5,original_records=15,programme_interfaces=14,models=70,mapped_rows=75,contact_mapping_rows=65,programme_only_mapping_rows=10,permutations=10000),execution_status='COMPLETE_SPATIAL_RECEIVER_PROGRAMME_AUDITED_NO_NOMINAL_OR_PRIMARY_BH_SIGNIFICANT_AXES',outcome='accepted_conditional_receiver_expression_results_not_axis_validation',accepted_baseline=old['baseline'],parent_version=old['latest_verified_handoff'],verified_evidence=evidence,verified_execution=execution,summary=summary,scientific_decision=policy,identity_decision='No changes',malignancy_decision='No changes',state_decision='No formal states assigned; continuous programmes remain separate',scientific_limits=['No nominal or BH-significant current fixed programme-exposure association; this is not proof of absent biological communication','Reused original discovery programme context is not ligand-specific target evidence or independent validation','Missing panel genes, author categories and eligible receiver coverage remain missing, not negative','One public author-mappable spatial source and unreviewed patient/core bridge; FOV/TMA/platform are not independent patients','Area-derived centroid neighbourhood is not membrane contact or absolute calibrated microns','Original subtype/name-marker identities remain distinct from formal states; broad author labels and background diagnostics have limits','Current family LR extension has0direct original9CMedge matches; noCMdriver demonstrated'],active_jobs=[],scientific_release=False,axis_validation_complete=False,next_bounded_task=next_task,files_to_avoid_modifying='CanonicalZarr172nodes9CM/originalglobalBH/alloriginalspatialinputs/allhistoricalandfailedruns')
M.mkdir(); (M/'previous_evidence_register_v15.json').write_bytes(reg.read_bytes()); (M/'previous_execution_checkpoint.json').write_bytes((ROOT/'current_execution_checkpoint.json').read_bytes()); write(M/'milestone_handoff.json',h)
write(reg,dict(schema='native_LR_program_evidence_register_v16',baseline=old['baseline'],latest_verified_handoff=str(M/'milestone_handoff.json'),parent_verified_handoff=old['latest_verified_handoff'],verified_evidence=evidence,active_jobs=[],scientific_decision=policy,scientific_release=False,axis_validation_complete=False,next_bounded_task=next_task))
write(ROOT/'current_execution_checkpoint.json',dict(verified_utc=t,latest_complete_handoff=str(M/'milestone_handoff.json'),verified_execution=execution,active_jobs=[],latest_bounded_stage_complete=True,spatial_receiver_programme_complete=True,scientific_release=False,axis_validation_complete=False))
f = AR/'control/handoffs/20261011_GSE243013_preflight_v1.json'; rt = read(f); (M/'previous_root_control_handoff.json').write_bytes(f.read_bytes()); rt.update(execution_status=h['execution_status'],active_jobs=[],last_verified_utc=t,scientific_release=False,latest_verified_handoff=str(M/'milestone_handoff.json'),latest_verified_evidence=evidence,next_task=next_task); rt.pop('current_execution_handoff',None);write(f,rt);write(AR/'control/handoffs/20261011_spatial_receiver_programme_milestone_v1.json',h)
subprocess.run(['python3',str(AR/'control/refresh_state.py')],check=True);conf = read(AR/'control/STATUS_CONFLICTS.json');write(M/'state_refresh_result.json',conf);assert not conf['conflicts']
print(json.dumps(dict(status=h['execution_status'],register='v16_evidence_not_annotation',registered_UTC=t,conflicts=0,handoff_sha256=sha(M/'milestone_handoff.json'),summary=summary)),flush=True)
