from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import os,json,hashlib
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
files={
'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json':'6a3c4311d89a9a02d55a19f2d249c960c508b5973368c06c2e5c43bb605edec9',
'all_source_native_UCell_projection_20261014_v2/inventory.json':'e45f75dbe2ee77c10f1aa44a059e1a24c974ac0011032ddde1980435e36f00c7',
'all_source_native_UCell_audit_aggregation_20261014_v2/inventory.json':'d2a3f6eb66abd8880a89a0a9bc5c36cec462926347122df90d8ad9f45b5e170c',
'original_CM_programme_pair_associations_20261015_v2/inventory.json':'798ffee075584f59dd5fd60f0c70790b0713cb0563dc7d0d95f54d7a5e06c298',
'state_after_DIALOGUE_frozen_projection_20261016_v1/inventory.json':'eed84449facbcdfa49fc64a6438304900e98f3bc5db6ddbd2c9fff18ce458bb5',
'state_DLG_LR_response_20261016_v1/inventory.json':'7565bf43f51fad95fa296bb16ce7996ac360c826557a577054927cde5225a4d1',
'state_DLG_orthogonal_TF_response_20261016_v1/inventory.json':'93ec71c1231e3816fc30dab17b4418174c9b2d4d1e7f0a278dfb1667f58ca335',
'state_DLG_orthogonal_TF_independent_audit_20261016_v1/inventory.json':'f74dd3e22a8a2c438b39db7877d113297b78b063d17a685f11b7696b73dbe99b',
'state_DLG_TF_joint_review_20261016_v1/inventory.json':'81025d5d10337bf59ee288998513a0d0e89d0ddd46129ae942d1a5db8dfa90dd',
'state_spatial_partial_panel_calibration_20261016_v2/inventory.json':'f4802973cf3ca061b35c0659389572ed513ae012d8e5957bf0bfeae684d515fd',
'state_spatial_partial_panel_independent_audit_20261016_v2/inventory.json':'3339daf356719ba15a187be2953f0090852cd692588ebeabca1492cb17c367cf',
'GSE273952_frozen_programme_projection_20261016_v1/inventory.json':'e2d6b55700712bc7896328fa1a5c944ee642ad6a448282f64fd052d924f169ce',
'GSE273952_frozen_programme_independent_audit_20261016_v1/inventory.json':'27d49af833bc87b2f07045a73114b75d92454e47059616d3d33cf6ff9a0622ed',
'GSE273952_spatial_geometry_audit_20261016_v1/inventory.json':'5cc1376909f440fda0e2b306d7eff063382db8bc7a63dbc63af3e335a142650d',
'GSE273952_spatial_RNA_descriptors_20261016_v1/inventory.json':'2c0cc821074f5a1f87c4f61554e6ef2d10daade630cd040ecb75a57235f63d93',
'GSE273952_spatial_RNA_descriptors_audit_20261016_v1/inventory.json':'e369443175d4ebdc7bfb37304b336e3de8a10bdf4d49485d64bd20e43489bdba'}
checks=[]
for rel,expected in files.items():
 p=V/rel;actual=sha(p);assert actual==expected,(rel,actual,expected);j=json.loads(p.read_text());checks.append(dict(relative_evidence=rel,sha256=actual,recorded_status=j.get('status','FROZEN_DICTIONARY'),technical_record_present=True))
d=json.loads((V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json').read_text());assert len(d['programmes'])==814 and not d['formal_states'];assert all(not p['formal_state_name_assigned'] for p in d['programmes']);a=json.loads((V/'GSE273952_frozen_programme_independent_audit_20261016_v1/inventory.json').read_text());b=json.loads((V/'GSE273952_spatial_RNA_descriptors_audit_20261016_v1/inventory.json').read_text());assert a['spots_checked']==71398 and a['manual_score_values']==57332594;assert b['Moran_records_checked']==46569 and b['LR_RNA_records_checked']==171
for name in ['GSE273952_frozen_programme_projection_20261016_v1','GSE273952_frozen_programme_independent_audit_20261016_v1','GSE273952_spatial_neighbourhood_preflight_20261016_v1','GSE273952_spatial_geometry_audit_20261016_v1','GSE273952_spatial_RNA_descriptors_20261016_v1','GSE273952_spatial_RNA_descriptors_audit_20261016_v1','GSE273952_spatial_public_20261016_v1','GSE273952_spatial_RNA_figures_20261016_v1']:
 j=json.loads((V/(name+'.execution.json')).read_text());assert j['exit_code']==0
O=V/'state_DIALOGUE_spatial_attempt_scope_acceptance_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);out=dict(status='PASS_FROZEN_STATE_DISCOVERY_NATIVE_PROJECTION_ORIGINAL_CM_DIALOGUE_LR_TF_AND_ACTUAL_SPATIAL_EXPLORATION_ATTEMPT_SCOPE',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),frozen_programmes=814,method_records=dict(Counter(p['method'] for p in d['programmes'])),evidence=checks,followup_jobs_all_terminal_zero=True,original172_and_20261010_CM9_remain_frozen=True,no_fresh_complete_canonical_matrix_rehash_claim=True,formal_new_state_names=False,validated_cross_cancer_communication_mechanism=False,scientific_unresolved=['exact172 endpoint spatial identity','full18major spatial composition reference','patient identity bridge19libraries','independent same-cancer validation of source-limited RNA candidate','protein receptor assembly and causal perturbation'],meaning='The requested exploratory attempt was executed and audited; negative/weak findings retained. This is not formal biological mechanism acceptance.');(O/'public_scope_acceptance.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
