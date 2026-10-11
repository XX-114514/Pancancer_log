from pathlib import Path
from datetime import datetime,timezone,timedelta
import os,json,hashlib,shutil
import pandas as pd
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';O=P/'spatial_receiver_programme_exposure_v1';M=P/'spatial_receiver_programme_milestone_v1';S=O/'public_checkpoint_v1'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
assert not S.exists()
reg=read(R/'CM_program_screen_evidence_register_v1.json');assert reg['schema']=='native_LR_program_evidence_register_v16' and reg['active_jobs']==[] and reg['latest_verified_handoff']==str(M/'milestone_handoff.json')
h=read(M/'milestone_handoff.json');v=read(O/'inventory.json');a=read(O/'independent_audit.json');c=read(O/'accepted_primary_consistency_checks.json');conf=read(M/'state_refresh_result.json');assert conf['conflicts']==[]
assert h['summary']['tests']==70 and h['summary']['estimable']==22 and h['summary']['nominal_p05']==h['summary']['primary_BH_q05']==0
S.mkdir();transforms={}
for n in ['all70_primary_programme_exposure.tsv','all75_record_mapping.tsv','all140_radius_sensitivity.tsv','all70_programme_name_sensitivity.tsv','all70_receptor_phenotype.tsv','native_gene_coverage.tsv','independent_all70_checks.tsv','accepted_all70_consistency_checks.tsv']:
    tab=pd.read_csv(O/n,sep='\t',float_precision='round_trip')
    removed=[]
    if 'covariate_names' in tab:
        tab=tab.drop(columns='covariate_names');removed=['covariate_names']
    tab.to_csv(S/n,sep='\t',index=False,na_rep='NA')
    transforms[n]=dict(private_source_sha256=sha(O/n),public_sha256=sha(S/n),rows=len(tab),omitted_columns=removed,reason='Omit exact FOV/author category levels; preserve nuisance column count and fixed schema in checkpoint' if removed else 'Aggregate rows only; no cell/FOV/clinical rows')
for n in ['independent_audit.json','accepted_primary_consistency_checks.json']:shutil.copy2(O/n,S/n)
execution={};timeline=[]
for phase in ['production','audit']:
    x=read(R/('family_LR_spatial_receiver_'+phase+'_v1.execution.json'));assert x['exit_code']==0
    execution[phase]={k:x[k] for k in ['status','started_UTC','finished_UTC','exit_code','script_sha256']}
    for kind in ['started','finished']:
        t=datetime.fromisoformat(x[kind+'_UTC']);timeline.append(dict(event_UTC=t.isoformat(),event_UTC8=t.astimezone(timezone(timedelta(hours=8))).isoformat(),stage=phase,event=kind,status=x['status'],evidence='checkpoint.json#execution.'+phase))
for timestamp,stage,event in [(read(O/'predeclared_plan.json')['created_UTC'],'plan','predeclared'),(c['started_UTC'],'acceptance_checks','started'),(c['finished_UTC'],'acceptance_checks','finished'),(h['verified_utc'],'register_v16','accepted_conditional_receiver_results')]:
    t=datetime.fromisoformat(timestamp);timeline.append(dict(event_UTC=t.isoformat(),event_UTC8=t.astimezone(timezone(timedelta(hours=8))).isoformat(),stage=stage,event=event,status='PREDECLARED' if stage=='plan' else 'COMPLETE',evidence='registration_receipt.json' if stage=='register_v16' else 'accepted_primary_consistency_checks.json' if stage=='acceptance_checks' else '../spatial_receiver_programme_plan_v1.json'))
pd.DataFrame(timeline).sort_values('event_UTC').to_csv(S/'timeline.tsv',sep='\t',index=False)
write(S/'registration_receipt.json',dict(schema='spatial_receiver_registration_receipt_v1',registered_UTC=h['verified_utc'],registered_UTC8=datetime.fromisoformat(h['verified_utc']).astimezone(timezone(timedelta(hours=8))).isoformat(),register=reg['schema'],handoff_sha256=sha(M/'milestone_handoff.json'),state_conflicts=0,active_jobs=0,scientific_release=False,axis_validation_complete=False))
write(S/'privacy_transforms.json',transforms)
cp=dict(schema='spatial_receiver_programme_public_checkpoint_v1',observed_UTC=datetime.now(timezone.utc).isoformat(),execution_status=h['execution_status'],summary=h['summary'],scope=h['scope'],execution=execution,plan_sha256=v['plan_sha256'],producer_inventory_sha256=sha(O/'inventory.json'),audit_sha256=sha(O/'independent_audit.json'),handoff_sha256=sha(M/'milestone_handoff.json'),producer_sha256=v['producer_sha256'],auditor_sha256=a['auditor_sha256'],code_commit='a7d065e97258222ff8255d4dcc00559a0fced613',private_native_input_sha256=v['input_sha256'],nuisance_schema=['intercept','standardized log1p native RNA','standardized log1p native nGenes','standardized logArea','standardized log1p allQC neighbour degree','standardized log1p sender-family degree','FOV dummy categories','original-author-label dummy categories','receiver sender-family-membership dummy category'],scientific_limits=h['scientific_limits'],next_bounded_task=h['next_bounded_task'],scientific_release=False,axis_validation_complete=False)
write(S/'checkpoint.json',cp);shutil.copytree(O/'figures_v1',S/'figures')
write(S/'publication_manifest.json',{str(p.relative_to(S)):sha(p) for p in S.rglob('*') if p.is_file()})
print(json.dumps(dict(status='EXPORTED_AGGREGATES_PENDING_LOCAL_HASH_AND_PIXEL_REVIEW',observed_UTC=cp['observed_UTC'],files=len(list(S.rglob('*'))),summary=h['summary'])),flush=True)
