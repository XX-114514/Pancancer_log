#!/usr/bin/env python3
"""Recover actual LOSO failures, reuse accepted inputs/fits, await original CNV."""
import argparse,hashlib,json,os,re,subprocess,time,traceback
from pathlib import Path
from datetime import datetime,timezone


def now():return datetime.now(timezone.utc).isoformat()


def save(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False)


def main(a):
    r=a.run;c=a.code_dir;s=a.scratch;m=r/'manifests';version=a.attempt_version
    if not re.fullmatch(r'v[0-9]+',version) or not re.fullmatch(r'discovery_v[0-9]+',a.discovery_attempt):raise ValueError('Invalid recovery attempt')
    save(m/('downstream_recovery_'+version+'_launch_private.json'),dict(started_utc=now(),pid=os.getpid(),
        code_git_sha=a.code_git_sha,original_native_CNV_pid=a.CNV_pid,
        source_contract_sha256=hashlib.sha256((r/'config/source_contract_v2_private.json').read_bytes()).hexdigest(),
        code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in c.iterdir() if p.suffix in ['.py','.R']}))
    env=os.environ.copy()
    for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:env[k]='4'
    def stage(name,command):
        save(m/(name+'_recovery_'+version+'_started_private.json'),dict(started_utc=now(),command=command))
        print(json.dumps(dict(phase=name,status='running',utc=now())),flush=True)
        with (r/'logs'/(name+'_recovery_'+version+'.log')).open('xb') as f:
            result=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env=env)
        save(m/(name+'_recovery_'+version+'_terminal_public.json'),dict(finished_utc=now(),exit_code=result.returncode))
        if result.returncode:raise RuntimeError(name+' failed; immutable recovery preserved')
    try:
        recount=m/'CM_input_recount_v1_public.json'
        if not recount.exists():
            stage('CM_input_recount',[a.python,'-u',str(c/'audit_immune_CM_inputs.py'),'--run',str(r),'--scratch',str(s)])
        audit=json.loads((m/'CM_input_recount_v1_public.json').read_text())
        if audit['n_atlas_cells']!=2670057:raise ValueError('Recovery input scope differs')
        if audit['status']!='PASS_SAVED_CM_INPUTS_RECOUNTED_FROM_ALL_CELL_METADATA':raise ValueError('CM recount not accepted')
        for proof in audit['scopes']:
            folder=r/'CM'/proof['scope']/'inputs'
            for name,key in [('frequency_private.tsv','input_frequency_sha256'),('sample_metadata_private.tsv','input_metadata_sha256')]:
                if hashlib.sha256((folder/name).read_bytes()).hexdigest()!=proof[key]:raise ValueError('Audited CM input changed')
        save(m/('CM_input_recount_reuse_'+version+'_public.json'),dict(observed_utc=now(),
            original_recount_sha256=hashlib.sha256(recount.read_bytes()).hexdigest(),all_input_file_checksums_unchanged=True))
        scopes=[]
        for name in ['CM_preparation_public.json','CM_material_preparation_public.json']:
            scopes+=json.loads((m/name).read_text())['scopes']
        reports=[]
        for scope in scopes:
            if not scope['ready']:continue
            name=scope['scope'];folder=r/'CM'/name
            prior=sorted([p for p in folder.glob('discovery_v*') if re.fullmatch(r'discovery_v[0-9]+',p.name)],key=lambda p:int(p.name.split('_v')[-1]),reverse=True)
            completed=next((p for p in prior if (p/'CM_public.json').exists()),None)
            if completed:
                proof=json.loads((completed/'CM_public.json').read_text())
                if proof['scope']!=name or proof['n_nodes']!=scope['n_nodes'] or proof['n_complete_denominator_units']!=scope['n_complete_denominator_units'] or proof['n_composition_blocks']!=scope['n_independent_blocks'] or proof['n_sources']!=scope['n_sources']:raise ValueError('Prior completed CM differs')
                files=['CM_public.json','NMF_discovery.rds','W_private.tsv','H_patient_block_private.tsv','rank_survey.tsv','global_node.tsv','global_edge.tsv','each_node.tsv','each_edge.tsv','CM_members.tsv','LOSO_scope_status.tsv']
                save(m/('CM_'+name+'_reuse_'+version+'_private.json'),dict(observed_utc=now(),attempt=completed.name,
                    unchanged_audited_input=True,file_sha256={f:hashlib.sha256((completed/f).read_bytes()).hexdigest() for f in files}))
                reports.append(dict(proof,discovery_attempt=completed.name))
                print(json.dumps(dict(phase='CM_scope_reused',scope=name,attempt=completed.name,utc=now())),flush=True)
                continue
            command=[a.rscript,str(c/'run_immune_CM.R'),str(folder),str(a.covarnet_R_dir),a.code_git_sha,a.discovery_attempt]
            old=next((p for p in prior if all((p/f).exists() for f in ['rank_survey.tsv','NMF_discovery.rds','W_private.tsv','H_patient_block_private.tsv'])),None)
            if old:command.append(str(old))
            stage('CM_'+name,command)
            proof=json.loads((folder/a.discovery_attempt/'CM_public.json').read_text())
            if proof['scope']!=name or proof['n_nodes']!=scope['n_nodes']:raise ValueError('CM output scope differs')
            reports.append(dict(proof,discovery_attempt=a.discovery_attempt))
        CM_index=m/('CM_recovery_'+version+'_public.json')
        save(CM_index,dict(observed_utc=now(),results=reports,
            all_ready_scopes_completed=len(reports)==sum(x['ready'] for x in scopes),
            actual_inputs_recounted=True,old_attempts_preserved=True,new_CM_is_not_fixed_P3_validation=True))
        gate=m/'native_CNV_v5_terminal_public.json'
        print(json.dumps(dict(phase='waiting_original_native_CNV',utc=now())),flush=True)
        while not gate.exists():
            if not Path('/proc/'+str(a.CNV_pid)).exists():raise RuntimeError('Original CNV stopped without terminal')
            time.sleep(20)
        if json.loads(gate.read_text())['exit_code']!=0:raise RuntimeError('Original CNV failed; library recovery required')
        proof=json.loads((m/'CNV_public.json').read_text())
        if proof['n_cells']!=2670057 or proof['n_libraries']!=400 or proof['n_method_failed']:
            raise ValueError('Original CNV artifact gate differs')
        stage('authority_pack',[a.python,'-u',str(c/'pack_immune_authority.py'),
            '--contract',str(r/'config/source_contract_v2_private.json'),'--run',str(r),'--scratch',str(s),
            '--quota-group',a.quota_group,'--CM-discovery-attempt',a.discovery_attempt,'--CM-discovery-index',str(CM_index)])
        proof=json.loads((m/'immune_authority_release_public.json').read_text())
        if not proof['all_checks_pass'] or not proof['technical_authority_released']:raise RuntimeError('Final archive gate failed')
        save(m/('downstream_recovery_'+version+'_terminal_public.json'),dict(finished_utc=now(),status='completed',
            n_cells=proof['n_cells'],n_features=proof['n_features'],technical_authority_released=True,
            n_CM_scopes=len(reports),independent_clinical_validation=False,fixed_P3_projection_completed=False))
    except Exception as exc:
        (r/'logs'/('downstream_recovery_'+version+'_traceback_private.txt')).write_text(traceback.format_exc())
        save(m/('downstream_recovery_'+version+'_terminal_public.json'),dict(finished_utc=now(),status='failed',failure=str(exc),history_preserved=True))
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['run','scratch','code-dir','covarnet-R-dir']:p.add_argument('--'+name,type=Path,required=True)
    for name in ['python','rscript','code-git-sha','quota-group']:p.add_argument('--'+name,required=True)
    p.add_argument('--CNV-pid',type=int,required=True)
    p.add_argument('--attempt-version',default='v1')
    p.add_argument('--discovery-attempt',default='discovery_v2')
    main(p.parse_args())
