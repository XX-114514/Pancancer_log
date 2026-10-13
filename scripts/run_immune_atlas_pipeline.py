#!/usr/bin/env python3
"""Detached, append-only stage driver for the authorized immune atlas run."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import hashlib,json,os,subprocess,time,traceback
from pathlib import Path


def now():return datetime.now(timezone.utc).isoformat()


def save(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False)


def main(args):
    run=args.run;scratch=args.scratch;code=args.code_dir
    lock=run/'manifests'/('production_pipeline_launch_'+args.attempt_version+'_private.json')
    save(lock,dict(started_utc=now(),pid=os.getpid(),code_git_sha=args.code_git_sha,
        script_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in code.iterdir() if p.suffix in ['.py','.R']},
        final_storage='original_project_volume',temporary_storage_is_not_authority=True))
    env=os.environ.copy()
    for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:env[k]='8'
    # Scanpy forwards settings.n_jobs=16 to PyNNDescent, which calls
    # numba.set_num_threads(16). The import-time maximum must permit that.
    env['NUMBA_NUM_THREADS']='16'
    env['MPLBACKEND']='Agg'
    def stage(name,command):
        save(run/'manifests'/(name+'_'+args.attempt_version+'_started_private.json'),dict(phase=name,started_utc=now(),command=command))
        print(json.dumps(dict(phase=name,status='running',utc=now())),flush=True)
        with (run/'logs'/(name+'_'+args.attempt_version+'.log')).open('xb') as f:
            result=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env=env)
        save(run/'manifests'/(name+'_'+args.attempt_version+'_terminal_public.json'),dict(phase=name,finished_utc=now(),exit_code=result.returncode,
             status='process_success_pending_artifact_gate' if result.returncode==0 else 'failed',authoritative_release=False))
        if result.returncode:raise RuntimeError(name+' failed; immutable attempt/log preserved')
        print(json.dumps(dict(phase=name,status='process_success',utc=now())),flush=True)
    try:
        terminal=scratch/'manifests/independent_audit_terminal.json'
        last=None
        while not terminal.exists():
            producer_terminal=scratch/'manifests/QC_DD_ambient_terminal.json'
            if producer_terminal.exists() and json.loads(producer_terminal.read_text())['exit_code']!=0:
                raise RuntimeError('Additional source QC failed before biological integration')
            if not Path('/proc/'+str(args.QC_pid)).exists():raise RuntimeError('Additional QC orchestrator stopped without audited terminal')
            count=len(list((scratch/'objects/prepared/GSE305078').glob('*/summary.json')))
            if count!=last:print(json.dumps(dict(phase='waiting_additional_QC',finished_libraries=count,utc=now())),flush=True);last=count
            time.sleep(20)
        if json.loads(terminal.read_text())['exit_code']!=0:raise RuntimeError('Additional QC independent audit failed')
        contract=run/'config/source_contract_v2_private.json'
        freeze=[args.python,'-u',str(code/'freeze_immune_cohort.py'),'--contract',str(run/'config/source_contract_private.json'),
             '--audit',str(scratch/'manifests/independent_QC_audit.json'),'--prepared',str(scratch/'objects/prepared'),'--scratch',str(scratch),'--run',str(run)]
        if args.resume_from=='integration':
            if not contract.exists():raise RuntimeError('Resume requires frozen source contract')
            frozen_contract=json.loads(contract.read_text())
            audit_hash=hashlib.sha256((scratch/'manifests/independent_QC_audit.json').read_bytes()).hexdigest()
            if frozen_contract['additional_QC_audit_sha256']!=audit_hash:raise RuntimeError('Resume QC audit SHA differs')
            if not args.resume_pilot_report:raise RuntimeError('Resume requires actual accepted pilot report')
            proof=json.loads(args.resume_pilot_report.read_text())
            save(run/'manifests'/('integration_resume_'+args.attempt_version+'_public.json'),dict(
                observed_utc=now(),frozen_contract_sha256=hashlib.sha256(contract.read_bytes()).hexdigest(),
                pilot_report_sha256=hashlib.sha256(args.resume_pilot_report.read_bytes()).hexdigest(),
                numerical_method_unchanged=True,NUMBA_max_threads=16,QC_and_cohort_freeze_not_repeated=True))
        else:
            if args.reuse_native_union:freeze+=['--reuse-native-union',str(args.reuse_native_union)]
            stage('cohort_freeze',freeze)
            pilot=run/('pilot_attempt_'+args.attempt_version);pilot.mkdir(mode=0o700)
            for d in ['manifests','tables']:(pilot/d).mkdir(mode=0o700)
            stage('biology_pilot_v2',[args.python,'-u',str(code/'immune_atlas_biology.py'),'--contract',str(contract),
                 '--run',str(pilot),'--scratch',str(scratch/('pilot_attempt_'+args.attempt_version)),'--panel-dir',str(run/'config'),'--phase','pilot'])
            proof=json.loads((pilot/'manifests/integration_pilot_public.json').read_text())
        if not proof['pilot_only'] or proof['n_shared_features']<2000:raise RuntimeError('Pilot artifact gate failed')
        pair=[m for m in proof['metrics'] if m.get('phase') in ['before','after']]
        if len(pair)==2 and pair[0]['known_marker_major_neighbor_purity'] is not None and pair[1]['known_marker_major_neighbor_purity']<pair[0]['known_marker_major_neighbor_purity']-.10:
            raise RuntimeError('Pilot marker-lineage purity loss exceeds 0.10; batch correction review required')
        common=[args.python,'-u',str(code/'immune_atlas_biology.py'),'--contract',str(contract),'--run',str(run),
                '--scratch',str(scratch),'--panel-dir',str(run/'config')]
        stage('integration',common+['--phase','integrate'])
        proof=json.loads((run/'manifests/integration_integrate_public.json').read_text())
        frozen=json.loads((run/'manifests/cohort_freeze_v2_public.json').read_text())
        if proof['pilot_only'] or proof['n_cells']!=frozen['n_model_cells']:raise RuntimeError('Full integration coverage gate failed')
        pair=[m for m in proof['metrics'] if m.get('phase') in ['before','after']]
        if len(pair)==2 and pair[0]['known_marker_major_neighbor_purity'] is not None and pair[1]['known_marker_major_neighbor_purity']<pair[0]['known_marker_major_neighbor_purity']-.10:
            raise RuntimeError('Full model marker-lineage purity loss exceeds 0.10; batch correction review required')
        stage('annotation',common+['--phase','annotate'])
        stage('fullnative_marker_refinement',[args.python,'-u',str(code/'refine_immune_native_markers.py'),
             '--contract',str(contract),'--run',str(run),'--scratch',str(scratch)])
        def cnv():
            stage('native_CNV',[args.python,'-u',str(code/'immune_native_cnv.py'),'--contract',str(contract),
                  '--run',str(run),'--scratch',str(scratch),'--core-dir',str(args.core_dir),'--gene-annotation',str(args.gene_annotation)])
        def cm():
            stage('CM_preparation',[args.python,'-u',str(code/'prepare_immune_CM.py'),'--run',str(run),'--scratch',str(scratch)])
            report=json.loads((run/'manifests/CM_preparation_public.json').read_text())
            completed=0
            for scope in report['scopes']:
                if scope['ready']:
                    stage('CM_'+scope['scope'],[args.rscript,str(code/'run_immune_CM.R'),str(run/'CM'/scope['scope']),str(args.covarnet_R_dir),args.code_git_sha])
                    if not (run/'CM'/scope['scope']/'discovery_v1/CM_public.json').exists():raise RuntimeError('CM output gate failed')
                    completed+=1
            if not completed:raise RuntimeError('No eligible CM scope reached discovery; composition evidence requires review')
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures=[pool.submit(cnv),pool.submit(cm)]
            for future in futures:future.result()
        stage('authority_pack_and_audit',[args.python,'-u',str(code/'pack_immune_authority.py'),'--contract',str(contract),
             '--run',str(run),'--scratch',str(scratch),'--quota-group',args.quota_group])
        release=json.loads((run/'manifests/immune_authority_release_public.json').read_text())
        if not release['technical_authority_released'] or not release['all_checks_pass']:raise RuntimeError('Final authority gate failed')
        save(run/'manifests'/('production_pipeline_terminal_'+args.attempt_version+'_public.json'),dict(finished_utc=now(),status='completed',
             n_cells=release['n_cells'],n_features=release['n_features'],technical_authority_released=True,
             independent_validation_frozen=False,clinical_eligibility_is_source_specific=True))
        print(json.dumps(dict(utc=now(),phase='pipeline_complete',n_cells=release['n_cells'])),flush=True)
    except Exception as exc:
        (run/'logs'/('production_pipeline_traceback_'+args.attempt_version+'_private.txt')).write_text(traceback.format_exc())
        save(run/'manifests'/('production_pipeline_terminal_'+args.attempt_version+'_public.json'),dict(finished_utc=now(),status='failed',
             failure=str(exc),authoritative_release=False,failed_history_preserved=True))
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['run','scratch','code-dir','core-dir','gene-annotation','covarnet-R-dir']:p.add_argument('--'+key,type=Path,required=True)
    for key in ['python','rscript','code-git-sha','quota-group']:p.add_argument('--'+key,required=True)
    p.add_argument('--QC-pid',type=int,required=True)
    p.add_argument('--attempt-version',default='v1',choices=['v1','v2','v3','v4'])
    p.add_argument('--reuse-native-union',type=Path)
    p.add_argument('--resume-from',choices=['integration'])
    p.add_argument('--resume-pilot-report',type=Path)
    main(p.parse_args())
