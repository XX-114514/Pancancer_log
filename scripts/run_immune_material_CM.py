#!/usr/bin/env python3
"""Additive material-stratified CM discovery; preserve frozen production code."""
import argparse
from datetime import datetime,timezone
import hashlib,json,os,subprocess,time,traceback
from pathlib import Path


def now():return datetime.now(timezone.utc).isoformat()


def save(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False)


def main(args):
    manifests=args.run/'manifests';code=args.code_dir
    save(manifests/'CM_material_launch_v1_private.json',dict(started_utc=now(),pid=os.getpid(),
        code_git_sha=args.code_git_sha,
        source_contract_sha256=hashlib.sha256((args.run/'config/source_contract_v2_private.json').read_bytes()).hexdigest(),
        script_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in code.iterdir() if p.suffix in ['.py','.R']}))
    env=os.environ.copy()
    for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:env[k]='4'
    def stage(name,command):
        print(json.dumps(dict(phase=name,status='running',utc=now())),flush=True)
        with (args.run/'logs'/(name+'_material_v1.log')).open('xb') as f:
            result=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env=env)
        save(manifests/(name+'_material_v1_terminal_public.json'),dict(phase=name,finished_utc=now(),exit_code=result.returncode))
        if result.returncode:raise RuntimeError(name+' failed; attempt preserved')
    try:
        gate=manifests/'fullnative_marker_refinement_v5_terminal_public.json'
        while not gate.exists():
            terminal=manifests/'production_pipeline_terminal_v5_public.json'
            if terminal.exists():raise RuntimeError('Production stopped before accepted marker output')
            time.sleep(20)
        if json.loads(gate.read_text())['exit_code']!=0:raise RuntimeError('Marker producer failed')
        report=json.loads((manifests/'marker_refinement_public.json').read_text())
        if report['n_cells']!=2670057 or not report['all_cell_native_expression_statistics']:
            raise RuntimeError('Full-cell marker prerequisite differs')
        stage('CM_material_preparation',[args.python,'-u',str(code/'prepare_immune_CM.py'),
            '--run',str(args.run),'--scratch',str(args.scratch),'--scope-mode','material'])
        scopes=json.loads((manifests/'CM_material_preparation_public.json').read_text())['scopes']
        completed=[]
        for scope in scopes:
            if not scope['ready']:continue
            name=scope['scope']
            stage('CM_'+name,[args.rscript,str(code/'run_immune_CM.R'),str(args.run/'CM'/name),str(args.covarnet_R_dir),args.code_git_sha])
            evidence=json.loads((args.run/'CM'/name/'discovery_v1/CM_public.json').read_text())
            if evidence['scope']!=name:raise RuntimeError('CM scope output differs')
            completed.append(evidence)
        save(manifests/'CM_material_terminal_v1_public.json',dict(finished_utc=now(),status='completed',
            n_prepared_scopes=len(scopes),n_completed_scopes=len(completed),results=completed,
            n_insufficient_scopes=sum(not s['ready'] for s in scopes),
            no_existing_atlas_or_CM_overwritten=True,material_and_capture_stratified=True,
            biological_or_clinical_validation=False))
    except Exception as exc:
        (args.run/'logs/CM_material_traceback_v1_private.txt').write_text(traceback.format_exc())
        save(manifests/'CM_material_terminal_v1_public.json',dict(finished_utc=now(),status='failed',failure=str(exc),failed_history_preserved=True))
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['run','scratch','code-dir','covarnet-R-dir']:p.add_argument('--'+name,type=Path,required=True)
    for name in ['python','rscript','code-git-sha']:p.add_argument('--'+name,required=True)
    main(p.parse_args())
