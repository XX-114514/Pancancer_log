#!/usr/bin/env python3
"""Resume archive serialization from accepted CNV and CM, without model refits."""
import argparse
import hashlib
import json
import os
import subprocess
import traceback
from datetime import datetime,timezone
from pathlib import Path


def now():return datetime.now(timezone.utc).isoformat()


def save(path,data):
    with path.open('x') as f:json.dump(data,f,indent=2)


def main(a):
    r=a.run;m=r/'manifests';c=a.code_dir
    save(m/'authority_metadata_recovery_v1_start_private.json',dict(started_utc=now(),pid=os.getpid(),
        code_git_sha=a.code_git_sha,script_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in c.glob('*.py')}))
    try:
        cnv=json.loads((m/'CNV_logspace_v2_public.json').read_text())
        assert cnv['all_cell_IDs_exact'] and cnv['all_sample_IDs_exact'] and cnv['n_method_failed']==0
        assert json.loads((m/'CM_recovery_v3_public.json').read_text())['all_ready_scopes_completed']
        assert not (r/'outputs/immune_validation_authority_v1.zarr').exists()
        commands=[('authority_metadata_recovery_v1',[a.python,'-u',str(c/'pack_immune_authority.py'),
            '--contract',str(r/'config/source_contract_v2_private.json'),'--run',str(r),'--scratch',str(a.scratch),
            '--quota-group',a.quota_group,'--CM-discovery-index',str(m/'CM_recovery_v3_public.json'),
            '--CNV-results',str(r/'tables/all_core_CNV_logspace_v2_private.tsv.gz'),
            '--candidate-name','immune_fullgene_authority_candidate_v2.zarr']),
            ('final_interface_audit_v1',[a.python,'-u',str(c/'audit_immune_final_interfaces.py'),
                '--run',str(r),'--scratch',str(a.scratch)])]
        env=os.environ.copy()
        for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:env[k]='4'
        for name,command in commands:
            save(m/(name+'_started_private.json'),dict(started_utc=now(),command=command))
            with (r/'logs'/(name+'.log')).open('xb') as f:
                result=subprocess.run(command,env=env,stdout=f,stderr=subprocess.STDOUT)
            save(m/(name+'_terminal_public.json'),dict(finished_utc=now(),exit_code=result.returncode))
            if result.returncode:raise RuntimeError(name+' failed; outputs preserved')
        assert json.loads((m/'immune_final_interface_audit_v1_public.json').read_text())['status']=='PASS_FINAL_ARCHIVE_INTERFACES'
        save(m/'authority_metadata_recovery_v1_pipeline_terminal_public.json',dict(finished_utc=now(),status='completed',
            technical_archive_and_interfaces_passed=True,biological_models_not_refit=True))
    except Exception as exc:
        (r/'logs/authority_metadata_recovery_v1_traceback_private.txt').write_text(traceback.format_exc())
        save(m/'authority_metadata_recovery_v1_pipeline_terminal_public.json',dict(finished_utc=now(),status='failed',
            error=str(exc),outputs_preserved=True))
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['run','scratch','code-dir']:p.add_argument('--'+name,type=Path,required=True)
    for name in ['python','code-git-sha','quota-group']:p.add_argument('--'+name,required=True)
    main(p.parse_args())
