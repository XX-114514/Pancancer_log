#!/usr/bin/env python3
"""Wait for the authorized run, merge once, independently audit, publish locally."""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def now():return datetime.now(timezone.utc).isoformat()


def save(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False)


def main():
    p=argparse.ArgumentParser();p.add_argument('--preprocess-root',required=True);p.add_argument('--collection',required=True)
    p.add_argument('--final-root',required=True);p.add_argument('--scratch-root',required=True);p.add_argument('--code-git-sha',required=True)
    p.add_argument('--attempt-version',default='')
    args=p.parse_args();root=Path(args.final_root);preprocess=Path(args.preprocess_root);scratch=Path(args.scratch_root)
    if args.attempt_version and not args.attempt_version.replace('_','').isalnum():raise ValueError('Unsafe attempt version')
    suffix=('_'+args.attempt_version) if args.attempt_version else ''
    for name in ['outputs','logs','manifests','scripts','diagnostics']:(root/name).mkdir(exist_ok=True)
    scratch.mkdir(mode=0o700,parents=True,exist_ok=True)
    lock=(root/'logs/assembly.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    components=[Path(__file__),Path(__file__).with_name('merge_public_preprocessed_zarr.py'),Path(__file__).with_name('audit_combined_public_zarr.py'),Path(__file__).with_name('public_gene_identity_contract.py')]
    save(root/'manifests'/('assembly_start'+suffix+'.json'),dict(started_utc=now(),pid=os.getpid(),code_git_sha=args.code_git_sha,
        components_sha256={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in components},
        temporary_volume='data1_personal_project_area',final_volume='data4_existing_Project_v5_area',
        source_final_AC=False,batch_correction_performed=False))
    sequence=0
    while not (preprocess/'manifests/execution_terminal.json').exists():
        start=json.loads((preprocess/'manifests/execution_start.json').read_text());pid=start['pid']
        try:os.kill(pid,0)
        except ProcessLookupError:raise RuntimeError('Preprocessing stopped without final audit/terminal')
        save(root/'manifests'/('assembly_wait'+suffix+'_%04d.json'%sequence),dict(utc=now(),phase='waiting_full_preprocessing_audit'))
        sequence+=1;time.sleep(20)
    terminal=json.loads((preprocess/'manifests/execution_terminal.json').read_text())
    if terminal['status']!='completed':raise RuntimeError('Preprocessing terminal requires repair')
    env=os.environ.copy()
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:env[key]='1'
    commands=[('merge',[sys.executable,'-u',str(components[1]),'--preprocess-root',str(preprocess),
                       '--collection',args.collection,'--final-root',str(root),'--scratch-root',str(scratch),'--code-git-sha',args.code_git_sha]),
              ('audit',[sys.executable,'-u',str(components[2]),'--final-root',str(root),'--report',str(root/'manifests/independent_final_audit.json')])]
    for phase,command in commands:
        with (root/'logs'/(phase+'.log')).open('x') as log:
            child=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT)
            while child.poll() is None:
                save(root/'manifests'/('assembly_phase'+suffix+'_%04d.json'%sequence),dict(utc=now(),phase=phase,pid=child.pid))
                sequence+=1;time.sleep(20)
        if child.returncode:
            save(root/'manifests/assembly_terminal.json',dict(completed_utc=now(),status='needs_repair',phase=phase,exit_code=child.returncode,
                source_final_AC=False,existing_attempt_preserved=True))
            return child.returncode
    final=json.loads((root/'manifests/independent_final_audit.json').read_text())
    if not final['all_checks_pass']:raise RuntimeError('Final audit did not pass')
    final['status']='completed_audited_preprocessed_Zarr_source_acceptance_pending'
    final['code_git_sha']=args.code_git_sha;final['completed_utc']=now()
    save(root/'manifests/assembly_terminal.json',final)
    readable=['# 已下载公共来源：统一前处理Zarr验收','',f"完成：{final['completed_utc']}（UTC）。",'',
        f"{final['n_sources']}来源、{final['n_libraries']}文库、{final['n_cells']:,}保留细胞、{final['n_genes']:,}基因，全量独立验收PASS。",'',
        '最终对象在原项目卷。临时卷用于R交换、分片及指纹工作区；原权威库不改。','',
        'counts是发布源计数；ambient_removed_counts仅保存实际采用的校正差值；分析计数=counts减该差值；X为每细胞1万总量/log1p。',
        '未测基因由measured_in_source标记，不作为真实零表达；技术、组织、富集及产品分支保留。',
        '作者原处理/脑谱系/高损失细胞受到采用门槛；不能称所有污染已清除或所有细胞是确认单胞。',
        '研究/原始文库别名复核仍待，正式新增AC保持0；这是完整技术合集，不是新生物学注释或批次校正模型。',
        '最终Git日志发布需脱敏审阅和普通push；不上传矩阵或患者行。','']
    with (root/'manifests/COMBINED_ZARR_ACCEPTANCE.md').open('x') as f:f.write('\n'.join(readable))
    return 0


if __name__=='__main__':raise SystemExit(main())
