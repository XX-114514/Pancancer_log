#!/usr/bin/env python3
"""Publish aggregate, time-stamped preprocessing status without patient rows."""
import argparse
from collections import Counter
from datetime import datetime,timezone,timedelta
import hashlib
import json
from pathlib import Path
import os

EXPECTED={
 'GSE294109':('肾透明细胞癌',25,137728),
 'GSE286827':('头颈癌',57,384628),
 'GSE292046':('胶质瘤疫苗：T富集snRNA',8,54156),
 'CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea':('CAR-T输注产品',59,416800),
 'CXG_b252b015_b488_4d5c_b16e_968c13e48a2c':('卵巢癌SPECTRUM',275,927205),
 'CXG_67b6b9ac_b841_42ce_b2d4_ed2ec12327bd':('大B细胞淋巴瘤',232,598289),
 'CXG_843c1459_e852_429c_9339_17b5760e7716':('儿童高级别胶质瘤',63,401253),
 'CXG_34f5307e_7b4d_4a48_b68f_2ba844c6414b':('转移性乳腺癌',64,394534),
 'GSE256326':('肾髓质癌',9,23880)}


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def collect(root):
    clock=datetime.now(timezone.utc);start=root/'manifests/execution_start.json'
    execution=json.loads(start.read_text()) if start.exists() else {}
    pid=execution.get('pid');alive=False
    if pid:
        try:
            os.kill(int(pid),0)
            alive=Path('/proc/%d/stat'%pid).read_text().split()[2]!='Z'
            command=Path('/proc/%d/cmdline'%pid).read_bytes()
            alive=alive and b'run_ambient_preparation.py' in command and str(root).encode() in command
        except (ProcessLookupError,FileNotFoundError):pass
    records=[]
    for dataset,(label,libraries,cells) in EXPECTED.items():
        directory=root/'objects'/dataset
        completed=[]
        for file in sorted(directory.glob('*/summary.json')):
            completed.append(json.loads(file.read_text()))
        success=[r for r in completed if r['status']=='success'];failed=[r for r in completed if r['status']=='failed']
        state='pending'
        if directory.exists():state='running' if alive else 'interrupted_or_terminal_review'
        final=directory/'dataset_summary.json'
        if final.exists():state=json.loads(final.read_text())['status']
        n_cor_selected=sum(r.get('n_output',0) for r in success if r.get('ambient_validation_pass',False))
        records.append(dict(dataset=dataset,name=label,status=state,expected_libraries=libraries,
            published_input_reference_cells=cells,n_libraries_completed=len(completed),
            n_libraries_success=len(success),n_libraries_failed=len(failed),
            n_output_cells=sum(r.get('n_output',0) for r in success),
            n_extra_doublets_removed=sum(r.get('n_additional_removed',0) for r in success),
            n_ambient_model_selected_cells=n_cor_selected,
            n_ambient_model_selected_libraries=sum(bool(r.get('ambient_validation_pass',False)) for r in success),
            n_ambient_unidentifiable_libraries=sum(r.get('ambient_status','').startswith('not_identifiable') for r in success),
            n_correction_review_libraries=sum(r.get('ambient_status','').startswith('DecontX') and not r.get('ambient_validation_pass',False) for r in success),
            n_doublet_cap_review_libraries=sum(bool(r.get('doublet_overcap_review',False)) for r in success),
            n_technical_preparation_eligible_cells=sum(r.get('n_preintegration_cell_eligible',0) for r in success),
            failure_types=dict(Counter(r.get('error','unknown').split('(')[0] for r in failed))))
    terminal=root/'manifests/execution_terminal.json';audit=root/'manifests/independent_preparation_audit.json'
    t=json.loads(terminal.read_text()) if terminal.exists() else None
    a=json.loads(audit.read_text()) if audit.exists() else None
    return dict(schema='ambient_preparation_public_status_v1',observed_utc=clock.isoformat(),
        observed_utc_plus8=clock.astimezone(timezone(timedelta(hours=8))).isoformat(),
        run_id=root.name,started_utc=execution.get('start_time'),code_git_sha=execution.get('code_git_sha'),
        status=t['status'] if t else 'running' if alive else 'not_confirmed_active',
        orchestrator_alive=alive,expected_sources=9,expected_libraries=792,input_scope_cells=3338473,
        input_scope_semantics='4 post-QC sources plus 5 pre-QC native-count sources; not retained/final accepted cells',
        n_completed_libraries=sum(r['n_libraries_completed'] for r in records),
        n_failed_libraries=sum(r['n_libraries_failed'] for r in records),
        n_output_cells=sum(r['n_output_cells'] for r in records),
        n_extra_doublets_removed=sum(r['n_extra_doublets_removed'] for r in records),
        n_technical_preparation_eligible_cells=sum(r['n_technical_preparation_eligible_cells'] for r in records),
        independent_audit_status='passed' if a and a['all_checks_pass'] else 'failed' if a else 'pending_full_run_completion',
        source_final_AC=False,final_accepted_new_cells=0,cross_study_integrated=False,
        records=records,execution_manifest_sha256=digest(start) if start.exists() else None,
        independent_audit_sha256=digest(audit) if audit.exists() else None,
        terminal_sha256=digest(terminal) if terminal.exists() else None)


def render(s):
    lines=['# 新增公共队列：双细胞、RNA污染与整合准备',
        '',f"取证：{s['observed_utc_plus8']}（UTC：{s['observed_utc']}）。运行版本：`{s['run_id']}`。",'',
        f"状态：**{s['status']}**。{s['n_completed_libraries']}/792文库已完成方法步骤，失败{s['n_failed_libraries']}；当前成功分片保留{s['n_output_cells']:,}细胞，额外移除{s['n_extra_doublets_removed']:,}个预测双细胞。",
        '',f"技术准备候选{s['n_technical_preparation_eligible_cells']:,}细胞；全量独立审计：{s['independent_audit_status']}。这不是来源最终AC，最终新增计数仍0，未做跨研究整合。",
        '', '333.8万为4个post-QC与5个pre-QC分支的输入范围，不代表最终通过数。', '',
        '| 队列 | 状态 | 文库完成/总数 | 成功分片细胞 | 额外双细胞移除 | 污染模型采用文库 | 不可辨识文库 | 校正待复核文库 |',
        '|---|---|---:|---:|---:|---:|---:|---:|']
    for r in s['records']:
        lines.append('| %s | %s | %d/%d | %s | %s | %d | %d | %d |'%(r['name'],r['status'],r['n_libraries_completed'],r['expected_libraries'],
            format(r['n_output_cells'],','),format(r['n_extra_doublets_removed'],','),r['n_ambient_model_selected_libraries'],
            r['n_ambient_unidentifiable_libraries'],r['n_correction_review_libraries']))
    lines+=['','## 如何解释','',
        '- 处理逐原始文库执行；原始counts保留，DecontX估计另列counts_decontx。',
        '- 无空液滴且单一谱系不能可靠识别污染，明确保留限制；不能把保留计数当成去污染成功。',
        '- 模型不收敛、marker明显损失、高UMI损失或双细胞超限进入复核。原始细胞不丢弃。',
        '- source_final_AC仍false：来源别名及原始文库与54个已有GSE的重叠门槛未闭合。',
        '- CAR-T产品、T富集snRNA、分选文库、正常组织对照需保持上下文，避免不适当的细胞组成比较。',
        '- TNBC文库来源及CTCL跨图谱来源仍隔离，未进入本792文库队列。', '',
        f"代码冻结：`{s['code_git_sha']}`。",'']
    return '\n'.join(lines)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run-root',required=True);p.add_argument('--output-json',required=True)
    p.add_argument('--output-md',required=True);args=p.parse_args()
    s=collect(Path(args.run_root))
    with open(args.output_json,'x') as f:json.dump(s,f,indent=2,ensure_ascii=False,allow_nan=False)
    with open(args.output_md,'x') as f:f.write(render(s))
    print(json.dumps({k:v for k,v in s.items() if k!='records'},ensure_ascii=False))
