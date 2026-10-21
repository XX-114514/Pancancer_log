"""Review fixed-case selection using accepted SOURCE aggregates only.

No individual inputs, new association tests, missingness weights or refits.
Different retention by outcome is a generalizability limitation, not proof
that any observed adjusted state association is an artifact.
"""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json


def main():
    p=argparse.ArgumentParser();p.add_argument('--task',type=Path,required=True);args=p.parse_args()
    task=args.task
    paths={
        'baseline':task/'private_complete_existing_clinical_metadata_v1.json',
        'capture':task/'native_measurement_v1/private_complete_source_coverage.json',
        'proof':task/'native_measurement_v1/private_source_input_proof.json',
        'main':task/'state_statistics_v1/private_main_same_case_coverage.json',
    }
    data={k:json.loads(v.read_text()) for k,v in paths.items()}
    assert data['baseline']['clinical_table_SHA256']==data['proof']['clinical_table_SHA256']
    baseline={r['cancer_type']:r for r in data['baseline']['all_source_scope']
              if r['role']=='EXPLICIT_ICI_AUTHOR_CASE'}
    capture={r['cancer_type']:r for r in data['capture']
             if r['QC']=='native_QC500_250_MT20' and r['parent']=='dual_regulatory_evidence'}
    assert set(baseline)=={'LUAD','LUSC'} and set(capture)=={'ALL','LUAD','LUSC'}
    rows=[]
    for cancer in ['ALL','LUAD','LUSC']:
        original=list(baseline.values()) if cancer=='ALL' else [baseline[cancer]]
        n1=sum(x['eligible_MPR_pCR'] for x in original);n0=sum(x['eligible_nonMPR'] for x in original)
        assert n1+n0==sum(x['new_state_clinical_metadata_eligible'] for x in original)
        c=capture[cancer]
        assert c['clinical_metadata_units']==n1+n0
        for group,denominator,retained in [
            ('MPR_pCR',n1,c['MPR_pCR_among_both']),('nonMPR',n0,c['nonMPR_among_both'])]:
            assert 0<=retained<=denominator and denominator>0
            rows.append(dict(source='GSE243013',cancer=cancer,stage='BOTH_PARENTS_GE20',outcome=group,
                original_clinical_eligible_units=denominator,retained_units=retained,
                excluded_units=denominator-retained,retained_fraction=retained/denominator))
    assert data['main']['clinical_and_both_parent_capture']==capture['ALL']['both_parents_atleast20']
    for group,field in [('MPR_pCR','response1'),('nonMPR','response0')]:
        original=next(x for x in rows if x['cancer']=='ALL' and x['outcome']==group)
        kept=data['main'][field];assert kept<=original['retained_units']
        n=original['original_clinical_eligible_units']
        rows.append(dict(source='GSE243013',cancer='ALL',stage='FINAL_AGENT_GROUP_GE10',outcome=group,
            original_clinical_eligible_units=n,retained_units=kept,excluded_units=n-kept,retained_fraction=kept/n))
    out=task/'selection_review_v1';out.mkdir(exist_ok=False)
    result=dict(status='COMPLETE_SOURCE_AGGREGATE_SELECTION_REVIEW_NO_NEW_INFERENCE',
        observed_UTC=datetime.now(timezone.utc).isoformat(),rows=rows,
        input_SHA256={k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in paths.items()},
        denominator='Original strict ICI and complete clinical/pathology author materials, not reviewed independent patients',
        no_new_P_CI_or_weight_or_state_refit=True,
        no_histotype_agent_denominators_inferred_from_point_sensitivity_subsets=True)
    (out/'private_source_capture_selection_review.json').write_text(json.dumps(result,indent=2)+'\n')
    lines=['# NSCLC：固定状态临床分析的来源级筛选轨迹', '',
           '只复用已接受的来源级覆盖，不恢复个体，不新增P/区间/选择权重或拟合。单位是作者材料，不是已审独立患者。', '',
           '| 癌种 | 阶段 | 病理标签 | 原临床合格 | 保留 | 排除 | 保留率 |', '|---|---|---|---|---|---|---|']
    for r in rows:
        lines.append(f"| {r['cancer']} | {r['stage']} | {r['outcome']} | {r['original_clinical_eligible_units']} | "
                     f"{r['retained_units']} | {r['excluded_units']} | {r['retained_fraction']:.1%} |")
    lines += ['', '## 能说明什么', '',
              '主分析是“临床字段完整、双方父群充分捕获、药物类别有覆盖”的条件人群。'
              '如果两类病理结局的保留率不同，其结果不能无条件推广到全部术后或全部治疗前病例。'
              '这是可核实的选择轨迹；不直接证明选择偏倚方向，也不证明已观察的状态关联是伪影。', '',
              '调整RNA/细胞数/父群比例有助于描述相同测量背景下的关系，不能重建未捕获群，也不能自动修复未知重复、取样或选择问题。'
              '术后CD8增殖与non-MPR关联仍可有研究价值，但残留抗原/组织组成是待检验解释，不能直接命名为促耐药细胞。', '',
              '分层敏感性在癌种内重新执行药物组覆盖门槛，不能拿其分母拼回主分析。'
              '本次完整表没有癌种×药物门槛后的主分母，保持未知。'
              '原六状态/九组合结果及校正家族保持，不以本检查重新选择队列/阈值。', '']
    (out/'REPORT_zh.md').write_text('\n'.join(lines))
    print(json.dumps(result))


if __name__=='__main__':main()
