"""Consolidate real evidence on server. Never print scientific payloads.

Report is for source-cohort interpretation, not a validated interaction claim.
"""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def block(frame):
    return '\n```text\n'+frame.to_string(index=False)+'\n```\n'


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1'
    out=crc/'private_scientific_review_v1'
    assert not out.exists()
    transfer=crc/'private_existing_source_transfer_v1'
    paths=dict(holdout=crc/'heldout_patient_results_v1.tsv',
               clinical=crc/'private_MMR_cohort_test_v1/all24_private_results.tsv',
               nuisance=crc/'private_program_interpretation_v1/all36_private_pathway_nuisance_diagnostics.tsv',
               source=transfer/'all8_private_source_consistency.tsv',
               coverage=transfer/'private_capture_and_gene_coverage.tsv')
    frames={k:pd.read_csv(p,sep='\t') for k,p in paths.items()}
    s=frames['source'];assert len(s)==8
    valid=s.p.notna();q=multipletests(s.p.fillna(1),method='fdr_bh')[1]
    assert np.allclose(s.loc[valid,'BH_eight_slot_family'],q[valid],atol=1e-12)
    assert s.loc[~valid,'BH_eight_slot_family'].isna().all()
    assert (s.loc[valid,'n_patients']>=10).all()
    assert s.loc[valid,'rho'].between(-1,1).all()
    nuisance=frames['nuisance'];assert len(nuisance)==36
    clinical=frames['clinical'];assert len(clinical)==24
    actions=[]
    for program in ['MCP1','MCP2']:
        same=s[s.program.eq(program)]
        eligible=same[same.p.notna()]
        supported=eligible[(eligible.rho>0)&(eligible.lower>0)&(eligible.BH_eight_slot_family<.05)]
        c=clinical[clinical.program.eq(program)&clinical.partition.eq('holdout')]
        clinical_support=c[c.p.notna()&(c.BH_four_slot_family<.05)&(c.lower*c.upper>0)]
        discovery=clinical[clinical.program.eq(program)&clinical.partition.eq('combined_exploratory')]
        actions.append(dict(program=program,transfer_estimable_sources=len(eligible),
                            positive_interval_BH_sources=len(supported),
                            heldout_MMR_supported_records=len(clinical_support),
                            combined_exploratory_MMR_BH_records=int((discovery.BH_four_slot_family<.05).sum()),
                            proposed_action='spatial_or_external_cohort_validation_priority'
                            if len(supported)>=2 or len(clinical_support)>0 else 'retain_exploratory_program_do_not_promote_to_mechanism',
                            criteria_are_posthoc_evidence_triage_not_new_hypothesis_tests=True,
                            biological_or_clinical_value_verified=False))
    sex=nuisance[nuisance.reference.eq('sex_linked_nuisance_check')&(nuisance.exploratory_BH<.05)]
    stress=nuisance[nuisance.reference.eq('immediate_early_stress_nuisance_check')&(nuisance.exploratory_BH<.05)]
    out.mkdir()
    body='# CRC固定细胞对：表达程序、MMR与来源转移\n\n'
    body+='此文件仅在授权服务器内部，含真实患者派生统计和模型级基因线索。未授权公开或导出。\n\n'
    body+='## 研究范围\n\n原PAIR015：CD8_T_cell_NKG7_9—Macrophage_IFI6_18，连接CM01/03。只代表这一接口，不能替代全部9CM。原细胞身份、编号、原Zarr及历史阳性/阴性均保留。\n\n'
    body+='## 1. 固定患者留出\n\n发现与留出只针对新MCP拟合；原研究及整个CM发现并非独立。MCP2若没有双方基因支持，应保留缺失，不把CCA分数当第二个基因程序。\n'
    body+=block(frames['holdout'])
    body+='## 2. 程序的功能与共同因素\n\n这是未加权签名重叠富集，背景为每种类型发现阶段候选基因；不是因果、蛋白结合或正式细胞状态。干扰素alpha/gamma重叠不能直接判定IFNG特异性。\n'
    body+=block(nuisance.sort_values('exploratory_BH'))
    body+='性别线索需原作者真实性别元数据；即时早期/应激可能是生物响应也可能包含制备因素。不能靠删掉不利基因后重做程序来恢复显著性。\n'
    body+='性别集合富集需审阅：'+str(not sex.empty)+'；即时早期/应激集合富集需审阅：'+str(not stress.empty)+'。\n\n'
    body+='## 3. MMR队列意义\n\nMMRd/MMRp是分子分型，不是治疗响应。患者组别不足、协变量缺失或共线记录为不可估计。病理T阶段按1–4的线性趋势处理，不称完整临床分期；左右部位非精细位置。无真实性别/年龄，残余混杂未消除。发现/合并只作为探索，不能替代留出。\n'
    body+=block(clinical)
    body+='## 4. 固定来源转移\n\n四个来源全部保留；基因不完整不填零；先按双方原细分各20细胞选可比库，再每位已审来源患者取第一个库，最多各50细胞。CRC发现归一化和质量系数原样转移，跨癌种/平台的校准外推仍有限。癌种限定的意义可以保留。\n'
    body+=block(frames['coverage'])+block(s)
    body+='来源均曾用于原CM探索。这里检验的是未用于CRC MCP拟合的来源转移，不是完全独立临床验证。两个来源支持同一MCP也不能直接说明某一LR主导。丰度/细胞数调整只给描述性partial rho，无新的显著性检验。\n\n'
    body+='## 5. 下一步证据分流\n'+block(pd.DataFrame(actions))
    body+='上述分流是后验研究优先级，不增加统计发现。真正空间验证必须核坐标、全基因/面板、细胞分割、患者/切片单位、程序与完整受体的可测性，并与切片留出和距离打乱比较。HALO图像特征与RNA不同，不直接输入全转录组算法。旧空间45/70负检验不原样重跑。\n\n'
    body+='## 6. 方法与原作者资料\n\n'
    body+='- CRC空间免疫枢纽：https://portals.broadinstitute.org/crc-immune-hubs/\n'
    body+='- DIALOGUE：https://pubmed.ncbi.nlm.nih.gov/35513526/\n'
    body+='- NCEM：https://www.nature.com/articles/s41587-022-01467-z\n'
    body+='- SpaTalk：https://www.nature.com/articles/s41467-022-32111-8\n'
    body+='\n这些原文用于方法设计或同研究解释，不代替当前真实验证。全文没有确认新的耐药/抑制模块或直接通讯机制。\n'
    (out/'CRC_MCP_MMR_SOURCE_REVIEW_PRIVATE.md').write_text(body,encoding='utf-8')
    (out/'private_next_action.json').write_text(json.dumps(actions,indent=2)+'\n')
    inventory=dict(status='PRIVATE_SCIENTIFIC_REVIEW_GENERATED_NOT_VALIDATED_INTERACTION',
                   finished_UTC=datetime.now(timezone.utc).isoformat(),
                   input_sha256={k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in paths.items()},
                   source_family_missingness_audit='PASS_ALL8_FIXED_SLOTS',
                   publication_scope_pending=True,outputs_remain_server=True)
    (out/'private_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    print('PRIVATE_SCIENTIFIC_REVIEW_AND_SOURCE_STRUCTURE_AUDIT_COMPLETE_NO_DATA_EXPORT')


if __name__=='__main__':
    main()
