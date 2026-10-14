"""Private evidence card integrating strict, partial, nuisance and MMR results."""
import json,os
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests


def table(d):return '\n```text\n'+d.to_string(index=False)+'\n```\n'


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1'
    out=crc/'private_followup_evidence_v2';assert not out.exists();out.mkdir()
    strict=pd.read_csv(crc/'private_existing_source_transfer_v1/all8_private_source_consistency.tsv',sep='\t')
    observed=pd.read_csv(crc/'private_existing_source_transfer_observed_part_v1/all8_private_source_consistency.tsv',sep='\t')
    calibration=pd.read_csv(crc/'private_observed_part_calibration_v1/all4_private_score_calibration.tsv',sep='\t')
    nuisance=pd.read_csv(crc/'private_partial_source_covariates_v1/all20_private_descriptive_source_diagnostics.tsv',sep='\t')
    hallmark=pd.read_csv(crc/'private_all50_hallmark_diagnostic_v1/all300_private_function_statistics.tsv',sep='\t')
    mmr=pd.read_csv(crc/'private_MMR_cohort_test_v1/all24_private_results.tsv',sep='\t')
    cm=pd.read_csv(root/'CRC_all9_frozen_CM_MMR_20261014_v2/all18_private_module_MMR_results.tsv',sep='\t')
    assert len(strict)==len(observed)==8 and len(calibration)==4 and len(nuisance)==20 and len(hallmark)==300
    q=multipletests(observed.p.fillna(1),method='fdr_bh')[1];valid=observed.p.notna()
    assert np.allclose(observed.loc[valid,'BH_eight_slot_family'],q[valid],atol=1e-12)
    assert np.allclose(hallmark.BH_all300,multipletests(hallmark.p,method='fdr_bh')[1],atol=1e-12)
    body='# PAIR015 / CM01–CM03：后续证据卡\n\n仅在服务器保存；汇总对话审阅已获准，公开统计/基因名单未获准。\n\n'
    body+='## 测量层次\n\n严格完整基因转移与可测部分敏感性是不同结果，不能覆盖历史。部分测量保留原权重/均值/SD/质量系数，不重拟合；缺基因是未观察贡献，不是RNA=0。CRC校准只证实测量近似，不是外部生物学验证。\n'+table(strict)+table(calibration)
    body+='## 来源协调\n\n四来源及两个程序槽全部保留；原CM研究已使用这些来源，转移只独立于CRC MCP拟合。不同癌种有限定意义，无信号和不可估计同时保留。\n'+table(observed)
    body+='## 丰度与质量\n\n以下是后验描述性诊断，无新增P/BH，不能称完整调整后显著。所有八个协变量联合调整不可估计。留一患者只是敏感性，不能替代外部验证。\n'+table(nuisance)
    body+='## 功能\n\n完整50Hallmark×2端×all/up/down，300检验统一BH。以下展示BH<0.1；正式解释需同时看方向。T端仅全签名富集的集合不能解释成单向激活。巨噬负向富集不代表因果性免疫抑制或耐药。\n'+table(hallmark[hallmark.BH_all300<.1].sort_values('BH_all300'))
    body+='## MMR及原全部CM\n\nMMR不是治疗响应；CRC MCP发现/合并与患者留出分别处理。原CM v1因分类不完整没有可估计效应，不是阴性；v2补全解剖采样分区保留冻结身份，实算结果如下。\n'+table(mmr)+table(cm)
    body+='## 决策\n\n当前仅作为来源限定的跨类型表达协调候选，未验证直接LR、空间接口或临床预测。下一优先核现有空间面板的冻结权重覆盖/患者与身份，或独立免疫治疗队列基线材料/真实结局和可捕获身份。不能只挑有利队列、重学基因恢复显著、强拼新270群与旧172身份、把共识rank当FDR或非零H当原CM复现。\n'
    (out/'PAIR015_MCP_FOLLOWUP_PRIVATE.md').write_text(body,encoding='utf-8')
    (out/'private_audit.json').write_text(json.dumps(dict(status='PASS_PARTIAL8SLOT_BH_AND_COMPLETE300_HALLMARK_BH_PRIVATE_CARD_WRITTEN',
      finished_UTC=datetime.now(timezone.utc).isoformat(),direct_interaction_or_clinical_value_verified=False),indent=2)+'\n')
    print('PRIVATE_FOLLOWUP_CARD_AND_FAMILY_AUDIT_COMPLETE_NO_SCIENTIFIC_FILE_EXPORT')


if __name__=='__main__':main()
