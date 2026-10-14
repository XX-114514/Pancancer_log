"""Independent all-nine-CM scalar/family audit and server-private report."""
import json,os
from pathlib import Path
from datetime import datetime,timezone

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out=root/'CRC_all9_frozen_CM_MMR_20261014_v1'
    results=pd.read_csv(out/'all18_private_module_MMR_results.tsv',sep='\t')
    coverage=pd.read_csv(out/'all9_private_capture_coverage.tsv',sep='\t')
    assert len(results)==18 and len(coverage)==9
    assert not results.duplicated(['cm','method']).any()
    replayed=0
    for _,group in results.groupby('method'):
        assert len(group)==9
        q=multipletests(group.p.fillna(1),method='fdr_bh')[1]
        valid=group.p.notna()
        assert np.allclose(group.loc[valid,'BH_nine_slot_family'],q[valid],atol=1e-12)
        assert group.loc[~valid,'BH_nine_slot_family'].isna().all()
    for row in results.itertuples():
        if not np.isfinite(row.p):continue
        selected=pd.read_csv(out/(row.cm+'_private_capture_selected_patients.tsv'),sep='\t')
        assert selected.analysis_block.is_unique and selected.sampleID.is_unique
        columns=['activity','MMRd']+(['right_site','T_stage','log_total_cells'] if row.method.endswith('HC3') else [])
        selected=selected.dropna(subset=columns);g=selected.MMRd.to_numpy().astype(bool);y=selected.activity.to_numpy()
        assert len(selected)==row.n_patients and g.sum()==row.n_MMRd and (~g).sum()==row.n_MMRp
        if row.method=='primary_rank':
            u=mannwhitneyu(y[g],y[~g]).statistic
            effect=2*u/(int(g.sum())*int((~g).sum()))-1
            assert abs(effect-row.effect)<1e-12
        else:
            X=np.column_stack((np.ones(len(selected)),g.astype(float),selected[['right_site','T_stage','log_total_cells']].to_numpy()))
            ref=sm.OLS(y,X).fit(cov_type='HC3',use_t=True)
            assert abs(ref.params[1]-row.effect)<1e-8
            assert abs(ref.pvalues[1]-row.p)<1e-8
        replayed+=1
    dest=out/'ALL9_CM_MMR_PRIVATE_REVIEW.md';assert not dest.exists()
    body='# 全9个原CM：CRC/MMR队列意义审阅\n\n本报告仅在服务器内部，统计公开及进入对话须另外明确授权。\n\n'
    body+='## 原输入与捕获范围\n\n原172身份、9CM与完整1290样本H保持。仅原MAJ158未分选原发CRC的已审来源患者。每模块先按原图核心权重、细胞数、大类和组织分区筛捕获机会，再选首个合格库；不同模块的可检验患者分母可能不同，不能直接比较最小P。\n'
    body+='\n```text\n'+coverage.to_string(index=False)+'\n```\n'
    body+='## 真实结果\n\n主患者秩效应方向MMRd−MMRp。另列部位/病理T趋势/总细胞数调整的HC3探索结果。完整9槽BH、2000患者bootstrap；不可估计缺失不补阴性。\n'
    body+='\n```text\n'+results.to_string(index=False)+'\n```\n'
    body+='## 可以与不能支持的结论\n\n这些是已用于发现的来源内生态分型关联，不是独立CM恢复、真实治疗响应预测、空间共定位或LR因果。总细胞数调整不等于消除全部采样、组成、性别/年龄或技术因素。病理T趋势不是完整临床分期。若差异仅在某个CRC模块成立，其癌种限定价值应保留，再找未用于发现的患者队列或空间区域作验证。若无支持，不改标签/阈值追显著；与当前完整状态程序和新免疫队列衔接。\n'
    body+='\n## 审计\n\n全部18键、2个9槽家族及缺失门槛已核；所有可估计主秩效应用SciPy U统计复算，调整系数/P用Statsmodels HC3复算。数值通过不等同生物学发现。\n'
    dest.write_text(body,encoding='utf-8')
    check=out/'private_independent_audit.json';assert not check.exists()
    check.write_text(json.dumps(dict(status='PASS_ALL18_KEYS_TWO9SLOT_FAMILIES_ALL_ESTIMATED_SCALAR_REPLAY',
      audited_UTC=datetime.now(timezone.utc).isoformat(),estimated_replayed=replayed,results_remain_server=True),indent=2)+'\n')
    print('ALL9_PRIVATE_CM_MMR_INDEPENDENT_AUDIT_AND_REPORT_COMPLETE_NO_RESULTS_EXPORTED')


if __name__=='__main__':
    main()
