# 自动通讯张量与完整留出临床家族

2026-10-18T19:32:32.042060+00:00 UTC / 2026-10-19T03:32:32.042060+08:00 UTC+8 效果前冻结；[前置资格](NSCLC-author-lineage-and-communication-factor-preflight-20261019.md)。

## 问题与解释

不人工挑选通讯轴，先从全合格亚类和完整可测LR资源学习共同表达变化，再检验固定因子在模型学习留出记录中的病理反应关联。这不是原2026年10月9CM复现，也不是完全未触碰的临床验证。

## name

2026-10-19 communication tensor programmes v1 GSE243013

## scores

log1p(min ligand component mean_logCP10K * min receptor component mean_logCP10K), all native cells per author unit; expression cooccurrence proxy, not full LIANA consensus scores or binding/flux

## scope

Frozen129discovery/87heldout and38eligible author types; full metadata-primary4520LR vocabulary, every position retained in support manifest. Global LRsupport: all ligand components detected>=10percent in at least one available sender AND all receptor components in at least one available receiver within same discovery context, >=ceil129/3 contexts. No heldout or outcome used. All measurable per-context values retained after global filter, including true zero/subthreshold detection.

## missing

Unit absent or fewer20cells: entire corresponding sender/receiver slice masked missing. Every retained LR native component required. Diagonal sender=receiver entries retained as within-population interfaces, not demonstrated single-cell autocrine.

## scaling

Only one global RMS over observed discovery tensor entries; no per-LR/clinical feature scaling. Exact same scalar for heldout. Noncontext factor columns each unit L2; absorb product of their norms and CP weights into context coefficients. Factor numbering decreasing discovery context column L2, deterministic tie order.

## fit

Official Tensorly0.9.0 non_negative_parafac, same core called by cell2cell0.8.0; random init seeds2026101920,2026101921,2026101922 for each rank2..6. Exactly400iterations/tol0, float32 GPU, TF32off, no outcome convergence/seed/rank selection. Retain every15model and choose minimum observed relative reconstruction error per rank. No single selected rank; all20factor clinical slots enter main BH/BY20.

## stability

Hungarian match each same-rank restart using mean cosine similarity of four normalized factor axes; report every matched similarity, median>=.8 descriptive stability criterion. No clinical test filtered by stability.

## projection

Fixed LR/sender/receiver factors. Per heldout context compute observed NNLS through Gram=(B.T@B)*(U_available.T@U_available)*(V_available.T@V_available), rhs exact scored tensor contraction. Cholesky normal equations converted to scipy NNLS, independent fixture equality required. >=3available author types, positive observed score norm, full Gram rank and sqrt(condition)<1e8. Unidentifiable remains missing. No relearning factors or rank on heldout.

## clinical

After all15fits and all20heldout projections saved: original pCR/MPR1 vs nonMPR0; postneoadjuvant surgery association only. Activity=log1p(context coefficient / median positive discovery coefficient), centered/scaled using discovery activity mean/SD only. Continuous activity is model readout; pathology response binary target regressor. Main covariates original8clinical + canonical agent fixed effects + log1p all immune cells + all-cell meanlogRNA/meanloggenes + sender/receiver factor-L1-weighted state fractions of all source immune cells. Complete nuisance column space retained and rank/span reported; no outcome/clinical P variable selection.

## clinical_gates

At least50complete author records, >=10each response, >=10residual df, target identifiable, nuisance and complete design stable/condition<1e8/leverage<.999999. Heldout87nonoutcome profile blocks eachsingleton; no reviewed independent-patient alias assertion. If future non-singleton block discovered, stop before clinical inference rather than sample-level randomization.

## statistics

Reuse previously independently accepted NSCLC cohort HC3 restricted-null Rademacher wild t,19999draws; two-sided add-one P. Fixed2000within cancer x canonical agent resamples, report valid; percentile CI only if>=1900valid, no resampling until desired result. All20main slots BH and BY with unestimable P=1 included. Fixed500seedoffsets defined in code/hash; no seed search. Full random vectors persisted privately for independent replay.

## sensitivity

All20unadjusted group contrasts, original8+agent clinical-only, clinical+RNA/count without factor abundance: point effect and asymptotic HC3 CI only, no extra significance selection. All20unadjusted two-sided Mann-Whitney U (Wilcoxon rank-sum) with tie/asymptotic correction, separate BH/BY20 as secondary only. Leave-one-author-record main effects; no claim independent-patient leaveout.

## limits

Source previously explored; heldout is separation of model learning, not untouched external validation. Immune-only labels are original author subtypes, not new unified formal states. Distinct from 2026 October9CM and 2026-10-18 RNA MCP. Original closed clinical/function families not reopened. No direct interaction, resistance mechanism or baseline prediction claim.

## resources

No downloads. Tensor CPU/GPU RAM only, no full tensor persisted. Max CPUavailable>8*tensorbytes+4GiB, GPUfree>6*tensorbytes+2GiB; otherwise explicit stop/no feature downselection or silent fallback. Output cap32MiB; group reserve500000000000+1073741824 maintained.

## 统计尺度

临床主模型为连续程序活动 ~ 原作者反应二元标签 + 完整控制空间。活动使用发现集系数的正值中位数作尺度，log1p后使用发现集均值/SD冻结标准化。正效应表示pCR/MPR记录较高，负效应表示nonMPR较高；不是病理反应风险比。主CI为固定2000次分层记录bootstrap的百分位CI，P为19999次HC3限制零假设wild-t；不是Wilcoxon的调整后版本。三个背景敏感性只报告效应及渐近HC3区间，Wilcoxon另校正20个次级问题。

保留87条均有原作者反应，87个非结局指纹均单条，但独立患者身份仍未审定。固定原8临床+药物+全来源数量/RNA质量的17列设计在因子拟合前已核满秩、条件数7.42。两端因子加权比例需在其无结局拟合后逐项核可识别；完整控制空间保留，不靠删协变量追求P值。

## 执行与验收

程序先写本方案/源码SHA，再等待全1840基因×11701单位原始复读接受。实际float32 CP和低内存Gram→Cholesky→NNLS在真实表达前另核合成边界；完整已接受软件门槛此前使用float64/direct NNLS。全部15模型、每起点误差/配对稳定性和全20投射先保存后才读临床结局；所有学得因子/病例行/效应私有。完整独立审核须检查资源支持、mask/尺度、各起点排序、固定投射、临床输入和全部随机流/BH/BY。当前启动是依赖等待，不是已完成真实拟合。
