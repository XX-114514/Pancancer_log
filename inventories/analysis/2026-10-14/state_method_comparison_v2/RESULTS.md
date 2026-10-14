# 完整cNMF与Spectra对照：真实结果及下一诊断

## 已完成范围

原600分解/6共识02:03:23 UTC真实exit0，原worker/child均退出；对照02:04:09—02:04:25 UTC exit0/退出，02:08:15独立审计完成。2major×3K全部90候选，9000重复对照、1300crossK、14355Spectra匹配全部保留；原fit/consensus全SHA前后不变。独立重算9000cosine/density、57420跨方法cos/rank统计、28710gene重叠、360source/node/复杂度统计、1300crossKgene/rank和90规则通过。K选择统计只核键/finite，没有独立重新拟合六套完整重建统计，不能宣称全部统计独立重演。 [审计](independent_audit.json)。

## 实际结果

- 90候选中68通过重复门槛，69在另一个K有稳定对应，57两者同时通过。它们跨K重复，**不是57独立状态**。
- 六个K轮廓系数约.730—.802，原密度保留98.35%—99.90%；不同细胞谱系的误差不能直接比较。全部K保留，未按LR/clinical p择K。
- T候选中可见细胞毒GZMA/GZMB/PRF1、FOXP3/IL2RA/CTLA4相关、S/G2M、热休克程序；Mac候选可见IFIT/MX1干扰素相关、SPP1/CTSL/VIM及增殖程序。只描述多基因内容，未赋新type或刺激特异性。
- Spectra八个重复稳定候选与cNMF的最佳gene-only匹配cosine最高约.568，其余多更低；当前两方法没有强同一性证据。它们优化目标/特征尺度不同，不能由不吻合单独认定某方法错误。
- 一些cNMF候选在多个K都很稳定，却主要是胶原/角蛋白/神经相关RNA；T和Mac两个粗身份内都出现。算法稳定不能保证它就是正确免疫状态。
- sourceR2很高的IFN或GPR34/C3程序可能反映癌种/微胶质背景，也可能包含实验差异；来源限定仍有价值，不能全部自动标artifact。

完整[90目录](all90_cnmf_candidates.tsv)、[全14355匹配](all14355_Spectra_cnmf_comparisons.tsv)、[六K诊断](all6_K_diagnostics.tsv)。该候选库无先验分解成功，生物状态验收仍未完成。

## 数据问题与后续动作

新[同细胞谱系RNA诊断](../state_lineage_coexpression_plan_v1.json)覆盖全部240533原生T/M细胞/5source，保留原编号。检查多marker是否在同一细胞出现、占完整RNA比例、是否同时表达本谱系核心，以及全90模块top5%usage细胞内情况。可能解释包括身份误分、混合细胞/双细胞、环境RNA及巨噬摄取RNA；当前没有区分它们的证据，不直接重注释或删细胞。此前通讯困难可能部分源于这些输入混合，但这是待检验解释。

官方pyUCell隔离环境已合成验证，尚未用真实数据赋状态。完成上述核查及字典冻结，再用全原生gene axis评分/来源投影，映射原9CM并补其他谱系，状态后DIALOGUE/空间。完整goal active，原172/9CM/Zarr及历史结果不改。
