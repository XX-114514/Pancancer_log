# 20261105_server_analysis_audit_v1

- 类型：既有分析的回顾取证、代码归档与交接；本次未运行新的生产分析。
- 取证：2026-11-05 00:49:42 UTC / 08:49:42 UTC+8起；事件时间沿用各原run收据，未记录的起止和退出码保持unknown。
- 基线：main `e6cf388ff5ca53c4e00cd4fcaa7559d72028a71e`。
- 范围：现有分析仓库、独立investigations、早期Xenium试点、指定PC_ROOT和labshare迁移收据；不扫描其他用户或项目。
- 环境：证据采集/文档生成Python标准库；既有scanpy环境Python3.10.19 / numpy1.26.4 / pandas2.3.3 / scipy1.11.4用于小型W尺度复核及合成测试；R基础parse检查。不是历史生产环境替代。
- 状态：科学梳理及轻量验证结果见[任务日志](../logs/codex/2026/2026-11-05_server-analysis-systematic-audit.md)；GitHub传输需独立核验。

## 输入与命令

```bash
python scripts/server_audit/collect_snapshot.py \
  --pc-root "$PC_ROOT" --workspace "$WORKSPACE" --source-repo "$CODE_ROOT" \
  --output inventories/server_analysis/2026-11-05

"${PC_ROOT}/envs/conda/scanpy/bin/python" -I \
  scripts/server_audit/check_nmf_scaling.py "$NMF_RUN" "$REFERENCE_W" "$NEW_SCALE_AUDIT_TSV"
```

两个采集入口均拒绝已有输出。实际运行在独立整理副本中；本次表格随后补录代码来源hash、迁移问题表汇总、配置及侧表存在情况，原/公开字节hash均登记。不会向原冻结运行回写结果。

## 交付与主要发现

- [系统分析报告](../docs/reports/20261105_server_analysis_audit/README.md)，覆盖数据版本、来源筛选、CNV、TME、去批次、CT75、V5/V7、blisa、NMF、候选注释及交接。
- [证据目录](../inventories/server_analysis/2026-11-05/)：主证据206项、补充69项、43个归档代码文件来源及对应公开hash、分析线/版本/状态表。
- 两个V7 rank survey和legacy分支已有完成收据；K9由override固定，不是rank survey发现。25µm C05现有摘要。
- 72来源恶性侧表均存在，未逐细胞重验；全量marker分流无完成收据。
- 新TME重整合日志到Harmony完成，最终验收缺；Scanorama已有图结果而缺scIB总分行。
- V7空间弱信号依赖W top10、来源/单位及身份识别限制，不构成CM空间验证。
- 新增NMF等价尺度敏感性复核：原余弦受尺度影响，规范后部分匹配改变；未重拟合或改变原生产结论。

## 验证与未验证

4个调度编排mock测试、3个合成收集测试通过；归档Python/shell/R语法通过。快照一致性、仓库validator、release validator和diff检查以任务日志最终实测状态为准。本次没有重hash大型矩阵、验证完整bin级p检验族、重跑身份准确率、联系旧服务器或验证Slurm最终退出；缺失不能用目录存在、输出时间或代码存在代替。

## 下一步

先验收TME run_v2，后冻结身份候选流程；NMF按规范尺度＋rank复核再选择独立验证；空间补测量/患者资格。详细优先级见系统报告和TODO顶部。保留历史失败与当前并发修改；不强推、不改历史。
