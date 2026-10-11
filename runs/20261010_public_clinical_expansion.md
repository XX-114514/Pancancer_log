# Run：新增公共临床单细胞收集

- ID：20261010_clinical_scrna_expansion_v1；范围：Project_v5
- 类型：public downloads + source audit + per-source preprocessing
- 运行代码commit：not_recorded；实际runner和manifest在外部run目录，本次只读取证
- 环境：服务器scanpy Python3.10；部分来源R/qs2支持，完整环境未在本次冻结
- 起始：2026-10-10，精确启动时间未统一核验；终态：尚未完成
- 本次证据：[快照与报告](../inventories/runtime/2026-10-11/REPORT.md)

## 范围与门槛

公共搜索/下载优先，排重基线4,676,787细胞/54GSE，包括未进入下游者。验收目标是QC去重后新增450万—550万；不是下载文件细胞数。人类scRNA/snRNA分开保留；空间、动物、细胞系及bulk不计。

当前11来源库存5,341,204源细胞；三份QC分片608,684，待最终来源验收；五份原始计数Zarr2,345,161细胞只做结构、计数及质量指标读回，不等于最终QC。头颈癌逐源计数转换进程仍在。

四个大型来源全矩阵计数/符号/barcode精确匹配基线均为0，但基因覆盖变化及别名可能隐藏重叠；来源和文库排重仍未关闭。TNBC/CTCL、CAR-T治疗产品及临床缺失边界详见报告。第二批及Zenodo未完成者不计验收。

## 命令与输出

此次取证命令、快照schema、逻辑路径与SHA见[RUNTIME_SYNC](../docs/RUNTIME_SYNC.md)及机器快照。原始下载证明、来源审计、accepted分片、计数Zarr及可读日志留在 `${PANCANCER_ROOT}/Project_v5/runs/20261010_clinical_scrna_expansion_v1`。

本记录不是完整科学重跑合同：原始worker参数、调度历史、资源峰值及每份来源环境还需在正式验收时锁定。跨研究整合未启动；没有以metadata、疑似重叠或未确认结局填补规模。

## 2026-10-11 18:56—19:02 UTC+8 复核节点

18:56:11取证完整源库存仍5,341,204；GSE286827新增QC384,628，四份合计993,312。头颈癌QC完成事件为03:42:39 UTC/11:42:39 UTC+8。19:02:09开始只读流式核对全部1,034,264原始账本行、993,312保留行、149文库和Zarr X/counts结构，全部检查通过；不重算生物过滤，不授予最终来源验收。

下载worker已退出；第二批8/14文件校验，Zenodo39%无完成证明。原queue_terminal与v6监控不能代表成功/运行。失败记录时间与SHA见[故障证据](../inventories/runtime/2026-10-11/download_failure_evidence_v1.json)。下一阶段准入见[报告](../inventories/runtime/2026-10-11/DOWNLOAD_READINESS.md)。本节点仅审计；没有恢复传输、启动QC新批次或跨研究整合。

复现账本核对：`python scripts/audit_qc_readiness.py --run "${RUN_ROOT}" --output "${NEW_AUDIT_JSON}"`。脚本使用stdlib，目标scanpy Python3.10已真实执行，新的输出使用排他创建。各来源ledger/threshold/proof SHA在汇总JSON；数值矩阵读回为既有producer凭证，本轮未重新加载全部表达矩阵。
