# Run：患者块证据状态覆盖表

- Run ID：20261012_patient_block_provenance_overlay_v1；父阶段v17完整基因准备。
- 预声明UTC 2026-10-12T02:38:40.719500+00:00 /UTC+8 2026-10-12T10:38:40.719500+08:00；状态PREDECLARED_NOT_STARTED。
- 范围1290行、121frozen/87author/194parent/888fallback；402文库对应175来源内患者块。
- 环境：实际执行将使用scanpy Python3.10.19；运行前commit以本合同普通push后的内容SHA为准，当前尚未派发。

## 输入、方法与验收

[冻结合同](../inventories/analysis/2026-10-12/patient_block_provenance_overlay_plan_v1.json)固定6输入SHA、producer/auditor SHA。全部旧列/行序/analysis_block必须不变，补有效证据字段；194旧状态逐键对应已审官方侧表。不把文库analysis_unit计患者，175来源内块不称完全跨GSE独立人群。私有1290行保留服务器，只发布来源级计数与SHA。

## 可复现执行结构

以PANCANCER_PREFLIGHT_ROOT提供根目录，固定解释器依次执行：

```bash
python build_patient_block_provenance_overlay_v1.py
python audit_patient_block_provenance_overlay_v1.py
```

真实成功终态及独立全行复核才验收；不重分配原块、改注释、拟合统计、发布轴。实际起止/exit尚未发生。

## 2026-10-12T10:48:39.301099+08:00 实际完成与独立验收

production02:42:52.462067—02:42:53.096537 UTC、audit02:42:53.096573—02:42:53.604247 UTC均exit0；执行前合同提交ad4ba3c90fe93f6c08035993296bd03eb0a2f4f6已云端独立核验。1290全行原字段/块未改、194官方父映射/87作者/121冻结状态均逐键复核，402文库/175来源内块。首次运输驱动换行SyntaxError未执行任何生产；修正仅本地driver，不变冻结科学代码。[SHA及独立审计](../inventories/analysis/2026-10-12/patient_block_provenance_overlay_v1/independent_audit.json)。仅metadata准备，正式状态/轴/临床机制未验证。
