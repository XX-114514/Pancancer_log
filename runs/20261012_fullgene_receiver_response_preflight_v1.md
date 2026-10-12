# Run：原编号细胞内响应可比性预检

- Run ID：20261012_fullgene_receiver_response_preflight_v1；父科学证据v17。
- 预声明UTC 2026-10-12T02:58:50.466041+00:00 /UTC+8 2026-10-12T10:58:50.466041+08:00；PREDECLARED_NOT_STARTED。
- scanpy Python3.10.19；原336记录/13目录节点/3已审来源/固定stratum。
- [合同及SHA](../inventories/analysis/2026-10-12/fullgene_receiver_response_preflight_plan_v1.json)；不做关联系数、p值、TF活性或正式状态。

## 执行与验收

提供PANCANCER_PREFLIGHT_ROOT后依次运行：

```bash
python fullgene_receiver_response_preflight_v1.py
python audit_fullgene_receiver_response_preflight_v1.py
```

全部336必须回映，来源未捕获原节点保留；旧病例块不重分配。九列RNA/count表达设计仅诊断，不是最终调整模型。独立audit以逐细胞原生cache与完整RNA分母重建L/R均值与计数，QR验秩及所有靶/块/映射。私有患者/样本行不导出。完成后记录真实起止/exit及SHA；当前未派发。

## 2026-10-12T11:15:24.098787+08:00 实际完成

生产03:01:16—03:01:54 UTC、独立audit03:01:54—03:02:16 UTC均exit0；执行合同commit7e7ef92187d618bbf601c937d23e06745fc13290预先云端核验。248接口/620原编号记录、336全回映、496设计原生cell+QR独立复算PASS。325计数/变异/基础rank通过，不代表完整模型、显著关系或轴验证。[结果及SHA](../inventories/analysis/2026-10-12/response_preflight_v1/independent_audit.json)。源未捕获/低覆盖/缺先验均保留，没有正式状态或canonical变动。
