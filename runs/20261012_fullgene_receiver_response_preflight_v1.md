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
