# 固定通讯接口的条件空间邻近检验

## 预声明阶段

创建时间：2026-10-11T11:43:34.447250+00:00 / 2026-10-11T19:43:34.447250+08:00。状态：predeclared_not_started；父证据v14，未改注释v20261010/172节点/9CM或全发现BH。

问题：在同一视野、相同细胞家族位置、RNA深度与局部细胞密度条件下，固定L+发送者和R+接收者是否比随机标签更常邻近？范围为全部9个唯一接口×5文件=45主检验，回映13记录×5=65行；不因低覆盖更换端点。

[完整固定计划](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_plan_v1.json)。原生空间v3缓存/输入已独立审计；新目录保留所有历史。新检验用作者粗T/B谱系+背景口径，缺作者/无B类别/缺基因记NA；RNA纯谱系和count≥2另作描述敏感性。

几何采用原生像素和面积推算的圆半径，距离≤2倍双方半径之和；1/3倍半径敏感性。这个定义是有限邻域，不能称真实膜接触或绝对微米标定。仅同FOV，排除自接触。条件零模型在FOV×双家族归属×RNA深度二档×全QC细胞局部密度二档内打乱联合L/R两位标签，保留两端数量及同细胞共表达。固定2000次、seed20261011；文件内FOV合计，主单侧富集p按完整45家族BH，NA仅校正内部p=1，输出仍NA。

≥20/FOV原门槛是可测性描述，本轮不以此删掉稀疏可测关系；报告null方差、可移动标签数、FOV集中度/留一FOV。患者映射未知，不能把FOV或TMA当独立患者。[原补充表结构审阅](../inventories/analysis/2026-10-11/spatial_core_supplement_review_v1.json)没有完整FOV—患者桥；原Excel/RDS/细胞行不上传。

## 执行与验收

生产代码、审计代码和正式时间后续记录；现在没有运行结果。验收包括全45/65映射、源SHA、原始几何独立计数、端点/重叠/条件null全量复算与合成边界检查。审计PASS不代表患者级验证、蛋白结合、下游反应或原9CM驱动。APP–NOTCH2原文主要支持同细胞蛋白关联，HLA-DRA是完整MHC-II代理，其余结合/糖链/细胞方向限制完整保留在计划；不能把统计邻近直接写成免疫抑制机制。

可复现入口使用环境根目录；正式代码提交后补齐调用。关键/探索节点按持续[任务日志](../logs/codex/2026/2026-10-11_milestone-log-and-cloud-sync.md)普通push并独立核验SHA。

## 代码冻结检查

生产/审计代码[SHA及合成验收](../inventories/analysis/2026-10-11/spatial_conditional_code_freeze_v1.json)保存；两份代码独立实现图与统计量，未import生产者进行审计。服务器2000次合成重放、重叠/自接触/缺失均PASS。真实执行尚未开始。

```bash
PANCANCER_PREFLIGHT_ROOT="${RUN_ROOT}" python scripts/spatial_conditional_neighborhood_v1.py
PANCANCER_PREFLIGHT_ROOT="${RUN_ROOT}" python scripts/audit_spatial_conditional_neighborhood_v1.py
```

原生细胞缓存、FOV/置换细表只留服务器；公开45接口汇总。预声明已普通push并独立核验，实际代码提交随后冻结。
