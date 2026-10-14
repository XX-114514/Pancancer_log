# 状态发现原生单细胞输入 v1

预声明UTC 2026-10-14T00:42:24.057206+00:00/UTC+8 2026-10-14T08:42:24.057206+08:00。固定5源/240533原cell/1832节点单位/212库/424major-unit，[合同](../inventories/analysis/2026-10-14/state_singlecell_input_plan_v1.json)。每源native fullgene counts/原身份/继承blocks及均衡训练H5AD，不赋状态。源轴交集16625，native RNA分母保存。

执行采用固定scanpy Python3.10，PANCANCER_PREFLIGHT_ROOT/ANNOTATION_ROOT/PANCANCER_STATE_INPUT_PLAN提供输入，`python assemble_state_singlecell_inputs_v1.py`。新目录拒绝覆盖；全部计数/positive/log聚合对照已审cache，独立H5AD读回与发现抽样/源投影检查后才训练。原current/core/source manifest及文件SHA固定，单源输出与完整inventory保存新目录。

状态：predeclared_not_started。代码语法通过，正式生产/读回待。完整goal的Spectra/cNMF/UCell/原9CM/DIALOGUE/空间均未完成；不将输入或模型健康作为最终完成。

## 2026-10-14T08:54:34.467285+08:00 实际故障和新读回合同

00:46:41—00:49:00 UTC parent exit1，五源和训练H5AD保留，library212误与sample-major424比较；未验收输入。新[v2审计](../inventories/analysis/2026-10-14/state_singlecell_readback_recovery_plan_v2.json)全面canonical raw/单位/ID/块/训练数据读回，不重复生成矩阵。隔离固定官方Spectra GPU 3epoch toy验证实际通过，正式生物训练尚未开始。
