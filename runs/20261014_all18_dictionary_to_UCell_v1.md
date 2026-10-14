# 全18字典到全来源UCell：受控顺序执行

等待正在运行的固定24个诊断单元全部完成，之后顺序执行：全18字典v4组装 → 独立重建全部签名验收 → 4532933细胞全来源原生UCell评分。等待上限8小时，各字典阶段上限1小时，评分上限24小时；异常停止并保留记录，不自动重跑模型。

v4在尚未组装前纠正v3审阅标记：自身谱系的角蛋白、胶原、mural或神经面板不能直接列为其他谱系RNA。依据每个allowed major分别给出自身身份面板和需审阅面板，shared程序也逐大类记录。6个自身/异谱系及shared边界检查通过；未改变基因签名、稳定门槛、模型或旧结果，v3保留。

字典验收从原cNMF系数/TPM与三种子Spectra权重独立重建全部签名；确认810个cNMF候选及全部合格Spectra候选、原18大类范围、原生支持和谱系标记。验收仅表示可开始测量，不表示生物状态已确认。完整评分使用已核容量方案，独立全量评分审计另行完成。

[字典计划](../inventories/analysis/2026-10-14/all_lineage_state_candidate_dictionary_plan_v4.json)、[评分计划](../inventories/analysis/2026-10-14/all_source_ucell_full_projection_plan_v1.json)、[容量预检](../inventories/analysis/2026-10-14/all_source_ucell_stream_preflight_v1/RESULTS.md)。

执行入口从PANCANCER_SCANPY_PYTHON读取已核对运行环境路径，服务器与公开代码一致。PANCANCER_PREFLIGHT_ROOT、PANCANCER_STATE_MODEL_OUTPUT_ROOT和PANCANCER_ANNOTATION_ROOT继承已运行任务。当前仅派发请求，真实进程和阶段另记，不把等待阶段称为评分已启动。
