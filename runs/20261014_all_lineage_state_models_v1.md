# 全18大类状态模型：固定范围与资源块

[输入已接受](../inventories/analysis/2026-10-14/all18_native_inputs_accepted_v2/RESULTS.md)。[Spectra冻结计划](../inventories/analysis/2026-10-14/all_lineage_spectra_model_plan_v1.json)、[cNMF冻结计划](../inventories/analysis/2026-10-14/all_lineage_cnmf_model_plan_v1.json)。

Spectra官方GPU固定commit，原身份明确进mask。6预先固定资源块（T/lymphoid、myeloid、B/mast/hematopoietic、vascular/stromal、epithelial/melanocytic、glial），各3seed；HVG总预算2000按major均分+全部可测process priors。GPU官方use_cell_types=False明确未实现；18类型p²邻接内存开销较大，分块不改算法或原CM。block_global仅指允许作用在块内，不表示所有major真的活跃；共享资格至少2major支持/变异，其他major结果保留，不以所有谱系都活跃为生物必要。单major Glial没有跨type共享资格。

cNMF同100766原major逐个拟合，K10/15/20×100，共5400runs/54共识；官方raw variance normalization，native fullRNA TPM。无先验、无LR/临床p选择，全部K/seed保留；跨K重复不视为独立状态。隔离固定runtime，CPU每job2thread，GPU串行模型；不改基础环境。

模型共享5GiB预算/500GiB留空，开始前检查500+5，各阶段500+1；Spectra实际features后12×scope×p²×4及输入/参数保守张量预算<=实际freeGPU70%，否则显式停止并保留记录。旧T/M实测cNMF385320KiB、Spectra81252KiB用作比例参考，5GiB属于预算非精确保证。模型成功仍需多seed/crossK/source/node/RNA/细胞支持审核和固定全数据投影。

入口train_all_lineage_spectra_v1.py / train_all_lineage_cnmf_v1.py，环境PANCANCER_PREFLIGHT_ROOT及PANCANCER_ALL_LINEAGE_SPECTRA_PLAN/PANCANCER_ALL_LINEAGE_CNMF_PLAN，使用同一既有模型venv。只有独立输入accepted manifest精确SHA后运行。完整goal不以模型fit完成。

## 真实启动

2026-10-14T04:46:53.817940+00:00 /UTC+8 2026-10-14T12:46:53.817940+08:00取证：04:44:30 UTC启动两个job、04:45:58四proc一致，首T/lymphoid/首B训练中，无完成seed或major。[快照](../inventories/analysis/2026-10-14/all18_state_model_start_v1/RESULTS.md)。
