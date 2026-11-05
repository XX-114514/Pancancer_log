# V7 panel repair v1 验证证据

正式回放：[replay_v2/manifest.json](replay_v2/manifest.json)。旧40群范围1073743细胞，生成4480条候选证据；T6最新12群范围156631细胞，生成1344条候选证据。范围相互重叠，不能加总。T6仅提供top10 DEG，未出现基因不解释为差异阴性。

- [40群全部证据](replay_v2/old40_all_panel_evidence.tsv.gz)
- [T6全部证据](replay_v2/T6_12_all_panel_evidence.tsv.gz)
- [T6选定实测marker](replay_v2/T6_selected_observed_markers.tsv)
- [服务器合成测试退出码与独立审计](server_smoke_v2.json)
- [测试/资源/代码核验](verification.json)
- [完整panel与文献](../../../../resources/annotation/v7_panel_repair_v1/README.md)

当前目录根部同名回放文件为首次聚合回放历史快照；最终正文引用replay_v2。新版Python证据计算最后仅增加说明性注释，最终文件哈希与实际回放时的代码哈希分别保留。它们不是新的全矩阵生物学分群，也没有更改细胞标签。

服务器合成测试2026-11-05 06:12:28–06:12:47 UTC / 14:12:28–14:12:47 UTC+8完成，fixture/workflow/verify均exit0。4个测试QC排除后356细胞重建，1227条表达统计及6条Treg同细胞证据独立验算通过。R的可选Harmony/UCell分支未在本轮测试运行；真实V7全体重新注释尚未执行。
